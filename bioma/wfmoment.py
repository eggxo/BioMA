"""Configuration-driven 2-D deme WFmoments workflow."""

from __future__ import annotations

import configparser
import hashlib
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Callable, Dict, Optional, Sequence

from .gf_frequency import InputError
from .provenance import attach_input_fingerprints
from .runtime import subprocess_environment


def _path(value: str, base: Path, label: str, required: bool = True) -> Optional[Path]:
    value = (value or "").strip()
    if not value:
        if required:
            raise InputError("Missing [inputs] {}".format(label))
        return None
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _bool(value: str, default: bool = False) -> bool:
    value = (value or "").strip().lower()
    return default if not value else value not in {"0", "false", "no", "off"}


def _python(value: str) -> str:
    value = (value or "").strip()
    if value:
        return value
    return sys.executable


class WFMomentConfig:
    def __init__(self, path: Path):
        self.config_path = path.expanduser().resolve()
        parser = configparser.ConfigParser(interpolation=None)
        try:
            with self.config_path.open("r", encoding="utf-8-sig") as handle:
                parser.read_file(handle)
        except (OSError, configparser.Error, UnicodeError) as exc:
            raise InputError("Cannot read WFmoments configuration: {}".format(exc))
        for section in ("inputs", "analysis", "parameters"):
            if not parser.has_section(section):
                raise InputError("WFmoments configuration requires [{}]".format(section))
        inputs, analysis, parameters = parser["inputs"], parser["analysis"], parser["parameters"]
        base = self.config_path.parent
        self.current_raster = _path(inputs.get("current_raster", ""), base, "current_raster")
        self.pi_file = _path(inputs.get("pi_file", ""), base, "pi_file")
        self.structure_file = _path(inputs.get("structure_file", ""), base, "structure_file", False)
        self.param_file = _path(inputs.get("param_file", ""), base, "param_file", False)
        self.area_file = _path(inputs.get("area_file", ""), base, "area_file", False)
        self.future_masks_json = _path(inputs.get("future_masks_json", ""), base, "future_masks_json", False)
        self.output_dir = _path(analysis.get("output_dir", ""), base, "output_dir")
        self.species = analysis.get("species", "my_species").strip()
        self.compute_python = _python(parameters.get("compute_python", ""))
        self.rscript = parameters.get("rscript", "Rscript").strip() or "Rscript"
        self.nx = int(parameters.get("nx", "20"))
        self.ny = int(parameters.get("ny", "20"))
        self.threshold = float(parameters.get("threshold", "0.25"))
        self.min_valid = int(parameters.get("min_valid", "0"))
        self.loss_mode = parameters.get("loss_mode", "edge").strip().lower()
        self.direction = parameters.get("direction", "east_to_west").strip()
        self.migration = parameters.get("migration", "25").strip()
        self.migration_grid = parameters.get("migration_grid", "0.1,0.3,1,3,10,25,50").strip()
        self.fst_metric = parameters.get("fst_metric", "hudson").strip().lower()
        self.theta = parameters.get("theta", "auto").strip()
        self.z_gdar = parameters.get("z_gdar", "").strip()
        self.theta_probe = float(parameters.get("theta_probe", "1e-4"))
        self.time3 = float(parameters.get("time3", "3"))
        self.time5 = float(parameters.get("time5", "5"))
        self.mu = float(parameters.get("mu", "3.75e-8"))
        self.midterm_generations = parameters.get("midterm_generations", "auto").strip()
        self.replicates = int(parameters.get("replicates", "1"))
        self.seed = int(parameters.get("seed", "12345"))
        self.include_equilibrium = _bool(parameters.get("include_equilibrium", "false"))
        self.plot_equilibrium = _bool(parameters.get("plot_equilibrium", "false"))
        self.plot_raw = _bool(parameters.get("plot_raw", "false"))
        self.plot_title = analysis.get("plot_title", "").strip()
        self._validate()

    def _validate(self) -> None:
        for label, path in (("current_raster", self.current_raster), ("pi_file", self.pi_file)):
            if path is None or not path.is_file():
                raise InputError("{} does not exist: {}".format(label, path))
        for label, path in (("structure_file", self.structure_file), ("param_file", self.param_file), ("area_file", self.area_file), ("future_masks_json", self.future_masks_json)):
            if path is not None and not path.is_file():
                raise InputError("{} does not exist: {}".format(label, path))
        if self.area_file is None and self.future_masks_json is None:
            raise InputError("Provide [inputs] area_file or future_masks_json")
        if not self.species:
            raise InputError("[analysis] species must not be empty")
        if self.nx < 1 or self.ny < 1 or self.replicates < 1:
            raise InputError("nx, ny and replicates must be positive")
        if not 0 <= self.threshold <= 1:
            raise InputError("threshold must be between 0 and 1")
        if self.min_valid < 0:
            raise InputError("min_valid must not be negative")
        if self.loss_mode not in {"edge", "random"}:
            raise InputError("loss_mode must be edge or random")
        if self.fst_metric not in {"hudson", "nei"}:
            raise InputError("fst_metric must be hudson or nei")
        if self.mu <= 0 or self.theta_probe <= 0:
            raise InputError("mu and theta_probe must be positive")
        if self.z_gdar:
            try:
                z_gdar = float(self.z_gdar)
            except ValueError:
                raise InputError("z_gdar must be numeric when supplied")
            if z_gdar <= 0:
                raise InputError("z_gdar must be positive when supplied")
        if self.loss_mode == "random" and self.replicates < 2:
            # One replicate is valid for a smoke test, but warn in the manifest.
            pass

    def payload(self) -> Dict[str, object]:
        return {
            "inputs": {
                "current_raster": str(self.current_raster),
                "pi_file": str(self.pi_file),
                "structure_file": str(self.structure_file) if self.structure_file else None,
                "param_file": str(self.param_file) if self.param_file else None,
                "area_file": str(self.area_file) if self.area_file else None,
                "future_masks_json": str(self.future_masks_json) if self.future_masks_json else None,
            },
            "analysis": {"output_dir": str(self.output_dir), "species": self.species, "plot_title": self.plot_title},
            "parameters": {
                "compute_python": self.compute_python,
                "rscript": self.rscript,
                "nx": self.nx,
                "ny": self.ny,
                "threshold": self.threshold,
                "min_valid": self.min_valid,
                "loss_mode": self.loss_mode,
                "direction": self.direction,
                "migration": self.migration,
                "migration_grid": self.migration_grid,
                "fst_metric": self.fst_metric,
                "theta": self.theta,
                "z_gdar": self.z_gdar or None,
                "theta_probe": self.theta_probe,
                "time3": self.time3,
                "time5": self.time5,
                "mu": self.mu,
                "midterm_generations": self.midterm_generations,
                "replicates": self.replicates,
                "seed": self.seed,
                "include_equilibrium": self.include_equilibrium,
                "plot_equilibrium": self.plot_equilibrium,
                "plot_raw": self.plot_raw,
            },
        }


def _run(command: Sequence[str], log_path: Path, cwd: Optional[Path] = None) -> None:
    result = subprocess.run(
        list(command), cwd=str(cwd) if cwd else None, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        env=subprocess_environment(command[0]),
    )
    log_path.write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("WFmoments command failed; see {}".format(log_path))


def run_wfmoment_workflow(
    config_path: Path,
    dry_run: bool = False,
    progress: Optional[Callable[[str], None]] = None,
) -> Dict[str, object]:
    config = WFMomentConfig(config_path)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    script_dir = Path(__file__).resolve().parent / "scripts"
    compute_script = script_dir / "wfmoment_compute.py"
    plot_script = script_dir / "wfmoment_plot.R"
    if not compute_script.is_file() or not plot_script.is_file():
        raise InputError("WFmoments scripts are missing from {}".format(script_dir))
    manifest: Dict[str, object] = {
        "module": "wfmoment-2d-deme",
        "status": "planned" if dry_run else "running",
        "config": config.payload(),
    }
    input_paths = {
        "current_raster": config.current_raster,
        "pi_file": config.pi_file,
        "structure_file": config.structure_file,
        "param_file": config.param_file,
        "area_file": config.area_file,
        "future_masks_json": config.future_masks_json,
        "compute_script": compute_script,
        "plot_script": plot_script,
    }
    # A future-mask mapping is itself a small input table, but its referenced
    # rasters are scientific inputs too.  Record them individually so changing
    # a raster cannot be hidden behind an unchanged JSON mapping.
    if config.future_masks_json:
        try:
            mapping = json.loads(config.future_masks_json.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise InputError("Cannot read future_masks_json: {}".format(error))
        if not isinstance(mapping, dict):
            raise InputError("future_masks_json must contain an object mapping scenarios to rasters")
        for index, raw_path in enumerate(mapping.values(), start=1):
            if not isinstance(raw_path, str) or not raw_path.strip():
                raise InputError("future_masks_json contains an invalid raster path")
            mask_path = Path(raw_path).expanduser()
            if not mask_path.is_absolute():
                mask_path = config.future_masks_json.parent / mask_path
            input_paths["future_mask_{}".format(index)] = mask_path.resolve()
    attach_input_fingerprints(manifest, input_paths)
    if dry_run:
        (config.output_dir / "wfmoment_dry_run.json").write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
        )
        return manifest

    started = time.time()
    compute_args = [
        config.compute_python, str(compute_script),
        "--current-raster", str(config.current_raster),
        "--pi-file", str(config.pi_file),
        "--species", config.species,
        "--nx", str(config.nx), "--ny", str(config.ny),
        "--threshold", str(config.threshold), "--min-valid", str(config.min_valid),
        "--loss-mode", config.loss_mode, "--direction", config.direction,
        "--migration", config.migration, "--migration-grid", config.migration_grid,
        "--fst-metric", config.fst_metric, "--theta", config.theta,
        "--theta-probe", str(config.theta_probe), "--time3", str(config.time3),
        "--time5", str(config.time5), "--mu", str(config.mu),
        "--midterm-generations", config.midterm_generations,
        "--replicates", str(config.replicates), "--seed", str(config.seed),
        "--outdir", str(config.output_dir),
    ]
    if config.z_gdar:
        compute_args += ["--z-gdar", config.z_gdar]
    if config.structure_file:
        compute_args += ["--structure-file", str(config.structure_file)]
    if config.param_file:
        compute_args += ["--param-file", str(config.param_file)]
    if config.area_file:
        compute_args += ["--area-file", str(config.area_file)]
    if config.future_masks_json:
        compute_args += ["--future-masks-json", str(config.future_masks_json)]
    if config.include_equilibrium:
        compute_args.append("--include-equilibrium")
    if progress:
        progress("WFmoments 2-D deme calculation ({}, {} loss)".format(config.species, config.loss_mode))
    _run(compute_args, config.output_dir / "wfmoment_compute.log")

    plot_args = [
        config.rscript, str(plot_script),
        str(config.output_dir / "curve_summary.tsv"),
        str(config.output_dir / "scenario_points.tsv"),
        str(config.output_dir),
        config.species,
        config.plot_title,
        "true" if config.plot_equilibrium else "false",
        "true" if config.plot_raw else "false",
    ]
    if progress:
        progress("WFmoments figure")
    _run(plot_args, config.output_dir / "wfmoment_plot.log")

    manifest.update(
        {
            "status": "complete",
            "elapsed_seconds": round(time.time() - started, 3),
            "outputs": {
                "output_dir": str(config.output_dir),
                "curve_summary": str(config.output_dir / "curve_summary.tsv"),
                "scenario_points": str(config.output_dir / "scenario_points.tsv"),
                "figure_png": str(config.output_dir / "Figure_wfmoment_2D_deme.png"),
                "figure_pdf": str(config.output_dir / "Figure_wfmoment_2D_deme.pdf"),
            },
            "config_sha256": hashlib.sha256(config.config_path.read_bytes()).hexdigest(),
        }
    )
    (config.output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    return manifest
