"""RONA calculation and plotting workflow.

The numerical calculation follows ``RONA.revise.R``.  This module adds a
configuration-driven wrapper, scenario discovery, validation, and manifests.
"""

from __future__ import annotations

import configparser
import hashlib
import json
import os
import re
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from .gf_frequency import InputError
from .provenance import attach_input_fingerprints
from .runtime import subprocess_environment


SCENARIO_RE = re.compile(r"^(\d{4}-\d{4})-(ssp\d+)-(.+)$", re.IGNORECASE)
BIO_COUNT = 19


@dataclass(frozen=True)
class RonaConfig:
    config_path: Path
    alt_frequency: Path
    unld_dir: Path
    environment: Path
    future_climate: Path
    mask: Optional[Path]
    output_dir: Path
    models: Tuple[str, ...]
    ssps: Tuple[str, ...]
    periods: Tuple[str, ...]
    expected_populations: Optional[int]
    rscript: str
    interpolate: bool
    grid_step: float
    grid_k: int

    def payload(self) -> Dict[str, object]:
        return {
            "inputs": {
                "alt_frequency": str(self.alt_frequency),
                "unld_dir": str(self.unld_dir),
                "environment": str(self.environment),
                "future_climate": str(self.future_climate),
                "mask": str(self.mask) if self.mask else None,
            },
            "analysis": {
                "output_dir": str(self.output_dir),
                "models": list(self.models),
                "ssps": list(self.ssps),
                "periods": list(self.periods),
                "expected_populations": self.expected_populations,
            },
            "parameters": {
                "rscript": self.rscript,
                "interpolate": self.interpolate,
                "grid_step": self.grid_step,
                "grid_k": self.grid_k,
            },
        }


@dataclass(frozen=True)
class RonaScenario:
    name: str
    directory: Path
    period: str
    ssp: str
    model: str


def _values(value: str) -> Tuple[str, ...]:
    return tuple(dict.fromkeys(x.strip() for x in value.split(",") if x.strip()))


def _path(value: str, base: Path, required: bool = True) -> Optional[Path]:
    value = value.strip()
    if not value:
        if required:
            raise InputError("Missing required RONA input path")
        return None
    result = Path(value).expanduser()
    if not result.is_absolute():
        result = base / result
    return result.resolve()


def _find_rscript(value: str) -> str:
    if value.strip():
        return value.strip()
    return "Rscript"


def load_rona_config(path: Path) -> RonaConfig:
    path = path.expanduser().resolve()
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as error:
        raise InputError("Cannot read RONA configuration: {}".format(error))
    if not parser.has_section("inputs") or not parser.has_section("analysis"):
        raise InputError("RONA configuration requires [inputs] and [analysis]")
    inputs, analysis = parser["inputs"], parser["analysis"]
    parameters = parser["parameters"] if parser.has_section("parameters") else {}
    base = path.parent
    models = _values(analysis.get("models", ""))
    ssps = tuple(x.lower().replace("ssp", "") for x in _values(analysis.get("ssps", "")))
    periods = _values(analysis.get("periods", ""))
    expected = analysis.get("expected_populations", "").strip()
    expected_populations = int(expected) if expected else None
    grid_step = float(parameters.get("grid_step", "0.1"))
    grid_k = int(parameters.get("grid_k", "15"))
    if grid_step <= 0 or grid_k < 3:
        raise InputError("grid_step must be positive and grid_k must be at least 3")
    mask = _path(inputs.get("mask", ""), base, required=False)
    config = RonaConfig(
        config_path=path,
        alt_frequency=_path(inputs.get("alt_frequency", ""), base),
        unld_dir=_path(inputs.get("unld_dir", ""), base),
        environment=_path(inputs.get("environment", ""), base),
        future_climate=_path(inputs.get("future_climate", ""), base),
        mask=mask,
        output_dir=_path(analysis.get("output_dir", ""), base),
        models=models,
        ssps=ssps,
        periods=periods,
        expected_populations=expected_populations,
        rscript=_find_rscript(parameters.get("rscript", "")),
        interpolate=parameters.get("interpolate", "true").strip().lower() not in {"0", "false", "no"},
        grid_step=grid_step,
        grid_k=grid_k,
    )
    for label, item in (("alt frequency", config.alt_frequency), ("environment", config.environment)):
        if not item.is_file():
            raise InputError("{} does not exist: {}".format(label, item))
    for label, item in (("unLD directory", config.unld_dir), ("future climate", config.future_climate)):
        if not item.is_dir():
            raise InputError("{} does not exist: {}".format(label, item))
    if config.mask is not None and not config.mask.is_file():
        raise InputError("mask does not exist: {}".format(config.mask))
    if expected_populations is not None and expected_populations < 1:
        raise InputError("expected_populations must be positive")
    return config


def discover_rona_scenarios(root: Path) -> Tuple[RonaScenario, ...]:
    found: List[RonaScenario] = []
    for directory in sorted(p for p in root.rglob("*") if p.is_dir()):
        match = SCENARIO_RE.match(directory.name.strip())
        if not match:
            continue
        bio_files = [directory / "{}bio{}.cut.tif".format(directory.name, i) for i in range(1, 20)]
        if not all(p.is_file() for p in bio_files):
            bio_files = list(directory.glob("*bio*.cut.tif"))
            if len({re.search(r"bio(\d+)\.cut\.tif$", p.name, re.I).group(1) for p in bio_files if re.search(r"bio(\d+)\.cut\.tif$", p.name, re.I)}) < 19:
                continue
        period, ssp, model = match.groups()
        found.append(RonaScenario(directory.name.strip(), directory, period, ssp.lower().replace("ssp", ""), model.strip()))
    unique = {item.name: item for item in found}
    return tuple(unique[name] for name in sorted(unique))


def _select(scenarios: Iterable[RonaScenario], config: RonaConfig) -> Tuple[RonaScenario, ...]:
    selected = []
    for item in scenarios:
        if config.models and item.model not in config.models:
            continue
        if config.ssps and item.ssp not in config.ssps:
            continue
        if config.periods and item.period not in config.periods:
            continue
        selected.append(item)
    if not selected:
        raise InputError("No RONA scenarios match the configured models, SSPs, and periods")
    return tuple(selected)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _ld_prune_records(
    unld_dir: Path, *, include_unique: bool = False
) -> Union[List[Dict[str, object]], Tuple[List[Dict[str, object]], int]]:
    """Validate and fingerprint the per-BIO LD-pruning lists.

    RONA does not derive LD from the environment table.  The R script reads
    one ``LD_BIO*.prune.in`` list for each BIO and intersects those IDs with
    the frequency matrix.  Checking the lists before launching R makes that
    scientific dependency explicit and gives the run manifest a content hash
    for every list.  When ``include_unique`` is true, also return the number
    of unique locus IDs in the union of all lists.
    """
    records: List[Dict[str, object]] = []
    unique_loci = set()
    for bio in range(1, BIO_COUNT + 1):
        expected = unld_dir / "LD_BIO{}.prune.in".format(bio)
        path = expected if expected.is_file() else None
        if path is None:
            # Be tolerant of case differences on filesystems copied from a
            # Windows workstation, while still rejecting ambiguous matches.
            matches = [
                item for item in unld_dir.iterdir()
                if item.is_file() and item.name.lower() == expected.name.lower()
            ]
            if len(matches) == 1:
                path = matches[0]
        if path is None:
            raise InputError("RONA is missing LD-pruning list: {}".format(expected))
        try:
            ids = [line.strip() for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines() if line.strip()]
        except OSError as error:
            raise InputError("Cannot read RONA LD-pruning list {}: {}".format(path, error))
        if not ids:
            raise InputError("RONA LD-pruning list is empty: {}".format(path))
        records.append(
            {
                "bio": bio,
                "path": str(path.resolve()),
                "sha256": _sha256(path),
                "loci": len(dict.fromkeys(ids)),
            }
        )
        unique_loci.update(ids)
    if include_unique:
        # Keep the public record shape stable while allowing the run manifest
        # to distinguish per-list counts from the union across BIO lists.
        return records, len(unique_loci)
    return records


def _scenario_bio_files(scenario: RonaScenario) -> List[Path]:
    """Return the 19 raster files consumed for one selected scenario.

    The R implementation searches the selected directory for one
    ``bioN.cut.tif`` file per BIO. Fingerprinting this exact set avoids
    recursively hashing unrelated files in the future-climate catalog.
    """
    by_bio: Dict[int, List[Path]] = {}
    for path in scenario.directory.glob("*"):
        if not path.is_file():
            continue
        match = re.search(r"bio(\d+)\.cut\.tif$", path.name, re.IGNORECASE)
        if not match:
            continue
        by_bio.setdefault(int(match.group(1)), []).append(path)
    selected: List[Path] = []
    for bio in range(1, BIO_COUNT + 1):
        matches = by_bio.get(bio, [])
        if len(matches) != 1:
            raise InputError(
                "RONA scenario {} must contain exactly one BIO{} raster".format(
                    scenario.name, bio
                )
            )
        selected.append(matches[0].resolve())
    return selected


def run_rona_workflow(config_path: Path, dry_run: bool = False, progress: Optional[Callable[[str], None]] = None) -> Dict[str, object]:
    config = load_rona_config(config_path)
    ld_records, unique_loci = _ld_prune_records(config.unld_dir, include_unique=True)
    all_scenarios = discover_rona_scenarios(config.future_climate)
    scenarios = _select(all_scenarios, config)
    script_dir = Path(__file__).resolve().parent / "scripts"
    compute_script = script_dir / "rona_compute.R"
    plot_script = script_dir / "rona_plot.R"
    if not compute_script.is_file() or not plot_script.is_file():
        raise InputError("RONA scripts are missing from {}".format(script_dir))
    manifest: Dict[str, object] = {
        "module": "rona-workflow",
        "status": "planned" if dry_run else "running",
        "config": config.payload(),
        "scenarios": [
            {
                "name": item.name,
                "directory": str(item.directory),
                "period": item.period,
                "ssp": item.ssp,
                "model": item.model,
            }
            for item in scenarios
        ],
        "counts": {"available_scenarios": len(all_scenarios), "selected_scenarios": len(scenarios), "populations": config.expected_populations},
        "inputs": {
            "ld_pruning": {
                "description": "Per-BIO PLINK prune.in lists consumed by rona_compute.R",
                "files": ld_records,
                "total_loci_listed": sum(int(row["loci"]) for row in ld_records),
                "total_unique_loci_listed": unique_loci,
            }
        },
    }
    input_paths: Dict[str, Optional[Path]] = {
        "config_path": config.config_path,
        "alt_frequency": config.alt_frequency,
        "environment": config.environment,
        "mask": config.mask,
        "compute_script": compute_script,
        "plot_script": plot_script,
    }
    # LD lists are consumed file-by-file; future rasters are consumed as a
    # selected scenario directory. Keep both scopes explicit while avoiding a
    # recursive hash of the catalog root, which may contain many unselected
    # scenarios.
    for row in ld_records:
        input_paths["ld_pruning.BIO{}".format(row["bio"])] = Path(str(row["path"]))
    for scenario in scenarios:
        # Validate the exact set consumed by the R code, then fingerprint the
        # selected directory once. Its member list records each raster hash;
        # unrelated scenarios elsewhere in the catalog are never traversed.
        _scenario_bio_files(scenario)
        input_paths["future_scenario__{}".format(scenario.name)] = scenario.directory
    attach_input_fingerprints(manifest, input_paths)
    config.output_dir.mkdir(parents=True, exist_ok=True)
    if dry_run:
        (config.output_dir / "rona_dry_run.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        return manifest
    started = time.time()
    if progress:
        progress("RONA compute: {} scenarios".format(len(scenarios)))
    args = [config.rscript, str(compute_script), str(config.alt_frequency), str(config.unld_dir), str(config.environment), str(config.future_climate), str(config.output_dir), ",".join(config.models), ",".join(config.ssps), ",".join(config.periods)]
    result = subprocess.run(args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(config.rscript))
    (config.output_dir / "rona_compute.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("RONA computation failed; see {}".format(config.output_dir / "rona_compute.log"))
    if progress:
        progress("RONA plotting")
    plot_args = [config.rscript, str(plot_script), str(config.output_dir), str(config.environment), str(config.mask or ""), ",".join(config.models), ",".join(config.ssps), ",".join(config.periods), str(config.grid_step), str(config.grid_k), "1" if config.interpolate else "0"]
    result = subprocess.run(plot_args, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(config.rscript))
    (config.output_dir / "rona_plot.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("RONA plotting failed; see {}".format(config.output_dir / "rona_plot.log"))
    manifest.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 3), "outputs": {"output_dir": str(config.output_dir), "compute_log": str(config.output_dir / "rona_compute.log"), "plot_log": str(config.output_dir / "rona_plot.log")}})
    manifest["manifest_sha256"] = _sha256(config.config_path)
    (config.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
