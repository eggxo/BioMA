"""Configurable MaxEnt niche modelling and future projection workflow."""

from __future__ import annotations

import configparser
import hashlib
import json
import re
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional

from .gf_frequency import InputError
from .provenance import attach_input_fingerprints
from .runtime import subprocess_environment


_BIO_RE = re.compile(r"bio[_-]?(\d+)", re.IGNORECASE)


def _split(value: str) -> List[str]:
    return [x.strip() for x in value.split(",") if x.strip()]


def _bool(value: str, default: bool = False) -> bool:
    if value is None:
        return default
    value = value.strip().lower()
    if value in {"1", "true", "yes", "y", "on"}:
        return True
    if value in {"0", "false", "no", "n", "off"}:
        return False
    raise InputError("Invalid boolean value: {}".format(value))


def _path(value: str, base: Path, label: str, required: bool = True) -> Optional[Path]:
    value = (value or "").strip()
    if not value:
        if required:
            raise InputError("Missing [inputs] {}".format(label))
        return None
    p = Path(value).expanduser()
    if not p.is_absolute():
        p = base / p
    return p.resolve()


def _bio_ids(directory: Path) -> List[int]:
    files = list(directory.glob("*.tif"))
    ids = []
    for f in files:
        hit = _BIO_RE.search(f.name)
        if hit:
            ids.append(int(hit.group(1)))
    return sorted(set(ids))


@dataclass(frozen=True)
class NicheConfig:
    config_path: Path
    work_dir: Path
    occurrence_csv: Path
    current_env_dir: Path
    future_root: Path
    mask_shp: Path
    output_dir: Path
    periods: List[str]
    scenarios: List[str]
    gcms: List[str]
    correlation_threshold: float
    correlation_method: str
    correlation_scope: str
    subset_sizes: List[int]
    max_models: int
    background_n: int
    cv_folds: int
    feature_classes: List[str]
    beta_multipliers: List[float]
    selection_metric: str
    ensemble_method: str
    binary_output: bool
    seed: int
    rscript: str
    plot_rscript: str
    maxent_jar: Optional[Path]

    def payload(self) -> Dict[str, object]:
        return {
            "inputs": {
                "work_dir": str(self.work_dir),
                "occurrence_csv": str(self.occurrence_csv),
                "current_env_dir": str(self.current_env_dir),
                "future_root": str(self.future_root),
                "mask_shp": str(self.mask_shp),
            },
            "variables": {
                "correlation_threshold": self.correlation_threshold,
                "correlation_method": self.correlation_method,
                "correlation_scope": self.correlation_scope,
                "subset_sizes": self.subset_sizes,
            },
            "tuning": {
                "max_models": self.max_models,
                "background_n": self.background_n,
                "cv_folds": self.cv_folds,
                "feature_classes": self.feature_classes,
                "beta_multipliers": self.beta_multipliers,
                "selection_metric": self.selection_metric,
            },
            "projection": {
                "periods": self.periods,
                "scenarios": self.scenarios,
                "gcms": self.gcms,
                "ensemble_method": self.ensemble_method,
                "binary_output": self.binary_output,
            },
            "parameters": {"seed": self.seed, "rscript": self.rscript, "plot_rscript": self.plot_rscript, "maxent_jar": str(self.maxent_jar) if self.maxent_jar else None},
            "output_dir": str(self.output_dir),
        }


def _discover_tokens(root: Path):
    periods, scenarios, gcms = set(), set(), set()
    for directory in root.iterdir():
        if not directory.is_dir():
            continue
        parts = directory.name.split("-", 2)
        if len(parts) == 3 and re.fullmatch(r"\d{4}-\d{4}", parts[0]) and parts[1].lower().startswith("ssp"):
            periods.add(parts[0])
            scenarios.add(parts[1])
            gcms.add(parts[2])
    return sorted(periods), sorted(scenarios), sorted(gcms)


def load_config(path: Path) -> NicheConfig:
    path = path.expanduser().resolve()
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as exc:
        raise InputError("Cannot read niche configuration: {}".format(exc))
    required = {"inputs", "variables", "tuning", "projection", "analysis"}
    missing = sorted(required.difference(parser.sections()))
    if missing:
        raise InputError("Niche configuration missing sections: {}".format(", ".join(missing)))

    base = path.parent
    i, v, t, p, a = (parser[s] for s in ("inputs", "variables", "tuning", "projection", "analysis"))
    work_dir = _path(i.get("work_dir", ""), base, "work_dir")
    occurrence = _path(i.get("occurrence_csv", ""), base, "occurrence_csv")
    current = _path(i.get("current_env_dir", ""), base, "current_env_dir")
    future_root = _path(i.get("future_root", ""), base, "future_root")
    mask = _path(i.get("mask_shp", ""), base, "mask_shp")
    periods = _split(p.get("periods", "auto"))
    scenarios = _split(p.get("scenarios", "auto"))
    gcms = _split(p.get("gcms", "auto"))
    discovered = _discover_tokens(future_root)
    if periods == ["auto"]:
        periods = discovered[0]
    if scenarios == ["auto"]:
        scenarios = discovered[1]
    if gcms == ["auto"]:
        gcms = discovered[2]
    threshold = float(v.get("correlation_threshold", "0.8"))
    method = v.get("correlation_method", "pearson").strip().lower()
    scope = v.get("correlation_scope", "occurrence").strip().lower()
    subset_sizes = [int(x) for x in _split(v.get("candidate_subset_sizes", "4,6,8"))]
    feature_classes = [x.upper() for x in _split(t.get("feature_classes", "L,LQ,LQH,LQHP"))]
    betas = [float(x) for x in _split(t.get("beta_multipliers", "0.5,1,2,3,4"))]
    cfg = NicheConfig(
        config_path=path,
        work_dir=work_dir,
        occurrence_csv=occurrence,
        current_env_dir=current,
        future_root=future_root,
        mask_shp=mask,
        output_dir=_path(a.get("output_dir", ""), base, "output_dir"),
        periods=periods,
        scenarios=scenarios,
        gcms=gcms,
        correlation_threshold=threshold,
        correlation_method=method,
        correlation_scope=scope,
        subset_sizes=subset_sizes,
        max_models=int(t.get("max_models", "500")),
        background_n=int(t.get("background_n", "10000")),
        cv_folds=int(t.get("cv_folds", "5")),
        feature_classes=feature_classes,
        beta_multipliers=betas,
        selection_metric=t.get("selection_metric", "auc").strip().lower(),
        ensemble_method=p.get("ensemble_method", "mean").strip().lower(),
        binary_output=_bool(p.get("binary_output", "false")),
        seed=int(t.get("seed", "123")),
        rscript=a.get("rscript", "Rscript").strip() or "Rscript",
        plot_rscript=a.get("plot_rscript", "Rscript").strip() or "Rscript",
        maxent_jar=_path(i.get("maxent_jar", ""), base, "maxent_jar", False),
    )
    for label, value in (("work_dir", cfg.work_dir), ("occurrence_csv", cfg.occurrence_csv), ("current_env_dir", cfg.current_env_dir), ("future_root", cfg.future_root), ("mask_shp", cfg.mask_shp)):
        if not value.exists():
            raise InputError("{} does not exist: {}".format(label, value))
    if cfg.maxent_jar and not cfg.maxent_jar.is_file():
        raise InputError("maxent_jar does not exist: {}".format(cfg.maxent_jar))
    if not periods or not scenarios or not gcms:
        raise InputError("periods, scenarios and gcms must not be empty")
    if not 0 < threshold < 1:
        raise InputError("correlation_threshold must be between 0 and 1")
    if method not in {"pearson", "spearman", "kendall"}:
        raise InputError("Unsupported correlation_method: {}".format(method))
    if scope not in {"occurrence", "background", "mask"}:
        raise InputError("Unsupported correlation_scope: {}".format(scope))
    if cfg.selection_metric not in {"auc", "tss"}:
        raise InputError("selection_metric must be auc or tss")
    if cfg.ensemble_method not in {"mean", "median"}:
        raise InputError("ensemble_method must be mean or median")
    if cfg.binary_output:
        raise InputError("binary_output must be false for the current continuous-suitability workflow")
    if cfg.max_models < 1 or cfg.background_n < 10 or cfg.cv_folds < 2 or not cfg.subset_sizes or any(x < 1 for x in cfg.subset_sizes):
        raise InputError("Invalid tuning parameters")
    if not feature_classes or not betas:
        raise InputError("feature_classes and beta_multipliers must not be empty")
    return cfg


def _future_dirs(cfg: NicheConfig) -> Dict[str, Path]:
    found = {}
    for period in cfg.periods:
        for scenario in cfg.scenarios:
            for gcm in cfg.gcms:
                key = "{}__{}__{}".format(period, scenario, gcm)
                directory = cfg.future_root / "{}-{}-{}".format(period, scenario, gcm)
                if not directory.is_dir():
                    raise InputError("Future climate directory missing: {}".format(directory))
                if not _bio_ids(directory):
                    raise InputError("No BIO tif files in: {}".format(directory))
                found[key] = directory
    return found


def _write_spec(cfg: NicheConfig, future_dirs: Dict[str, Path]) -> Path:
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    spec = cfg.output_dir / "niche_config.tsv"
    rows = {"work_dir": str(cfg.work_dir), "occurrence_csv": str(cfg.occurrence_csv), "current_env_dir": str(cfg.current_env_dir), "future_root": str(cfg.future_root), "mask_shp": str(cfg.mask_shp), "output_dir": str(cfg.output_dir), "periods": ",".join(cfg.periods), "scenarios": ",".join(cfg.scenarios), "gcms": ",".join(cfg.gcms), "correlation_threshold": str(cfg.correlation_threshold), "correlation_method": cfg.correlation_method, "correlation_scope": cfg.correlation_scope, "subset_sizes": ",".join(map(str, cfg.subset_sizes)), "max_models": str(cfg.max_models), "background_n": str(cfg.background_n), "cv_folds": str(cfg.cv_folds), "feature_classes": ",".join(cfg.feature_classes), "beta_multipliers": ",".join(map(str, cfg.beta_multipliers)), "selection_metric": cfg.selection_metric, "ensemble_method": cfg.ensemble_method, "binary_output": str(cfg.binary_output).lower(), "seed": str(cfg.seed), "maxent_jar": str(cfg.maxent_jar or "")}
    with spec.open("w", encoding="utf-8", newline="") as handle:
        for key, value in rows.items():
            handle.write("{}\t{}\n".format(key, value.replace("\t", " ")))
    return spec


def run_niche_workflow(config_path: Path, dry_run: bool = False, progress: Optional[Callable[[str], None]] = None) -> Dict[str, object]:
    cfg = load_config(config_path)
    current_ids = _bio_ids(cfg.current_env_dir)
    if not current_ids:
        raise InputError("No BIO tif files in current_env_dir: {}".format(cfg.current_env_dir))
    future_dirs = _future_dirs(cfg)
    script_dir = Path(__file__).resolve().parent / "scripts"
    pipeline = script_dir / "niche_pipeline.R"
    plotter = script_dir / "niche_plot.R"
    if not pipeline.is_file() or not plotter.is_file():
        raise InputError("Niche scripts are missing from {}".format(script_dir))
    jar_input = None
    if cfg.maxent_jar:
        jar_input = {
            "path": str(cfg.maxent_jar),
            "sha256": hashlib.sha256(cfg.maxent_jar.read_bytes()).hexdigest(),
        }
    manifest: Dict[str, object] = {"module": "niche-workflow", "status": "planned" if dry_run else "running", "config": cfg.payload(), "inputs": {"maxent_jar": jar_input}, "discovery": {"current_bio": current_ids, "future_directories": {k: str(v) for k, v in future_dirs.items()}}}
    input_paths = {
        "occurrence_csv": cfg.occurrence_csv,
        "current_env_dir": cfg.current_env_dir,
        "future_root": cfg.future_root,
        "mask_shp": cfg.mask_shp,
        "maxent_jar": cfg.maxent_jar,
        "pipeline_script": pipeline,
        "plot_script": plotter,
    }
    input_paths.update({"future_scenario__{}".format(key): path for key, path in future_dirs.items()})
    attach_input_fingerprints(manifest, input_paths)
    if dry_run:
        cfg.output_dir.mkdir(parents=True, exist_ok=True)
        (cfg.output_dir / "niche_dry_run.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
        return manifest

    started = time.time()
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    spec = _write_spec(cfg, future_dirs)
    if progress:
        progress("Tuning MaxEnt variables and feature parameters")
    result = subprocess.run([cfg.rscript, str(pipeline), str(spec)], cwd=str(cfg.work_dir), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.rscript))
    (cfg.output_dir / "niche_pipeline.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("MaxEnt niche pipeline failed; see {}".format(cfg.output_dir / "niche_pipeline.log"))
    if progress:
        progress("Drawing continuous suitability and maladaptation maps")
    result = subprocess.run([cfg.plot_rscript, str(plotter), str(cfg.output_dir), str(cfg.occurrence_csv), str(cfg.mask_shp)], cwd=str(cfg.work_dir), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.plot_rscript))
    (cfg.output_dir / "niche_plot.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise InputError("Niche plotting failed; see {}".format(cfg.output_dir / "niche_plot.log"))
    manifest.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 3), "outputs": {"output_dir": str(cfg.output_dir), "tuning_metrics": str(cfg.output_dir / "tuning_metrics.csv"), "selected_model": str(cfg.output_dir / "selected_model.tsv"), "rasters": str(cfg.output_dir / "rasters"), "figures": str(cfg.output_dir / "figures")}, "config_sha256": hashlib.sha256(cfg.config_path.read_bytes()).hexdigest()})
    (cfg.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
