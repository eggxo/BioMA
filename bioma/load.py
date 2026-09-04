"""Derived-VCF genetic-load calculation, RF tuning, prediction and maps."""

from __future__ import annotations

import configparser
import csv
import hashlib
import json
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from .gf_frequency import InputError
from .runtime import subprocess_environment


@dataclass(frozen=True)
class LoadConfig:
    config_path: Path
    vcf_dir: Path
    population_dir: Path
    predictors: Path
    future_dir: Path
    mask: Optional[Path]
    output_dir: Path
    n_splits: int
    cv_v: int
    cv_repeats: int
    grid_size: int
    workers: int
    scheme: str
    rscript: str
    calc_python: str
    seed: int
    calc_script: Path
    rf_script: Path
    predict_script: Path
    future_pattern: str
    expected_future_files: int

    def payload(self):
        return {
            "inputs": {"vcf_dir": str(self.vcf_dir), "population_dir": str(self.population_dir), "predictors": str(self.predictors), "future_dir": str(self.future_dir), "mask": str(self.mask) if self.mask else None},
            "parameters": {"n_splits": self.n_splits, "cv_v": self.cv_v, "cv_repeats": self.cv_repeats, "grid_size": self.grid_size, "workers": self.workers, "scheme": self.scheme, "rscript": self.rscript, "calc_python": self.calc_python, "seed": self.seed, "calc_script": str(self.calc_script), "rf_script": str(self.rf_script), "predict_script": str(self.predict_script), "future_pattern": self.future_pattern, "expected_future_files": self.expected_future_files},
            "output_dir": str(self.output_dir),
        }


def _path(value: str, base: Path, label: str, required=True):
    value = value.strip()
    if not value:
        if required:
            raise InputError("Missing [inputs] {}".format(label))
        return None
    p = Path(value).expanduser()
    if not p.is_absolute(): p = base / p
    return p.resolve()


def _script_path(value: str, base: Path, bundled_name: str) -> Path:
    value = (value or "").strip()
    if value:
        path = Path(value).expanduser()
        if not path.is_absolute():
            path = base / path
        return path.resolve()
    return (Path(__file__).resolve().parent / "scripts" / bundled_name).resolve()


def load_config(path: Path) -> LoadConfig:
    path = path.expanduser().resolve(); parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as h: parser.read_file(h)
    except (OSError, configparser.Error, UnicodeError) as e:
        raise InputError("Cannot read load configuration: {}".format(e))
    if not parser.has_section("inputs") or not parser.has_section("analysis"):
        raise InputError("Load configuration requires [inputs] and [analysis]")
    i, a = parser["inputs"], parser["analysis"]; p = parser["parameters"] if parser.has_section("parameters") else {}; base = path.parent
    cfg = LoadConfig(
        path,
        _path(i.get("vcf_dir", ""), base, "vcf_dir"),
        _path(i.get("population_dir", ""), base, "population_dir"),
        _path(i.get("predictors", ""), base, "predictors"),
        _path(i.get("future_dir", ""), base, "future_dir"),
        _path(i.get("mask", ""), base, "mask", False),
        _path(a.get("output_dir", ""), base, "output_dir"),
        int(p.get("n_splits", "20")),
        int(p.get("cv_v", "5")),
        int(p.get("cv_repeats", "2")),
        int(p.get("grid_size", "30")),
        int(p.get("workers", "4")),
        p.get("scheme", "svd").strip().lower(),
        p.get("rscript", "Rscript").strip() or "Rscript",
        p.get("calc_python", sys.executable).strip(),
        int(p.get("seed", "1234")),
        _script_path(p.get("calc_script", ""), base, "calc_loadM_loadD_from_3vcf.py"),
        _script_path(p.get("rf_script", ""), base, "rf_one_target_svd.R"),
        _script_path(p.get("predict_script", ""), base, "predict_future_one_target.R"),
        p.get("future_pattern", "future_env_ssp*_mean.csv").strip(),
        int(p.get("expected_future_files", "0")),
    )
    for label, q in (("vcf_dir", cfg.vcf_dir), ("population_dir", cfg.population_dir), ("predictors", cfg.predictors), ("future_dir", cfg.future_dir)):
        if not q.exists(): raise InputError("{} does not exist: {}".format(label, q))
    if cfg.mask and not cfg.mask.is_file(): raise InputError("mask does not exist: {}".format(cfg.mask))
    for label, script in (("calc_script", cfg.calc_script), ("rf_script", cfg.rf_script), ("predict_script", cfg.predict_script)):
        if not script.is_file(): raise InputError("{} does not exist: {}".format(label, script))
    if cfg.n_splits < 1 or cfg.cv_v < 2 or cfg.cv_repeats < 1 or cfg.grid_size < 1 or cfg.workers < 1:
        raise InputError("n_splits, cv_v, cv_repeats, grid_size and workers must be positive; cv_v >= 2")
    if cfg.scheme != "svd":
        raise InputError("scheme must be svd in this BioMA release")
    if not cfg.future_pattern:
        raise InputError("future_pattern must not be empty")
    if cfg.expected_future_files < 0:
        raise InputError("expected_future_files must be zero or positive")
    return cfg


def _find_vcf(root: Path, pattern: str) -> Path:
    found = sorted(root.glob(pattern))
    if len(found) != 1: raise InputError("Expected one VCF matching {} under {}, found {}".format(pattern, root, len(found)))
    return found[0]


def _count_vcf_records(path: Path) -> int:
    """Count variant records without leaking a file descriptor."""
    count = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line and not line.startswith("#"):
                count += 1
    return count


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_calc(cfg: LoadConfig, syn: Path, deleterious: Path, nonsyn: Path, prefix: Path):
    prefix.parent.mkdir(parents=True, exist_ok=True)
    cmd = [cfg.calc_python, str(cfg.calc_script), str(syn), str(deleterious), str(nonsyn), str(cfg.population_dir), str(prefix)]
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.calc_python))
    (prefix.parent / (prefix.name + ".log")).write_text(result.stdout, encoding="utf-8")
    if result.returncode: raise InputError("load calculation failed; see {}.log".format(prefix))
    return prefix.with_name(prefix.name + ".population_load.tsv")


def _merge_training(base_csv: Path, strict_tsv: Path, relax_tsv: Path, output: Path):
    loads = {}
    for label, path in (("strict", strict_tsv), ("relax", relax_tsv)):
        with path.open(encoding="utf-8") as h:
            for row in csv.DictReader(h, delimiter="\t"):
                loads.setdefault(row["pop"], {})["mean_loadM_" + label] = row["mean_loadM"]
                loads[row["pop"]]["mean_loadD_" + label] = row["mean_loadD"]
    with base_csv.open(newline="", encoding="utf-8-sig") as src:
        rows = list(csv.DictReader(src)); fields = list(rows[0].keys()) if rows else []
    for key in ("mean_loadM_strict", "mean_loadD_strict", "mean_loadM_relax", "mean_loadD_relax"):
        if key not in fields: fields.append(key)
    for row in rows:
        if row.get("pop") not in loads: raise InputError("predictor table population missing from load table: {}".format(row.get("pop")))
        row.update(loads[row["pop"]])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8") as h:
        w = csv.DictWriter(h, fieldnames=fields); w.writeheader(); w.writerows(rows)
    return len(rows)


def run_load_workflow(config_path: Path, dry_run=False, progress: Optional[Callable[[str], None]] = None) -> Dict[str, object]:
    cfg = load_config(config_path)
    manifest = {"module": "load-workflow", "status": "planned" if dry_run else "running", "config": cfg.payload()}
    syn = _find_vcf(cfg.vcf_dir, "strict.synonymous.vcf")
    nonsyn = _find_vcf(cfg.vcf_dir, "strict.nonsynonymous.vcf")
    strict_del = _find_vcf(cfg.vcf_dir, "strict.deleterious_plain.vcf")
    relax_del = _find_vcf(cfg.vcf_dir, "relaxed.deleterious_with_warning.vcf")
    site_counts = {
        "strict_synonymous": _count_vcf_records(syn),
        "strict_nonsynonymous": _count_vcf_records(nonsyn),
        "strict_deleterious": _count_vcf_records(strict_del),
        "relaxed_deleterious": _count_vcf_records(relax_del),
    }
    manifest["inputs"] = {
        "vcfs": {"synonymous": str(syn), "nonsynonymous": str(nonsyn), "strict_deleterious": str(strict_del), "relaxed_deleterious": str(relax_del)},
        "site_counts": site_counts,
        "scripts": {
            "calculator": {"path": str(cfg.calc_script), "sha256": _sha256(cfg.calc_script)},
            "rf_tuning": {"path": str(cfg.rf_script), "sha256": _sha256(cfg.rf_script)},
            "future_prediction": {"path": str(cfg.predict_script), "sha256": _sha256(cfg.predict_script)},
        },
    }
    if dry_run:
        cfg.output_dir.mkdir(parents=True, exist_ok=True); (cfg.output_dir / "load_dry_run.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8"); return manifest
    started = time.time(); loads_dir = cfg.output_dir / "loads"; training_csv = cfg.output_dir / "pop_geo_niche_predictors_from_TSS.csv"
    if progress: progress("Calculating strict and relaxed load tables")
    strict_tsv = _run_calc(cfg, syn, strict_del, nonsyn, loads_dir / "load_strict")
    relax_tsv = _run_calc(cfg, syn, relax_del, nonsyn, loads_dir / "load_relax")
    _merge_training(cfg.predictors, strict_tsv, relax_tsv, training_csv)
    rf_script = cfg.rf_script
    targets = ("mean_loadM_relax", "mean_loadD_relax")
    rf_dirs = {}
    for target in targets:
        if progress: progress("Tuning random forest for {}".format(target))
        cmd = [
            cfg.rscript,
            str(rf_script),
            target,
            str(cfg.output_dir),
            str(cfg.n_splits),
            str(cfg.cv_v),
            str(cfg.cv_repeats),
            str(cfg.workers),
            str(cfg.grid_size),
            str(cfg.seed),
        ]
        result = subprocess.run(cmd, cwd=str(cfg.output_dir), text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.rscript))
        (cfg.output_dir / (target + ".rf.log")).write_text(result.stdout, encoding="utf-8")
        if result.returncode: raise InputError("RF tuning failed for {}; see {}.rf.log".format(target, target))
        rf_dirs[target] = cfg.output_dir / ("RF_" + re.sub(r"[^A-Za-z0-9]+", "_", target))
    if progress: progress("Predicting configured future scenarios and drawing figures")
    plot_script = Path(__file__).resolve().parent / "scripts" / "load_plot.R"
    pred_script = cfg.predict_script
    future_files = sorted(cfg.future_dir.glob(cfg.future_pattern))
    if not future_files:
        raise InputError("No future files matching {} under {}".format(cfg.future_pattern, cfg.future_dir))
    if cfg.expected_future_files and len(future_files) != cfg.expected_future_files:
        raise InputError("Expected {} future files matching {} under {}, found {}".format(cfg.expected_future_files, cfg.future_pattern, cfg.future_dir, len(future_files)))
    pred_dir = cfg.output_dir / "future_predictions"; pred_dir.mkdir(parents=True, exist_ok=True)
    for target in targets:
        metrics_path = rf_dirs[target] / ("RF_metrics_by_split_" + target + ".csv")
        best = subprocess.run([cfg.rscript, "-e", "x<-read.csv(commandArgs(TRUE)[1]); x<-x[order(-x$cv_rsq_mean,-x$test_rsq),]; cat(x$split_id[1],x$mtry[1],x$min_n[1],x$trees[1])", str(metrics_path)], text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=subprocess_environment(cfg.rscript))
        if best.returncode: raise InputError("Cannot select best RF parameters for {}".format(target))
        split_id, mtry, min_n, trees = best.stdout.strip().split()
        # The RF script already renders the reference-style RMSE tuning panel
        # for every split. Promote the selected split's panel to a stable,
        # target-specific S24 filename instead of rebuilding a combined plot
        # from only the selected parameter rows.
        split_plot_dir = rf_dirs[target] / "plots_tuning_each_split"
        split_plot_stem = "Figure_Tuning_{}_split{:03d}_facet_SVD".format(target, int(split_id))
        tuning_png = split_plot_dir / (split_plot_stem + ".png")
        tuning_pdf = split_plot_dir / (split_plot_stem + ".pdf")
        if not tuning_png.is_file() or not tuning_pdf.is_file():
            raise InputError("RF tuning plot missing for {} split {}".format(target, split_id))
        shutil.copy2(tuning_png, cfg.output_dir / ("Figure_S24_{}_tuning.png".format("loadM" if target == "mean_loadM_relax" else "loadD")))
        shutil.copy2(tuning_pdf, cfg.output_dir / ("Figure_S24_{}_tuning.pdf".format("loadM" if target == "mean_loadM_relax" else "loadD")))
        for future in future_files:
            out = pred_dir / (future.stem + "." + target + ".csv")
            cmd = [cfg.rscript, str(pred_script), target, split_id, mtry, min_n, trees, str(training_csv), str(future), str(out), str(cfg.seed), str(cfg.workers)]
            result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.rscript))
            if result.returncode: raise InputError("future prediction failed for {} {}".format(target, future.name))
    # Avoid leaving the old combined S24 artifact in a reused output directory.
    for stale in (cfg.output_dir / "Figure_S24_load_tuning.png", cfg.output_dir / "Figure_S24_load_tuning.pdf"):
        stale.unlink(missing_ok=True)
    cmd = [cfg.rscript, str(plot_script), str(cfg.output_dir), str(cfg.mask or ""), ",".join(targets)]
    result = subprocess.run(cmd, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, env=subprocess_environment(cfg.rscript))
    (cfg.output_dir / "load_plot.log").write_text(result.stdout, encoding="utf-8")
    if result.returncode: raise InputError("load plotting failed; see load_plot.log")
    manifest.update({"status": "complete", "elapsed_seconds": round(time.time() - started, 3), "outputs": {"output_dir": str(cfg.output_dir), "training_csv": str(training_csv), "future_predictions": str(pred_dir), "rf_targets": list(targets)}})
    manifest["config_sha256"] = hashlib.sha256(cfg.config_path.read_bytes()).hexdigest()
    (cfg.output_dir / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    return manifest
