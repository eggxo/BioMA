"""Compute local, forward, and reverse Gradient Forest offsets."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Sequence, TextIO, Tuple

from . import __version__
from .gf_frequency import InputError
from .runtime import subprocess_environment


def parse_forward_radii(value: str) -> List[float]:
    radii: List[float] = []
    for token in value.split(","):
        token = token.strip().lower()
        if not token:
            continue
        if token in ("inf", "infinity", "unlimited", "all"):
            radius = math.inf
        else:
            try:
                radius = float(token)
            except ValueError:
                raise InputError("Invalid forward radius: {}".format(token))
            if not math.isfinite(radius) or radius <= 0:
                raise InputError("Forward radii must be positive kilometres or 'inf'")
        if radius not in radii:
            radii.append(radius)
    if not radii:
        raise InputError("At least one forward radius is required")
    return radii


def _open_text_auto(path: Path) -> TextIO:
    with path.open("rb") as raw:
        magic = raw.read(2)
    if magic == b"\x1f\x8b":
        return gzip.open(str(path), "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _profile_climate_table(path: Path) -> Dict[str, object]:
    if not path.is_file():
        raise InputError("Climate table does not exist: {}".format(path))
    coordinate_digest = hashlib.sha256()
    rows = 0
    missing_rows = 0
    with _open_text_auto(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header is None or len(header) < 3:
            raise InputError("Climate table has no usable header: {}".format(path))
        if header[0:2] != ["lon", "lat"]:
            raise InputError("Climate table must begin with lon and lat: {}".format(path))
        predictors = header[2:]
        if len(predictors) != len(set(predictors)):
            raise InputError("Climate table contains duplicate predictors: {}".format(path))
        for line_number, row in enumerate(reader, start=2):
            if len(row) != len(header):
                raise InputError("Malformed climate row {} in {}".format(line_number, path))
            try:
                lon = float(row[0])
                lat = float(row[1])
            except ValueError:
                raise InputError("Non-numeric coordinates at row {} in {}".format(line_number, path))
            if not (-180 <= lon <= 180 and -90 <= lat <= 90):
                raise InputError("Invalid coordinates at row {} in {}".format(line_number, path))
            coordinate_digest.update((row[0] + "\t" + row[1] + "\n").encode("utf-8"))
            row_missing = False
            for value in row[2:]:
                if value.strip().upper() in ("", "NA", "NAN", "."):
                    row_missing = True
                else:
                    try:
                        number = float(value)
                    except ValueError:
                        raise InputError("Non-numeric climate value at row {} in {}".format(line_number, path))
                    if not math.isfinite(number):
                        row_missing = True
            missing_rows += int(row_missing)
            rows += 1
    if rows == 0:
        raise InputError("Climate table contains no rows: {}".format(path))
    return {
        "rows": rows,
        "missing_rows": missing_rows,
        "predictors": predictors,
        "coordinate_sha256": coordinate_digest.hexdigest(),
    }


def _resolve_rscript(requested: Optional[str]) -> str:
    candidates = []
    if requested:
        candidates.append(requested)
    discovered = shutil.which("Rscript")
    if discovered:
        candidates.append(discovered)
    checked = set()
    expression = (
        "suppressPackageStartupMessages(library(gradientForest));"
        "suppressPackageStartupMessages(library(FNN));"
        "suppressPackageStartupMessages(library(data.table));"
        "suppressPackageStartupMessages(library(geosphere))"
    )
    for candidate in candidates:
        if candidate in checked:
            continue
        checked.add(candidate)
        try:
            result = subprocess.run(
                [candidate, "--vanilla", "-e", expression],
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                env=subprocess_environment(candidate),
            )
        except FileNotFoundError:
            continue
        if result.returncode == 0:
            return candidate
    raise InputError(
        "No Rscript with gradientForest, FNN, data.table, and geosphere was found"
    )


def _read_key_values(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["key", "value"]:
            raise InputError("Malformed GF offset summary output")
        for row in reader:
            if len(row) == 2:
                values[row[0]] = row[1]
    return values


def compute_gf_offsets(
    model: Path,
    climate_dir: Path,
    scenario: str,
    output_dir: Path,
    forward_radii_km: Sequence[float] = (100.0, 250.0, 500.0, 1000.0, math.inf),
    initial_knn_k: int = 64,
    batch_size: int = 10000,
    verify_sample: int = 10,
    tie_tolerance: float = 1e-12,
    rscript: Optional[str] = None,
) -> Dict[str, object]:
    model = model.expanduser().resolve()
    if model.is_dir():
        model = model / "all_gfmod.data"
    climate_dir = climate_dir.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if not model.is_file():
        raise InputError("Gradient Forest model does not exist: {}".format(model))
    if not climate_dir.is_dir():
        raise InputError("Climate preparation directory does not exist: {}".format(climate_dir))
    if Path(scenario).name != scenario or scenario in ("", ".", ".."):
        raise InputError("Scenario must be one directory-safe scenario name")
    if output_dir.exists():
        raise InputError("Output directory already exists; choose a new directory: {}".format(output_dir))
    if initial_knn_k < 1:
        raise InputError("initial_knn_k must be at least 1")
    if batch_size < 1:
        raise InputError("batch_size must be at least 1")
    if verify_sample < 0:
        raise InputError("verify_sample cannot be negative")
    if tie_tolerance < 0:
        raise InputError("tie_tolerance cannot be negative")
    radii = list(forward_radii_km)
    if not radii or any((not math.isinf(value) and value <= 0) for value in radii):
        raise InputError("Forward radii must be positive kilometres or infinity")

    current_path = climate_dir / "current_background.tsv.gz"
    future_local_path = climate_dir / "future_local" / (scenario + ".tsv.gz")
    future_search_path = climate_dir / "future_search" / (scenario + ".tsv.gz")
    current_profile = _profile_climate_table(current_path)
    local_profile = _profile_climate_table(future_local_path)
    search_profile = _profile_climate_table(future_search_path)
    if current_profile["predictors"] != local_profile["predictors"]:
        raise InputError("Current and future-local predictor columns differ")
    if current_profile["predictors"] != search_profile["predictors"]:
        raise InputError("Current and future-search predictor columns differ")
    if current_profile["rows"] != local_profile["rows"]:
        raise InputError("Current and future-local row counts differ")
    if current_profile["coordinate_sha256"] != local_profile["coordinate_sha256"]:
        raise InputError("Current and future-local coordinates are not identical and ordered")
    if current_profile["missing_rows"]:
        raise InputError("Current background contains incomplete climate rows")

    resolved_rscript = _resolve_rscript(rscript)
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent))
    )
    try:
        helper = Path(__file__).resolve().parent / "scripts" / "compute_gf_offsets.R"
        if not helper.is_file():
            raise InputError("Bundled GF offset helper is missing: {}".format(helper))
        radii_argument = ",".join("Inf" if math.isinf(value) else "{:.17g}".format(value) for value in radii)
        command = [
            resolved_rscript,
            "--vanilla",
            str(helper),
            str(model),
            str(current_path),
            str(future_local_path),
            str(future_search_path),
            str(temporary_dir),
            scenario,
            radii_argument,
            str(initial_knn_k),
            str(batch_size),
            str(verify_sample),
            "{:.17g}".format(tie_tolerance),
        ]
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=subprocess_environment(resolved_rscript),
        )
        (temporary_dir / "offset_log.txt").write_text(
            result.stdout + ("\nSTDERR\n" + result.stderr if result.stderr else ""),
            encoding="utf-8",
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip() or "unknown R error"
            raise InputError("Gradient Forest offset computation failed: {}".format(detail))
        summary = _read_key_values(temporary_dir / ".r_summary.tsv")
        (temporary_dir / ".r_summary.tsv").unlink()

        warnings: List[Tuple[str, str, str]] = []
        local_missing = int(summary["local_missing"])
        reverse_missing = int(summary["reverse_missing"])
        forward_unmatched = int(summary["forward_unmatched"])
        if local_missing:
            warnings.append(
                (
                    "LOCAL_OFFSET_NODATA",
                    scenario,
                    "{} current cells have no complete future-local climate and receive NA".format(local_missing),
                )
            )
        if reverse_missing:
            warnings.append(
                (
                    "REVERSE_OFFSET_NODATA",
                    scenario,
                    "{} future-search cells have incomplete climate and receive NA".format(reverse_missing),
                )
            )
        if forward_unmatched:
            warnings.append(
                (
                    "FORWARD_NO_CANDIDATE",
                    scenario,
                    "{} source-radius combinations have no valid destination".format(forward_unmatched),
                )
            )
        with (temporary_dir / "warnings.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["warning_code", "scope", "message"])
            writer.writerows(warnings)

        output_files = {}
        for name in (
            "local_offset.tsv.gz",
            "forward_offset.tsv.gz",
            "reverse_offset.tsv.gz",
            "all_offsets.tsv.gz",
            "current_transformed.tsv.gz",
            "future_local_transformed.tsv.gz",
            "future_search_transformed.tsv.gz",
            "scenario_summary.tsv",
        ):
            path = temporary_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                raise InputError("Expected offset output was not created: {}".format(name))
            output_files[name] = {
                "size_bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }

        manifest: Dict[str, object] = {
            "module": "gf-offset",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": "completed_with_warning" if warnings else "completed_validated",
            "scenario": scenario,
            "inputs": {
                "model": str(model),
                "model_sha256": _sha256(model),
                "climate_dir": str(climate_dir),
                "current_background": str(current_path),
                "current_background_sha256": _sha256(current_path),
                "future_local": str(future_local_path),
                "future_local_sha256": _sha256(future_local_path),
                "future_search": str(future_search_path),
                "future_search_sha256": _sha256(future_search_path),
            },
            "parameters": {
                "forward_radii_km": ["unlimited" if math.isinf(value) else value for value in radii],
                "initial_knn_k": initial_knn_k,
                "batch_size": batch_size,
                "verify_sample": verify_sample,
                "tie_tolerance": tie_tolerance,
                "tie_break": "minimum geographic distance, then minimum input row",
                "distance": "Euclidean distance in Gradient Forest transformed BIO space",
            },
            "counts": {
                "current_cells": int(summary["current_cells"]),
                "future_local_cells": int(summary["future_local_cells"]),
                "future_search_cells": int(summary["future_search_cells"]),
                "future_search_valid_cells": int(summary["future_search_valid_cells"]),
                "combined_rows": int(summary["combined_rows"]),
                "local_missing": local_missing,
                "reverse_missing": reverse_missing,
                "forward_unmatched": forward_unmatched,
                "verification_cases": int(summary["verification_cases"]),
                "warnings": len(warnings),
            },
            "algorithm": {
                "nearest_neighbor": "exact FNN kd-tree with expanding k",
                "maximum_k_used": int(summary["maximum_k_used"]),
                "verification": "brute-force sample comparison",
            },
            "outputs": output_files,
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "rscript": resolved_rscript,
                "r": summary.get("r_version"),
                "gradientForest": summary.get("gradientforest_version"),
                "FNN": summary.get("fnn_version"),
                "elapsed_seconds": float(summary["elapsed_seconds"]),
            },
        }
        with (temporary_dir / "run_manifest.json").open("w", encoding="utf-8") as handle:
            json.dump(manifest, handle, indent=2, ensure_ascii=False)
            handle.write("\n")
        os.replace(str(temporary_dir), str(output_dir))
        return manifest
    except Exception:
        shutil.rmtree(str(temporary_dir), ignore_errors=True)
        raise
