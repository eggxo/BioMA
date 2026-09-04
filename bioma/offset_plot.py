"""Create publication-ready maps and relationship plots from final offset tables."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import shutil
import subprocess
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Optional, Sequence

from . import __version__
from .gf_frequency import InputError
from .runtime import subprocess_environment


REQUIRED_COLUMNS = {
    "scenario", "period", "ssp", "model", "radius_km", "lon", "lat",
    "local_offset", "forward_offset", "reverse_offset",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _resolve_input(path: Path) -> Path:
    path = path.expanduser().resolve()
    if path.is_dir():
        path = path / "all_offsets.tsv.gz"
    if not path.is_file():
        raise InputError("Final offset table does not exist: {}".format(path))
    return path


def _validate_header(path: Path) -> None:
    import gzip

    with path.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
    opener = gzip.open if compressed else open
    with opener(str(path), "rt", encoding="utf-8", newline="") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if header is None or not REQUIRED_COLUMNS.issubset(header):
        missing = sorted(REQUIRED_COLUMNS - set(header or []))
        raise InputError("Offset table is missing columns: {}".format(", ".join(missing)))


def _r_subprocess_env(rscript: str) -> Dict[str, str]:
    """Isolate a conda R runtime from incompatible libraries in the login shell."""
    return subprocess_environment(rscript)


def _resolve_rscript(requested: Optional[str]) -> str:
    candidates = []
    if requested:
        candidates.append(requested)
    discovered = shutil.which("Rscript")
    if discovered:
        candidates.append(discovered)
    project_environment = Path(__file__).resolve().parents[1] / "envs" / "bioma-plot" / "bin" / "Rscript"
    candidates.extend(
        [
            str(project_environment),
        ]
    )
    expression = ";".join(
        "suppressPackageStartupMessages(library({}))".format(package)
        for package in ("data.table", "ggplot2", "sf", "cowplot", "ragg")
    )
    checked = set()
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
                env=_r_subprocess_env(candidate),
            )
        except FileNotFoundError:
            continue
        if result.returncode == 0:
            return candidate
    raise InputError("No Rscript with data.table, ggplot2, sf, cowplot, and ragg was found")


def _read_key_values(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        if next(reader, None) != ["key", "value"]:
            raise InputError("Malformed plot summary output")
        for row in reader:
            if len(row) == 2:
                values[row[0]] = row[1]
    return values


def plot_offsets(
    inputs: Sequence[Path],
    output_dir: Path,
    radius: str = "unlimited",
    models: Sequence[str] = (),
    period: Optional[str] = None,
    ssp: Optional[str] = None,
    populations: Optional[Path] = None,
    boundary: Optional[Path] = None,
    minimum_models: int = 2,
    rscript: Optional[str] = None,
) -> Dict[str, object]:
    if not inputs:
        raise InputError("At least one final offset table is required")
    if minimum_models < 2:
        raise InputError("minimum_models must be at least 2 for ensemble plotting")
    models = tuple(dict.fromkeys(models))
    if models and len(models) < minimum_models:
        raise InputError(
            "Ensemble plotting requires at least {} distinct models".format(minimum_models)
        )
    resolved_inputs = [_resolve_input(path) for path in inputs]
    for path in resolved_inputs:
        _validate_header(path)
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise InputError("Output directory already exists; choose a new directory: {}".format(output_dir))
    if populations is not None:
        populations = populations.expanduser().resolve()
        if not populations.is_file():
            raise InputError("Population coordinate file does not exist: {}".format(populations))
    if boundary is not None:
        boundary = boundary.expanduser().resolve()
        if not boundary.is_file():
            raise InputError("Boundary file does not exist: {}".format(boundary))
    resolved_rscript = _resolve_rscript(rscript)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent))
    )
    try:
        helper = Path(__file__).resolve().parent / "scripts" / "plot_gf_offsets.R"
        if not helper.is_file():
            raise InputError("Bundled offset plotting helper is missing: {}".format(helper))
        command = [
            resolved_rscript, "--vanilla", str(helper), str(temporary_dir), radius,
            ",".join(models), period or "", ssp or "",
            str(populations) if populations else "", str(boundary) if boundary else "",
            str(minimum_models),
        ] + [str(path) for path in resolved_inputs]
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=_r_subprocess_env(resolved_rscript),
        )
        (temporary_dir / "plot_log.txt").write_text(
            result.stdout + ("\nSTDERR\n" + result.stderr if result.stderr else ""),
            encoding="utf-8",
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip() or "unknown R error"
            raise InputError("Offset plotting failed: {}".format(detail))
        summary = _read_key_values(temporary_dir / ".plot_summary.tsv")
        (temporary_dir / ".plot_summary.tsv").unlink()

        output_files = {}
        for name in (
            "ensemble_mean_offsets.tsv.gz",
            "all_offsets_maps.png", "all_offsets_maps.pdf",
            "offset_relationships.png", "offset_relationships.pdf",
            "plot_data.tsv.gz",
        ):
            path = temporary_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                raise InputError("Expected plot output was not created: {}".format(name))
            output_files[name] = {"size_bytes": path.stat().st_size, "sha256": _sha256(path)}

        manifest: Dict[str, object] = {
            "module": "offset-plot",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": "completed_validated",
            "inputs": [
                {"path": str(path), "sha256": _sha256(path)} for path in resolved_inputs
            ],
            "parameters": {
                "radius": radius,
                "models": list(models),
                "minimum_models": minimum_models,
                "period": period,
                "ssp": ssp,
                "populations": str(populations) if populations else None,
                "boundary": str(boundary) if boundary else None,
                "ensemble_aggregation": "unweighted arithmetic mean by coordinate; all selected models required",
                "map_scale": "min-max normalized within the selected ensemble",
                "rgb_channels": {"red": "local", "green": "forward", "blue": "reverse"},
            },
            "counts": {
                "input_rows": int(summary["input_rows"]),
                "selected_rows": int(summary["selected_rows"]),
                "plot_cells": int(summary["plot_cells"]),
                "models": int(summary["models"]),
                "full_local_coverage_cells": int(summary["full_local_coverage_cells"]),
                "full_forward_coverage_cells": int(summary["full_forward_coverage_cells"]),
                "full_reverse_coverage_cells": int(summary["full_reverse_coverage_cells"]),
            },
            "selection": {
                "period": summary["period"],
                "ssp": summary["ssp"],
                "radius": summary["radius"],
                "models": summary["model_names"].split(","),
            },
            "outputs": output_files,
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "rscript": resolved_rscript,
                "r": summary.get("r_version"),
                "ggplot2": summary.get("ggplot2_version"),
                "sf": summary.get("sf_version"),
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
