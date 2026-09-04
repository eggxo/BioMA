"""Plot ensemble-mean forward offset across migration-distance limits."""

from __future__ import annotations

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
from .offset_plot import (
    _read_key_values,
    _resolve_input,
    _resolve_rscript,
    _r_subprocess_env,
    _sha256,
    _validate_header,
)


def _normalize_ssp(value: str) -> str:
    normalized = value.strip().lower()
    if normalized.isdigit():
        normalized = "ssp" + normalized
    if not normalized.startswith("ssp") or not normalized[3:].isdigit():
        raise InputError("Invalid SSP name: {}".format(value))
    return normalized


def plot_forward_distance(
    inputs: Sequence[Path],
    output_dir: Path,
    models: Sequence[str] = (),
    period: Optional[str] = None,
    ssps: Sequence[str] = (),
    minimum_models: int = 2,
    rscript: Optional[str] = None,
) -> Dict[str, object]:
    if not inputs:
        raise InputError("At least one final offset table is required")
    if minimum_models < 2:
        raise InputError("minimum_models must be at least 2 for ensemble plotting")
    models = tuple(dict.fromkeys(value.strip() for value in models if value.strip()))
    if models and len(models) < minimum_models:
        raise InputError(
            "Forward-distance plotting requires at least {} distinct models".format(
                minimum_models
            )
        )
    normalized_ssps = tuple(dict.fromkeys(_normalize_ssp(value) for value in ssps))

    resolved_inputs = [_resolve_input(path) for path in inputs]
    for path in resolved_inputs:
        _validate_header(path)
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise InputError(
            "Output directory already exists; choose a new directory: {}".format(output_dir)
        )
    resolved_rscript = _resolve_rscript(rscript)

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent))
    )
    try:
        helper = Path(__file__).resolve().parent / "scripts" / "plot_forward_by_distance.R"
        if not helper.is_file():
            raise InputError("Bundled forward-distance plotting helper is missing: {}".format(helper))
        command = [
            resolved_rscript,
            "--vanilla",
            str(helper),
            str(temporary_dir),
            ",".join(models),
            period or "",
            ",".join(normalized_ssps),
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
            raise InputError("Forward-distance plotting failed: {}".format(detail))
        summary = _read_key_values(temporary_dir / ".plot_summary.tsv")
        (temporary_dir / ".plot_summary.tsv").unlink()

        output_files = {}
        for name in (
            "forward_offset_ensemble_grid.tsv.gz",
            "forward_offset_by_distance.tsv",
            "forward_offset_by_distance.png",
            "forward_offset_by_distance.pdf",
        ):
            path = temporary_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                raise InputError("Expected forward-distance output was not created: {}".format(name))
            output_files[name] = {"size_bytes": path.stat().st_size, "sha256": _sha256(path)}

        manifest: Dict[str, object] = {
            "module": "forward-distance-plot",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": "completed_validated",
            "inputs": [
                {"path": str(path), "sha256": _sha256(path)} for path in resolved_inputs
            ],
            "parameters": {
                "models": list(models),
                "minimum_models": minimum_models,
                "period": period,
                "ssps": list(normalized_ssps),
                "model_aggregation": (
                    "unweighted arithmetic mean by SSP, radius, and coordinate; "
                    "all selected models required"
                ),
                "spatial_summary": "median with 25th-75th percentile interval",
                "x_axis": "ordered migration-distance categories; Unlimited is last",
            },
            "counts": {
                "input_rows": int(summary["input_rows"]),
                "selected_rows": int(summary["selected_rows"]),
                "ensemble_grid_rows": int(summary["ensemble_grid_rows"]),
                "valid_ensemble_grid_rows": int(summary["valid_ensemble_grid_rows"]),
                "summary_rows": int(summary["summary_rows"]),
                "models": int(summary["models"]),
                "ssps": int(summary["ssps"]),
                "radii": int(summary["radii"]),
            },
            "selection": {
                "period": summary["period"],
                "models": summary["model_names"].split(","),
                "ssps": summary["ssp_names"].split(","),
                "radii": summary["radius_names"].split(","),
            },
            "outputs": output_files,
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "rscript": resolved_rscript,
                "r": summary.get("r_version"),
                "ggplot2": summary.get("ggplot2_version"),
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
