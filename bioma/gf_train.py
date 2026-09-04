"""Train an auditable Gradient Forest model from standardized BioMA tables."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import platform
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from . import __version__
from .gf_frequency import InputError
from .runtime import subprocess_environment


DEFAULT_PREDICTORS = tuple("bio{}".format(i) for i in range(1, 20))


@dataclass
class TrainingInputs:
    population_ids: List[str]
    predictors: List[str]
    responses: List[str]
    environment: Dict[str, List[float]]
    frequencies: Dict[str, List[float]]
    frequency_row_indexes: Dict[str, int]
    environment_row_indexes: Dict[str, int]


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _read_table(path: Path, label: str) -> Tuple[List[str], List[List[str]]]:
    if not path.is_file():
        raise InputError("{} file does not exist: {}".format(label, path))
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = list(reader)
    if not rows or not rows[0]:
        raise InputError("{} table is empty".format(label))
    header = rows[0]
    if len(header) != len(set(header)):
        raise InputError("{} table contains duplicate column names".format(label))
    for line_number, row in enumerate(rows[1:], start=2):
        if len(row) != len(header):
            raise InputError(
                "{} line {} has {} columns; expected {}".format(
                    label, line_number, len(row), len(header)
                )
            )
    if len(rows) == 1:
        raise InputError("{} table contains no data rows".format(label))
    return header, rows[1:]


def _numeric(value: str, label: str) -> float:
    if value.strip().upper() in ("", "NA", "NAN", "."):
        raise InputError("Missing numeric value: {}".format(label))
    try:
        result = float(value)
    except ValueError:
        raise InputError("Non-numeric value: {}".format(label))
    if not math.isfinite(result):
        raise InputError("Non-finite numeric value: {}".format(label))
    return result


def read_training_inputs(
    frequency_path: Path,
    environment_path: Path,
    predictors: Sequence[str] = DEFAULT_PREDICTORS,
    expected_sites: Optional[int] = None,
) -> TrainingInputs:
    """Validate inputs and align both tables by explicit population ID."""
    frequency_header, frequency_rows = _read_table(frequency_path, "Frequency")
    environment_header, environment_rows = _read_table(environment_path, "Environment")
    if not frequency_header or frequency_header[0] != "population_id":
        raise InputError("Frequency table first column must be population_id")
    if "population_id" not in environment_header:
        raise InputError("Environment table must contain population_id")
    if not predictors:
        raise InputError("At least one predictor is required")
    if len(predictors) != len(set(predictors)):
        raise InputError("Predictor names must be unique")
    missing_predictors = [name for name in predictors if name not in environment_header]
    if missing_predictors:
        raise InputError(
            "Environment table is missing predictors: {}".format(
                ", ".join(missing_predictors)
            )
        )

    responses = frequency_header[1:]
    if not responses:
        raise InputError("Frequency table contains no response loci")
    if expected_sites is not None and len(responses) != expected_sites:
        raise InputError(
            "Frequency table contains {} loci; expected {}".format(
                len(responses), expected_sites
            )
        )

    frequency_data: Dict[str, List[float]] = {}
    frequency_row_indexes: Dict[str, int] = {}
    for row_index, row in enumerate(frequency_rows, start=1):
        population_id = row[0].strip()
        if not population_id:
            raise InputError("Empty population_id in frequency table row {}".format(row_index))
        if population_id in frequency_data:
            raise InputError("Duplicate frequency population_id: {}".format(population_id))
        values = [
            _numeric(value, "frequency {}/{}".format(population_id, response))
            for response, value in zip(responses, row[1:])
        ]
        for response, value in zip(responses, values):
            if not 0.0 <= value <= 1.0:
                raise InputError(
                    "Frequency outside [0,1] for population {} locus {}".format(
                        population_id, response
                    )
                )
        frequency_data[population_id] = values
        frequency_row_indexes[population_id] = row_index

    population_column = environment_header.index("population_id")
    predictor_columns = [environment_header.index(name) for name in predictors]
    environment_data: Dict[str, List[float]] = {}
    environment_row_indexes: Dict[str, int] = {}
    environment_order: List[str] = []
    for row_index, row in enumerate(environment_rows, start=1):
        population_id = row[population_column].strip()
        if not population_id:
            raise InputError("Empty population_id in environment row {}".format(row_index))
        if population_id in environment_data:
            raise InputError("Duplicate environment population_id: {}".format(population_id))
        values = [
            _numeric(row[column], "environment {}/{}".format(population_id, predictor))
            for predictor, column in zip(predictors, predictor_columns)
        ]
        environment_data[population_id] = values
        environment_row_indexes[population_id] = row_index
        environment_order.append(population_id)

    frequency_ids = set(frequency_data)
    environment_ids = set(environment_data)
    if frequency_ids != environment_ids:
        only_frequency = sorted(frequency_ids - environment_ids)
        only_environment = sorted(environment_ids - frequency_ids)
        details = []
        if only_frequency:
            details.append("frequency-only={}".format(",".join(only_frequency[:10])))
        if only_environment:
            details.append("environment-only={}".format(",".join(only_environment[:10])))
        raise InputError("Population sets do not match: {}".format("; ".join(details)))

    if len(environment_order) < 5:
        raise InputError("Gradient Forest training requires at least five populations")
    for predictor_index, predictor in enumerate(predictors):
        values = [environment_data[population_id][predictor_index] for population_id in environment_order]
        if min(values) == max(values):
            raise InputError("Predictor has zero variance: {}".format(predictor))

    return TrainingInputs(
        population_ids=environment_order,
        predictors=list(predictors),
        responses=responses,
        environment=environment_data,
        frequencies=frequency_data,
        frequency_row_indexes=frequency_row_indexes,
        environment_row_indexes=environment_row_indexes,
    )


def _resolve_rscript(requested: Optional[str]) -> str:
    candidates = []
    if requested:
        candidates.append(requested)
    discovered = shutil.which("Rscript")
    if discovered:
        candidates.append(discovered)
    checked = set()
    for candidate in candidates:
        if candidate in checked:
            continue
        checked.add(candidate)
        try:
            result = subprocess.run(
                [
                    candidate,
                    "--vanilla",
                    "-e",
                    "suppressPackageStartupMessages(library(gradientForest))",
                ],
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
        "No Rscript with the 'gradientForest' package was found; provide one with --rscript"
    )


def _read_key_values(path: Path) -> Dict[str, str]:
    values: Dict[str, str] = {}
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader, None)
        if header != ["key", "value"]:
            raise InputError("Malformed Gradient Forest summary output")
        for row in reader:
            if len(row) == 2:
                values[row[0]] = row[1]
    return values


def train_gradient_forest(
    frequency_path: Path,
    environment_path: Path,
    output_dir: Path,
    predictors: Sequence[str] = DEFAULT_PREDICTORS,
    ntree: int = 500,
    nbin: int = 1001,
    corr_threshold: float = 0.5,
    max_level: Optional[float] = None,
    seed: int = 1,
    expected_sites: Optional[int] = None,
    rscript: Optional[str] = None,
) -> Dict[str, object]:
    frequency_path = frequency_path.expanduser().resolve()
    environment_path = environment_path.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise InputError("Output directory already exists; choose a new directory: {}".format(output_dir))
    if ntree < 1:
        raise InputError("ntree must be at least 1")
    if nbin < 2:
        raise InputError("nbin must be at least 2")
    if not 0.0 <= corr_threshold <= 1.0:
        raise InputError("corr_threshold must be between 0 and 1")

    inputs = read_training_inputs(
        frequency_path,
        environment_path,
        predictors=predictors,
        expected_sites=expected_sites,
    )
    if max_level is None:
        max_level = math.log(0.368 * len(inputs.population_ids) / 2.0, 2)
    if max_level <= 0:
        raise InputError("max_level must be positive")
    resolved_rscript = _resolve_rscript(rscript)

    warnings: List[Tuple[str, str, str]] = []
    ratio = len(inputs.population_ids) / float(len(inputs.predictors))
    if ratio < 2.0:
        warnings.append(
            (
                "LOW_POPULATION_TO_PREDICTOR_RATIO",
                "training",
                "{} populations for {} predictors (ratio {:.3f}); model uncertainty may be high".format(
                    len(inputs.population_ids), len(inputs.predictors), ratio
                ),
            )
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent))
    )
    try:
        training_matrix = temporary_dir / ".training_matrix.tsv"
        _write_tsv(
            training_matrix,
            ["population_id"] + inputs.predictors + inputs.responses,
            (
                [population_id]
                + ["{:.17g}".format(value) for value in inputs.environment[population_id]]
                + ["{:.17g}".format(value) for value in inputs.frequencies[population_id]]
                for population_id in inputs.population_ids
            ),
        )
        _write_tsv(
            temporary_dir / "training_alignment.tsv",
            [
                "training_row",
                "population_id",
                "environment_input_row",
                "frequency_input_row",
            ],
            (
                [
                    index,
                    population_id,
                    inputs.environment_row_indexes[population_id],
                    inputs.frequency_row_indexes[population_id],
                ]
                for index, population_id in enumerate(inputs.population_ids, start=1)
            ),
        )

        helper_script = Path(__file__).resolve().parent / "scripts" / "train_gradient_forest.R"
        if not helper_script.is_file():
            raise InputError("Bundled GF training helper is missing: {}".format(helper_script))
        command = [
            resolved_rscript,
            "--vanilla",
            str(helper_script),
            str(training_matrix),
            str(temporary_dir),
            str(ntree),
            str(nbin),
            "{:.17g}".format(corr_threshold),
            "{:.17g}".format(max_level),
            str(seed),
            ",".join(inputs.predictors),
        ]
        result = subprocess.run(
            command,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            env=subprocess_environment(resolved_rscript),
        )
        (temporary_dir / "training_log.txt").write_text(
            result.stdout + ("\nSTDERR\n" + result.stderr if result.stderr else ""),
            encoding="utf-8",
        )
        if result.returncode != 0:
            detail = result.stderr.strip() or result.stdout.strip() or "unknown R error"
            raise InputError("Gradient Forest training failed: {}".format(detail))

        summary = _read_key_values(temporary_dir / ".r_summary.tsv")
        positive_responses = int(summary["positive_responses"])
        if positive_responses == 0:
            warnings.append(
                (
                    "NO_POSITIVE_RESPONSE_R2",
                    "training",
                    "No response locus had positive out-of-bag R-squared",
                )
            )
        elif positive_responses / float(len(inputs.responses)) < 0.5:
            warnings.append(
                (
                    "LOW_POSITIVE_RESPONSE_FRACTION",
                    "training",
                    "Only {} of {} loci had positive out-of-bag R-squared".format(
                        positive_responses, len(inputs.responses)
                    ),
                )
            )
        _write_tsv(temporary_dir / "warnings.tsv", ["warning_code", "scope", "message"], warnings)

        for private_file in (training_matrix, temporary_dir / ".r_summary.tsv"):
            private_file.unlink()

        model_path = temporary_dir / "all_gfmod.data"
        if not model_path.is_file() or model_path.stat().st_size == 0:
            raise InputError("Gradient Forest model file was not created")
        manifest: Dict[str, object] = {
            "module": "gf-train",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": "completed_with_warning" if warnings else "completed_validated",
            "inputs": {
                "frequencies": str(frequency_path),
                "frequencies_sha256": _sha256(frequency_path),
                "environment": str(environment_path),
                "environment_sha256": _sha256(environment_path),
            },
            "parameters": {
                "predictors": list(inputs.predictors),
                "ntree": ntree,
                "nbin": nbin,
                "corr_threshold": corr_threshold,
                "compact": True,
                "max_level": max_level,
                "max_level_formula": "log2(0.368 * n_populations / 2)" if max_level == math.log(0.368 * len(inputs.population_ids) / 2.0, 2) else None,
                "seed": seed,
            },
            "counts": {
                "populations": len(inputs.population_ids),
                "predictors": len(inputs.predictors),
                "response_loci": len(inputs.responses),
                "positive_oob_r2_loci": positive_responses,
                "warnings": len(warnings),
            },
            "model": {
                "file": "all_gfmod.data",
                "size_bytes": model_path.stat().st_size,
                "sha256": _sha256(model_path),
            },
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
                "rscript": resolved_rscript,
                "r": summary.get("r_version"),
                "gradientForest": summary.get("gradientforest_version"),
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
