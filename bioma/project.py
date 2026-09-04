"""Project-level orchestration for the BioMA workflow modules."""

from __future__ import annotations

import configparser
import csv
import hashlib
import html
import json
import os
import platform
import shutil
import sys
import time
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

from . import __version__
from .gf_frequency import InputError


MODULE_ORDER: Tuple[str, ...] = (
    "gf",
    "rona",
    "mar",
    "load",
    "niche",
    "wfmoment",
    "vulnerability",
)

MODULE_DIRECTORIES: Mapping[str, str] = {
    "gf": "01_gf",
    "rona": "02_rona",
    "mar": "03_mar",
    "load": "04_load",
    "niche": "05_niche",
    "wfmoment": "06_wfmoment",
    "vulnerability": "07_vulnerability",
}

VULNERABILITY_DEPENDENCIES: Mapping[str, str] = {
    "gf_offsets": "gf",
    "rona_ensemble": "rona",
    "load_predictors": "load",
    "niche_raster": "niche",
}

FALSE_VALUES = {"", "0", "false", "no", "off", "disabled", "none"}
TRUE_VALUES = {"1", "true", "yes", "on", "enabled"}

PROJECT_INPUT_KEYS = {
    "adaptive_vcf",
    "adaptive_frequency",
    "whole_genome_vcf",
    "load_vcf_dir",
    "population_samples",
    "population_environment",
    "population_dir",
    "rona_unld_dir",
    "climate_current",
    "climate_future",
    "species_mask",
    "load_predictors",
    "load_future_dir",
    "mar_lonlat",
    "occurrence_csv",
    "maxent_jar",
    "wfmoment_current_raster",
    "pi_file",
    "structure_file",
    "species_parameters",
    "area_file",
    "future_masks_json",
}

PROJECT_INPUT_TARGETS = {
    "adaptive_vcf": (("gf", "inputs", "vcf"),),
    "population_environment": (("gf", "inputs", "coordinates"), ("rona", "inputs", "environment")),
    "climate_current": (("gf", "inputs", "present_climate"), ("niche", "inputs", "current_env_dir")),
    "climate_future": (("gf", "inputs", "future_climate"), ("rona", "inputs", "future_climate"), ("niche", "inputs", "future_root")),
    "species_mask": (("gf", "inputs", "current_mask"), ("rona", "inputs", "mask"), ("niche", "inputs", "mask_shp"), ("load", "inputs", "mask")),
    "whole_genome_vcf": (("mar", "inputs", "vcf"),),
    "load_vcf_dir": (("load", "inputs", "vcf_dir"),),
    "population_dir": (("load", "inputs", "population_dir"),),
    "rona_unld_dir": (("rona", "inputs", "unld_dir"),),
    "load_predictors": (("load", "inputs", "predictors"),),
    "load_future_dir": (("load", "inputs", "future_dir"),),
    "mar_lonlat": (("mar", "inputs", "lonlat"),),
    "occurrence_csv": (("niche", "inputs", "occurrence_csv"),),
    "maxent_jar": (("niche", "inputs", "maxent_jar"),),
    "wfmoment_current_raster": (("wfmoment", "inputs", "current_raster"),),
    "pi_file": (("wfmoment", "inputs", "pi_file"),),
    "structure_file": (("wfmoment", "inputs", "structure_file"),),
    "species_parameters": (("wfmoment", "inputs", "param_file"),),
    "area_file": (("wfmoment", "inputs", "area_file"),),
    "future_masks_json": (("wfmoment", "inputs", "future_masks_json"),),
}


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _sha256_file(path: Path) -> str:
    """Return a SHA-256 digest without loading the file into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _canonical_json_bytes(value: object) -> bytes:
    """Encode a manifest record deterministically for aggregate hashing."""
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _shapefile_components(path: Path) -> List[Path]:
    """Return a shapefile and all same-stem sidecars in stable order.

    A shapefile is a logical dataset spread across several files (normally
    ``.shp``, ``.shx``, ``.dbf`` and ``.prj``).  Hashing only the ``.shp``
    would allow a changed attribute table or projection to go unnoticed.
    ``Path.name`` matching deliberately also catches less common sidecars such
    as ``.cpg``, ``.qix`` and ``.shp.xml``.
    """
    if path.suffix.lower() != ".shp":
        return [path]
    # ``foo.shp`` has sidecars such as ``foo.shx`` and ``foo.dbf``; matching
    # the full filename (``foo.shp.*``) would miss those normal companions.
    stem = path.stem.casefold()
    sidecar_suffixes = {
        ".shp",
        ".shx",
        ".dbf",
        ".prj",
        ".cpg",
        ".qpj",
        ".qix",
        ".fix",
        ".sbn",
        ".sbx",
        ".fbn",
        ".fbx",
        ".ain",
        ".aih",
        ".atx",
        ".ixs",
        ".mxs",
        ".xml",
    }
    candidates = [
        item
        for item in path.parent.iterdir()
        if item.is_file()
        and (
            item == path
            or item.name.casefold() == path.name.casefold()
            or (
                item.stem.casefold() == stem
                and item.suffix.casefold() in sidecar_suffixes
            )
            or item.name.casefold() == stem + ".shp.xml"
        )
    ]
    # Include the primary path even when a platform presents a case variant.
    if path not in candidates:
        candidates.append(path)
    return sorted({item.resolve() for item in candidates}, key=lambda item: item.name.casefold())


def _iter_directory_entries(root: Path) -> Iterable[Dict[str, object]]:
    """Yield stable, content-addressed entries for a directory tree.

    Symlinks are recorded as links rather than followed, preventing cycles and
    making a manifest independent of files outside the declared input root.
    Regular files are hashed in streaming chunks.  Relative POSIX paths are
    used in the aggregate so moving a directory does not alter its content
    digest.
    """
    root = root.resolve()

    def visit(directory: Path, relative: str) -> Iterable[Dict[str, object]]:
        try:
            entries = sorted(directory.iterdir(), key=lambda item: item.name.casefold())
        except OSError as error:
            raise InputError("Cannot read project input directory {}: {}".format(directory, error))
        for item in entries:
            rel = "{}/{}".format(relative, item.name) if relative else item.name
            try:
                if item.is_symlink():
                    yield {
                        "relative_path": rel.replace(os.sep, "/"),
                        "kind": "symlink",
                        "target": os.readlink(str(item)).replace(os.sep, "/"),
                    }
                elif item.is_dir():
                    # Empty directories are represented as entries so the
                    # manifest remains a faithful description of the input.
                    children = iter(visit(item, rel))
                    try:
                        first_child = next(children)
                    except StopIteration:
                        yield {"relative_path": rel.replace(os.sep, "/"), "kind": "directory"}
                    else:
                        yield first_child
                        yield from children
                elif item.is_file():
                    stat = item.stat()
                    yield {
                        "relative_path": rel.replace(os.sep, "/"),
                        "kind": "file",
                        "size_bytes": stat.st_size,
                        "sha256": _sha256_file(item),
                    }
                else:
                    # Special files are not valid scientific inputs, but
                    # recording their type makes an accidental replacement
                    # visible without attempting to read from a device.
                    yield {"relative_path": rel.replace(os.sep, "/"), "kind": "special"}
            except OSError as error:
                raise InputError("Cannot inspect project input {}: {}".format(item, error))

    yield from visit(root, "")


def _aggregate_manifest_hash(entries: Iterable[Mapping[str, object]]) -> str:
    digest = hashlib.sha256()
    for entry in entries:
        digest.update(_canonical_json_bytes(entry))
        digest.update(b"\n")
    return digest.hexdigest()


def _fingerprint_path(path: Path) -> Dict[str, object]:
    """Describe a file or directory and return a stable content fingerprint.

    The returned directory record includes per-file hashes for the independent
    provenance manifest.  Callers can use :func:`_compact_fingerprint` when
    storing a small summary in ``resolved_project.json``.
    """
    path = path.expanduser().resolve()
    if path.is_dir():
        entries = list(_iter_directory_entries(path))
        # ``entries`` is sorted by traversal, but sort again to guard against
        # platform-specific directory ordering and future walker changes.
        entries.sort(
            key=lambda item: (
                str(item.get("relative_path", "")).casefold(),
                str(item.get("relative_path", "")),
            )
        )
        total_size = sum(int(item.get("size_bytes", 0) or 0) for item in entries)
        return {
            "path": str(path),
            "kind": "directory",
            "sha256": _aggregate_manifest_hash(entries),
            "entry_count": len(entries),
            "total_size_bytes": total_size,
            "files": entries,
        }
    if path.is_file():
        try:
            components = _shapefile_components(path)
        except OSError as error:
            raise InputError("Cannot inspect project input {}: {}".format(path, error))
        component_rows: List[Dict[str, object]] = []
        for component in components:
            try:
                stat = component.stat()
                component_rows.append(
                    {
                        "relative_path": component.name,
                        "kind": "file",
                        "size_bytes": stat.st_size,
                        "sha256": _sha256_file(component),
                    }
                )
            except OSError as error:
                raise InputError("Cannot read project input {}: {}".format(component, error))
        component_rows.sort(key=lambda item: str(item["relative_path"]).casefold())
        aggregate = _aggregate_manifest_hash(component_rows)
        primary = next(
            (
                row
                for row in component_rows
                if Path(str(row["relative_path"])).name.casefold() == path.name.casefold()
            ),
            component_rows[0],
        )
        result: Dict[str, object] = {
            "path": str(path),
            "kind": "file",
            "size_bytes": int(primary["size_bytes"]),
            "sha256": aggregate if len(component_rows) > 1 else str(primary["sha256"]),
            "primary_sha256": str(primary["sha256"]),
            "component_count": len(component_rows),
        }
        if len(component_rows) > 1:
            result["components"] = component_rows
        return result
    return {"path": str(path), "kind": "missing"}


def _compact_fingerprint(record: Mapping[str, object]) -> Dict[str, object]:
    """Drop per-entry detail while retaining the content identity summary."""
    return {
        str(key): value
        for key, value in record.items()
        if key not in {"files", "components"}
    }


def _fingerprint_paths(paths: Mapping[str, Path]) -> Dict[str, Dict[str, object]]:
    """Fingerprint paths with a cache so shared rasters are hashed once."""
    cache: Dict[str, Dict[str, object]] = {}
    result: Dict[str, Dict[str, object]] = {}
    for key, path in sorted(paths.items()):
        resolved = str(path.expanduser().resolve())
        if resolved not in cache:
            cache[resolved] = _fingerprint_path(Path(resolved))
        result[key] = cache[resolved]
    return result


def _write_input_manifest(path: Path, fingerprints: Mapping[str, Mapping[str, object]]) -> None:
    """Write the detailed input provenance manifest atomically.

    ``resolved_project.json`` intentionally stores only compact records.  This
    sidecar retains every directory entry and shapefile component so a reader
    can audit exactly which files contributed to a run.
    """
    payload = {
        "schema_version": 1,
        "hash_algorithm": "sha256",
        "directory_hash": "sha256(canonical-json-lines)",
        "inputs": {key: fingerprints[key] for key in sorted(fingerprints)},
    }
    _write_json_atomic(path, payload)


def _check_existing_generated_inputs(
    existing: Mapping[str, object],
    project_dir: Path,
) -> None:
    """Reject silently edited generated compatibility files before writes."""
    records = existing.get("resolved_input_fingerprints")
    paths = existing.get("resolved_inputs")
    if not isinstance(records, Mapping) or not isinstance(paths, Mapping):
        return
    generated_keys = {
        "gf_samples",
        "generated_population_dir",
        "generated_load_predictors",
        "generated_occurrence",
        "generated_mar_lonlat",
    }
    for key in generated_keys.intersection(records.keys(), paths.keys()):
        raw_path = paths.get(key)
        expected = records.get(key)
        if not isinstance(raw_path, str) or not isinstance(expected, Mapping):
            continue
        path = Path(raw_path).expanduser().resolve()
        # Older metadata may contain a relative generated path.  Resolve it
        # relative to the metadata directory when it is not found as-is.
        if not path.exists() and not Path(raw_path).is_absolute():
            path = (project_dir / raw_path).resolve()
        current = _compact_fingerprint(_fingerprint_path(path))
        expected_compact = _compact_fingerprint(expected)
        if current != expected_compact:
            raise InputError(
                "Generated project input changed outside BioMA: {}. "
                "Use --overwrite to rebuild compatibility files.".format(path)
            )


def _check_existing_input_manifest(existing: Mapping[str, object], project_dir: Path) -> None:
    """Verify the provenance sidecar itself was not edited or truncated."""
    expected = existing.get("input_manifest_sha256")
    if not isinstance(expected, str) or not expected:
        return
    raw_path = existing.get("input_manifest")
    manifest_path = Path(raw_path).expanduser() if isinstance(raw_path, str) and raw_path else project_dir / "input_manifest.json"
    if not manifest_path.is_absolute():
        manifest_path = project_dir / manifest_path
    manifest_path = manifest_path.resolve()
    if not manifest_path.is_file():
        raise InputError("Recorded project input manifest is missing: {}".format(manifest_path))
    actual = _sha256_file(manifest_path)
    if actual != expected:
        raise InputError(
            "Project input manifest changed outside BioMA: {}. "
            "Use --overwrite to regenerate provenance metadata.".format(manifest_path)
        )


def _configured_input_paths(config: ProjectConfig) -> Dict[str, Path]:
    """Collect external paths referenced by module ``[inputs]`` sections.

    Unified ``[inputs]`` values are handled separately.  This additional scan
    covers legacy standalone INI files so changing a VCF, raster directory, or
    MaxEnt jar referenced there also invalidates the project signature.
    Missing paths are retained as ``missing`` fingerprints later, allowing a
    newly-created file to be detected on the next invocation.
    """
    paths: Dict[str, Path] = {}
    overridden_targets = {
        (module, section, option)
        for input_key, targets in PROJECT_INPUT_TARGETS.items()
        if input_key in config.inputs
        for module, section, option in targets
    }
    # Canonical tables materialize a few compatibility options whose names do
    # not appear verbatim in PROJECT_INPUT_TARGETS.
    if "population_samples" in config.inputs:
        overridden_targets.add(("gf", "inputs", "samples"))
        overridden_targets.add(("load", "inputs", "population_dir"))
        overridden_targets.add(("mar", "inputs", "lonlat"))
    if "population_environment" in config.inputs:
        overridden_targets.add(("load", "inputs", "predictors"))
        overridden_targets.add(("niche", "inputs", "occurrence_csv"))
    for module, source in config.module_configs.items():
        parser = _read_ini(source)
        if not parser.has_section("inputs"):
            continue
        for key, raw in parser.items("inputs"):
            # work_dir is an execution scratch/output root, not a scientific
            # input.  Recursing through it would hash logs and prior results
            # and make otherwise identical projects non-reproducible.
            if key.lower() in {"work_dir"}:
                continue
            if (module, "inputs", key) in overridden_targets:
                continue
            # In project mode RONA consumes the GF-generated ALT frequency
            # table whenever an adaptive VCF and GF are selected; a legacy
            # path in the standalone RONA INI is therefore not an input.
            if (
                module == "rona"
                and key.lower() == "alt_frequency"
                and (
                    "adaptive_frequency" in config.inputs
                    or ("adaptive_vcf" in config.inputs and "gf" in config.module_configs)
                )
            ):
                continue
            value = raw.strip()
            if not value or value.lower() in FALSE_VALUES or value.lower() == "auto":
                continue
            candidate = Path(value).expanduser()
            if not candidate.is_absolute():
                candidate = source.parent / candidate
            paths["module.{}.inputs.{}".format(module, key)] = candidate.resolve()
    return paths


def _json_sha256(value: object) -> str:
    encoded = json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    return _sha256_bytes(encoded.encode("utf-8"))


def _write_json_atomic(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name("{}.{}.tmp".format(path.name, os.getpid()))
    # Stream JSON to disk so a directory manifest is not duplicated as one
    # enormous intermediate string in memory.
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    temporary.replace(path)


def _resolve_path(value: str, base: Path, label: str) -> Path:
    value = (value or "").strip()
    if not value:
        raise InputError("Missing {}".format(label))
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _split(value: str) -> Tuple[str, ...]:
    return tuple(part.strip() for part in (value or "").split(",") if part.strip())


def _normalise_ssp(value: str) -> str:
    normalised = value.strip().lower()
    return normalised[3:] if normalised.startswith("ssp") else normalised


def _bool(value: str, default: bool) -> bool:
    normalised = (value or "").strip().lower()
    if not normalised:
        return default
    if normalised in TRUE_VALUES:
        return True
    if normalised in FALSE_VALUES:
        return False
    raise InputError("Expected a boolean value, found '{}'".format(value))


@dataclass(frozen=True)
class ProjectConfig:
    config_path: Path
    name: str
    output_dir: Path
    resume: bool
    stop_on_error: bool
    report: bool
    module_configs: Mapping[str, Path]
    inputs: Mapping[str, Path]
    species: str
    models: Tuple[str, ...]
    ssps: Tuple[str, ...]
    periods: Tuple[str, ...]
    seed: Optional[int]
    rscript: str
    compute_python: str
    integration_period: str
    integration_ssp: str
    ensemble_method: str
    # Optional vulnerability controls.  Empty values deliberately mean
    # "leave the module profile unchanged" so a project can inherit a
    # standalone vulnerability INI while still allowing one-command overrides.
    rona_variables: str = ""
    rona_summary: str = ""
    rona_weights: str = ""
    vulnerability_group_order: str = ""
    vulnerability_exclude_groups: str = ""

    def payload(self) -> Dict[str, object]:
        return {
            "project": {
                "name": self.name,
                "output_dir": str(self.output_dir),
                "resume": self.resume,
                "stop_on_error": self.stop_on_error,
                "report": self.report,
            },
            "modules": {key: str(path) for key, path in self.module_configs.items()},
            "inputs": {key: str(path) for key, path in self.inputs.items()},
            "shared": {
                "species": self.species or None,
                "models": list(self.models),
                "ssps": list(self.ssps),
                "periods": list(self.periods),
                "seed": self.seed,
                "rscript": self.rscript or None,
                "compute_python": self.compute_python or None,
            },
            "integration": {
                "period": self.integration_period or None,
                "ssp": self.integration_ssp or None,
                "ensemble_method": self.ensemble_method,
                "rona_variables": self.rona_variables or None,
                "rona_summary": self.rona_summary or None,
                "rona_weights": self.rona_weights or None,
                "group_order": self.vulnerability_group_order or None,
                "exclude_groups": self.vulnerability_exclude_groups or None,
            },
        }


def load_project_config(path: Path) -> ProjectConfig:
    config_path = path.expanduser().resolve()
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with config_path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as error:
        raise InputError("Cannot read project configuration: {}".format(error))
    for section in ("project", "modules"):
        if not parser.has_section(section):
            raise InputError("Project configuration requires [{}]".format(section))

    base = config_path.parent
    project = parser["project"]
    shared = parser["shared"] if parser.has_section("shared") else parser["DEFAULT"]
    integration = parser["integration"] if parser.has_section("integration") else parser["DEFAULT"]
    name = project.get("name", config_path.stem).strip()
    if not name:
        raise InputError("[project] name must not be empty")
    output_dir = _resolve_path(project.get("output_dir", ""), base, "[project] output_dir")

    module_configs: Dict[str, Path] = {}
    unknown = sorted(set(parser["modules"].keys()).difference(MODULE_ORDER))
    if unknown:
        raise InputError("Unknown project modules: {}".format(", ".join(unknown)))
    for key in MODULE_ORDER:
        raw = parser["modules"].get(key, "").strip()
        if raw.lower() in FALSE_VALUES:
            continue
        if raw.lower() in TRUE_VALUES:
            raise InputError("[modules] {} must be a configuration path, not '{}'".format(key, raw))
        module_path = _resolve_path(raw, base, "[modules] {}".format(key))
        if not module_path.is_file():
            raise InputError("Module configuration does not exist: {}".format(module_path))
        module_configs[key] = module_path
    if not module_configs:
        raise InputError("Enable at least one module with a configuration path in [modules]")

    project_inputs: Dict[str, Path] = {}
    if parser.has_section("inputs"):
        unknown_inputs = sorted(set(parser["inputs"].keys()).difference(PROJECT_INPUT_KEYS))
        if unknown_inputs:
            raise InputError("Unknown project inputs: {}".format(", ".join(unknown_inputs)))
        for key, raw in parser["inputs"].items():
            value = raw.strip()
            if not value or value.lower() in FALSE_VALUES:
                continue
            if value.lower() == "auto":
                if key != "adaptive_frequency":
                    raise InputError("[inputs] {} does not accept auto".format(key))
                continue
            project_inputs[key] = _resolve_path(value, base, "[inputs] {}".format(key))

    seed_value = shared.get("seed", "").strip()
    try:
        seed = int(seed_value) if seed_value else None
    except ValueError:
        raise InputError("[shared] seed must be an integer")

    models = _split(shared.get("models", ""))
    ssps = tuple(_normalise_ssp(value) for value in _split(shared.get("ssps", "")))
    periods = _split(shared.get("periods", ""))
    integration_period = integration.get("period", "").strip() or (periods[0] if len(periods) == 1 else "")
    integration_ssp = _normalise_ssp(integration.get("ssp", ""))
    if not integration_ssp and len(ssps) == 1:
        integration_ssp = ssps[0]
    ensemble_method = integration.get("ensemble_method", "mean").strip().lower()
    if ensemble_method not in {"mean", "median"}:
        raise InputError("[integration] ensemble_method must be mean or median")
    rona_variables = integration.get("rona_variables", "").strip()
    rona_summary = integration.get("rona_summary", "").strip().lower()
    if rona_summary == "weighted_mean":
        rona_summary = "weighted"
    if rona_summary and rona_summary not in {"mean", "max", "weighted"}:
        raise InputError("[integration] rona_summary must be mean, max, or weighted")
    rona_weights = integration.get("rona_weights", "").strip()
    vulnerability_group_order = integration.get("group_order", "").strip()
    vulnerability_exclude_groups = integration.get("exclude_groups", "").strip()

    return ProjectConfig(
        config_path=config_path,
        name=name,
        output_dir=output_dir,
        resume=_bool(project.get("resume", "true"), True),
        stop_on_error=_bool(project.get("stop_on_error", "true"), True),
        report=_bool(project.get("report", "true"), True),
        module_configs=module_configs,
        inputs=project_inputs,
        species=shared.get("species", "").strip(),
        models=models,
        ssps=ssps,
        periods=periods,
        seed=seed,
        rscript=shared.get("rscript", "").strip(),
        compute_python=shared.get("compute_python", "").strip(),
        integration_period=integration_period,
        integration_ssp=integration_ssp,
        ensemble_method=ensemble_method,
        rona_variables=rona_variables,
        rona_summary=rona_summary,
        rona_weights=rona_weights,
        vulnerability_group_order=vulnerability_group_order,
        vulnerability_exclude_groups=vulnerability_exclude_groups,
    )


def _read_ini(path: Path) -> configparser.ConfigParser:
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (OSError, configparser.Error, UnicodeError) as error:
        raise InputError("Cannot read module configuration {}: {}".format(path, error))
    return parser


def _read_project_table(path: Path, label: str) -> List[Dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            first = handle.readline()
            if not first:
                raise InputError("{} is empty: {}".format(label, path))
            delimiter = "\t" if "\t" in first else ","
            handle.seek(0)
            rows = list(csv.DictReader(handle, delimiter=delimiter))
    except (OSError, UnicodeError, csv.Error) as error:
        raise InputError("Cannot read {} {}: {}".format(label, path, error))
    if not rows:
        raise InputError("{} contains no data rows: {}".format(label, path))
    return [{str(key).strip(): (value or "").strip() for key, value in row.items()} for row in rows]


def _column(row: Mapping[str, str], *names: str) -> str:
    lower = {key.lower(): value for key, value in row.items()}
    for name in names:
        value = lower.get(name.lower(), "").strip()
        if value:
            return value
    return ""


def _write_project_table(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]], delimiter: str = "\t") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _materialize_project_inputs(config: ProjectConfig, project_dir: Path, selected: Sequence[str]) -> Dict[str, Path]:
    """Create small compatibility tables from the canonical project tables."""
    resolved = dict(config.inputs)
    generated = project_dir / "generated_inputs"
    samples_path = config.inputs.get("population_samples")
    environment_path = config.inputs.get("population_environment")
    sample_rows: Optional[List[Dict[str, str]]] = None
    environment_rows: Optional[List[Dict[str, str]]] = None

    if samples_path:
        sample_rows = _read_project_table(samples_path, "population_samples")
        for row in sample_rows:
            if not _column(row, "sample_id") or not _column(row, "population_id"):
                raise InputError("population_samples must contain sample_id and population_id")
        gf_samples = generated / "gf_samples.tsv"
        _write_project_table(
            gf_samples,
            ("sample_id", "population_id"),
            ((_column(row, "sample_id"), _column(row, "population_id")) for row in sample_rows),
        )
        resolved["gf_samples"] = gf_samples.resolve()

        if "load" in selected and "population_dir" not in resolved:
            population_dir = generated / "load_population_lists"
            population_dir.mkdir(parents=True, exist_ok=True)
            grouped: Dict[str, List[str]] = {}
            for row in sample_rows:
                grouped.setdefault(_column(row, "population_id"), []).append(_column(row, "sample_id"))
            for population, members in grouped.items():
                if not population or Path(population).name != population or population in {".", ".."}:
                    raise InputError("Unsafe population_id for generated load keep-list: {}".format(population))
                population_file = population_dir / population
                population_file.write_text("".join(member + "\n" for member in members), encoding="utf-8")
            resolved["generated_population_dir"] = population_dir.resolve()

    if environment_path:
        environment_rows = _read_project_table(environment_path, "population_environment")
        if "load" in selected and "load_predictors" not in resolved:
            predictor_rows = []
            for row in environment_rows:
                population = _column(row, "ID", "population_id", "pop")
                lon = _column(row, "lon", "longitude")
                lat = _column(row, "lat", "latitude")
                group = _column(row, "group", "cluster", "pop") or "all"
                bios = [_column(row, "bio{}".format(index)) for index in range(1, 20)]
                if not population or not lon or not lat or any(value == "" for value in bios):
                    raise InputError("population_environment needs ID/population_id, lon, lat and bio1-bio19 for load auto-generation")
                predictor_rows.append([population, lat, lon, group] + bios)
            predictors = generated / "load_predictors.csv"
            _write_project_table(predictors, ("pop", "lat", "lon", "cluster") + tuple("bio{}".format(i) for i in range(1, 20)), predictor_rows, delimiter=",")
            resolved["generated_load_predictors"] = predictors.resolve()

        if "niche" in selected and "occurrence_csv" not in resolved:
            species = config.species or config.name
            points = []
            seen = set()
            for row in environment_rows:
                lon = _column(row, "lon", "longitude")
                lat = _column(row, "lat", "latitude")
                if not lon or not lat:
                    raise InputError("population_environment needs lon and lat for occurrence auto-generation")
                key = (lon, lat)
                if key not in seen:
                    points.append([species, lon, lat])
                    seen.add(key)
            occurrence = generated / "occurrence.csv"
            _write_project_table(occurrence, ("species", "lon", "lat"), points, delimiter=",")
            resolved["generated_occurrence"] = occurrence.resolve()

    if "mar" in selected and "mar_lonlat" not in resolved and sample_rows:
        environment_rows = environment_rows or (_read_project_table(environment_path, "population_environment") if environment_path else [])
        coordinates: Dict[str, Tuple[str, str]] = {}
        for row in environment_rows:
            population = _column(row, "ID", "population_id")
            lon = _column(row, "lon", "longitude")
            lat = _column(row, "lat", "latitude")
            if population and lon and lat:
                coordinates[population] = (lon, lat)
        mar_rows = []
        for row in sample_rows:
            sample = _column(row, "sample_id")
            population = _column(row, "population_id")
            lon = _column(row, "lon", "longitude") or (coordinates.get(population, ("", ""))[0])
            lat = _column(row, "lat", "latitude") or (coordinates.get(population, ("", ""))[1])
            if not lon or not lat:
                raise InputError("MAR lonlat auto-generation needs lon/lat in population_samples or population_environment")
            mar_rows.append([sample, lon, lat])
        mar_lonlat = generated / "mar_lonlat.tsv"
        _write_project_table(mar_lonlat, ("ID", "LONGITUDE", "LATITUDE"), mar_rows)
        resolved["generated_mar_lonlat"] = mar_lonlat.resolve()
    return resolved


def _validate_project_input_paths(config: ProjectConfig) -> None:
    directory_keys = {
        "load_vcf_dir",
        "population_dir",
        "rona_unld_dir",
        "climate_current",
        "climate_future",
        "load_future_dir",
    }
    for key, path in config.inputs.items():
        expected = "directory" if key in directory_keys else "file"
        if expected == "directory" and not path.is_dir():
            raise InputError("Project input {} must be an existing directory: {}".format(key, path))
        if expected == "file" and not path.is_file():
            raise InputError("Project input {} must be an existing file: {}".format(key, path))


def _input_contract_rows(
    config: ProjectConfig,
    resolved_inputs: Mapping[str, Path],
    fingerprints: Optional[Mapping[str, Mapping[str, object]]] = None,
) -> List[List[str]]:
    descriptions = {
        "adaptive_vcf": "Adaptive loci VCF for GF/RONA",
        "adaptive_frequency": "Adaptive loci population frequency table for RONA",
        "whole_genome_vcf": "Beagle-filtered whole-genome VCF for MAR",
        "load_vcf_dir": "Derived + SIFT annotated load VCF directory",
        "population_samples": "Canonical sample/population table",
        "population_environment": "Canonical population environment table",
        "rona_unld_dir": "RONA LD-pruning files",
        "climate_current": "Shared current BIO rasters",
        "climate_future": "Shared future BIO rasters",
        "species_mask": "Shared species mask",
        "wfmoment_current_raster": "Current binary habitat raster",
        "pi_file": "Observed pi table",
        "structure_file": "Population structure/FST table",
        "species_parameters": "Species parameter table",
        "area_file": "Area scenario table",
        "future_masks_json": "Future binary masks",
        "gf_samples": "Generated GF sample/population table",
        "generated_population_dir": "Generated load keep-list directory",
        "generated_load_predictors": "Generated load predictor table",
        "generated_occurrence": "Generated MaxEnt occurrence table",
        "generated_mar_lonlat": "Generated MAR sample coordinates",
    }
    fingerprints = fingerprints or {}
    rows: List[List[str]] = []
    for key, path in config.inputs.items():
        record = fingerprints.get(key, {})
        rows.append(
            [
                key,
                str(path),
                descriptions.get(key, "Project input"),
                "user",
                str(record.get("kind", "")),
                str(record.get("size_bytes", record.get("total_size_bytes", ""))),
                str(record.get("sha256", "")),
                str(record.get("entry_count", record.get("component_count", ""))),
            ]
        )
    for key, path in resolved_inputs.items():
        if key not in config.inputs:
            record = fingerprints.get(key, {})
            rows.append(
                [
                    key,
                    str(path),
                    descriptions.get(key, "Generated compatibility input"),
                    "generated",
                    str(record.get("kind", "")),
                    str(record.get("size_bytes", record.get("total_size_bytes", ""))),
                    str(record.get("sha256", "")),
                    str(record.get("entry_count", record.get("component_count", ""))),
                ]
            )
    # Standalone module configurations and their referenced external files are
    # included as provenance rows as well.  They are not canonical project
    # inputs, so label them separately instead of hiding them from readers.
    seen = set(config.inputs).union(resolved_inputs)
    for key in sorted(set(fingerprints).difference(seen)):
        record = fingerprints.get(key, {})
        path = record.get("path", "")
        if not path:
            continue
        origin = "configuration" if key in {"project.configuration"} or key.endswith(".configuration") else "module-config"
        rows.append(
            [
                key,
                str(path),
                "Configuration or module-referenced input",
                origin,
                str(record.get("kind", "")),
                str(record.get("size_bytes", record.get("total_size_bytes", ""))),
                str(record.get("sha256", "")),
                str(record.get("entry_count", record.get("component_count", ""))),
            ]
        )
    return rows


def _ensure_section(parser: configparser.ConfigParser, name: str) -> None:
    if not parser.has_section(name):
        parser.add_section(name)


def _normalise_input_paths(parser: configparser.ConfigParser, source: Path) -> None:
    if not parser.has_section("inputs"):
        return
    for key, raw in list(parser.items("inputs")):
        value = raw.strip()
        if not value or value.lower() == "auto":
            continue
        path = Path(value).expanduser()
        if not path.is_absolute():
            parser.set("inputs", key, str((source.parent / path).resolve()))


def _set_if(value: str, parser: configparser.ConfigParser, section: str, key: str) -> None:
    if not value:
        return
    _ensure_section(parser, section)
    parser.set(section, key, value)


def _apply_shared(parser: configparser.ConfigParser, module: str, config: ProjectConfig) -> None:
    if config.species and module == "wfmoment":
        _set_if(config.species, parser, "analysis", "species")
    if config.models:
        value = ",".join(config.models)
        if module in {"gf", "rona"}:
            _set_if(value, parser, "analysis", "models")
        elif module == "niche":
            _set_if(value, parser, "projection", "gcms")
    if config.ssps:
        plain = ",".join(config.ssps)
        if module in {"gf", "rona"}:
            _set_if(plain, parser, "analysis", "ssps")
        elif module == "niche":
            _set_if(",".join("ssp{}".format(value) for value in config.ssps), parser, "projection", "scenarios")
    if config.periods:
        value = ",".join(config.periods)
        if module in {"gf", "rona"}:
            _set_if(value, parser, "analysis", "periods")
        elif module == "niche":
            _set_if(value, parser, "projection", "periods")
    if config.seed is not None:
        section_key = {
            "gf": ("parameters", "seed"),
            "mar": ("parameters", "randseed"),
            "load": ("parameters", "seed"),
            "niche": ("tuning", "seed"),
            "wfmoment": ("parameters", "seed"),
        }.get(module)
        if section_key:
            _set_if(str(config.seed), parser, section_key[0], section_key[1])
    if config.rscript:
        section = "analysis" if module == "niche" else "parameters"
        _set_if(config.rscript, parser, section, "rscript")
        if module == "niche":
            _set_if(config.rscript, parser, "analysis", "plot_rscript")
    if config.compute_python:
        if module == "wfmoment":
            _set_if(config.compute_python, parser, "parameters", "compute_python")
        elif module == "load":
            _set_if(config.compute_python, parser, "parameters", "calc_python")
    if module == "vulnerability":
        # These are intentionally opt-in project overrides.  A blank project
        # value must not erase an explicit standalone module profile, while
        # values such as ``auto`` or ``none`` provide an explicit reset.
        _set_if(config.rona_variables, parser, "parameters", "rona_variables")
        _set_if(config.rona_summary, parser, "parameters", "rona_summary")
        _set_if(config.rona_weights, parser, "parameters", "rona_weights")
        _set_if(config.vulnerability_group_order, parser, "parameters", "group_order")
        _set_if(config.vulnerability_exclude_groups, parser, "parameters", "exclude_groups")


def _auto_vulnerability_inputs(parser: configparser.ConfigParser, config: ProjectConfig) -> List[str]:
    if not parser.has_section("inputs"):
        return []
    auto_keys = [key for key, value in parser.items("inputs") if value.strip().lower() == "auto"]
    if not auto_keys:
        return []
    if not config.integration_period or not config.integration_ssp:
        raise InputError(
            "Automatic vulnerability inputs require [integration] period and ssp "
            "(or exactly one shared period and SSP)"
        )
    period = config.integration_period
    ssp = config.integration_ssp
    scenario = "ssp{}".format(ssp)
    paths = {
        "gf_offsets": config.output_dir / "01_gf" / "05_offset_plots" / period / scenario / "ensemble_mean_offsets.tsv.gz",
        "rona_ensemble": config.output_dir / "02_rona" / "ensemble_mean" / "{}_{}_RONA_ensemble_mean.tsv".format(period, scenario),
        "load_predictors": config.output_dir / "04_load" / "pop_geo_niche_predictors_from_TSS.csv",
        "niche_raster": config.output_dir / "05_niche" / "rasters" / "maladaptation_{}_{}_ensemble_{}.tif".format(period, scenario, config.ensemble_method),
    }
    for key in auto_keys:
        if key not in paths:
            raise InputError("No automatic vulnerability mapping is defined for [inputs] {}".format(key))
        dependency = VULNERABILITY_DEPENDENCIES[key]
        if dependency not in config.module_configs and not paths[key].is_file():
            raise InputError("Automatic vulnerability input {} requires the {} module".format(key, dependency))
        parser.set("inputs", key, str(paths[key]))
    _ensure_section(parser, "analysis")
    if not parser.get("analysis", "scenario_label", fallback="").strip():
        parser.set("analysis", "scenario_label", "{} SSP{}".format(period, ssp))
    return auto_keys


def _apply_project_inputs(
    parser: configparser.ConfigParser,
    module: str,
    config: ProjectConfig,
    resolved_inputs: Mapping[str, Path],
    selected: Sequence[str],
) -> None:
    def set_target(section: str, key: str, input_key: str) -> None:
        value = resolved_inputs.get(input_key)
        if value is not None:
            _ensure_section(parser, section)
            parser.set(section, key, str(value))

    for input_key, targets in PROJECT_INPUT_TARGETS.items():
        for target_module, section, key in targets:
            if target_module == module:
                set_target(section, key, input_key)

    if module == "gf":
        set_target("inputs", "samples", "gf_samples")
    elif module == "load":
        set_target("inputs", "population_dir", "generated_population_dir")
        set_target("inputs", "predictors", "generated_load_predictors")
    elif module == "mar":
        set_target("inputs", "lonlat", "generated_mar_lonlat")
    elif module == "niche":
        set_target("inputs", "occurrence_csv", "generated_occurrence")

    if module == "rona" and parser.has_section("inputs"):
        raw_frequency = parser.get("inputs", "alt_frequency", fallback="").strip().lower()
        explicit_frequency = resolved_inputs.get("adaptive_frequency")
        if explicit_frequency is not None:
            parser.set("inputs", "alt_frequency", str(explicit_frequency))
        elif raw_frequency == "auto" or ("adaptive_vcf" in resolved_inputs and "gf" in selected):
            if "gf" not in selected:
                raise InputError("RONA adaptive frequency auto-wiring requires the GF module")
            parser.set("inputs", "alt_frequency", str(config.output_dir / MODULE_DIRECTORIES["gf"] / "01_frequency" / "population_alt_frequency.tsv"))


def _effective_config(
    module: str,
    source: Path,
    config: ProjectConfig,
    destination: Path,
    resolved_inputs: Mapping[str, Path],
    selected: Sequence[str],
) -> Tuple[Path, List[str]]:
    parser = _read_ini(source)
    _normalise_input_paths(parser, source)
    _ensure_section(parser, "analysis")
    parser.set("analysis", "output_dir", str(config.output_dir / MODULE_DIRECTORIES[module]))
    _apply_shared(parser, module, config)
    _apply_project_inputs(parser, module, config, resolved_inputs, selected)

    # User-supplied script overrides remain relative to the source module
    # configuration, even though the effective INI is stored under results.
    if module == "load" and parser.has_section("parameters"):
        for key in ("calc_script", "rf_script", "predict_script"):
            value = parser.get("parameters", key, fallback="").strip()
            if value:
                path = Path(value).expanduser()
                if not path.is_absolute():
                    parser.set("parameters", key, str((source.parent / path).resolve()))

    auto_inputs = _auto_vulnerability_inputs(parser, config) if module == "vulnerability" else []
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        parser.write(handle)
    return destination, auto_inputs


def _preflight_config(source: Path, output_dir: Path, destination: Path) -> Path:
    parser = _read_ini(source)
    _ensure_section(parser, "analysis")
    parser.set("analysis", "output_dir", str(output_dir))
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        parser.write(handle)
    return destination


def _runner(module: str) -> Callable[..., Dict[str, object]]:
    if module == "gf":
        from .gf_workflow import run_gf_workflow
        return run_gf_workflow
    if module == "rona":
        from .rona import run_rona_workflow
        return run_rona_workflow
    if module == "mar":
        from .mar import run_mar_workflow
        return run_mar_workflow
    if module == "load":
        from .load import run_load_workflow
        return run_load_workflow
    if module == "niche":
        from .niche import run_niche_workflow
        return run_niche_workflow
    if module == "wfmoment":
        from .wfmoment import run_wfmoment_workflow
        return run_wfmoment_workflow
    if module == "vulnerability":
        from .vulnerability import run_vulnerability_workflow
        return run_vulnerability_workflow
    raise InputError("Unknown BioMA module: {}".format(module))


def _manifest_complete(manifest: Mapping[str, object]) -> bool:
    return manifest.get("status") == "complete" or manifest.get("engineering_status") == "completed"


def _manifest_config_hash(manifest: Mapping[str, object]) -> Optional[str]:
    direct = manifest.get("config_sha256") or manifest.get("manifest_sha256")
    if isinstance(direct, str) and direct:
        return direct
    inputs = manifest.get("inputs")
    if isinstance(inputs, Mapping):
        value = inputs.get("configuration_sha256")
        if isinstance(value, str) and value:
            return value
    return None


def _existing_module(output_dir: Path, effective_config: Path) -> Tuple[str, Optional[Dict[str, object]]]:
    if not output_dir.exists():
        return "empty", None
    if not output_dir.is_dir():
        return "invalid", None
    if not any(output_dir.iterdir()):
        return "empty", None
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.is_file():
        return "incomplete", None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as error:
        raise InputError("Cannot read existing module manifest {}: {}".format(manifest_path, error))
    if not _manifest_complete(manifest):
        return "incomplete", manifest
    recorded_hash = _manifest_config_hash(manifest)
    if recorded_hash and recorded_hash != _sha256_file(effective_config):
        return "changed", manifest
    if not recorded_hash:
        return "unverified", manifest
    return "complete", manifest


def _archive_module(output_root: Path, module_dir: Path) -> Optional[Path]:
    if not module_dir.exists():
        return None
    root = output_root.resolve()
    target = module_dir.resolve()
    allowed_names = set(MODULE_DIRECTORIES.values()).union({"00_project"})
    if target.parent != root or target.name not in allowed_names:
        raise InputError("Refusing to archive unexpected path: {}".format(target))
    stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
    backup_root = root / "_bioma_backups" / stamp
    backup_root.mkdir(parents=True, exist_ok=True)
    destination = backup_root / target.name
    counter = 1
    while destination.exists():
        counter += 1
        destination = backup_root / "{}-{}".format(target.name, counter)
    shutil.move(str(target), str(destination))
    return destination


def _selection(config: ProjectConfig, only_modules: Sequence[str]) -> List[str]:
    requested = list(dict.fromkeys(only_modules))
    unknown = sorted(set(requested).difference(MODULE_ORDER))
    if unknown:
        raise InputError("Unknown --module value(s): {}".format(", ".join(unknown)))
    selected = [key for key in MODULE_ORDER if key in config.module_configs]
    if requested:
        missing = [key for key in requested if key not in config.module_configs]
        if missing:
            raise InputError("Requested modules are disabled in project.ini: {}".format(", ".join(missing)))
        selected = [key for key in selected if key in requested]
    return selected


def _write_summary(path: Path, records: Iterable[Mapping[str, object]]) -> None:
    fields = ("order", "module", "status", "output_dir", "elapsed_seconds", "message")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for record in records:
            writer.writerow(record)


def _write_report(path: Path, manifest: Mapping[str, object]) -> None:
    records = manifest.get("modules", [])
    if not records and manifest.get("status") == "planned":
        preflight = manifest.get("preflight", {})
        selection = manifest.get("selection", [])
        project = manifest.get("project", {})
        output_root = Path(str(project.get("output_dir", ""))) if isinstance(project, Mapping) else Path()
        if isinstance(preflight, Mapping) and isinstance(selection, list):
            planned_records = []
            for index, module in enumerate(selection, start=1):
                detail = preflight.get(module, {})
                if not isinstance(detail, Mapping):
                    detail = {}
                waiting_for = detail.get("waiting_for", [])
                message = "Waiting for: {}".format(", ".join(str(value) for value in waiting_for)) if waiting_for else ""
                planned_records.append(
                    {
                        "order": index,
                        "module": module,
                        "status": detail.get("status", "unknown"),
                        "output_dir": str(output_root / MODULE_DIRECTORIES.get(module, module)),
                        "message": message,
                    }
                )
            records = planned_records
    rows = []
    if isinstance(records, list):
        for record in records:
            if not isinstance(record, Mapping):
                continue
            rows.append(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td><code>{}</code></td><td>{}</td></tr>".format(
                    html.escape(str(record.get("order", ""))),
                    html.escape(str(record.get("module", ""))),
                    '<span class="{}">{}</span>'.format(
                        html.escape(str(record.get("status", ""))),
                        html.escape(str(record.get("status", ""))),
                    ),
                    html.escape(str(record.get("output_dir", ""))),
                    html.escape(str(record.get("message", ""))),
                )
            )
    shared = manifest.get("project", {})
    document = """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>{title}</title><style>
body{{font-family:system-ui,sans-serif;margin:2rem auto;max-width:1100px;padding:0 1rem;color:#202124}}
h1{{font-size:1.8rem}} table{{width:100%;border-collapse:collapse;margin-top:1.5rem}}
th,td{{border-bottom:1px solid #d8dadd;padding:.65rem;text-align:left;vertical-align:top}}
th{{background:#f3f4f5}} code{{font-size:.85rem;overflow-wrap:anywhere}}
.meta{{color:#555}} .completed,.reused,.validated{{color:#176b35;font-weight:600}}
.deferred{{color:#8a5a00;font-weight:600}}
.failed,.blocked{{color:#a12622;font-weight:600}}
</style></head><body><h1>{title}</h1>
<p class="meta">BioMA {version} | status: <strong>{status}</strong> | generated: {generated}</p>
<p>Project output: <code>{output}</code></p>
<table><thead><tr><th>#</th><th>Module</th><th>Status</th><th>Output</th><th>Message</th></tr></thead>
<tbody>{rows}</tbody></table>
<p class="meta">Machine-readable provenance is stored in <code>00_project/run_manifest.json</code> and
<code>00_project/project_summary.tsv</code>.</p></body></html>""".format(
        title=html.escape(str(shared.get("name", "BioMA project report")) if isinstance(shared, Mapping) else "BioMA project report"),
        version=html.escape(__version__),
        status=html.escape(str(manifest.get("status", "unknown"))),
        generated=html.escape(str(manifest.get("completed_at_utc") or manifest.get("updated_at_utc") or "")),
        output=html.escape(str(shared.get("output_dir", "")) if isinstance(shared, Mapping) else ""),
        rows="".join(rows),
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(document, encoding="utf-8")


def run_project_workflow(
    config_path: Path,
    dry_run: bool = False,
    only_modules: Sequence[str] = (),
    resume: Optional[bool] = None,
    overwrite: bool = False,
    progress: Optional[Callable[[str], None]] = None,
) -> Dict[str, object]:
    config = load_project_config(config_path)
    selected = _selection(config, only_modules)
    resume_enabled = config.resume if resume is None else resume
    announce = progress or (lambda message: None)
    output_root = config.output_dir
    if output_root.exists() and not output_root.is_dir():
        raise InputError("Project output path is not a directory: {}".format(output_root))
    output_root.mkdir(parents=True, exist_ok=True)
    project_dir = output_root / "00_project"
    if project_dir.exists() and not project_dir.is_dir():
        if not overwrite:
            raise InputError("Project metadata path is not a directory: {}".format(project_dir))
        archived_project = _archive_module(output_root, project_dir)
        announce("Archived invalid project metadata to {}".format(archived_project))
    effective_dir = project_dir / "effective_configs"
    source_dir = project_dir / "source_configs"
    _validate_project_input_paths(config)
    sources: Dict[str, Dict[str, str]] = {}

    # Resolve and hash all declared inputs before creating any generated
    # compatibility files.  This ordering is important: a changed input must
    # fail fast without mutating the previous project's metadata or generated
    # tables.  Directory/file content is represented by a stable fingerprint,
    # rather than only by its path.
    for module in MODULE_ORDER:
        if module not in config.module_configs:
            continue
        source = config.module_configs[module]
        source_hash = _sha256_file(source)
        sources[module] = {"path": str(source), "sha256": source_hash}

    configured_paths: Dict[str, Path] = dict(config.inputs)
    # Configuration files are user-controlled inputs too.  Keeping them in the
    # same manifest makes the provenance sidecar self-contained; their hashes
    # are also retained in ``sources`` for backwards compatibility.
    configured_paths["project.configuration"] = config.config_path
    configured_paths.update(
        {
            "module.{}.configuration".format(module): source
            for module, source in config.module_configs.items()
        }
    )
    configured_paths.update(_configured_input_paths(config))
    input_fingerprints = _fingerprint_paths(configured_paths)
    compact_input_fingerprints = {
        key: _compact_fingerprint(value) for key, value in input_fingerprints.items()
    }
    signature_payload = {
        "input_hash_schema": 1,
        "project_config_sha256": _sha256_file(config.config_path),
        "configured_modules": [key for key in MODULE_ORDER if key in config.module_configs],
        "selected_modules": list(selected),
        "source_configs": sources,
        "input_fingerprints": compact_input_fingerprints,
        "bioma_version": __version__,
    }
    input_signature = _json_sha256(signature_payload)
    resolved_path = project_dir / "resolved_project.json"
    if resolved_path.is_file():
        try:
            existing = json.loads(resolved_path.read_text(encoding="utf-8"))
        except (OSError, ValueError) as error:
            raise InputError("Cannot read existing project metadata {}: {}".format(resolved_path, error))
        recorded_input_signature = existing.get("input_signature_sha256") or existing.get("project_signature_sha256")
        if recorded_input_signature != input_signature and not overwrite:
            raise InputError(
                "Project configuration changed for {}. Use a new output_dir or explicit --overwrite; "
                "old module directories will be backed up when overwrite is used.".format(output_root)
            )
        if not overwrite:
            _check_existing_input_manifest(existing, project_dir)
            _check_existing_generated_inputs(existing, project_dir)
    if overwrite and project_dir.is_dir() and any(project_dir.iterdir()):
        archived_project = _archive_module(output_root, project_dir)
        announce("Archived project metadata to {}".format(archived_project))

    # The content check above has passed; generated compatibility files and
    # effective configurations can now be written safely.
    resolved_inputs = _materialize_project_inputs(config, project_dir, selected)
    resolved_input_fingerprints: Dict[str, Dict[str, object]] = {}
    for key, path in resolved_inputs.items():
        if key in input_fingerprints:
            resolved_input_fingerprints[key] = input_fingerprints[key]
        else:
            resolved_input_fingerprints[key] = _fingerprint_path(path)

    # The final project signature includes generated compatibility files as
    # well.  The pre-materialization ``input_signature`` above remains the
    # immutable gate used to protect an existing project from accidental
    # writes; this second digest is the complete post-resolution identity.
    compact_resolved_fingerprints = {
        key: _compact_fingerprint(value) for key, value in resolved_input_fingerprints.items()
    }
    final_signature_payload = {
        "input_signature_sha256": input_signature,
        "resolved_input_fingerprints": compact_resolved_fingerprints,
        "bioma_version": __version__,
    }
    signature = _json_sha256(final_signature_payload)

    effective_configs: Dict[str, Path] = {}
    auto_inputs: Dict[str, List[str]] = {}

    for module in MODULE_ORDER:
        if module not in config.module_configs:
            continue
        source = config.module_configs[module]
        source_hash = sources[module]["sha256"]
        source_copy = source_dir / "{}.{}.ini".format(module, source_hash[:12])
        source_copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(source_copy))
        effective, auto = _effective_config(
            module,
            source,
            config,
            effective_dir / "{}.ini".format(module),
            resolved_inputs,
            selected,
        )
        effective_configs[module] = effective
        auto_inputs[module] = auto
        sources[module] = {
            "path": str(source),
            "sha256": source_hash,
            "snapshot": str(source_copy),
            "effective_config": str(effective),
            "effective_sha256": _sha256_file(effective),
        }

    # Keep the detailed per-file records in a separate sidecar.  The compact
    # summaries in resolved_project.json and input_contract.tsv remain small
    # even when a climate directory contains many large rasters.
    manifest_fingerprints: Dict[str, Dict[str, object]] = {}
    for key, value in input_fingerprints.items():
        manifest_fingerprints["project.{}".format(key)] = value
    for key, value in resolved_input_fingerprints.items():
        manifest_fingerprints["resolved.{}".format(key)] = value
    input_manifest_path = project_dir / "input_manifest.json"
    _write_input_manifest(input_manifest_path, manifest_fingerprints)
    input_manifest_sha256 = _sha256_file(input_manifest_path)
    _write_json_atomic(
        resolved_path,
        {
            "project_signature_sha256": signature,
            "input_signature_sha256": input_signature,
            "signature_payload": {
                "input_hash_schema": signature_payload["input_hash_schema"],
                "project_config_sha256": signature_payload["project_config_sha256"],
                "configured_modules": signature_payload["configured_modules"],
                "selected_modules": signature_payload["selected_modules"],
                "input_fingerprints": compact_input_fingerprints,
                "resolved_input_fingerprints": compact_resolved_fingerprints,
            },
            "config": config.payload(),
            "resolved_inputs": {key: str(value) for key, value in resolved_inputs.items()},
            "input_fingerprints": compact_input_fingerprints,
            "resolved_input_fingerprints": {
                key: value for key, value in compact_resolved_fingerprints.items()
            },
            "input_manifest": str(input_manifest_path),
            "input_manifest_sha256": input_manifest_sha256,
            "sources": sources,
            "updated_at_utc": _utc_now(),
        },
    )
    _write_project_table(
        project_dir / "input_contract.tsv",
        ("key", "path", "description", "origin", "kind", "size_bytes", "sha256", "entries_or_components"),
        _input_contract_rows(
            config,
            resolved_inputs,
            {**input_fingerprints, **resolved_input_fingerprints},
        ),
    )

    preflight: Dict[str, Dict[str, object]] = {}
    announce("Preflight {} module(s)".format(len(selected)))
    for module in selected:
        if module == "rona":
            effective_rona = _read_ini(effective_configs[module])
            auto_frequency = Path(effective_rona.get("inputs", "alt_frequency", fallback="")).expanduser()
            expected_frequency = output_root / MODULE_DIRECTORIES["gf"] / "01_frequency" / "population_alt_frequency.tsv"
            if auto_frequency == expected_frequency and not auto_frequency.is_file():
                if "gf" not in selected:
                    raise InputError("Project preflight failed for RONA: automatic adaptive frequency requires selected module gf")
                preflight[module] = {"status": "deferred", "waiting_for": ["gf"], "auto_input": "alt_frequency"}
                continue
        if module == "vulnerability" and auto_inputs[module]:
            parser = _read_ini(effective_configs[module])
            missing_auto: List[Tuple[str, Path]] = []
            for key, value in parser.items("inputs"):
                input_path = Path(value).expanduser()
                if key in auto_inputs[module]:
                    if not input_path.is_file():
                        missing_auto.append((key, input_path))
                elif not input_path.is_file():
                    raise InputError(
                        "Project preflight failed for vulnerability: explicit input {} does not exist: {}".format(
                            key, input_path
                        )
                    )
            unavailable = [
                (key, VULNERABILITY_DEPENDENCIES[key], path)
                for key, path in missing_auto
                if VULNERABILITY_DEPENDENCIES[key] not in selected
            ]
            if unavailable:
                details = "; ".join(
                    "{} requires selected module {} or existing file {}".format(key, dependency, path)
                    for key, dependency, path in unavailable
                )
                raise InputError("Project preflight failed for vulnerability: {}".format(details))
            if missing_auto:
                preflight[module] = {
                    "status": "deferred",
                    "auto_inputs": auto_inputs[module],
                    "waiting_for": sorted({VULNERABILITY_DEPENDENCIES[key] for key, _ in missing_auto}),
                }
                continue
        try:
            preflight_path = _preflight_config(
                effective_configs[module],
                project_dir / "preflight" / module,
                project_dir / "preflight_configs" / "{}.ini".format(module),
            )
            planned = _runner(module)(preflight_path, dry_run=True, progress=None)
            preflight[module] = {"status": "validated", "plan": planned}
        except Exception as error:
            preflight[module] = {"status": "failed", "error": str(error)}
            raise InputError("Project preflight failed for {}: {}".format(module, error))

    base_manifest: Dict[str, object] = {
        "module": "bioma-project",
        "bioma_version": __version__,
        "status": "planned" if dry_run else "running",
        "updated_at_utc": _utc_now(),
        "project_signature_sha256": signature,
        "input_signature_sha256": input_signature,
        "input_fingerprints": compact_input_fingerprints,
        "input_manifest": str(input_manifest_path),
        "input_manifest_sha256": input_manifest_sha256,
        "project": config.payload()["project"],
        "selection": selected,
        "shared": config.payload()["shared"],
        "integration": config.payload()["integration"],
        "sources": sources,
        "preflight": preflight,
        "runtime": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
        },
        "modules": [],
    }
    state_path = project_dir / "project_state.json"
    _write_json_atomic(state_path, base_manifest)
    if dry_run:
        _write_json_atomic(project_dir / "project_plan.json", base_manifest)
        _write_summary(project_dir / "project_summary.tsv", [])
        if config.report:
            _write_report(output_root / "report" / "index.html", base_manifest)
        return base_manifest

    started = time.monotonic()
    records: List[Dict[str, object]] = []
    failures = 0
    for index, module in enumerate(selected, start=1):
        module_dir = output_root / MODULE_DIRECTORIES[module]
        if module == "vulnerability" and auto_inputs[module]:
            failed_modules = {
                str(record["module"])
                for record in records
                if record.get("status") in {"failed", "blocked"}
            }
            blocked_by = sorted(
                {
                    VULNERABILITY_DEPENDENCIES[key]
                    for key in auto_inputs[module]
                    if VULNERABILITY_DEPENDENCIES[key] in failed_modules
                }
            )
            if blocked_by:
                record = {
                    "order": index,
                    "module": module,
                    "status": "blocked",
                    "output_dir": str(module_dir),
                    "elapsed_seconds": 0.0,
                    "message": "Required upstream module(s) failed: {}".format(", ".join(blocked_by)),
                }
                records.append(record)
                base_manifest["modules"] = records
                base_manifest["status"] = "failed"
                base_manifest["updated_at_utc"] = _utc_now()
                _write_json_atomic(state_path, base_manifest)
                announce("Block {}: {}".format(module, record["message"]))
                continue
        status, existing_manifest = _existing_module(module_dir, effective_configs[module])
        if overwrite and status != "empty":
            archived = _archive_module(output_root, module_dir)
            announce("Archived {} to {}".format(module, archived))
            status, existing_manifest = "empty", None
        if resume_enabled and status == "complete" and existing_manifest is not None:
            record = {
                "order": index,
                "module": module,
                "status": "reused",
                "output_dir": str(module_dir),
                "elapsed_seconds": 0.0,
                "message": "Completed result and configuration hash matched",
                "manifest": str(module_dir / "run_manifest.json"),
            }
            records.append(record)
            announce("Reuse {}".format(module))
            base_manifest["modules"] = records
            base_manifest["updated_at_utc"] = _utc_now()
            _write_json_atomic(state_path, base_manifest)
            continue
        if status != "empty":
            raise InputError(
                "Existing {} result is {} and cannot be safely reused: {}. "
                "Use --overwrite to move it into _bioma_backups before rerunning.".format(module, status, module_dir)
            )

        announce("Run {} ({}/{})".format(module, index, len(selected)))
        module_started = time.monotonic()
        try:
            result = _runner(module)(
                effective_configs[module],
                dry_run=False,
                progress=lambda message, module=module: announce("{}: {}".format(module, message)),
            )
            record = {
                "order": index,
                "module": module,
                "status": "completed",
                "output_dir": str(module_dir),
                "elapsed_seconds": round(time.monotonic() - module_started, 3),
                "message": "",
                "manifest": str(module_dir / "run_manifest.json"),
                "result": result,
            }
        except Exception as error:
            failures += 1
            record = {
                "order": index,
                "module": module,
                "status": "failed",
                "output_dir": str(module_dir),
                "elapsed_seconds": round(time.monotonic() - module_started, 3),
                "message": str(error),
            }
            records.append(record)
            base_manifest["modules"] = records
            base_manifest["status"] = "failed"
            base_manifest["updated_at_utc"] = _utc_now()
            _write_json_atomic(state_path, base_manifest)
            if config.stop_on_error:
                _write_summary(project_dir / "project_summary.tsv", records)
                if config.report:
                    _write_report(output_root / "report" / "index.html", base_manifest)
                raise
            continue
        records.append(record)
        base_manifest["modules"] = records
        base_manifest["updated_at_utc"] = _utc_now()
        _write_json_atomic(state_path, base_manifest)

    base_manifest.update(
        {
            "status": "complete" if failures == 0 else "complete_with_failures",
            "completed_at_utc": _utc_now(),
            "elapsed_seconds": round(time.monotonic() - started, 3),
            "counts": {
                "selected_modules": len(selected),
                "completed": sum(record["status"] == "completed" for record in records),
                "reused": sum(record["status"] == "reused" for record in records),
                "failed": failures,
                "blocked": sum(record["status"] == "blocked" for record in records),
            },
            "modules": records,
        }
    )
    manifest_path = project_dir / "run_manifest.json"
    _write_json_atomic(manifest_path, base_manifest)
    _write_json_atomic(state_path, base_manifest)
    _write_summary(project_dir / "project_summary.tsv", records)
    if config.report:
        _write_report(output_root / "report" / "index.html", base_manifest)
    announce("Complete {}".format(output_root))
    return base_manifest
