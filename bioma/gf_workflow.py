"""Run the complete Gradient Forest workflow from one INI configuration."""

from __future__ import annotations

import configparser
import hashlib
import json
import os
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence, Tuple

from . import __version__
from .climate_prepare import (
    Scenario,
    discover_bioclim_rasters,
    discover_scenarios,
    prepare_climate_inputs,
)
from .forward_distance_plot import plot_forward_distance
from .gf_frequency import InputError, build_population_frequency
from .gf_offset import compute_gf_offsets, parse_forward_radii
from .gf_train import DEFAULT_PREDICTORS, train_gradient_forest
from .offset_plot import plot_offsets


@dataclass(frozen=True)
class WorkflowConfig:
    config_path: Path
    vcf: Path
    samples: Optional[Path]
    sample_groups_dir: Optional[Path]
    coordinates: Path
    present_climate: Path
    future_climate: Path
    current_mask: Path
    future_mask: Path
    output_dir: Path
    models: Tuple[str, ...]
    ssps: Tuple[str, ...]
    periods: Tuple[str, ...]
    predictors: Tuple[str, ...]
    expected_sites: Optional[int]
    min_population_samples: int
    warn_population_samples: int
    ntree: int
    nbin: int
    corr_threshold: float
    max_level: Optional[float]
    seed: int
    forward_radii_km: Tuple[float, ...]
    initial_knn_k: int
    batch_size: int
    verify_sample: int
    tie_tolerance: float
    map_radius: str
    minimum_models: int
    supplied_bio_tolerance: float
    rscript: Optional[str]
    gdalinfo: str
    gdallocationinfo: str
    ogrinfo: str

    def payload(self) -> Dict[str, object]:
        return {
            "inputs": {
                "vcf": str(self.vcf),
                "samples": str(self.samples) if self.samples else None,
                "sample_groups_dir": str(self.sample_groups_dir) if self.sample_groups_dir else None,
                "coordinates": str(self.coordinates),
                "present_climate": str(self.present_climate),
                "future_climate": str(self.future_climate),
                "current_mask": str(self.current_mask),
                "future_mask": str(self.future_mask),
            },
            "analysis": {
                "output_dir": str(self.output_dir),
                "models": list(self.models),
                "ssps": list(self.ssps),
                "periods": list(self.periods),
                "expected_sites": self.expected_sites,
            },
            "parameters": {
                "predictors": list(self.predictors),
                "min_population_samples": self.min_population_samples,
                "warn_population_samples": self.warn_population_samples,
                "ntree": self.ntree,
                "nbin": self.nbin,
                "corr_threshold": self.corr_threshold,
                "max_level": self.max_level,
                "seed": self.seed,
                "forward_radii_km": [
                    "unlimited" if value == float("inf") else value
                    for value in self.forward_radii_km
                ],
                "initial_knn_k": self.initial_knn_k,
                "batch_size": self.batch_size,
                "verify_sample": self.verify_sample,
                "tie_tolerance": self.tie_tolerance,
                "map_radius": self.map_radius,
                "minimum_models": self.minimum_models,
                "supplied_bio_tolerance": self.supplied_bio_tolerance,
                "rscript": self.rscript,
                "gdalinfo": self.gdalinfo,
                "gdallocationinfo": self.gdallocationinfo,
                "ogrinfo": self.ogrinfo,
            },
        }


@dataclass(frozen=True)
class WorkflowPlan:
    scenarios: Tuple[Scenario, ...]
    models: Tuple[str, ...]
    ssps: Tuple[str, ...]
    periods: Tuple[str, ...]

    def payload(self) -> Dict[str, object]:
        return {
            "models": list(self.models),
            "ssps": ["ssp{}".format(value) for value in self.ssps],
            "periods": list(self.periods),
            "scenarios": [
                {
                    "name": item.name,
                    "model": item.model,
                    "ssp": "ssp{}".format(item.ssp),
                    "period": item.period,
                }
                for item in self.scenarios
            ],
        }


def _csv_values(value: str) -> Tuple[str, ...]:
    return tuple(dict.fromkeys(part.strip() for part in value.split(",") if part.strip()))


def _required(section: configparser.SectionProxy, name: str) -> str:
    value = section.get(name, "").strip()
    if not value:
        raise InputError("Workflow configuration is missing [{}] {}".format(section.name, name))
    return value


def _resolved_path(value: str, base: Path) -> Path:
    path = Path(value).expanduser()
    if not path.is_absolute():
        path = base / path
    return path.resolve()


def _optional_int(section: configparser.SectionProxy, name: str) -> Optional[int]:
    value = section.get(name, "").strip()
    return int(value) if value else None


def _optional_float(section: configparser.SectionProxy, name: str) -> Optional[float]:
    value = section.get(name, "").strip()
    return float(value) if value else None


def load_workflow_config(path: Path) -> WorkflowConfig:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise InputError("Workflow configuration does not exist: {}".format(path))
    parser = configparser.ConfigParser(interpolation=None)
    try:
        with path.open("r", encoding="utf-8-sig") as handle:
            parser.read_file(handle)
    except (configparser.Error, UnicodeError) as error:
        raise InputError("Cannot read workflow configuration: {}".format(error))
    for section_name in ("inputs", "analysis"):
        if not parser.has_section(section_name):
            raise InputError("Workflow configuration is missing [{}]".format(section_name))
    inputs = parser["inputs"]
    analysis = parser["analysis"]
    parameters = parser["parameters"] if parser.has_section("parameters") else parser["DEFAULT"]
    base = path.parent

    samples_value = inputs.get("samples", "").strip()
    groups_value = inputs.get("sample_groups_dir", "").strip()
    if bool(samples_value) == bool(groups_value):
        raise InputError(
            "Set exactly one of [inputs] samples or sample_groups_dir"
        )
    current_mask = _resolved_path(_required(inputs, "current_mask"), base)
    future_mask_value = inputs.get("future_mask", "").strip()
    future_mask = _resolved_path(future_mask_value, base) if future_mask_value else current_mask
    expected_sites = _optional_int(analysis, "expected_sites")
    if expected_sites is not None and expected_sites < 1:
        raise InputError("[analysis] expected_sites must be positive")
    minimum_models = parameters.getint("minimum_models", fallback=2)
    if minimum_models < 2:
        raise InputError("[parameters] minimum_models must be at least 2")
    predictors = _csv_values(parameters.get("predictors", ",".join(DEFAULT_PREDICTORS)))
    if not predictors:
        raise InputError("[parameters] predictors cannot be empty")

    config = WorkflowConfig(
        config_path=path,
        vcf=_resolved_path(_required(inputs, "vcf"), base),
        samples=_resolved_path(samples_value, base) if samples_value else None,
        sample_groups_dir=_resolved_path(groups_value, base) if groups_value else None,
        coordinates=_resolved_path(_required(inputs, "coordinates"), base),
        present_climate=_resolved_path(_required(inputs, "present_climate"), base),
        future_climate=_resolved_path(_required(inputs, "future_climate"), base),
        current_mask=current_mask,
        future_mask=future_mask,
        output_dir=_resolved_path(_required(analysis, "output_dir"), base),
        models=_csv_values(analysis.get("models", "")),
        ssps=tuple(
            value.lower().replace("ssp", "")
            for value in _csv_values(analysis.get("ssps", ""))
        ),
        periods=_csv_values(analysis.get("periods", "")),
        predictors=predictors,
        expected_sites=expected_sites,
        min_population_samples=parameters.getint("min_population_samples", fallback=3),
        warn_population_samples=parameters.getint("warn_population_samples", fallback=5),
        ntree=parameters.getint("ntree", fallback=500),
        nbin=parameters.getint("nbin", fallback=1001),
        corr_threshold=parameters.getfloat("corr_threshold", fallback=0.5),
        max_level=_optional_float(parameters, "max_level"),
        seed=parameters.getint("seed", fallback=1),
        forward_radii_km=tuple(
            parse_forward_radii(
                parameters.get("forward_radii", "100,250,500,1000,inf")
            )
        ),
        initial_knn_k=parameters.getint("initial_knn_k", fallback=64),
        batch_size=parameters.getint("batch_size", fallback=10000),
        verify_sample=parameters.getint("verify_sample", fallback=10),
        tie_tolerance=parameters.getfloat("tie_tolerance", fallback=1e-12),
        map_radius=parameters.get("map_radius", "unlimited").strip() or "unlimited",
        minimum_models=minimum_models,
        supplied_bio_tolerance=parameters.getfloat("supplied_bio_tolerance", fallback=1e-6),
        rscript=parameters.get("rscript", "").strip() or None,
        gdalinfo=parameters.get("gdalinfo", "gdalinfo").strip() or "gdalinfo",
        gdallocationinfo=parameters.get("gdallocationinfo", "gdallocationinfo").strip() or "gdallocationinfo",
        ogrinfo=parameters.get("ogrinfo", "ogrinfo").strip() or "ogrinfo",
    )
    required_files = [
        ("VCF", config.vcf),
        ("population coordinates", config.coordinates),
        ("current mask", config.current_mask),
        ("future mask", config.future_mask),
    ]
    if config.samples is not None:
        required_files.append(("sample table", config.samples))
    for label, input_path in required_files:
        if not input_path.is_file():
            raise InputError("{} does not exist: {}".format(label, input_path))
    if config.sample_groups_dir is not None and not config.sample_groups_dir.is_dir():
        raise InputError("Sample-group directory does not exist: {}".format(config.sample_groups_dir))
    return config


def build_workflow_plan(config: WorkflowConfig) -> WorkflowPlan:
    scenarios = discover_scenarios(
        config.future_climate,
        models=config.models,
        ssps=config.ssps,
        periods=config.periods,
    )
    models = config.models or tuple(sorted({item.model for item in scenarios}))
    ssps = config.ssps or tuple(sorted({item.ssp for item in scenarios}, key=int))
    periods = config.periods or tuple(sorted({item.period for item in scenarios}))
    if len(models) < config.minimum_models:
        raise InputError(
            "Complete GF workflow requires at least {} climate models; found {}".format(
                config.minimum_models, len(models)
            )
        )
    scenario_by_key = {(item.model, item.ssp, item.period): item for item in scenarios}
    missing = [
        (model, ssp, period)
        for period in periods
        for ssp in ssps
        for model in models
        if (model, ssp, period) not in scenario_by_key
    ]
    if missing:
        preview = ", ".join("{}/{}/{}".format(*item) for item in missing[:10])
        raise InputError("Workflow requires a complete model x SSP x period grid; missing {}".format(preview))
    ordered = tuple(
        scenario_by_key[(model, ssp, period)]
        for period in periods
        for ssp in ssps
        for model in models
    )
    return WorkflowPlan(ordered, tuple(models), tuple(ssps), tuple(periods))


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_sha256(value: object) -> str:
    encoded = json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _source_snapshot(config: WorkflowConfig, plan: WorkflowPlan) -> Dict[str, object]:
    files = {
        "vcf": config.vcf,
        "coordinates": config.coordinates,
    }
    if config.samples is not None:
        files["samples"] = config.samples
    file_rows = {
        label: {"path": str(path), "sha256": _sha256_file(path)}
        for label, path in files.items()
    }
    mask_rows: Dict[str, List[Dict[str, object]]] = {}
    for label, mask in (("current_mask", config.current_mask), ("future_mask", config.future_mask)):
        components = sorted(path for path in mask.parent.glob(mask.stem + ".*") if path.is_file())
        mask_rows[label] = [
            {"path": str(path.resolve()), "sha256": _sha256_file(path)} for path in components
        ]
    raster_paths = list(discover_bioclim_rasters(config.present_climate))
    for scenario in plan.scenarios:
        raster_paths.extend(scenario.rasters)
    raster_rows = []
    for path in dict.fromkeys(raster_paths):
        stat = path.stat()
        raster_rows.append(
            {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
        )
    sample_group_rows = []
    if config.sample_groups_dir is not None:
        for path in sorted(
            item
            for item in config.sample_groups_dir.iterdir()
            if item.is_file() and not item.name.startswith(".")
        ):
            sample_group_rows.append(
                {"path": str(path.resolve()), "sha256": _sha256_file(path)}
            )
        if not sample_group_rows:
            raise InputError("Sample-group directory contains no non-hidden files")
        _sample_design_text(config.sample_groups_dir)
    return {
        "files": file_rows,
        "sample_groups": sample_group_rows,
        "masks": mask_rows,
        "rasters": raster_rows,
    }


def _sample_design_text(group_dir: Path) -> str:
    rows = ["sample_id\tpopulation_id"]
    seen_samples: Dict[str, str] = {}
    group_files = sorted(
        path for path in group_dir.iterdir() if path.is_file() and not path.name.startswith(".")
    )
    if not group_files:
        raise InputError("Sample-group directory contains no non-hidden files: {}".format(group_dir))
    for path in group_files:
        population_id = path.name.strip()
        if not population_id:
            raise InputError("A sample-group file has an empty population name")
        sample_count = 0
        for line_number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), start=1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 1:
                raise InputError(
                    "Sample-group row must contain exactly one sample ID: {}:{}".format(
                        path, line_number
                    )
                )
            sample_id = fields[0]
            if sample_id in seen_samples:
                raise InputError(
                    "Sample {} occurs in both {} and {}".format(
                        sample_id, seen_samples[sample_id], population_id
                    )
                )
            seen_samples[sample_id] = population_id
            rows.append("{}\t{}".format(sample_id, population_id))
            sample_count += 1
        if sample_count == 0:
            raise InputError("Sample-group file is empty: {}".format(path))
    return "\n".join(rows) + "\n"


def _materialize_sample_design(group_dir: Path, output_path: Path) -> Path:
    text = _sample_design_text(group_dir)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        if output_path.read_text(encoding="utf-8") != text:
            raise InputError("Existing derived sample table does not match sample-group inputs")
    else:
        temporary = output_path.with_name(".{}.tmp".format(output_path.name))
        temporary.write_text(text, encoding="utf-8")
        os.replace(str(temporary), str(output_path))
    return output_path


def _write_json_atomic(path: Path, value: object) -> None:
    temporary = path.with_name(".{}.tmp".format(path.name))
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write("\n")
    os.replace(str(temporary), str(path))


def _completed_stage(
    output_dir: Path,
    module: str,
    required_files: Sequence[str],
) -> Optional[Dict[str, object]]:
    if not output_dir.exists():
        return None
    if not output_dir.is_dir():
        raise InputError("Workflow stage path is not a directory: {}".format(output_dir))
    manifest_path = output_dir / "run_manifest.json"
    if not manifest_path.is_file():
        raise InputError("Existing workflow stage is incomplete (no manifest): {}".format(output_dir))
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except (ValueError, UnicodeError) as error:
        raise InputError("Cannot read existing stage manifest {}: {}".format(manifest_path, error))
    if manifest.get("module") != module or manifest.get("engineering_status") != "completed":
        raise InputError("Existing workflow stage is not a completed {} run: {}".format(module, output_dir))
    missing = [name for name in required_files if not (output_dir / name).is_file()]
    if missing:
        raise InputError(
            "Existing {} stage is missing outputs: {}".format(module, ", ".join(missing))
        )
    return manifest


def _stage_record(
    key: str,
    output_dir: Path,
    manifest: Dict[str, object],
    reused: bool,
) -> Dict[str, object]:
    manifest_path = output_dir / "run_manifest.json"
    return {
        "key": key,
        "module": manifest.get("module"),
        "status": "reused" if reused else "completed",
        "scientific_status": manifest.get("scientific_status"),
        "output_dir": str(output_dir.resolve()),
        "manifest": str(manifest_path.resolve()),
        "manifest_sha256": _sha256_file(manifest_path),
    }


def run_gf_workflow(
    config_path: Path,
    dry_run: bool = False,
    progress: Optional[Callable[[str], None]] = None,
) -> Dict[str, object]:
    config = load_workflow_config(config_path)
    plan = build_workflow_plan(config)
    config_payload = config.payload()
    source_snapshot = _source_snapshot(config, plan)
    signature_payload = {"config": config_payload, "sources": source_snapshot}
    workflow_signature = _json_sha256(signature_payload)
    planned_stages = 3 + len(plan.scenarios) + len(plan.periods) * len(plan.ssps) + len(plan.periods)
    if dry_run:
        return {
            "module": "gf-workflow",
            "bioma_version": __version__,
            "engineering_status": "planned",
            "workflow_signature_sha256": workflow_signature,
            "config": config_payload,
            "selection": plan.payload(),
            "counts": {"scenarios": len(plan.scenarios), "planned_stages": planned_stages},
        }

    announce = progress or (lambda message: None)
    output_root = config.output_dir
    if output_root.exists() and not output_root.is_dir():
        raise InputError("Workflow output path is not a directory: {}".format(output_root))
    output_root.mkdir(parents=True, exist_ok=True)
    resolved_config_path = output_root / "resolved_config.json"
    if resolved_config_path.exists():
        existing = json.loads(resolved_config_path.read_text(encoding="utf-8"))
        if existing.get("workflow_signature_sha256") != workflow_signature:
            raise InputError(
                "Workflow inputs or configuration changed for {}; choose a new output_dir".format(output_root)
            )
    else:
        unrelated = [path for path in output_root.iterdir() if path.name != resolved_config_path.name]
        if unrelated:
            raise InputError(
                "Output directory is not an initialized BioMA workflow directory: {}".format(output_root)
            )
        _write_json_atomic(
            resolved_config_path,
            {
                "workflow_signature_sha256": workflow_signature,
                "config": config_payload,
                "source_snapshot": source_snapshot,
            },
        )

    started = time.monotonic()
    state_path = output_root / "workflow_state.json"
    state: Dict[str, object] = {
        "module": "gf-workflow",
        "bioma_version": __version__,
        "workflow_signature_sha256": workflow_signature,
        "status": "running",
        "updated_at_utc": datetime.now(timezone.utc).isoformat(),
        "steps": {},
    }
    _write_json_atomic(state_path, state)
    stage_records: List[Dict[str, object]] = []
    sample_path = config.samples
    if sample_path is None:
        if config.sample_groups_dir is None:
            raise InputError("No sample design input was configured")
        sample_path = _materialize_sample_design(
            config.sample_groups_dir,
            output_root / "00_inputs" / "sample_design.tsv",
        )

    def run_stage(
        key: str,
        module: str,
        output_dir: Path,
        required_files: Sequence[str],
        function: Callable[[], Dict[str, object]],
    ) -> Dict[str, object]:
        existing_manifest = _completed_stage(output_dir, module, required_files)
        reused = existing_manifest is not None
        announce("Reuse {}".format(key) if reused else "Run {}".format(key))
        manifest = existing_manifest if existing_manifest is not None else function()
        record = _stage_record(key, output_dir, manifest, reused)
        stage_records.append(record)
        steps = state["steps"]
        if isinstance(steps, dict):
            steps[key] = record
        state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        _write_json_atomic(state_path, state)
        return manifest

    frequency_dir = output_root / "01_frequency"
    climate_dir = output_root / "02_climate"
    model_dir = output_root / "03_model"
    offsets_root = output_root / "04_offsets"
    maps_root = output_root / "05_offset_plots"
    distance_root = output_root / "06_forward_distance_plots"
    try:
        run_stage(
            "01_frequency",
            "gf-frequency",
            frequency_dir,
            ("population_alt_frequency.tsv", "variant_manifest.tsv"),
            lambda: build_population_frequency(
                vcf_path=config.vcf,
                samples_path=sample_path,
                output_dir=frequency_dir,
                min_population_samples=config.min_population_samples,
                warn_population_samples=config.warn_population_samples,
                expected_sites=config.expected_sites,
            ),
        )
        run_stage(
            "02_climate",
            "climate-prepare",
            climate_dir,
            ("population_environment.tsv", "current_background.tsv.gz", "future_manifest.tsv"),
            lambda: prepare_climate_inputs(
                present_dir=config.present_climate,
                future_root=config.future_climate,
                coordinates_path=config.coordinates,
                current_mask=config.current_mask,
                future_mask=config.future_mask,
                output_dir=climate_dir,
                models=plan.models,
                ssps=plan.ssps,
                periods=plan.periods,
                validate_only=False,
                rscript=config.rscript,
                gdalinfo=config.gdalinfo,
                gdallocationinfo=config.gdallocationinfo,
                ogrinfo=config.ogrinfo,
                supplied_bio_tolerance=config.supplied_bio_tolerance,
            ),
        )
        run_stage(
            "03_model",
            "gf-train",
            model_dir,
            ("all_gfmod.data", "training_alignment.tsv"),
            lambda: train_gradient_forest(
                frequency_path=frequency_dir / "population_alt_frequency.tsv",
                environment_path=climate_dir / "population_environment.tsv",
                output_dir=model_dir,
                predictors=config.predictors,
                ntree=config.ntree,
                nbin=config.nbin,
                corr_threshold=config.corr_threshold,
                max_level=config.max_level,
                seed=config.seed,
                expected_sites=config.expected_sites,
                rscript=config.rscript,
            ),
        )

        scenario_outputs: Dict[Tuple[str, str, str], Path] = {}
        for scenario in plan.scenarios:
            scenario_dir = offsets_root / scenario.name
            key = "04_offset/{}".format(scenario.name)
            run_stage(
                key,
                "gf-offset",
                scenario_dir,
                ("all_offsets.tsv.gz", "scenario_summary.tsv"),
                lambda scenario=scenario, scenario_dir=scenario_dir: compute_gf_offsets(
                    model=model_dir,
                    climate_dir=climate_dir,
                    scenario=scenario.name,
                    output_dir=scenario_dir,
                    forward_radii_km=config.forward_radii_km,
                    initial_knn_k=config.initial_knn_k,
                    batch_size=config.batch_size,
                    verify_sample=config.verify_sample,
                    tie_tolerance=config.tie_tolerance,
                    rscript=config.rscript,
                ),
            )
            scenario_outputs[(scenario.model, scenario.ssp, scenario.period)] = scenario_dir

        for period in plan.periods:
            for ssp in plan.ssps:
                inputs = [scenario_outputs[(model, ssp, period)] for model in plan.models]
                plot_dir = maps_root / period / "ssp{}".format(ssp)
                key = "05_offset_plot/{}/ssp{}".format(period, ssp)
                run_stage(
                    key,
                    "offset-plot",
                    plot_dir,
                    ("ensemble_mean_offsets.tsv.gz", "all_offsets_maps.png", "offset_relationships.png"),
                    lambda inputs=inputs, plot_dir=plot_dir, period=period, ssp=ssp: plot_offsets(
                        inputs=inputs,
                        output_dir=plot_dir,
                        radius=config.map_radius,
                        models=plan.models,
                        period=period,
                        ssp="ssp{}".format(ssp),
                        populations=climate_dir / "population_environment.tsv",
                        boundary=config.current_mask,
                        minimum_models=config.minimum_models,
                        rscript=config.rscript,
                    ),
                )

        for period in plan.periods:
            inputs = [
                scenario_outputs[(model, ssp, period)]
                for ssp in plan.ssps
                for model in plan.models
            ]
            plot_dir = distance_root / period
            key = "06_forward_distance_plot/{}".format(period)
            run_stage(
                key,
                "forward-distance-plot",
                plot_dir,
                ("forward_offset_by_distance.tsv", "forward_offset_by_distance.png"),
                lambda inputs=inputs, plot_dir=plot_dir, period=period: plot_forward_distance(
                    inputs=inputs,
                    output_dir=plot_dir,
                    models=plan.models,
                    period=period,
                    ssps=["ssp{}".format(value) for value in plan.ssps],
                    minimum_models=config.minimum_models,
                    rscript=config.rscript,
                ),
            )
    except Exception as error:
        state["status"] = "failed"
        state["last_error"] = str(error)
        state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
        _write_json_atomic(state_path, state)
        raise

    scientific_statuses = {record.get("scientific_status") for record in stage_records}
    workflow_manifest: Dict[str, object] = {
        "module": "gf-workflow",
        "bioma_version": __version__,
        "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        "engineering_status": "completed",
        "scientific_status": (
            "completed_with_warning"
            if "completed_with_warning" in scientific_statuses
            else "completed_validated"
        ),
        "inputs": {
            "configuration": str(config.config_path),
            "configuration_sha256": _sha256_file(config.config_path),
            "workflow_signature_sha256": workflow_signature,
        },
        "selection": plan.payload(),
        "counts": {
            "models": len(plan.models),
            "ssps": len(plan.ssps),
            "periods": len(plan.periods),
            "scenarios": len(plan.scenarios),
            "stages": len(stage_records),
            "reused_stages": sum(record["status"] == "reused" for record in stage_records),
        },
        "steps": stage_records,
        "outputs": {
            "frequency": str(frequency_dir),
            "climate": str(climate_dir),
            "model": str(model_dir),
            "offsets": str(offsets_root),
            "offset_plots": str(maps_root),
            "forward_distance_plots": str(distance_root),
        },
        "runtime": {"elapsed_seconds": time.monotonic() - started},
    }
    _write_json_atomic(output_root / "run_manifest.json", workflow_manifest)
    state["status"] = "completed"
    state["updated_at_utc"] = datetime.now(timezone.utc).isoformat()
    _write_json_atomic(state_path, state)
    announce("Complete {}".format(output_root))
    return workflow_manifest
