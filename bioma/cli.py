"""BioMA command-line interface."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional, Sequence

from . import __version__
from .climate_prepare import prepare_climate_inputs
from .gf_frequency import InputError, build_population_frequency
from .gf_offset import compute_gf_offsets, parse_forward_radii
from .gf_train import DEFAULT_PREDICTORS, train_gradient_forest
from .gf_workflow import run_gf_workflow
from .forward_distance_plot import plot_forward_distance
from .offset_plot import plot_offsets
from .rona import run_rona_workflow
from .mar import run_mar_workflow
from .load import run_load_workflow
from .vulnerability import run_vulnerability_workflow
from .niche import run_niche_workflow
from .wfmoment import run_wfmoment_workflow
from .project import MODULE_ORDER, run_project_workflow
from .doctor import format_doctor_report, run_doctor


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="bioma",
        description="Auditable climate-vulnerability workflow modules",
    )
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    project = subparsers.add_parser(
        "project",
        aliases=["project-run"],
        help="Run selected BioMA modules as one resumable project",
    )
    project.add_argument("config", type=Path)
    project.add_argument("--dry-run", action="store_true")
    project.add_argument("--module", action="append", choices=MODULE_ORDER, default=[])
    project.add_argument("--overwrite", action="store_true")
    project.add_argument("--no-resume", action="store_true")

    workflow = subparsers.add_parser(
        "run",
        aliases=["gf-run"],
        help="Run the complete GF workflow from one INI configuration",
    )
    workflow.add_argument("config", type=Path)
    workflow.add_argument("--dry-run", action="store_true")

    frequency = subparsers.add_parser(
        "gf-frequency",
        help="Build a population ALT-frequency matrix for Gradient Forest",
    )
    frequency.add_argument("--vcf", required=True, type=Path)
    frequency.add_argument("--samples", required=True, type=Path)
    frequency.add_argument("--outdir", required=True, type=Path)
    frequency.add_argument("--min-pop-samples", type=int, default=3)
    frequency.add_argument("--warn-pop-samples", type=int, default=5)
    frequency.add_argument("--expected-sites", type=int)

    climate = subparsers.add_parser(
        "climate-prepare",
        help="Validate and standardize present/future BIO1-BIO19 climate data",
    )
    climate.add_argument("--present-dir", required=True, type=Path)
    climate.add_argument("--future-root", required=True, type=Path)
    climate.add_argument("--coordinates", required=True, type=Path)
    climate.add_argument("--current-mask", required=True, type=Path)
    climate.add_argument("--future-mask", type=Path)
    climate.add_argument("--outdir", required=True, type=Path)
    climate.add_argument("--model", action="append", default=[])
    climate.add_argument("--ssp", action="append", default=[])
    climate.add_argument("--period", action="append", default=[])
    climate.add_argument("--validate-only", action="store_true")
    climate.add_argument("--rscript")
    climate.add_argument("--gdalinfo", default="gdalinfo")
    climate.add_argument("--gdallocationinfo", default="gdallocationinfo")
    climate.add_argument("--ogrinfo", default="ogrinfo")
    climate.add_argument("--supplied-bio-tolerance", type=float, default=1e-6)

    train = subparsers.add_parser(
        "gf-train",
        help="Train a reproducible Gradient Forest model",
    )
    train.add_argument("--frequencies", required=True, type=Path)
    train.add_argument("--environment", required=True, type=Path)
    train.add_argument("--outdir", required=True, type=Path)
    train.add_argument("--predictors", default=",".join(DEFAULT_PREDICTORS))
    train.add_argument("--ntree", type=int, default=500)
    train.add_argument("--nbin", type=int, default=1001)
    train.add_argument("--corr-threshold", type=float, default=0.5)
    train.add_argument("--max-level", type=float)
    train.add_argument("--seed", type=int, default=1)
    train.add_argument("--expected-sites", type=int)
    train.add_argument("--rscript")

    offset = subparsers.add_parser(
        "gf-offset",
        help="Compute local, forward, and reverse Gradient Forest offsets",
    )
    offset.add_argument("--model", required=True, type=Path)
    offset.add_argument("--climate", required=True, type=Path)
    offset.add_argument("--scenario", required=True)
    offset.add_argument("--outdir", required=True, type=Path)
    offset.add_argument("--forward-radii", default="100,250,500,1000,inf")
    offset.add_argument("--initial-knn-k", type=int, default=64)
    offset.add_argument("--batch-size", type=int, default=10000)
    offset.add_argument("--verify-sample", type=int, default=10)
    offset.add_argument("--tie-tolerance", type=float, default=1e-12)
    offset.add_argument("--rscript")

    plot = subparsers.add_parser(
        "offset-plot",
        help="Plot one or more final GF offset tables",
    )
    plot.add_argument("--input", action="append", required=True, type=Path)
    plot.add_argument("--outdir", required=True, type=Path)
    plot.add_argument("--radius", default="unlimited")
    plot.add_argument("--models", default="")
    plot.add_argument("--period")
    plot.add_argument("--ssp")
    plot.add_argument("--populations", type=Path)
    plot.add_argument("--boundary", type=Path)
    plot.add_argument("--minimum-models", type=int, default=2)
    plot.add_argument("--rscript")

    distance_plot = subparsers.add_parser(
        "forward-distance-plot",
        help="Plot ensemble-mean forward offset across migration-distance limits",
    )
    distance_plot.add_argument("--input", action="append", required=True, type=Path)
    distance_plot.add_argument("--outdir", required=True, type=Path)
    distance_plot.add_argument("--models", default="")
    distance_plot.add_argument("--period")
    distance_plot.add_argument("--ssp", action="append", default=[])
    distance_plot.add_argument("--minimum-models", type=int, default=2)
    distance_plot.add_argument("--rscript")

    rona = subparsers.add_parser(
        "rona",
        aliases=["rona-run"],
        help="Run the RONA calculation and plotting workflow from one INI configuration",
    )
    rona.add_argument("config", type=Path)
    rona.add_argument("--dry-run", action="store_true")

    mar = subparsers.add_parser(
        "mar",
        aliases=["mar-run"],
        help="Run the MAR population-genomics workflow from one INI configuration",
    )
    mar.add_argument("config", type=Path)
    mar.add_argument("--dry-run", action="store_true")

    load = subparsers.add_parser(
        "load",
        aliases=["load-run"],
        help="Calculate and predict genetic loadM/loadD from derived VCFs",
    )
    load.add_argument("config", type=Path)
    load.add_argument("--dry-run", action="store_true")

    niche = subparsers.add_parser(
        "niche",
        aliases=["niche-run"],
        help="Tune MaxEnt variables and project continuous future suitability",
    )
    niche.add_argument("config", type=Path)
    niche.add_argument("--dry-run", action="store_true")

    vulnerability = subparsers.add_parser(
        "vulnerability",
        aliases=["vulnerability-run"],
        help="Combine GF, RONA, niche, loadM and loadD population vulnerability values",
    )
    vulnerability.add_argument("config", type=Path)
    vulnerability.add_argument("--dry-run", action="store_true")

    wfmoment = subparsers.add_parser(
        "wfmoment",
        aliases=["wfmoment-run"],
        help="Run the 2-D deme WFmoments genetic-diversity workflow",
    )
    wfmoment.add_argument("config", type=Path)
    wfmoment.add_argument("--dry-run", action="store_true")

    doctor = subparsers.add_parser(
        "doctor",
        aliases=["doctor-run"],
        help="Check BioMA runtimes, external tools, and scientific packages",
    )
    doctor.add_argument(
        "config",
        nargs="?",
        type=Path,
        help="Optional project or module INI used to scope required checks",
    )
    doctor.add_argument("--config", dest="config_option", type=Path, help=argparse.SUPPRESS)
    doctor.add_argument("--project", dest="project", type=Path, help="Alias for the project INI path")
    doctor.add_argument("--json", action="store_true", help="Emit a machine-readable JSON report")
    doctor.add_argument("--format", choices=("text", "json"), default="text", help="Report format (default: text)")
    doctor.add_argument("--strict", action="store_true", help="Treat optional checks as failures")
    doctor.add_argument("--python", "--python-executable", dest="python_executable", help="Python executable to probe")
    doctor.add_argument("--rscript", help="Rscript executable to probe")
    doctor.add_argument("--gdalinfo", help="gdalinfo executable to probe")
    doctor.add_argument("--gdallocationinfo", help="gdallocationinfo executable to probe")
    doctor.add_argument("--ogrinfo", help="ogrinfo executable to probe")
    doctor.add_argument("--java", help="Java executable to probe")
    doctor.add_argument("--maxent-jar", "--maxent", dest="maxent_jar", type=Path, help="Path to maxent.jar")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = _parser()
    args = parser.parse_args(argv)
    try:
        if args.command in ("project", "project-run"):
            manifest = run_project_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                only_modules=args.module,
                resume=False if args.no_resume else None,
                overwrite=args.overwrite,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0 if manifest.get("status") in {"planned", "complete"} else 1
        if args.command in ("run", "gf-run"):
            manifest = run_gf_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command in ("rona", "rona-run"):
            manifest = run_rona_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("mar", "mar-run"):
            manifest = run_mar_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("load", "load-run"):
            manifest = run_load_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("niche", "niche-run"):
            manifest = run_niche_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("vulnerability", "vulnerability-run"):
            manifest = run_vulnerability_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("wfmoment", "wfmoment-run"):
            manifest = run_wfmoment_workflow(
                config_path=args.config,
                dry_run=args.dry_run,
                progress=lambda message: print("[BioMA] {}".format(message), file=sys.stderr, flush=True),
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False, default=str))
            return 0
        if args.command in ("doctor", "doctor-run"):
            config_candidates = [value for value in (args.config, args.config_option, args.project) if value is not None]
            if len(config_candidates) > 1:
                parser.error("doctor accepts only one configuration path")
            report = run_doctor(
                config_path=config_candidates[0] if config_candidates else None,
                python=args.python_executable,
                rscript=args.rscript,
                gdalinfo=args.gdalinfo,
                gdallocationinfo=args.gdallocationinfo,
                ogrinfo=args.ogrinfo,
                java=args.java,
                maxent_jar=args.maxent_jar,
                strict=args.strict,
            )
            if args.json or args.format == "json":
                print(json.dumps(report, indent=2, ensure_ascii=False))
            else:
                print(format_doctor_report(report))
            return 0 if report.get("ok") else 1
        if args.command == "gf-frequency":
            manifest = build_population_frequency(
                vcf_path=args.vcf,
                samples_path=args.samples,
                output_dir=args.outdir,
                min_population_samples=args.min_pop_samples,
                warn_population_samples=args.warn_pop_samples,
                expected_sites=args.expected_sites,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command == "gf-offset":
            manifest = compute_gf_offsets(
                model=args.model,
                climate_dir=args.climate,
                scenario=args.scenario,
                output_dir=args.outdir,
                forward_radii_km=parse_forward_radii(args.forward_radii),
                initial_knn_k=args.initial_knn_k,
                batch_size=args.batch_size,
                verify_sample=args.verify_sample,
                tie_tolerance=args.tie_tolerance,
                rscript=args.rscript,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command == "offset-plot":
            models = [value.strip() for value in args.models.split(",") if value.strip()]
            manifest = plot_offsets(
                inputs=args.input,
                output_dir=args.outdir,
                radius=args.radius,
                models=models,
                period=args.period,
                ssp=args.ssp,
                populations=args.populations,
                boundary=args.boundary,
                minimum_models=args.minimum_models,
                rscript=args.rscript,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command == "forward-distance-plot":
            models = [value.strip() for value in args.models.split(",") if value.strip()]
            manifest = plot_forward_distance(
                inputs=args.input,
                output_dir=args.outdir,
                models=models,
                period=args.period,
                ssps=args.ssp,
                minimum_models=args.minimum_models,
                rscript=args.rscript,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command == "gf-train":
            predictors = [value.strip() for value in args.predictors.split(",") if value.strip()]
            manifest = train_gradient_forest(
                frequency_path=args.frequencies,
                environment_path=args.environment,
                output_dir=args.outdir,
                predictors=predictors,
                ntree=args.ntree,
                nbin=args.nbin,
                corr_threshold=args.corr_threshold,
                max_level=args.max_level,
                seed=args.seed,
                expected_sites=args.expected_sites,
                rscript=args.rscript,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
        if args.command == "climate-prepare":
            manifest = prepare_climate_inputs(
                present_dir=args.present_dir,
                future_root=args.future_root,
                coordinates_path=args.coordinates,
                current_mask=args.current_mask,
                future_mask=args.future_mask,
                output_dir=args.outdir,
                models=args.model,
                ssps=args.ssp,
                periods=args.period,
                validate_only=args.validate_only,
                rscript=args.rscript,
                gdalinfo=args.gdalinfo,
                gdallocationinfo=args.gdallocationinfo,
                ogrinfo=args.ogrinfo,
                supplied_bio_tolerance=args.supplied_bio_tolerance,
            )
            print(json.dumps(manifest, indent=2, ensure_ascii=False))
            return 0
    except InputError as error:
        print("BioMA input error: {}".format(error), file=sys.stderr)
        return 2
    except Exception as error:
        print("BioMA runtime error: {}".format(error), file=sys.stderr)
        return 1
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
