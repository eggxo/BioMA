"""Build population ALT-allele-frequency matrices for Gradient Forest."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import os
import platform
import re
import shutil
import tempfile
from collections import OrderedDict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, TextIO, Tuple

from . import __version__


class InputError(RuntimeError):
    """Raised when an input violates the GF frequency data contract."""


@dataclass
class SampleDesign:
    sample_to_population: Dict[str, str]
    population_to_samples: "OrderedDict[str, List[str]]"


@dataclass
class VariantRecord:
    index: int
    variant_id: str
    chrom: str
    pos: str
    ref: str
    alt: str


def _open_text_auto(path: Path) -> TextIO:
    """Open plain or gzip text by file signature, not filename extension."""
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


def read_sample_design(path: Path) -> SampleDesign:
    if not path.is_file():
        raise InputError("Sample metadata file does not exist: {}".format(path))

    sample_to_population: Dict[str, str] = {}
    population_to_samples: "OrderedDict[str, List[str]]" = OrderedDict()
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "population_id"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise InputError(
                "Sample metadata must be tab-separated and contain columns: "
                "sample_id, population_id"
            )
        for line_number, row in enumerate(reader, start=2):
            sample_id = (row.get("sample_id") or "").strip()
            population_id = (row.get("population_id") or "").strip()
            if not sample_id or not population_id:
                raise InputError(
                    "Empty sample_id or population_id at metadata line {}".format(
                        line_number
                    )
                )
            if sample_id in sample_to_population:
                raise InputError(
                    "Sample occurs more than once in metadata: {}".format(sample_id)
                )
            sample_to_population[sample_id] = population_id
            population_to_samples.setdefault(population_id, []).append(sample_id)

    if not sample_to_population:
        raise InputError("Sample metadata contains no samples")
    return SampleDesign(sample_to_population, population_to_samples)


def _split_gt(gt: str) -> List[str]:
    return re.split(r"[|/]", gt)


def _format_frequency(value: Optional[float]) -> str:
    if value is None:
        return "NA"
    return "{:.12g}".format(value)


def _write_tsv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def build_population_frequency(
    vcf_path: Path,
    samples_path: Path,
    output_dir: Path,
    min_population_samples: int = 3,
    warn_population_samples: int = 5,
    expected_sites: Optional[int] = None,
) -> Dict[str, object]:
    """Compute population ALT frequencies and write an auditable result bundle."""
    vcf_path = vcf_path.expanduser().resolve()
    samples_path = samples_path.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()

    if not vcf_path.is_file():
        raise InputError("VCF file does not exist: {}".format(vcf_path))
    if min_population_samples < 1:
        raise InputError("min_population_samples must be at least 1")
    if warn_population_samples < min_population_samples:
        raise InputError(
            "warn_population_samples cannot be smaller than min_population_samples"
        )
    if output_dir.exists():
        raise InputError(
            "Output directory already exists; choose a new directory: {}".format(
                output_dir
            )
        )

    design = read_sample_design(samples_path)
    for population_id, members in design.population_to_samples.items():
        if len(members) < min_population_samples:
            raise InputError(
                "Population {} has {} samples; minimum is {}".format(
                    population_id, len(members), min_population_samples
                )
            )

    population_ids = list(design.population_to_samples.keys())
    frequencies: Dict[str, List[Optional[float]]] = {
        population_id: [] for population_id in population_ids
    }
    min_called_chromosomes: Dict[str, Optional[int]] = {
        population_id: None for population_id in population_ids
    }
    sum_called_chromosomes: Dict[str, int] = {
        population_id: 0 for population_id in population_ids
    }
    complete_genotypes: Dict[str, int] = {
        population_id: 0 for population_id in population_ids
    }
    missing_frequency_sites: Dict[str, int] = {
        population_id: 0 for population_id in population_ids
    }

    variants: List[VariantRecord] = []
    warnings: List[Tuple[str, str, str]] = []
    vcf_samples: Optional[List[str]] = None
    population_indexes: Dict[str, List[int]] = {}
    seen_variant_ids = set()

    with _open_text_auto(vcf_path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 10:
                    raise InputError("VCF contains no genotype samples")
                vcf_samples = fields[9:]
                if len(vcf_samples) != len(set(vcf_samples)):
                    raise InputError("VCF header contains duplicate sample IDs")
                sample_index = {sample_id: i for i, sample_id in enumerate(vcf_samples)}
                missing = [
                    sample_id
                    for sample_id in design.sample_to_population
                    if sample_id not in sample_index
                ]
                if missing:
                    preview = ", ".join(missing[:10])
                    raise InputError(
                        "{} metadata samples are absent from the VCF: {}".format(
                            len(missing), preview
                        )
                    )
                extra = [
                    sample_id
                    for sample_id in vcf_samples
                    if sample_id not in design.sample_to_population
                ]
                if extra:
                    warnings.append(
                        (
                            "VCF_EXTRA_SAMPLES",
                            "run",
                            "{} VCF samples are not assigned in metadata and were ignored".format(
                                len(extra)
                            ),
                        )
                    )
                population_indexes = {
                    population_id: [sample_index[s] for s in members]
                    for population_id, members in design.population_to_samples.items()
                }
                continue
            if line.startswith("#"):
                continue
            if vcf_samples is None:
                raise InputError("VCF #CHROM header was not found before data lines")

            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9 + len(vcf_samples):
                raise InputError(
                    "VCF line {} has {} columns; expected {}".format(
                        line_number, len(fields), 9 + len(vcf_samples)
                    )
                )
            chrom, pos, raw_id, ref, alt = fields[0:5]
            alts = alt.split(",")
            if len(alts) != 1 or alts[0] == ".":
                raise InputError(
                    "GF frequency input must be biallelic; line {} has ALT={}".format(
                        line_number, alt
                    )
                )
            variant_id = raw_id if raw_id not in ("", ".") else "{}:{}:{}:{}".format(
                chrom, pos, ref, alt
            )
            if variant_id in seen_variant_ids:
                raise InputError("Duplicate variant ID: {}".format(variant_id))
            seen_variant_ids.add(variant_id)

            format_fields = fields[8].split(":")
            if "GT" not in format_fields:
                raise InputError(
                    "VCF line {} FORMAT does not contain GT".format(line_number)
                )
            gt_index = format_fields.index("GT")
            genotype_fields = fields[9:]

            variants.append(
                VariantRecord(len(variants) + 1, variant_id, chrom, pos, ref, alt)
            )
            for population_id in population_ids:
                alt_count = 0
                called_chromosomes = 0
                for sample_position in population_indexes[population_id]:
                    sample_parts = genotype_fields[sample_position].split(":")
                    gt = sample_parts[gt_index] if gt_index < len(sample_parts) else "."
                    alleles = _split_gt(gt)
                    called = [allele for allele in alleles if allele != "."]
                    if len(called) == len(alleles) and called:
                        complete_genotypes[population_id] += 1
                    for allele in called:
                        if allele not in ("0", "1"):
                            raise InputError(
                                "Unexpected allele index {} at variant {}".format(
                                    allele, variant_id
                                )
                            )
                        called_chromosomes += 1
                        alt_count += int(allele == "1")

                if called_chromosomes == 0:
                    frequency = None
                    missing_frequency_sites[population_id] += 1
                else:
                    frequency = alt_count / float(called_chromosomes)
                frequencies[population_id].append(frequency)
                sum_called_chromosomes[population_id] += called_chromosomes
                previous_min = min_called_chromosomes[population_id]
                if previous_min is None or called_chromosomes < previous_min:
                    min_called_chromosomes[population_id] = called_chromosomes

    if vcf_samples is None:
        raise InputError("VCF #CHROM header was not found")
    if not variants:
        raise InputError("VCF contains no variants")
    if expected_sites is not None and len(variants) != expected_sites:
        raise InputError(
            "VCF contains {} variants; expected {}".format(len(variants), expected_sites)
        )

    for population_id, members in design.population_to_samples.items():
        if len(members) < warn_population_samples:
            warnings.append(
                (
                    "LOW_POPULATION_SAMPLE_SIZE",
                    population_id,
                    "Population has {} samples; allele-frequency uncertainty may be high".format(
                        len(members)
                    ),
                )
            )
        if missing_frequency_sites[population_id]:
            warnings.append(
                (
                    "MISSING_POPULATION_FREQUENCY",
                    population_id,
                    "{} variants have no called alleles in this population".format(
                        missing_frequency_sites[population_id]
                    ),
                )
            )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=".{}-tmp-".format(output_dir.name), dir=str(output_dir.parent))
    )
    try:
        frequency_rows = (
            [population_id]
            + [_format_frequency(value) for value in frequencies[population_id]]
            for population_id in population_ids
        )
        _write_tsv(
            temporary_dir / "population_alt_frequency.tsv",
            ["population_id"] + [variant.variant_id for variant in variants],
            frequency_rows,
        )

        _write_tsv(
            temporary_dir / "variant_manifest.tsv",
            ["variant_index", "variant_id", "chrom", "pos", "ref", "alt"],
            (
                [v.index, v.variant_id, v.chrom, v.pos, v.ref, v.alt]
                for v in variants
            ),
        )

        qc_rows = []
        genotype_denominator = len(variants)
        for population_id in population_ids:
            n_samples = len(design.population_to_samples[population_id])
            total_genotypes = n_samples * genotype_denominator
            missing_rate = 1.0 - (
                complete_genotypes[population_id] / float(total_genotypes)
            )
            qc_rows.append(
                [
                    population_id,
                    n_samples,
                    len(variants),
                    min_called_chromosomes[population_id],
                    "{:.6f}".format(
                        sum_called_chromosomes[population_id] / float(len(variants))
                    ),
                    missing_frequency_sites[population_id],
                    "{:.8f}".format(missing_rate),
                    "warning"
                    if n_samples < warn_population_samples
                    or missing_frequency_sites[population_id]
                    else "pass",
                ]
            )
        _write_tsv(
            temporary_dir / "population_frequency_qc.tsv",
            [
                "population_id",
                "n_samples",
                "n_variants",
                "min_called_chromosomes",
                "mean_called_chromosomes",
                "n_missing_frequency_variants",
                "missing_genotype_rate",
                "qc_status",
            ],
            qc_rows,
        )

        _write_tsv(
            temporary_dir / "warnings.tsv",
            ["warning_code", "scope", "message"],
            warnings,
        )

        manifest: Dict[str, object] = {
            "module": "gf-frequency",
            "bioma_version": __version__,
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
            "engineering_status": "completed",
            "scientific_status": (
                "completed_with_warning" if warnings else "completed_validated"
            ),
            "inputs": {
                "vcf": str(vcf_path),
                "vcf_sha256": _sha256(vcf_path),
                "samples": str(samples_path),
                "samples_sha256": _sha256(samples_path),
            },
            "parameters": {
                "min_population_samples": min_population_samples,
                "warn_population_samples": warn_population_samples,
                "expected_sites": expected_sites,
            },
            "counts": {
                "vcf_samples": len(vcf_samples),
                "assigned_samples": len(design.sample_to_population),
                "populations": len(population_ids),
                "variants": len(variants),
                "warnings": len(warnings),
            },
            "runtime": {
                "python": platform.python_version(),
                "platform": platform.platform(),
            },
        }
        with (temporary_dir / "run_manifest.json").open(
            "w", encoding="utf-8"
        ) as handle:
            json.dump(manifest, handle, indent=2, ensure_ascii=False)
            handle.write("\n")

        os.replace(str(temporary_dir), str(output_dir))
        return manifest
    except Exception:
        shutil.rmtree(str(temporary_dir), ignore_errors=True)
        raise
