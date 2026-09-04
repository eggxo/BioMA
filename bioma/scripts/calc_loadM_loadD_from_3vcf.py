#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
import gzip
import math
import statistics
from collections import defaultdict

# ------------------------------------------------------------
# Usage
# ------------------------------------------------------------
USAGE = """
Usage:
  python calc_loadM_loadD_from_3vcf.py \
    <syn_vcf> <del_vcf> <nonsyn_vcf> <pop_dir> <out_prefix>

Example:
  python calc_loadM_loadD_from_3vcf.py \
    neutral_SYN.vcf.gz \
    deleterious_MIS_LOF.vcf.gz \
    nonsyn_all.vcf.gz \
    /path/to/population_lists \
    pade_load

Notes:
  1) 默认这3个VCF里 ALT 就是 derived allele。
  2) 默认都是 biallelic SNP VCF。
  3) 群体目录 pop_dir 下每个文件名视为群体名，文件内容每行一个样本名。
  4) 计算逻辑保持不变：
       Ps = syn derived allele frequency per individual
       Pn = nonsyn derived allele frequency per individual
       Pd = deleterious derived allele frequency per individual
       loadM = Pn / (Pn + Ps)
       loadD = Pd / (Pn + Ps)
"""

if len(sys.argv) != 6:
    sys.stderr.write(USAGE + "\n")
    sys.exit(1)

syn_vcf = sys.argv[1]
del_vcf = sys.argv[2]
nonsyn_vcf = sys.argv[3]
pop_dir = sys.argv[4]
out_prefix = sys.argv[5]


# ------------------------------------------------------------
# Helpers
# ------------------------------------------------------------
def open_text_auto(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gt_to_dosage(sample_field):
    """
    sample_field like:
      0/0:...
      0|1:...
      ./.:...
      .|.:...
    Return:
      0, 1, 2 or None if missing / non-diploid / unreadable
    """
    gt = sample_field.split(":", 1)[0]

    if gt in {"./.", ".|.", ".", "./", ".|"}:
        return None

    sep = "/" if "/" in gt else "|" if "|" in gt else None
    if sep is None:
        return None

    alleles = gt.split(sep)
    if len(alleles) != 2:
        return None

    if "." in alleles:
        return None

    dosage = 0
    for a in alleles:
        if a == "0":
            dosage += 0
        elif a == "1":
            dosage += 1
        else:
            # multi-allelic or malformed GT, skip safely
            return None
    return dosage


def read_pop_file_auto(fp):
    last_err = None
    for enc in ("utf-8-sig", "gb18030"):
        try:
            with open(fp, "r", encoding=enc) as f:
                return [line.strip() for line in f if line.strip()]
        except UnicodeDecodeError as e:
            last_err = e
    raise UnicodeError(f"Cannot decode file: {fp}") from last_err


def load_populations(pop_dir):
    sample2pop = {}
    pop2samples = defaultdict(list)

    if not os.path.isdir(pop_dir):
        raise FileNotFoundError(f"Population directory not found: {pop_dir}")

    bad_files = []

    for fn in sorted(os.listdir(pop_dir)):
        if fn.startswith("."):
            continue

        fp = os.path.join(pop_dir, fn)
        if not os.path.isfile(fp):
            continue

        try:
            lines = read_pop_file_auto(fp)
        except Exception as e:
            bad_files.append((fp, str(e)))
            continue

        pop = fn
        for line in lines:
            s = line.split()[0]
            if not s:
                continue
            if s in sample2pop:
                raise ValueError(f"Sample found in multiple population files: {s}")
            sample2pop[s] = pop
            pop2samples[pop].append(s)

    if bad_files:
        sys.stderr.write("\nWARNING: these files in pop_dir could not be decoded and were skipped:\n")
        for fp, msg in bad_files:
            sys.stderr.write(f"  {fp}\t{msg}\n")
        sys.stderr.write("\n")

    if len(sample2pop) == 0:
        raise ValueError(f"No valid samples were read from pop_dir: {pop_dir}")

    return sample2pop, pop2samples


def read_vcf_header_samples(vcf_path):
    with open_text_auto(vcf_path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                return parts[9:]
    raise ValueError(f"No #CHROM header found in {vcf_path}")


def mean_or_na(vals):
    vals2 = [x for x in vals if x is not None and not math.isnan(x)]
    if len(vals2) == 0:
        return "NA"
    return sum(vals2) / len(vals2)


def sd_or_na(vals):
    vals2 = [x for x in vals if x is not None and not math.isnan(x)]
    if len(vals2) <= 1:
        return "NA"
    return statistics.stdev(vals2)


def fmt(x):
    if x == "NA" or x is None:
        return "NA"
    if isinstance(x, float):
        if math.isnan(x):
            return "NA"
        return f"{x:.8f}"
    return str(x)


# ------------------------------------------------------------
# Read populations
# ------------------------------------------------------------
sample2pop, pop2samples = load_populations(pop_dir)


# ------------------------------------------------------------
# Read header/sample order
# Use syn_vcf as the canonical sample order
# ------------------------------------------------------------
samples = read_vcf_header_samples(syn_vcf)
sample_set = set(samples)

all_pop_samples = set(sample2pop.keys())
vcf_only = sorted(sample_set - all_pop_samples)
pop_only = sorted(all_pop_samples - sample_set)
matched = sorted(sample_set & all_pop_samples)

with open(f"{out_prefix}.sample_check.tsv", "w") as out:
    out.write("category\tsample\n")
    for s in matched:
        out.write(f"matched\t{s}\n")
    for s in vcf_only:
        out.write(f"vcf_only\t{s}\n")
    for s in pop_only:
        out.write(f"pop_only\t{s}\n")

calc_samples = samples


# ------------------------------------------------------------
# Data structures
# ------------------------------------------------------------
syn_sum = defaultdict(int)
syn_n = defaultdict(int)

nonsyn_sum = defaultdict(int)
nonsyn_n = defaultdict(int)

del_sum = defaultdict(int)
del_n = defaultdict(int)

syn_sites = set()
nonsyn_sites = set()
del_sites = set()


# ------------------------------------------------------------
# Core parser
# ------------------------------------------------------------
def process_vcf(vcf_path, target, site_set=None):
    """
    target:
      - "syn"
      - "nonsyn"
      - "deleterious"
    site_set:
      stores unique site keys in this file
    """
    seen_in_this_file = set()

    with open_text_auto(vcf_path) as f:
        for line in f:
            if line.startswith("##"):
                continue

            if line.startswith("#CHROM"):
                parts = line.rstrip("\n").split("\t")
                local_samples = parts[9:]
                if local_samples != samples:
                    raise ValueError(
                        f"Sample order mismatch between {vcf_path} and canonical sample list"
                    )
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue

            chrom, pos, _vid, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]

            # defensively skip multi-allelic or non-SNP lines
            if "," in alt:
                continue
            if len(ref) != 1 or len(alt) != 1:
                continue

            key = (chrom, pos, ref, alt)
            if key in seen_in_this_file:
                continue
            seen_in_this_file.add(key)

            if site_set is not None:
                site_set.add(key)

            sample_fields = parts[9:]
            for s, sf in zip(samples, sample_fields):
                dosage = parse_gt_to_dosage(sf)
                if dosage is None:
                    continue

                if target == "syn":
                    syn_sum[s] += dosage
                    syn_n[s] += 1
                elif target == "nonsyn":
                    nonsyn_sum[s] += dosage
                    nonsyn_n[s] += 1
                elif target == "deleterious":
                    del_sum[s] += dosage
                    del_n[s] += 1
                else:
                    raise ValueError(f"Unknown target: {target}")


# ------------------------------------------------------------
# Process files
# ------------------------------------------------------------
process_vcf(syn_vcf, target="syn", site_set=syn_sites)
process_vcf(nonsyn_vcf, target="nonsyn", site_set=nonsyn_sites)
process_vcf(del_vcf, target="deleterious", site_set=del_sites)


# ------------------------------------------------------------
# Overlap / sanity report
# ------------------------------------------------------------
with open(f"{out_prefix}.site_overlap.tsv", "w") as out:
    out.write("metric\tvalue\n")
    out.write(f"syn_sites\t{len(syn_sites)}\n")
    out.write(f"nonsyn_sites\t{len(nonsyn_sites)}\n")
    out.write(f"deleterious_sites\t{len(del_sites)}\n")
    out.write(f"del_intersect_nonsyn\t{len(del_sites & nonsyn_sites)}\n")
    out.write(f"del_not_in_nonsyn\t{len(del_sites - nonsyn_sites)}\n")
    out.write(f"nonsyn_not_in_del\t{len(nonsyn_sites - del_sites)}\n")


# ------------------------------------------------------------
# Individual results
# ------------------------------------------------------------
individual_rows = []

for s in calc_samples:
    pop = sample2pop.get(s, "NA")

    Ps = (syn_sum[s] / (2.0 * syn_n[s])) if syn_n[s] > 0 else float("nan")
    Pn = (nonsyn_sum[s] / (2.0 * nonsyn_n[s])) if nonsyn_n[s] > 0 else float("nan")
    Pd = (del_sum[s] / (2.0 * del_n[s])) if del_n[s] > 0 else float("nan")

    loadM = float("nan")
    if not math.isnan(Pn) and not math.isnan(Ps) and (Pn + Ps) > 0:
        loadM = Pn / (Pn + Ps)

    loadD = float("nan")
    if not math.isnan(Pd) and not math.isnan(Ps) and not math.isnan(Pn) and (Pn + Ps) > 0:
        loadD = Pd / (Pn + Ps)

    row = {
        "sample": s,
        "pop": pop,
        "syn_called_sites": syn_n[s],
        "nonsyn_called_sites": nonsyn_n[s],
        "deleterious_called_sites": del_n[s],
        "syn_derived_sum": syn_sum[s],
        "nonsyn_derived_sum": nonsyn_sum[s],
        "deleterious_derived_sum": del_sum[s],
        "Ps": Ps,
        "Pn": Pn,
        "Pd": Pd,
        "loadM": loadM,
        "loadD": loadD,
    }
    individual_rows.append(row)

with open(f"{out_prefix}.individual_load.tsv", "w") as out:
    header = [
        "sample", "pop",
        "syn_called_sites", "nonsyn_called_sites", "deleterious_called_sites",
        "syn_derived_sum", "nonsyn_derived_sum", "deleterious_derived_sum",
        "Ps", "Pn", "Pd", "loadM", "loadD"
    ]
    out.write("\t".join(header) + "\n")
    for r in individual_rows:
        out.write("\t".join(fmt(r[h]) for h in header) + "\n")


# ------------------------------------------------------------
# Population summary
# ------------------------------------------------------------
pop_to_rows = defaultdict(list)
for r in individual_rows:
    if r["pop"] != "NA":
        pop_to_rows[r["pop"]].append(r)

with open(f"{out_prefix}.population_load.tsv", "w") as out:
    header = [
        "pop", "n_samples",
        "mean_Ps", "sd_Ps",
        "mean_Pn", "sd_Pn",
        "mean_Pd", "sd_Pd",
        "mean_loadM", "sd_loadM",
        "mean_loadD", "sd_loadD"
    ]
    out.write("\t".join(header) + "\n")

    for pop in sorted(pop_to_rows.keys()):
        rows = pop_to_rows[pop]
        Ps_vals = [r["Ps"] for r in rows]
        Pn_vals = [r["Pn"] for r in rows]
        Pd_vals = [r["Pd"] for r in rows]
        loadM_vals = [r["loadM"] for r in rows]
        loadD_vals = [r["loadD"] for r in rows]

        out.write("\t".join([
            pop,
            str(len(rows)),
            fmt(mean_or_na(Ps_vals)),
            fmt(sd_or_na(Ps_vals)),
            fmt(mean_or_na(Pn_vals)),
            fmt(sd_or_na(Pn_vals)),
            fmt(mean_or_na(Pd_vals)),
            fmt(sd_or_na(Pd_vals)),
            fmt(mean_or_na(loadM_vals)),
            fmt(sd_or_na(loadM_vals)),
            fmt(mean_or_na(loadD_vals)),
            fmt(sd_or_na(loadD_vals)),
        ]) + "\n")

print("Done.")
print(f"Output: {out_prefix}.individual_load.tsv")
print(f"Output: {out_prefix}.population_load.tsv")
print(f"Output: {out_prefix}.site_overlap.tsv")
print(f"Output: {out_prefix}.sample_check.tsv")
