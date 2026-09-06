# BioMA RONA module

Complete parameter, method, and result interpretation is available in the
[English user guide](docs/USER_GUIDE.en.md#5-rona) and
[中文手册](docs/USER_GUIDE.zh-CN.md#5-rona).

The RONA module follows `RONA.revise.R` and is exposed as one command. It
calculates per-BIO RONA using linear allele-frequency/environment regressions,
weights loci by R-squared, calculates row-wise standard errors, and writes
population tables plus population plots and masked GAM maps.

## Configure

Copy `workflow.rona.example.ini` and edit the input paths, model/SSP/period
selection, and output directory. The `future_climate` directory may contain
scenario directories directly or below one additional directory level.

## Validate

```bash
bin/bioma rona workflow.rona.ini --dry-run
```

## Run

```bash
bin/bioma rona workflow.rona.ini
```

Outputs are written below the configured directory:

- `weighted/`: R-squared weighted RONA values for BIO1-BIO19
- `SE/`: standard errors
- `weighted_SE/`: formatted `RONA +/- SE` values
- `ensemble_mean/`: coordinate-level arithmetic means used for every plot
- `maps/`: one PDF and PNG per BIO x period x SSP
- `boxplots/`: one PDF and PNG per BIO
- `run_manifest.json`, `rona_compute.log`, `rona_plot.log`

The workflow accepts `--dry-run`, validates all selected scenarios before
calculation, and leaves the original source data unchanged.

The `unld_dir` input is a scientific dependency, not a plotting option. It
must contain `LD_BIO1.prune.in` through `LD_BIO19.prune.in`, normally produced
by PLINK from the adaptive-site genotype VCF. During a run, each list is read
by `rona_compute.R`; BioMA first requires at least one frequency-column match,
then writes its locus count, overlap count, and SHA-256 digest to the manifest
under `inputs.ld_pruning`. The environment table must contain `bio1`-`bio19`
and share at least three population IDs with the frequency table;
`expected_populations` optionally enforces their exact count. This table
supplies BIO predictor values and cannot replace genotype-based LD lists.

For plotting, every selected model must have a valid value at the coordinate.
The arithmetic mean is computed first and only that ensemble mean is passed to
the GAM interpolation. The default grid is 0.1 degrees with `grid_k = 15`.
