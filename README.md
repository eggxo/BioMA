# BioMA

[![CI](https://github.com/eggxo/BioMA/actions/workflows/ci.yml/badge.svg)](https://github.com/eggxo/BioMA/actions/workflows/ci.yml)
[![License: GPL-3.0-only](https://img.shields.io/badge/license-GPL--3.0--only-blue.svg)](LICENSE)

BioMA is an auditable command-line workflow for multidimensional climate-vulnerability analyses. It currently integrates seven modules while keeping every module independently runnable: Gradient Forest (GF), RONA, MAR, genetic loadM/loadD, MaxEnt niche modelling, WFmoments genetic-diversity loss, and multidimensional vulnerability synthesis.

BioMA is released under the GNU General Public License, version 3 only
(GPL-3.0-only). The full license text is included in `LICENSE`. The MaxEnt jar
and the two research R packages remain subject to their own distribution terms.

## Quick start

BioMA is validated on Linux. The shortest installation path uses the locked
Conda environment:

```bash
conda-lock install --name bioma conda-lock.yml
conda activate bioma
python -m pip install --no-deps .
bioma --version
bioma doctor
```

Run the server-independent seven-module smoke project before supplying real
data:

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bin/bioma project /tmp/bioma-demo/project.ini --dry-run
```

The demo validates configuration, input wiring, provenance generation, and
module contracts. It does not replace a scientific run with licensed MaxEnt
and research R packages; those dependencies are checked by `bioma doctor` and
documented in [INSTALL.md](INSTALL.md).

For a complete analysis, copy [workflow.project.example.ini](workflow.project.example.ini),
replace the paths in `[inputs]`, select the modules in `[modules]`, and run:

```bash
bin/bioma project project.ini --dry-run
bin/bioma project project.ini
```

The project command is resumable and writes configuration snapshots, input
hashes, module status, a TSV summary, and an HTML report. It never edits the
declared source VCFs, tables, rasters, or masks.

## Whole project: one command

Each scientific module has its own INI file. `workflow.project.example.ini`
connects them without hiding their scientific parameters. Edit the paths once,
validate the complete project, and run it:

```bash
cp workflow.project.example.ini project.ini
bin/bioma project project.ini --dry-run
bin/bioma project project.ini
```

The standard result layout is:

```text
results/
  00_project/       resolved configuration, checksums, state, and manifest
  01_gf/
  02_rona/
  03_mar/
  04_load/
  05_niche/
  06_wfmoment/
  07_vulnerability/
  report/index.html
```

The project writes `00_project/input_contract.tsv`, which records every
canonical user input and every compatibility file generated for a module. It
also writes `00_project/input_manifest.json`: each declared file is recorded
with a streaming SHA-256 digest, and each directory has a stable recursive
manifest (relative paths, sizes, and per-file digests). A shapefile input is
treated as one logical dataset, so same-stem `.dbf`, `.shx`, `.prj`, `.cpg`,
and other sidecar files are included automatically. The project signature uses
these content fingerprints, not paths alone; changing an input is rejected
before generated compatibility files are rewritten unless `--overwrite` is
given.
If the project provenance sidecar is missing while module results remain,
BioMA also refuses to adopt those results; use a new output directory or
explicitly pass `--overwrite` after confirming the inputs.

Standalone module commands apply the same content-level provenance contract:
their `run_manifest.json` includes `input_fingerprints` and an aggregate
`input_signature_sha256`. Directory inputs are recorded recursively, and
shapefile sidecars are included as part of the logical dataset. This makes a
standalone run auditable even when it is not launched through `bioma project`.

Shared species, climate-model, SSP, period, and seed values can be set once in
`[shared]`. The vulnerability configuration may use `auto` for GF, RONA, load,
and niche inputs; BioMA then connects the canonical upstream products for the
selected integration period and SSP. Missing upstream dependencies are caught
before calculations, and vulnerability is marked `blocked` if a required
upstream module fails in a continue-on-error run.

The vulnerability panel is species-agnostic by default. Set
`rona_variables` to a comma-separated list of RONA fields (or `auto` to use
all `*_RONA` columns), choose `rona_summary = mean`, `max`, or `weighted`, and
provide `rona_weights` for a weighted summary. Leave `group_order` and
`exclude_groups` blank to infer the group order and include all groups; set
them explicitly in a species profile when a particular ordering or exclusion
is scientifically justified.

In project mode the same five controls may be placed in `[integration]` of the
single `project.ini`; non-empty values override the vulnerability module
profile for that run. The adaptive-site ALT frequency is generated once from
the adaptive VCF and is shared by GF and RONA; a separate RONA `alt.frq` is a
legacy standalone-input option only.

Completed module results are reused only when their saved effective
configuration hash matches. Changed or incomplete outputs are rejected. To
rerun intentionally, use `--overwrite`; BioMA moves the old module and project
metadata directories under `_bioma_backups/` before starting. Use one or more
`--module` options to run an enabled subset:

```bash
bin/bioma project project.ini --module rona --module mar
```

The project command writes a machine-readable manifest, a TSV module summary,
and a static HTML report. A run containing failed or blocked modules returns a
nonzero command status so schedulers do not mistake partial completion for
success.

Installation requirements, including the two research R packages and the
MaxEnt jar requirement, are documented in `INSTALL.md`.

For the validated Linux platform, `conda-lock.yml` pins the resolved Conda and
pip artifacts (including checksums). Use `conda-lock install` for a reproducible
environment; `environment.yml` remains the editable source specification.

Before a long run, inspect the local toolchain with the read-only doctor:

```bash
bin/bioma doctor
bin/bioma doctor workflow.project.ini
bin/bioma doctor workflow.project.ini --json > doctor.json
```

With a project INI, checks needed by enabled modules are marked as required;
other scientific packages are reported as optional. Use `--strict` when a
complete installation is required. The doctor probes Python, R/Rscript, GDAL
utilities, Java, MaxEnt's `maxent.jar`, the `gradientForest` and `mar` R
packages, and the Python `wfmoments` stack without changing any input or result
file. Module-specific executable paths are read from the project/module INI;
options such as `--rscript` and `--maxent-jar` override them.

The complete per-module input inventory and the next-stage input consolidation
recommendations are in `INPUTS.md`.

The module dependency graph and the boundary between shared inputs and genomic
data streams are in `ARCHITECTURE.md`, with editable Graphviz source and PNG/SVG
renders included in the distribution.

For a simplified project setup, fill the canonical `[inputs]` section in
`workflow.project.example.ini`. BioMA keeps the three genomic data streams
separate: adaptive-site VCFs for GF/RONA, the Beagle-filtered whole-genome VCF
for MAR, and derived/SIFT-annotated VCFs for loadM/loadD. The climate raster
catalog, population environment table, and species mask are shared and wired
automatically.

## Public demo and tests

The repository includes a server-independent fixture generator. It creates all
seven module configurations, 19 BIO scenario placeholders, LD lists, and the
three genomic input streams under a disposable directory:

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bin/bioma project /tmp/bioma-demo/project.ini --dry-run
```

Use `--valid-rasters` in the locked Conda environment when a small readable
GeoTIFF set is needed. The automated suite runs one minimum contract test per
module and never requires the production server paths. R/MaxEnt numerical tests
are separate because `maxent.jar` and the research R packages have their own
licenses and installation channels.

The repository also includes `tests/data/gf_frequency_50/`, a fixed 50-site
GF fixture with 175 samples and 25 populations. It contains the input VCF,
sample table, expected population ALT-frequency output, and a fixture manifest
for a fast regression check:

```bash
python -m unittest tests.test_gf_frequency -v
```

Both fixtures are included in the source distribution and contain no private
server paths or production data. Additional species-specific test datasets can
be added under `tests/data/<name>/` with a manifest, license/permission note,
and a deterministic expected-output check.

## Complete GF workflow: one command

Copy `workflow.example.ini`, edit only the paths and the requested model/SSP/period lists, then run:

```bash
bin/bioma run workflow.ini
```

This single command runs population allele frequencies, climate preparation, GF training, local/forward/reverse offset for every requested scenario, ensemble maps for every period × SSP, and the forward-offset migration-distance figure for every period. Outputs are organized as `01_frequency` through `06_forward_distance_plots` under the configured output directory.

For sample grouping, the configuration accepts either one two-column `samples` table or the historical `sample_groups_dir` layout containing one non-hidden keep-list file per population. In the latter case BioMA creates and records the two-column sample table automatically.

The workflow requires a complete model × SSP × period grid and at least two climate models. It records a snapshot of all source inputs. Re-running the same command safely reuses completed stages; a changed input or configuration is rejected for the existing output directory so stale results cannot be mixed with new results. Use `--dry-run` to validate inputs and list all planned scenarios without calculating them:

```bash
bin/bioma run workflow.ini --dry-run
```

The individual commands below remain available for testing or rerunning a single module.

## GF population-frequency module

Input sample metadata must be a tab-separated file with one sample per row:

```text
sample_id	population_id
sample-A	population-1
sample-B	population-1
```

Run the module from the source directory:

```bash
bin/bioma gf-frequency \
  --vcf adaptive_sites.vcf \
  --samples samples.tsv \
  --outdir results/gf_frequency \
  --min-pop-samples 3 \
  --warn-pop-samples 5
```

It can also be installed as a user-level Linux command:

```bash
python3 -m pip install --user .
bioma --version
```

The reader detects gzip compression from file content rather than the filename extension. Inputs may therefore be plain VCF or genuinely gzip-compressed VCF, but every variant must be biallelic and contain a `GT` field.

Outputs:

- `population_alt_frequency.tsv`: GF-ready population × variant matrix.
- `population_frequency_qc.tsv`: sample size, called chromosomes, missingness, and QC status by population.
- `variant_manifest.tsv`: ordered variant identity and alleles.
- `warnings.tsv`: explicit scientific or data warnings.
- `run_manifest.json`: input checksums, parameters, versions, counts, and run status.

An existing output directory is never overwritten.

## Climate preparation module

The module discovers BIO1-BIO19 from filenames, validates raster grids, extracts present climate at population coordinates, and materializes masked present/future background tables. Source TIFF files are only read and are never copied or modified.

```bash
bin/bioma climate-prepare \
  --present-dir "climate/wc2.1_2.5m_bio present" \
  --future-root climate/select3 \
  --coordinates GF_ENV_19bio.txt \
  --current-mask current_distribution.shp \
  --future-mask future_search_domain.shp \
  --model BCC-CSM2-MR \
  --ssp 245 \
  --period 2061-2080 \
  --outdir results/climate_prepare
```

Repeat `--model`, `--ssp`, and `--period` to request a Cartesian set of scenarios. If none are supplied, all discoverable scenarios are included. `--future-mask` defaults to `--current-mask`. Use `--validate-only` to create manifests and the population environment table without materializing background tables.

Outputs include `population_environment.tsv`, `current_background.tsv.gz`, one table per scenario under `future_local/` and `future_search/`, `future_manifest.tsv`, `raster_manifest.tsv`, `warnings.tsv`, and `run_manifest.json`.

The current mask defines the present/background domain used by reverse offset. The future mask defines the candidate destination domain used by forward offset. Local offset always uses future climate at the current-mask cells. A global run is intentionally not an implicit default.

## Gradient Forest training module

Train a model directly from the standardized outputs of modules 1 and 2:

```bash
bin/bioma gf-train \
  --frequencies results/gf_frequency/population_alt_frequency.tsv \
  --environment results/climate_prepare/population_environment.tsv \
  --outdir results/gf_train
```

Use `--expected-sites` only when an independent site-count check is part of
your study design; otherwise BioMA infers the count from the frequency table.

Production defaults reproduce the historical 19-BIO setup: 500 trees, 1001 bins, correlation threshold 0.5, compact model storage, and `maxLevel = log2(0.368 * n_populations / 2)`. BioMA additionally fixes and records the random seed. Population rows are joined by `population_id`; row-number deletion is never used.

Outputs include the downstream-compatible `all_gfmod.data`, predictor importance, response out-of-bag performance, transformed training environments, population alignment, warnings, R session details, training log, and a checksum-bearing run manifest.

## Gradient Forest offset module

Compute all three offset directions for one prepared future scenario:

```bash
bin/bioma gf-offset \
  --model results/gf_train \
  --climate results/climate_prepare \
  --scenario 2061-2080-ssp245-MPI-ESM1-2-HR \
  --forward-radii 100,250,500,1000,inf \
  --outdir results/gf_offset_mpi_ssp245
```

Local offset compares present and future climate at the same current-mask cell. Forward offset finds the most similar future-search cell for every current cell within each requested migration radius. Reverse offset finds the most similar current-background cell for every future-search cell. Distances are Euclidean in the 19-dimensional GF-transformed space.

The implementation uses an exact expanding-k nearest-neighbour search. A match is finalized only after all possible equal-GF-distance candidates have been considered; ties are resolved by geographic distance and then input row. A deterministic brute-force sample check is run before results are accepted.

Outputs include compressed local, forward, and reverse tables, transformed climate tables, scenario summaries, warnings, runtime/session information, input and output checksums, and algorithm-verification details.

`all_offsets.tsv.gz` is the canonical final table. Its grain is one coordinate × one forward radius, with local, forward, and reverse offsets plus destination coordinates, geographic distances, bearings, and status fields in one row. Separate component files remain available for auditing.

Create maps and relationship plots from an ensemble of model tables:

```bash
bin/bioma offset-plot \
  --input results/offset_bcc \
  --input results/offset_mpi \
  --input results/offset_ukesm \
  --models BCC-CSM2-MR,MPI-ESM1-2-HR,UKESM1-0-LL \
  --period 2061-2080 \
  --ssp ssp245 \
  --radius unlimited \
  --populations results/climate_prepare/population_environment.tsv \
  --boundary current_distribution.shp \
  --outdir results/offset_plot
```

The plotting module requires at least two distinct climate models. It first validates that every model has the same coordinate grid, then computes an unweighted arithmetic mean for local, forward, and reverse offset at each coordinate. A mean is reported only when that metric is valid in every selected model. Duplicate model-coordinate rows are rejected so no model can receive extra weight. All normalization and plotting happen only after this averaging step.

Outputs include the canonical `ensemble_mean_offsets.tsv.gz`, an RGB composite (mean local = red, mean forward = green, mean reverse = blue), three normalized ensemble-mean maps, three pairwise relationship plots, and the exact `plot_data.tsv.gz` used by every panel. PNG and PDF figures are both written.

A standalone plotting environment specification is provided as `environment-plot.yml`. It can be created with:

```bash
source ~/.bashrc
mamba env create -p ./envs/bioma-plot -f environment-plot.yml
mamba activate ./envs/bioma-plot
bin/bioma offset-plot --help
```

When BioMA is run from its source directory it automatically discovers this project-local plotting environment. `--rscript` remains available for users who keep the required R packages in another environment.

### Forward offset by migration-distance limit

Create the migration-distance response figure from the SSPs that were actually calculated:

```bash
bin/bioma forward-distance-plot \
  --input results/offset_bcc_ssp245 \
  --input results/offset_mpi_ssp245 \
  --input results/offset_ukesm_ssp245 \
  --input results/offset_bcc_ssp585 \
  --input results/offset_mpi_ssp585 \
  --input results/offset_ukesm_ssp585 \
  --models BCC-CSM2-MR,MPI-ESM1-2-HR,UKESM1-0-LL \
  --period 2081-2100 \
  --ssp 245 \
  --ssp 585 \
  --minimum-models 3 \
  --outdir results/forward_distance_plot
```

If `--ssp` is omitted, all SSPs present in the supplied tables are plotted. For each SSP, radius, and coordinate, BioMA first computes the unweighted climate-model mean and requires every selected model to be valid. It then plots the median of those coordinate-level ensemble means; the ribbon is their interquartile range. Distance limits are displayed as ordered categories with `Unlimited` last; no artificial axis-break mark is drawn.

Outputs include `forward_offset_ensemble_grid.tsv.gz`, `forward_offset_by_distance.tsv`, publication-ready PNG/PDF figures, and a checksum-bearing run manifest.

## RONA workflow: one command

The RONA module follows `RONA.revise.R`: it reads the population allele-frequency
table and BIO1-BIO19 environment, performs the per-locus regressions, computes
R-squared weighted RONA and standard errors, then writes population plots and
masked GAM maps.

Copy `workflow.rona.example.ini`, edit the paths and requested scenario lists,
then validate and run:

```bash
bin/bioma rona workflow.rona.ini --dry-run
bin/bioma rona workflow.rona.ini
```

Standalone RONA configurations currently require `alt_frequency`, `unld_dir`,
`environment`, and `future_climate`; `mask` is optional for maps. In the
project-level workflow, `alt_frequency` is generated automatically from the
adaptive-site VCF and the same frequency table is passed to GF and RONA. The
`unld_dir` files remain a legacy/implementation input until an explicit
genotype-based LD-generation rule is defined. The future-climate root may
contain scenario directories directly or below a nested directory. Outputs are
written under the configured directory as `weighted/`, `SE/`, `weighted_SE/`,
`ensemble_mean/`, `maps/`, `boxplots/`, and the compute/plot logs plus
`run_manifest.json`.

For every plot, all selected climate models must be present and valid at a
coordinate. BioMA computes the arithmetic model mean first, writes it to
`ensemble_mean/`, and uses only that mean for GAM interpolation and plotting.
The default interpolation settings are `grid_step = 0.1` degrees and
`grid_k = 15`.

The R environment must provide `terra`, `plotrix`, `raster`, and the plotting
dependencies used by the module. Put its `Rscript` on `PATH` or set the
`rscript` field explicitly. Source files are never modified.

## Regression test

The bundled fixture contains 50 adaptive variants, 175 samples, and 25
populations. Its expected frequencies are derived from a small deterministic
reference matrix and are independent of any production site count.

```bash
python3 -m unittest discover -s tests -v
# Requires R, data.table, ggplot2, and raster
Rscript tests/test_vulnerability_r.R
```

## WFmoments 2-D deme workflow

The WFmoments module projects the loss of nucleotide diversity after habitat
loss on a two-dimensional deme grid. It starts from the observed species-level
pi value, calibrates the neutral moment model, removes demes by an edge order
or by reproducible random permutations, and evaluates the immediate, 3-
generation, and 5-generation projections. The GDAR short-term curve is kept
as the independent area-only reference `A_remaining ^ z_gdar`; it is not
substituted with the WF immediate curve.

Copy `workflow.wfmoment.example.ini`, edit the paths and analysis settings,
then validate or run the module with one command:

```bash
bin/bioma wfmoment workflow.wfmoment.ini --dry-run
bin/bioma wfmoment workflow.wfmoment.ini
```

Required inputs are a binary current-habitat raster and a species row in the
`pi_file` table (`species`, `pi_obs`). Supply either `area_file` (with
`Scenario` and `Area_km2`) or `future_masks_json` for scenario points. A
`param_file` with `z_gdar` and a `structure_file` with `fst_global_est` are
optional when `z_gdar` and migration are supplied explicitly. Relative paths
are resolved relative to the configuration file.

Important parameters include:

- `nx`, `ny`, `threshold`, and `min_valid`: the 2-D deme grid and occupancy rule.
- `loss_mode = edge` or `random`; `direction` controls the edge order
  (`east_to_west`, `west_to_east`, `north_to_south`, or `south_to_north`).
- `migration` and `theta`: numeric values or `auto`; auto migration uses the
  supplied global FST and the configured `migration_grid`, while auto theta
  uses the configured low-theta probe to match the observed pi scale and
  records the resulting model baseline for audit.
- `replicates` and `seed`: random-loss replicate count and reproducibility.
- `time3`, `time5`, and `mu`: projection and demographic settings.
  `midterm_generations = auto` is recorded from pi and mu for planning/audit;
  the current figure intentionally plots the configured 3- and 5-generation
  trajectories and does not claim to be the paper's Ne/2 mid-term run.
- `include_equilibrium`, `plot_equilibrium`, and `plot_raw`: optional
  equilibrium series and raw-point overlays.

Outputs are written below `output_dir`:

- `deme_table.tsv`, `deme_order.tsv`, and `metadata.json/.tsv` document the
  grid, removal order, calibrated parameters, and input checksums.
- `curve_replicates.tsv` contains every replicate and removal step;
  `curve_summary.tsv` contains the arithmetic replicate means used for plotting;
  `scenario_points.tsv` contains interpolated values for each area scenario.
- `Figure_wfmoment_2D_deme.png` and `.pdf` are the publication-style figure.
  `Figure_wfmoment_2D_deme_curve_data.tsv` preserves the finite raw curve
  values, whereas `Figure_wfmoment_2D_deme_curve_display.tsv` contains the
  interpolated/monotone display line (including explicit NA breaks). The raw
  `curve_summary.tsv` and `curve_replicates.tsv` files are never rewritten by
  the plotting step. For the default edge display, the continuous WF line
  uses states with at least five extant demes, then adds a display-only,
  monotone terminal bridge to `(100%, 0%)` when the raw table contains a
  complete-loss (`occupied_demes=0`) endpoint. The bridge is shaped as
  `1 - u^1.5`, so it reaches the complete-loss boundary smoothly without
  joining the one-deme or zero-deme sentinel as an observed state. The
  synthetic bridge is marked in the display table and its settings are recorded
  in `Figure_wfmoment_2D_deme_diagnostics.tsv`. If the input is a partial curve
  without a complete-loss endpoint, `--tail-mode auto` leaves the line at the
  observed range; it does not invent a zero at a zoomed `x-max`.

  Use `--tail-mode none` for a strict truncation at the last retained deme;
  `--tail-mode smooth` enables the bridge explicitly, while `--tail-mode zero`
  and `--tail-mode linear` remain legacy display-only choices. `--tail-end`
  changes the bridge endpoint (the default is `x-max`) and `--tail-exponent`
  changes its curvature power. `--min-display-demes` and
  `--tail-cutoff-demes` can make the retained connected core more conservative.
  `--include-zero-deme-endpoint` is an intentional QA override and restores the
  raw zero-deme row to the display. `Figure_wfmoment_2D_deme_point_data.tsv`
  retains all scenario points, including points excluded from the default
  period-matched display.

The default display matches the reference convention: GDAR is blue, WF 3
generations is orange dashed, and WF 5 generations is green; WF points are
shown for their matching periods (`2061-2080` for 3 generations and
`2081-2100` for 5 generations). In a connected edge-contraction trajectory,
WF 3 generations is expected to remain above WF 5 generations. Under
`loss_mode=random`, the species-wide pi statistic also contains between-deme
divergence, so fragmentation can legitimately produce local WF 5 > WF 3
values; those values remain in the raw output and are not globally reordered.
Use the plotting script's `--all-wf-points` option when all WF scenario points
are needed. A zero-occupied-deme endpoint is retained in the raw tables but is
not treated as an observed continuous state by the default display.

The production smoke test using a real species input completed with the edge
configuration (`20 x 20`, threshold `0.25`, one replicate) and
produced PNG/PDF figures plus all audit tables. A second random-loss test
(`10 x 10`, two replicates) also completed successfully. The small tests are
workflow checks; use production grid/replicate settings for final estimates.

## MAR workflow: one command

The MAR module wraps `MAR.all.R` through `bin/bioma mar`. `maxsnps = auto`
counts non-header VCF records and passes that value to `MARPIPELINE`, avoiding a
hard-coded SNP limit. The direction scheme defaults to `random` and can be
changed to `eastwest`, `westeast`, `northsouth`, or `southnorth`.

```bash
bin/bioma mar workflow.mar.ini --dry-run
bin/bioma mar workflow.mar.ini
```
