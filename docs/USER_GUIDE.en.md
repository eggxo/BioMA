# BioMA User and Methods Guide (English)

This manual applies to BioMA 0.8.x. It provides a consistent description of
project mode and all seven scientific modules: purpose, computational
principle, inputs, every configurable field, command examples, principal
outputs, and interpretation limits. The configuration parsers and pinned
dependencies define what BioMA actually computes. The publications listed at
the end provide methodological context; they do not constitute independent
biological validation of BioMA results.

[中文手册](USER_GUIDE.zh-CN.md) | [Installation](../INSTALL.md) |
[Input inventory](../INPUTS.md) | [Architecture](../ARCHITECTURE.md)

## 1. Interpretive boundaries

- GF, RONA, and MaxEnt extrapolate relationships fitted in current geographic
  samples to future environments. A high value indicates greater modeled
  change or mismatch, not extinction probability, percentage fitness loss, or
  a causal effect.
- loadM/loadD are derived-allele ratio proxies defined by the bundled BioMA
  script. They are not directly measured fitness and are not classical
  selection-coefficient-weighted genetic load.
- MAR describes a static diversity-area power law. WFmoments describes the
  dynamics of neutral nucleotide diversity after habitat loss. They answer
  different questions.
- The vulnerability panel rescales five continuous metrics to one to three
  squares relative to groups in the current project. It is a ranking and
  visualization tool, not a new risk index with an absolute cross-species
  interpretation.
- Raw scales from different modules are not directly comparable. GCM
  aggregation is an unweighted arithmetic mean or the configured median, not
  a probability-weighted climate-model estimate.

## 2. Install, diagnose, and run the demo

```bash
conda-lock install --mamba --name bioma conda-lock.yml
mamba activate bioma
bin/install-bioma-r-deps.sh
python -m pip install --no-deps .
bioma --version
bioma doctor
```

After configuring a project, run `bioma doctor project.ini --strict` for a
strict check of enabled modules. Without a configuration, `--strict` also
counts unconfigured optional components as failures.

Generate the server-independent demo and validate all seven module contracts:

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bioma project /tmp/bioma-demo/project.ini --dry-run
```

`--dry-run` validates inputs, parameters, scenarios, and dependencies but does
not perform the scientific calculations. Start a real run with:

```bash
bioma project project.ini
```

Relative paths are resolved from the directory containing their INI file.
Use `true` or `false` for Boolean settings. Source VCFs, tables, rasters, and
masks are read-only; manifests record content-level SHA-256 fingerprints.

## 3. Project mode: connect every module once

Copy `workflow.project.example.ini`, fill shared inputs in `[inputs]`, and
enable modules in `[modules]`:

```bash
cp workflow.project.example.ini project.ini
bioma doctor project.ini
bioma project project.ini --dry-run
bioma project project.ini
```

Command-line options:

| Option | Meaning |
| --- | --- |
| `config` | Project INI path. |
| `--dry-run` | Perform preflight and write effective configuration snapshots only. |
| `--module NAME` | Run only a named enabled module; repeat as needed. Names are `gf`, `rona`, `mar`, `load`, `niche`, `wfmoment`, and `vulnerability`. |
| `--overwrite` | When configuration or input content changed, move old module and project metadata under `_bioma_backups/` and rerun. |
| `--no-resume` | Ignore `[project] resume` for this invocation. |

### 3.1 `[project]`

| Field | Default | Meaning |
| --- | --- | --- |
| `name` | INI stem | Project name stored in reports and manifests. |
| `output_dir` | required | Project result root. |
| `resume` | `true` | Reuse completed results only when the effective configuration and input-content signatures match exactly. |
| `stop_on_error` | `true` | Stop after the first failed module. `false` lets independent modules continue and marks dependent modules as blocked. |
| `report` | `true` | Write `report/index.html`. |

### 3.2 `[modules]`

The values of `gf`, `rona`, `mar`, `load`, `niche`, `wfmoment`, and
`vulnerability` are module INI paths. Disable a module with `false`, `no`,
`off`, `disabled`, `none`, `0`, or an empty value. Execution follows this
fixed order. Vulnerability depends on GF, RONA, load, and niche products.

### 3.3 Canonical `[inputs]`

| Field | Purpose and contract |
| --- | --- |
| `adaptive_vcf` | Biallelic adaptive-locus VCF shared by GF/RONA, with `GT`. GF generates population ALT frequencies. |
| `adaptive_frequency` | Legacy compatibility input. Do not set it in new projects; use only for an existing frequency table when GF is not run. |
| `whole_genome_vcf` | Whole-genome, imputed or no-missing VCF for MAR. Do not substitute the adaptive or load VCF. |
| `load_vcf_dir` | Directory containing the four derived-polarized, SIFT-annotated load VCFs. |
| `population_samples` | One row per individual; at least `sample_id,population_id`, optionally `lon,lat`. It can generate GF sample metadata, MAR coordinates, load keep lists, and MaxEnt occurrences. |
| `population_environment` | One row per population with ID, coordinates, optional group, and `bio1`-`bio19`. Shared by GF and RONA. |
| `rona_unld_dir` | The 19 files `LD_BIO1.prune.in` through `LD_BIO19.prune.in`. These are scientific locus-selection inputs to RONA. |
| `climate_current` | Current BIO TIFF directory shared by GF and niche. |
| `climate_future` | Future scenario root shared by GF, RONA, and niche; each module still enforces its naming contract. |
| `species_mask` | Species-range shapefile. Same-stem `.shx/.dbf/.prj/.cpg` and other sidecars are included in the content fingerprint. |
| `population_dir` | Optional load keep-list directory; generated from `population_samples` when blank. |
| `load_predictors` | Optional current population predictor CSV; a compatible table is generated from canonical population data when blank. |
| `load_future_dir` | Directory of future BIO CSVs for load. Version 0.8.x does not yet derive these directly from TIFFs. |
| `mar_lonlat` | Optional per-sample `ID,LONGITUDE,LATITUDE`; generated from canonical individual/population data when blank. |
| `occurrence_csv` | Optional `species,lon,lat`; deduplicated occurrences are generated from canonical coordinates when blank. |
| `maxent_jar` | Legally obtained `maxent.jar`, not distributed by BioMA. An installed `dismo` copy may also supply it. |
| `wfmoment_current_raster` | Current binary habitat/suitability raster, preferably in an equal-area CRS. |
| `pi_file` | At least `species,pi_obs`. |
| `structure_file` | At least `species,fst_global_est`; required for `migration=auto`. |
| `species_parameters` | At least `species,z_gdar`; it may hold other species parameters. |
| `area_file` | `Scenario,Area_km2`, including `Scenario=current`; mutually exclusive with `future_masks_json`. |
| `future_masks_json` | JSON object mapping scenario names to future binary masks; mutually exclusive with `area_file`. |

### 3.4 `[shared]`

| Field | Default | Meaning |
| --- | --- | --- |
| `species` | empty | Overrides the species row key in modules such as WFmoments. |
| `models` | empty | Comma-separated GCMs overriding GF/RONA; an empty value preserves module settings or discovery. |
| `ssps` | empty | Comma-separated SSPs, written as `245` or `ssp245`. |
| `periods` | empty | Comma-separated periods such as `2061-2080,2081-2100`. |
| `seed` | empty | A nonempty value overrides modules that expose a random seed. |
| `rscript` | empty | A nonempty value overrides module Rscript settings. Leave blank when modules need different environments. |
| `compute_python` | empty | A nonempty value overrides the WFmoments Python interpreter. |

### 3.5 `[integration]`

| Field | Default | Meaning |
| --- | --- | --- |
| `period` | conditionally required | Period used for automatic vulnerability wiring. It can be inferred when `[shared] periods` has one value. |
| `ssp` | conditionally required | SSP used for automatic vulnerability wiring. |
| `ensemble_method` | `mean` | `mean` or `median`, used in the auto-wired niche product name. |
| `rona_variables` | empty | Vulnerability override; `auto` discovers all `*_RONA` columns. |
| `rona_summary` | empty | Vulnerability override: `mean`, `max`, or `weighted`. |
| `rona_weights` | empty | Comma-separated `field:weight` pairs for weighted RONA. |
| `group_order` | empty | Explicit group order; empty preserves first-seen order. |
| `exclude_groups` | empty | Groups to exclude. BioMA has no species-specific default exclusions. |

### 3.6 `[doctor]`

`python`, `rscript`, `gdalinfo`, `gdallocationinfo`, `ogrinfo`, `java`, and
`maxent_jar` are optional runtime path overrides. Empty values are discovered
from PATH, module INIs, and installed `dismo`. `bioma doctor project.ini
--json` emits a machine-readable report. `--strict` makes missing optional
components fail the check as well.

## 4. Gradient Forest (GF)

### 4.1 Purpose, principle, and metrics

GF first calculates ALT allele frequencies at adaptive loci for every
population. It then fits `gradientForest` with locus frequencies as responses
and environmental variables as predictors. Important random-forest splits are
accumulated along each environmental gradient to produce nonlinear predictor
transformations. BioMA calculates Euclidean distances in this transformed
space:

- `local offset` is the transformed distance between current and future
  climate at the same cell. It represents the multilocus compositional change
  required to remain in place.
- `forward offset` finds the future destination with minimum transformed
  distance from each current cell, subject to a migration radius. A high value
  means migration still finds no close future environment.
- `reverse offset` finds the closest current source for each future cell. A
  high value means that future habitat lacks a similar current genetic source.

These values are model-space distances, not allele-frequency percentages,
migration probabilities, or fitness. At each coordinate, BioMA requires a
valid value from every selected GCM before taking an unweighted mean.

### 4.2 Inputs and analysis fields

| Field | Default | Meaning |
| --- | --- | --- |
| `[inputs] vcf` | required | Biallelic adaptive-site VCF/VCF.GZ with `GT`. ALT is the counted allele; BioMA does not infer derived state. |
| `[inputs] samples` | exactly one of two | TSV with `sample_id,population_id`. |
| `[inputs] sample_groups_dir` | exactly one of two | One file per population, one sample ID per line; filename is the population name. |
| `[inputs] coordinates` | required | Whitespace-delimited table accepting `ID`/`population_id`, `lon`/`longitude`, `lat`/`latitude`, optional `pop`/`group`, and optional supplied BIO values. |
| `[inputs] present_climate` | required | Current BIO1-BIO19 TIFF directory. |
| `[inputs] future_climate` | required | Complete period x SSP x GCM scenario grid. |
| `[inputs] current_mask` | required | Current distribution/search mask. |
| `[inputs] future_mask` | `current_mask` | Forward/reverse search mask; a larger mask permits novel future habitat. |
| `[analysis] output_dir` | required | Module result directory. |
| `[analysis] models` | discovered | Selected comma-separated GCMs. |
| `[analysis] ssps` | discovered | Selected comma-separated SSPs. |
| `[analysis] periods` | discovered | Selected comma-separated periods. |
| `[analysis] expected_sites` | empty | Optional independent locus-count QC. Empty infers the count from the VCF. Do not hard-code a fixture count. |

### 4.3 Every `[parameters]` field

| Field | Code default | Meaning and constraint |
| --- | --- | --- |
| `predictors` | `bio1,...,bio19` | GF environmental columns. Each must exist in population and raster tables and vary among populations. |
| `min_population_samples` | `3` | Minimum samples per population; fewer is an error. |
| `warn_population_samples` | `5` | Low-sample warning threshold; at least the minimum threshold. |
| `ntree` | `500` | Trees per response forest, at least 1. Larger values are generally more stable and slower. |
| `nbin` | `1001` | Bins used to accumulate split importance along gradients, at least 2. |
| `corr_threshold` | `0.5` | `gradientForest` response-correlation threshold in [0,1]. It is not the MaxEnt BIO collinearity threshold. |
| `max_level` | automatic | Maximum tree level. Automatic is `log2(0.368*n_populations/2)`; an explicit value must be positive. |
| `seed` | `1` | GF random seed. |
| `forward_radii` | `100,250,500,1000,inf` | Maximum great-circle search distances in km; `inf`/`unlimited` removes the limit. |
| `initial_knn_k` | `64` | Initial transformed-space nearest-neighbor candidate count. The algorithm expands it until the exact answer is resolved. |
| `batch_size` | `10000` | Offset query batch ceiling. It affects memory and speed, not the intended result. |
| `verify_sample` | `10` | Query cells checked against brute-force search; `0` disables this QA. |
| `tie_tolerance` | `1e-12` | Transformed-distance tie tolerance. Ties use geographic distance and then stable row order. |
| `map_radius` | `unlimited` | Forward radius shown in the combined map; it must match a `forward_radii` output label. |
| `minimum_models` | `2` (template `3`) | Minimum distinct GCMs for a complete grid and ensemble, at least 2. |
| `supplied_bio_tolerance` | `1e-6` | Allowed difference between supplied population BIO values and values re-extracted from rasters. Exceedance warns; raster values remain authoritative. |
| `rscript` | discovered | Rscript containing `gradientForest`, `FNN`, `data.table`, and plotting packages. |
| `gdalinfo` | `gdalinfo` | Raster geometry checker. |
| `gdallocationinfo` | `gdallocationinfo` | Population-coordinate raster extractor. |
| `ogrinfo` | `ogrinfo` | Mask validator. |

### 4.4 Run and outputs

```bash
bioma run workflow.example.ini --dry-run
bioma run workflow.example.ini
```

`01_frequency/population_alt_frequency.tsv` is shared with RONA.
`03_model/all_gfmod.data` is the model; `predictor_importance.tsv` and
`response_performance.tsv` report predictor importance and locus OOB R-squared.
`04_offsets/*/all_offsets.tsv.gz` contains all offsets, matched coordinates,
distances, and bearings. `05_offset_plots/*/ensemble_mean_offsets.tsv.gz` is
the default vulnerability input. `06_forward_distance_plots/` reports
migration-radius sensitivity. Inspect positive-OOB-R-squared coverage, NoData,
`no_candidate`, and full GCM coverage before interpreting maps.

## 5. RONA

### 5.1 Purpose, principle, and metric

RONA (risk of non-adaptedness) fits a univariate current-population regression
`p = alpha + beta*E` for every BIO and locus. Given current frequency `p_obs`
and future environment `E_future`, BioMA calculates
`abs(alpha + beta*E_future - p_obs)`, then takes an R-squared-weighted mean
over loci. High RONA means a larger modeled mean allele-frequency shift is
required; low RONA means a smaller required shift. It does not model selection
strength, gene flow, drift, linkage, or the realized rate of change.

Each BIO uses its own `LD_BION.prune.in` and intersects those IDs with the
frequency columns. Version 0.8.x does not run PLINK automatically, and an
environment table cannot replace these lists. `SE` is the ordinary standard
error `sd/sqrt(n)` of finite locus-level absolute differences. It is not
between-GCM uncertainty or a regression prediction interval.

### 5.2 Every configuration field

| Field | Default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] alt_frequency` | required | Population x adaptive-locus frequency table. The first column may be `Pop`; project mode uses GF ALT frequencies. |
| `[inputs] unld_dir` | required | Nineteen nonempty `LD_BIO1.prune.in` through `LD_BIO19.prune.in`. |
| `[inputs] environment` | required | `ID,pop,lon,lat,bio1,...,bio19`; `ID` matches frequency populations. |
| `[inputs] future_climate` | required | Directories named `<period>-ssp<code>-<model>`, each with exactly one `bioN.cut.tif`. |
| `[inputs] mask` | empty | Optional plotting mask; it does not change population RONA. |
| `[analysis] output_dir` | required | Result directory. |
| `[analysis] models` | empty | GCM filter; empty retains all discovered models. |
| `[analysis] ssps` | empty | SSP filter, with or without `ssp`. |
| `[analysis] periods` | empty | Period filter. |
| `[analysis] expected_populations` | empty | Optional expected overlap count. Dry and full runs require the actual intersection of environment `ID` and the frequency-table first column to equal this positive integer. |
| `[parameters] rscript` | `Rscript` | Interpreter with `data.table`, `terra`, `plotrix`, `raster`, `mgcv`, and plotting dependencies. |
| `[parameters] interpolate` | `true` | Produce GAM spatial surfaces; it does not change population tables. |
| `[parameters] grid_step` | `0.1` | Longitude/latitude interpolation step, greater than 0. A finer grid increases size and time, not observation density. |
| `[parameters] grid_k` | `15` | `mgcv::gam` 2-D smooth basis dimension, at least 3. Too high may overfit; too low may oversmooth. |

```bash
bioma rona workflow.rona.example.ini --dry-run
bioma rona workflow.rona.example.ini
```

`weighted/` stores per-GCM R-squared-weighted RONA; `SE/` stores locus standard
errors; `weighted_SE/` stores display strings; `ensemble_mean/` stores
coordinate-wise arithmetic means only where every selected GCM is valid;
`maps/` and `boxplots/` hold figures. Report valid loci per BIO, regression
R-squared distributions, the LD rule, and GCM count. Do not interpret RONA as
survival probability.
Before execution, BioMA requires all `bio1`-`bio19` environment columns, at
least three populations shared by the two tables, and at least one frequency
locus matching each LD list.

## 6. MAR

### 6.1 Purpose, principle, and metrics

The mutations-area relationship fits genetic diversity after spatial sampling
or simulated cell extinction to `S = c*A^z`. BioMA uses pinned `mar 0.2.0` to
construct the spatial genomics object and extinction trajectories, then fits
`M` (sites at which ALT is observed), `E` (endemic sites found only in the
retained area), `thetaW` (Watterson's theta), and `thetaPi` (nucleotide
diversity estimate). Normalized remaining diversity is
`D_remaining = A_remaining^z`. For `0<A_remaining<1`, a larger positive `z`
implies more predicted diversity loss for the same area loss.

MAR is a static relationship inferred from current spatial samples and is
primarily an immediate/short-term proxy after area reduction. It does not
model drift over future generations. Sampling density, spatial coverage, VCF
ascertainment, and extinction scheme all affect `z`.

### 6.2 Every configuration field

| Field | Default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] vcf` | required | Biallelic whole-genome VCF, with sample order matching coordinates; imputed/no-missing data are recommended. |
| `[inputs] lonlat` | required | Exactly unique `ID,LON/LONGITUDE,LAT/LATITUDE`, one row per sample in genotype order. |
| `[inputs] scenario_file` | empty | Optional table. Columns `scenario,A_remaining,geom_job_id` add points to official curves. |
| `[analysis] output_dir` | required | MAR work and result directory. |
| `[analysis] name` | `bioma_mar` | Base name for package objects. |
| `[analysis] geom_id` | `7` | Selects `geom_job_id` from `scenario_file`; it does not affect fitting. |
| `[parameters] maxsnps` | `auto` | `auto/all/empty` counts all VCF records. A positive integer is a cap; `mar` randomly downsamples above it using `randseed`. |
| `[parameters] scheme` | `random` | `random`, `inwards`, `outwards`, `northsouth`, or `southnorth`. Centered schemes use the most densely sampled cell; pole schemes model directional loss. |
| `[parameters] nrep` | `10` | Spatial sampling/extinction replicates, at least 1. |
| `[parameters] xfrac` | `0.01` | Range/cell fraction advanced per step, in `(0,1]`; the package processes at least one cell when fewer than 100 cells exist. |
| `[parameters] quorum` | `true` | Ask MAR sampling windows to contain samples where possible; it does not alter extinction cell deletion. |
| `[parameters] randseed` | `123` | Seed for locus subsampling and spatial replicates. |
| `[parameters] marsteps` | `data,gm,sfs,mar,ext,plot` | `mar` package stages. BioMA official curves read `extdflist`, so normal full runs must retain `data,gm,ext`; skip prerequisites only with known compatible RDA objects. |
| `[parameters] rscript` | `Rscript` | Interpreter with `mar 0.2.0`, `SeqArray`, `sars`, and plotting dependencies. |

```bash
bioma mar workflow.mar.example.ini --dry-run
bioma mar workflow.mar.example.ini
```

`MAR_official_fit_params.tsv` reports `c` and `z` for every metric;
`MAR_display_curves.tsv` holds normalized curves;
`MAR_scenario_points.tsv` holds optional scenario points; and
`MAR_habitat_loss_curves.png/.pdf` is the primary figure. RDA/GDS objects,
logs, and the manifest support audit. Report VCF filtering, `scheme`,
replicates, and the definition of area.

## 7. loadM/loadD

### 7.1 Purpose, principle, and metrics

ALT must already represent the derived allele. For each diploid individual,
BioMA maps genotype to 0/1/2 derived dosage and calculates mean derived allele
frequency per callable site class: `Ps` (synonymous), `Pn` (nonsynonymous),
and `Pd` (deleterious). The current implementation defines:

```text
loadM = Pn / (Pn + Ps)
loadD = Pd / (Pn + Ps)
```

loadM is nonsynonymous derived burden relative to the nonsynonymous plus
synonymous background. loadD is the SIFT-deleterious burden using the same
denominator. These names denote BioMA/source-script proxies, not realized
fitness load. Incorrect polarization, unequal callable categories, and SIFT
classification error affect them.

The module then uses current `bio1`-`bio19`. Standardization and SVD are fitted
on training data only, followed by a `ranger` random forest for relaxed loadM
and relaxed loadD. Each outer split is 80/20; repeated v-fold CV within its
training set selects `mtry`, `min_n`, and `trees` by RMSE. One split is chosen
for future prediction by `cv_rsq_mean` first and `test_rsq` second. Future
values are statistical extrapolations, not evolutionary simulations of future
mutation accumulation.

### 7.2 Every configuration field

| Field | Code default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] vcf_dir` | required | Contains `strict.synonymous.vcf`, `strict.nonsynonymous.vcf`, `strict.deleterious_plain.vcf`, and `relaxed.deleterious_with_warning.vcf`; gzip content is accepted. |
| `[inputs] population_dir` | required | One keep-list per population, one sample ID per line. |
| `[inputs] predictors` | required | Current CSV with at least `pop,bio1,...,bio19`; mapping also uses `lon,lat` and optional `cluster`. |
| `[inputs] future_dir` | required | Future predictor CSV directory; every selected file has `bio1`-`bio19`. |
| `[inputs] mask` | empty | Optional plotting mask; it does not affect load calculation or RF fitting. |
| `[analysis] output_dir` | required | Result directory. |
| `[parameters] n_splits` | `20` (template `1000`) | Independent 80/20 outer splits, at least 1. Cost grows approximately linearly. |
| `[parameters] cv_v` | `5` | CV folds within each training set, at least 2. |
| `[parameters] cv_repeats` | `2` (template `10`) | Repetitions of v-fold CV, at least 1. |
| `[parameters] grid_size` | `30` (template `50`) | Latin-hypercube parameter combinations per split, at least 1. |
| `[parameters] workers` | `4` (template `16`) | Parallel workers, at least 1 and no more than scheduler allocation. |
| `[parameters] scheme` | `svd` | Version 0.8.x accepts only `svd`. |
| `[parameters] seed` | `1234` | Seed for outer splits, CV, grids, and future model reconstruction. |
| `[parameters] rscript` | `Rscript` | Interpreter with `tidymodels`, `ranger`, `future`, and plotting packages. |
| `[parameters] calc_python` | current BioMA Python | Python used by the load calculator. |
| `[parameters] calc_script` | bundled | Advanced override; replacement content is recorded as a scientific SHA-256 input. |
| `[parameters] rf_script` | bundled | Advanced RF implementation override. Use only a reviewed, interface-compatible script. |
| `[parameters] predict_script` | bundled | Advanced future-prediction implementation override. |
| `[parameters] future_pattern` | `future_env_ssp*_mean.csv` | Glob selecting files in `future_dir`. |
| `[parameters] expected_future_files` | `0` | A positive value enforces the count; `0` accepts any nonempty set. |

```bash
bioma load workflow.load.example.ini --dry-run
bioma load workflow.load.example.ini
```

`loads/` stores individual/population strict and relaxed load plus sample and
site-overlap QC. `pop_geo_niche_predictors_from_TSS.csv` is a merged copy;
source predictors are unchanged. Two `RF_*` directories store every split's
CV, test, importance, and model outputs. `Figure_S24_loadM_tuning.*` and
`Figure_S24_loadD_tuning.*` remain separate. `future_predictions/`, `maps/`,
and `Figure_S25_load_maps.*` hold future results. Manuscript CV/test R-squared
must come from production settings and report dispersion, sample size, and
spatial extrapolation limits, not smoke-test values.

## 8. MaxEnt niche model

### 8.1 Purpose, principle, and metrics

MaxEnt estimates a maximum-entropy distribution constrained by environmental
features at presence locations relative to background. BioMA first calculates
BIO correlations over the user-selected scope and iteratively removes highly
correlated variables. It then combines candidate variable subsets, feature
classes, and beta multipliers, selecting by k-fold presence/background AUC or
TSS. The final model is refitted with all occurrences and background and
projected to each GCM.

Output is continuous cloglog suitability, not a directly observed occurrence
probability. `maladaptation = current suitability - future suitability`:
positive values mean predicted suitability loss, negative values predicted
improvement, and zero no change. Version 0.8.x does not threshold outputs.

### 8.2 Every configuration field

| Field | Default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] work_dir` | existing path required | Compatibility/work path label. The bundled script currently stores formal intermediates in `output_dir`. |
| `[inputs] occurrence_csv` | required | At least `lon,lat`; at least five unique in-mask records with complete BIO values. `species/cluster` may label plots. |
| `[inputs] current_env_dir` | required | Current BIO TIFFs with uniquely identifiable BIO numbers. |
| `[inputs] future_root` | required | Directories named `<period>-<scenario>-<gcm>`. |
| `[inputs] mask_shp` | required | Background calibration, clipping, and plotting extent. |
| `[inputs] maxent_jar` | empty | User-supplied jar; empty uses `dismo/java/maxent.jar`. At least one source must exist. |
| `[variables] correlation_threshold` | `0.8` | Remove variables with absolute correlation above this value, strictly between 0 and 1. |
| `[variables] correlation_method` | `pearson` | `pearson`, `spearman`, or `kendall`. |
| `[variables] correlation_scope` | `occurrence` | Calculate correlations over `occurrence`, `background`, or up to 10,000 sampled `mask` cells. |
| `[variables] candidate_subset_sizes` | `4,6,8` | Positive candidate subset sizes after correlation filtering; values above remaining variable count are truncated. |
| `[tuning] max_models` | `500` | Maximum variable-set x feature x beta candidates, at least 1. Excess combinations are seeded samples. |
| `[tuning] background_n` | `10000` | Target random-background count, at least 10, excluding occurrences and requiring complete BIO values. |
| `[tuning] cv_folds` | `5` | Presence/background CV folds, at least 2 and capped at occurrence count. |
| `[tuning] feature_classes` | `L,LQ,LQH,LQHP` | Combinations of `L/Q/H/P/T`: linear, quadratic, hinge, product, and threshold. |
| `[tuning] beta_multipliers` | `0.5,1,2,3,4` | Regularization multiplier candidates; larger values generally smooth the response. Positive values are recommended. |
| `[tuning] selection_metric` | `auc` | `auc` or `tss`. Both measure presence/background discrimination, not calibration directly. |
| `[tuning] seed` | `123` | Seed for background points, variable subsets, grid truncation, and folds. |
| `[projection] periods` | `auto` | Period list; `auto` discovers directory tokens. |
| `[projection] scenarios` | `auto` | For example `ssp245,ssp585`; `auto` discovers. |
| `[projection] gcms` | `auto` | GCM list; `auto` discovers. |
| `[projection] ensemble_method` | `mean` | Pixelwise GCM `mean` or `median`; one GCM is used directly. |
| `[projection] binary_output` | `false` | Must be `false` in version 0.8.x; no thresholded habitat is created. |
| `[analysis] output_dir` | required | Result directory. |
| `[analysis] rscript` | `Rscript` | MaxEnt tuning and projection environment. |
| `[analysis] plot_rscript` | `Rscript` | Plotting environment, which may differ from compute R. |

```bash
bioma niche workflow.niche.example.ini --dry-run
bioma niche workflow.niche.example.ini
```

`tables/correlation_matrix.csv`, `variables_after_correlation.txt`, and
`selected_variables.txt` record variable filtering. `tuning_metrics.csv` and
`selected_model.tsv` record model selection. `models/` stores fold and final
fits; `rasters/` stores current, per-GCM, ensemble, and maladaptation rasters;
`figures/` stores tuning, current, future, and maladaptation plots. Examine
occurrence bias, background definition, extrapolation/clamping, AUC/TSS
dispersion, and GCM agreement. Suitability change is not genetic
maladaptation.

## 9. WFmoments 2-D deme

### 9.1 Purpose, principle, and metrics

The module aggregates a current binary habitat into `nx * ny` demes and uses
a four-neighbor, 2-D stepping-stone Wright-Fisher model. The pinned
`wfmoments` numerical implementation uses a closed ordinary-differential-
equation system for the first two moments of the joint allele-frequency
distribution to calculate species-wide neutral nucleotide diversity `pi`
after demes are removed. The current zero-loss state is normalized to 100%.
Outputs include immediate, configured `time3`, configured `time5`, and
optional equilibrium trajectories.

The independent GDAR reference is `A_remaining^z_gdar`; it does not replace
the WF immediate trajectory. In connected edge contraction, an earlier curve
is normally above a later curve. Under random fragmentation, species-wide pi
contains between-deme divergence, so local nonmonotonic behavior is not
automatically an error.

### 9.2 Every configuration field

| Field | Default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] current_raster` | required | Current binary mask; positive is suitable/occupied. Equal-area CRS is recommended. |
| `[inputs] pi_file` | required | `species,pi_obs`, with positive `pi_obs`. |
| `[inputs] structure_file` | empty | `species,fst_global_est`; required for automatic migration calibration. |
| `[inputs] param_file` | empty | `species,z_gdar`; alternatively set `z_gdar` under parameters. |
| `[inputs] area_file` | one of two | `Scenario,Area_km2`, including `current`; future area ratios are clipped to 0-1. |
| `[inputs] future_masks_json` | one of two | JSON mapping scenario names to binary TIFFs with the same shape as current. |
| `[analysis] species` | `my_species` | Exact row key in species tables. |
| `[analysis] output_dir` | required | Result directory. |
| `[analysis] plot_title` | empty | Optional figure title. |
| `[parameters] compute_python` | discovered | Python with `numpy,pandas,rasterio,wfmoments`. |
| `[parameters] rscript` | `Rscript` | Plotting interpreter with `ggplot2`. |
| `[parameters] nx` | `20` | East-west deme count, from 1 through raster column count. |
| `[parameters] ny` | `20` | North-south deme count, from 1 through raster row count. |
| `[parameters] threshold` | `0.25` | Fraction of valid positive pixels required to occupy a block, in [0,1]. |
| `[parameters] min_valid` | `0` | Minimum valid pixels in a block, a nonnegative integer. |
| `[parameters] loss_mode` | `edge` | `edge` for contiguous range contraction; `random` for fragmentation sensitivity. |
| `[parameters] direction` | `east_to_west` | Edge deletion direction: `east_to_west`, `west_to_east`, `north_to_south`, `south_to_north`, or aliases `eastwest/e2w/westeast/w2e/northsouth/n2s/southnorth/s2n`. Ignored by random mode. |
| `[parameters] migration` | `25` | Positive model migration rate between adjacent demes, or `auto`; it is not directly an observed individual proportion. |
| `[parameters] migration_grid` | `0.1,0.3,1,3,10,25,50` | Positive candidates for `migration=auto`. The winner minimizes model east-west FST difference from `fst_global_est`. |
| `[parameters] fst_metric` | `hudson` | `hudson` or `nei` FST for calibration. |
| `[parameters] theta` | `auto` | Positive value or `auto`; automatic scaling matches current equilibrium pi to `pi_obs`. |
| `[parameters] z_gdar` | empty | Explicit positive GDAR exponent, overriding `param_file`. |
| `[parameters] theta_probe` | `1e-4` | Positive initial theta for migration/theta calibration. |
| `[parameters] time3` | `3` | Evolution time for the first future WF curve, normally interpreted as 3 generations. |
| `[parameters] time5` | `5` | Evolution time for the second future WF curve, normally interpreted as 5 generations. |
| `[parameters] mu` | `3.75e-8` | Positive mutation rate used for `Ne = pi_obs/(4*mu)` and the recorded automatic midterm. |
| `[parameters] midterm_generations` | `auto` | `auto` records `Ne/2`. The current figure still computes only `time3/time5`; this field does not replace them. |
| `[parameters] replicates` | `1` | Deletion-order replicates. Production random analyses should exceed 1; identical edge order is deterministic. |
| `[parameters] seed` | `12345` | Random deletion-order seed. |
| `[parameters] include_equilibrium` | `false` | Calculate a new equilibrium for every retained-deme set; computationally expensive. |
| `[parameters] plot_equilibrium` | `false` | Plot equilibrium values, meaningful only if calculated. |
| `[parameters] plot_raw` | `false` | Overlay raw computed points. |

```bash
bioma wfmoment workflow.wfmoment.example.ini --dry-run
bioma wfmoment workflow.wfmoment.example.ini
```

`deme_table.tsv` and `deme_order.tsv` record aggregation and deletion order.
`curve_replicates.tsv` is the immutable per-replicate calculation and
`curve_summary.tsv` its arithmetic mean. `migration_scan.tsv` exists only for
automatic calibration. `scenario_points.tsv` interpolates area scenarios;
`metadata.*` records calibration error and settings. The publication figure's
continuous tail may include an explicitly marked display-only smooth bridge,
but `Figure_wfmoment_2D_deme_curve_data.tsv`, raw summaries, and replicates are
not changed. The zero-deme point at 100% habitat loss is a boundary condition,
not an ordinary observed biological state.

## 10. Multidimensional vulnerability

### 10.1 Purpose, principle, and metric

The module combines continuous niche maladaptation, GF offset, a RONA summary,
relaxed loadM, and relaxed loadD at population coordinates. It first takes the
median of populations in each group, then rescales each metric relative to all
groups in the current analysis:

`score = round(1 + 2*(value-min)/(max-min))`, clipped to 1-3.

If all finite values are equal, the score is 2; if all are missing, it is NA.
High values are assumed to mean higher vulnerability, so replacement fields
must have the same direction. Squares are not significance levels,
probabilities, weights, or counts of independent evidence.

### 10.2 Every configuration field

| Field | Default | Meaning and constraint |
| --- | --- | --- |
| `[inputs] gf_offsets` | required/project `auto` | GF ensemble table with `lon,lat` and `gf_field`. |
| `[inputs] rona_ensemble` | required/project `auto` | Table with `ID` or `pop` and selected `*_RONA` fields. |
| `[inputs] load_predictors` | required/project `auto` | `pop,lon,lat`, grouping field, and loadM/loadD fields. |
| `[inputs] niche_raster` | required/project `auto` | Continuous maladaptation raster. |
| `[analysis] output_dir` | required | Result directory. |
| `[analysis] scenario_label` | `2061-2080 SSP245` | Figure title only. |
| `[parameters] gf_field` | `local_offset_mean` | GF numeric column. Forward/reverse means may be selected if present and documented. |
| `[parameters] rona_variables` | `auto` | `auto/all` discovers every `*_RONA`, or give an explicit comma-separated list. |
| `[parameters] rona_summary` | `mean` | Reduce BIO RONA to one column with `mean`, `max`, or `weighted`; missing values are ignored per population. |
| `[parameters] rona_weights` | empty | For example `BIO1_RONA:0.3,BIO12_RONA:0.7`; nonnegative and renormalized over available positive weights. |
| `[parameters] rona_bio3_field` | empty, deprecated | Legacy compatibility only. Use `rona_variables`; BioMA does not bind BIO3 by default. |
| `[parameters] rona_bio15_field` | empty, deprecated | Legacy compatibility only. BioMA does not bind BIO15 by default. |
| `[parameters] loadM_field` | `mean_loadM_relax` | loadM column. |
| `[parameters] loadD_field` | `mean_loadD_relax` | loadD column. |
| `[parameters] group_column` | `cluster` | Group column in the load predictor table. |
| `[parameters] group_order` | empty | Explicit order. Empty/`auto` preserves first-seen order; unlisted observed groups are appended with a warning. |
| `[parameters] exclude_groups` | empty | Exclusion list. Empty/`none` excludes nothing. Admixed and similar rules belong only in species profiles. |
| `[parameters] nearest_tolerance` | `0.25` | Maximum longitude/latitude Euclidean distance for nearest GF-grid matching, greater than 0; not km. |
| `[parameters] rscript` | `Rscript` | Interpreter with `data.table,ggplot2,raster`. |

```bash
bioma vulnerability workflow.vulnerability.example.ini --dry-run
bioma vulnerability workflow.vulnerability.example.ini
```

`vulnerability_population_values.tsv` stores unscaled population values and
GF match distance. `vulnerability_rona_selection.tsv` fixes auto-discovered
fields and weights. `vulnerability_group_summary.tsv` stores group medians and
1-3 scores. `Figure_multidimensional_vulnerability.png/.pdf` is the panel.
Always publish the raw table with the panel and avoid treating metrics derived
from the same inputs as fully independent evidence.

## 11. Provenance, quality control, and reproducibility

Each standalone `run_manifest.json` records effective configuration, input
content fingerprints, script fingerprints, versions, scenarios, and outputs.
Project mode additionally writes `00_project/input_contract.tsv`,
`input_manifest.json`, effective INIs, state, and an HTML report. Directories
are fingerprinted recursively using stable relative paths. Shapefile sidecars
form one logical dataset. Changing content at the same path changes the
signature.

For publication, retain the actual INIs, manifests, software version/tag,
Conda lock, all QC and model-performance tables, GCM list, random seeds,
locus/sample filtering, and raw uninterpolated results. Spatial interpolation
is visualization, not added observation, and smoothed pixels should not
replace population values in statistical tests.

## 12. Method references

1. Ellis N, Smith SJ, Pitcher CR. Gradient forests: calculating importance gradients on physical predictors. *Ecology* (2012). https://doi.org/10.1890/11-0252.1
2. Fitzpatrick MC, Keller SR. Ecological genomics meets community-level modelling of biodiversity. *Ecology Letters* (2015). https://doi.org/10.1111/ele.12376
3. Rellstab C et al. Signatures of local adaptation in candidate genes of oaks with respect to present and future climatic conditions. *Molecular Ecology* (2016). https://doi.org/10.1111/mec.13889
4. Rellstab C, Dauphin B, Exposito-Alonso M. Prospects and limitations of genomic offset in conservation management. *Evolutionary Applications* (2021). https://doi.org/10.1111/eva.13205
5. Pina-Martins F et al. New insights into adaptation and population structure of cork oak using genotyping by sequencing. *Global Change Biology* (2019). https://doi.org/10.1111/gcb.14497
6. Sang Y et al. Genomic insights into local adaptation and future climate-induced vulnerability of a keystone forest tree in East Asia. *Nature Communications* (2022). https://doi.org/10.1038/s41467-022-34206-8
7. Exposito-Alonso M et al. Genetic diversity loss in the Anthropocene. *Science* (2022). https://doi.org/10.1126/science.abn5642
8. Lin M et al. marApp: An R package and web portal to calculate mutations- and genetic diversity-area relationship for conservation. bioRxiv (2025). https://doi.org/10.1101/2025.09.09.675155
9. Breiman L. Random Forests. *Machine Learning* (2001). https://doi.org/10.1023/A:1010933404324
10. Vaser R et al. SIFT missense predictions for genomes. *Nature Protocols* (2016). https://doi.org/10.1038/nprot.2015.123
11. Phillips SJ, Anderson RP, Schapire RE. Maximum entropy modeling of species geographic distributions. *Ecological Modelling* (2006). https://doi.org/10.1016/j.ecolmodel.2005.03.026
12. Elith J et al. A statistical explanation of MaxEnt for ecologists. *Diversity and Distributions* (2011). https://doi.org/10.1111/j.1472-4642.2010.00725.x
13. Mualim KS et al. Large future genetic diversity losses are predicted from conservation indicators even with habitat protection. *PNAS* (2026). https://doi.org/10.1073/pnas.2514371123

These references explain origins and limitations. A BioMA analysis should also
cite the external packages and exact versions used. `CITATION.cff` will be
added after author and manuscript metadata are finalized.
