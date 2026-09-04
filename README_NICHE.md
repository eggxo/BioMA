# BioMA niche module

The niche module runs a MaxEnt-only species-distribution workflow with
automatic variable selection, parameter tuning, and configurable future
projections. It keeps suitability continuous; binary habitat maps are not
created in the current version.

## One-command use

Edit `workflow.niche.example.ini` and validate the paths and scenario grid:

```bash
bin/bioma niche workflow.niche.example.ini --dry-run
bin/bioma niche workflow.niche.example.ini
```

For a first portable test, generate the public demo and use its
`configs/niche.ini` (the smoke run stops at configuration validation):

```bash
python tests/data/demo/make_demo.py --outdir /tmp/bioma-demo
bin/bioma niche /tmp/bioma-demo/configs/niche.ini --dry-run
```

The historical `workflow.niche.test.ini` profile is kept only in a checkout
for site-specific regression runs and is not included in release archives.
The production template uses a larger tuning grid and 10,000 background
points.

## Variable selection and tuning

The correlation threshold is user-defined (default `0.8`) and is used only as
a pre-filter. The correlation scope can be `occurrence`, `background`, or
`mask`. After filtering, the module generates candidate variable subsets and
combines them with the configured MaxEnt feature classes and beta multipliers.
Cross-validated AUC (or TSS when selected) chooses the final model. The full
candidate table and selected variables are saved for reproducibility.

Candidate model count is capped by `max_models`; this prevents an uncontrolled
combinatorial explosion when many BIO variables are available. Set
`candidate_subset_sizes` and the feature/beta lists to control the search.

The MaxEnt calculation uses the environment containing `dismo` and
`maxent.jar`; map rendering can use a separate `plot_rscript` environment.

## Projections

Periods, SSP scenarios, and GCM names are read from the configuration. The
module expects future directories named `<period>-<scenario>-<gcm>` and checks
that every requested directory contains uniquely identifiable BIO tif files.
Each GCM is projected separately, then the configured `mean` or `median`
ensemble is calculated for each period and scenario. Maladaptation is the
continuous `current suitability - future suitability`.
When one GCM is configured, that projection is used directly as the ensemble;
there is no artificial second model or change to the suitability values.

## Outputs

Under `output_dir`:

- `tables/correlation_matrix.csv` and `variables_after_correlation.txt`;
- `tuning_metrics.csv` and `selected_model.tsv`;
- `models/` with candidate fold fits and the final `best_model`;
- `rasters/` with current, per-GCM, ensemble suitability, and maladaptation;
- `tables/maladaptation_points_*.csv` and `projection_manifest.csv`;
- `figures/` with current suitability, future suitability maps, continuous
  maladaptation maps, combined map panels, and a tuning plot;
- `niche_pipeline.log`, `niche_plot.log`, and `run_manifest.json`.

The run manifest records content fingerprints for the occurrence table,
climate directories, mask components, MaxEnt jar, and the bundled scripts.

The module uses the species mask for raster clipping. Country or nine-dash
boundaries are not required for calculation and can be added to the plotting
script later as optional display layers.

## Multidimensional vulnerability panel

After GF, RONA, niche, and loadM/loadD have been run, the saved population-level
tables can be combined using the companion vulnerability command:

```bash
bin/bioma vulnerability workflow.vulnerability.example.ini --dry-run
bin/bioma vulnerability workflow.vulnerability.example.ini
```

The panel follows `多维度不适应性进度图.R`: rows are population groups and each
metric is represented by three filled squares. The module writes raw values for
every population, group medians, and 1--3 scores. RONA variables are selected
with `rona_variables` (a comma-separated list such as
`BIO3_RONA,BIO15_RONA`, or `auto` to discover every `*_RONA` column). The
selected values are reduced to one panel metric with `rona_summary = mean`,
`max`, or `weighted`; the latter uses `rona_weights = BIO3_RONA:0.5,BIO15_RONA:0.5`.
Selected per-BIO values remain in the raw and group-summary tables for audit.
Missing values are ignored; weighted summaries renormalize the positive weights
available for each population, and remain missing when no selected value is
available.

`group_order` and `exclude_groups` are both optional. Leave them blank (or use
`auto`/`none`) to retain groups in first-seen input order and include every
group. An explicit order is honoured and groups present in the data but absent
from that list are appended with a warning. This keeps species-specific group
names out of the software defaults; the test profile contains its own explicit
values only to reproduce the historical example.

The panel metrics are continuous niche maladaptation, GF local offset, the
configured RONA summary, relaxed loadM, and relaxed loadD. Higher values receive
more filled squares. GF grid values are matched to each population by nearest
coordinate within `nearest_tolerance`.

The bundled `tests/test_vulnerability_r.R` exercises automatic BIO discovery,
mean/max/weighted summaries, automatic group ordering, and explicit ordering
with a small synthetic raster fixture.
