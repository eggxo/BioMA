# BioMA loadM/loadD module

This module wraps the derived-allele VCF load calculation, random-forest
hyperparameter tuning, future prediction, and spatial interpolation used for
the genetic-load figures. It is exposed as one BioMA command.

The VCF calculator, RF tuning script, and future-prediction script are bundled
with BioMA; no site-specific script path is required. Advanced users can
replace any of them with the reviewed `calc_script`, `rf_script`, and
`predict_script` configuration keys. Their paths and SHA-256 checksums are
recorded in `run_manifest.json`.

## Input assumptions

`vcf_dir` must contain exactly these four files:

- `strict.synonymous.vcf`
- `strict.nonsynonymous.vcf`
- `strict.deleterious_plain.vcf`
- `relaxed.deleterious_with_warning.vcf`

The VCFs are expected to be the already polarised/derived-allele files. The
module does not infer ancestral or derived states. Population names must match
the `pop` column in the predictor table and the population coordinate/niche
files used by the supplied load calculator.

## Configure and run

Copy `workflow.load.example.ini`, edit paths if needed, and validate it:

```bash
bin/bioma load workflow.load.example.ini --dry-run
```

Run the complete workflow with:

```bash
bin/bioma load workflow.load.example.ini
```

For a first test, use `workflow.load.test.ini` (two outer splits, one CV
repeat). Production settings in the example use more splits and workers; the
full run is computationally expensive because each target performs a 50-point
Latin-hypercube random-forest search within every split.

## Parameters

- `n_splits`, `cv_v`, `cv_repeats`, `workers`: RF resampling and parallelism.
- `grid_size`: number of Latin-hypercube RF parameter combinations evaluated
  in every outer split.
- `seed`: reproducibility seed used for outer splitting, tuning grids,
  cross-validation, and future prediction.
- `scheme`: currently must be `svd`; unsupported values are rejected rather
  than silently ignored.
- `rscript`: R runtime, discovered from `PATH` unless explicitly configured.
- `calc_python`: Python interpreter for
  `calc_loadM_loadD_from_3vcf.py` (defaults to the BioMA interpreter).

The RF targets are `mean_loadM_relax` and `mean_loadD_relax`, matching the
future plotting script. Strict and relaxed population load tables are both
calculated and written; the strict/relaxed columns are merged into a copied
predictor table, so the source table is never modified.

## Outputs

Under `output_dir`:

- `loads/`: strict and relaxed `population_load.tsv`, logs, and calculator
  diagnostics;
- `pop_geo_niche_predictors_from_TSS.csv`: predictor table with four load
  columns;
- `RF_mean_loadM_relax/` and `RF_mean_loadD_relax/`: tuning metrics, selected
  splits, model objects, current checks, and plots;
- `future_predictions/`: one prediction table per target and SSP/period;
- `Figure_S24_loadM_tuning.png/.pdf` and
  `Figure_S24_loadD_tuning.png/.pdf`: reference-style RF tuning panels
  (mean RMSE across CV folds versus `mtry`, `min_n`, and `trees`), with the
  selected split's optimum marked by a red diamond;
- `maps/` and `Figure_S25_load_maps.png/.pdf`: GAM-interpolated species-mask
  maps. Interpolation uses the arithmetic population/model input values and a
  0.05-degree grid;
- `run_manifest.json`, `load_plot.log`, and target-specific RF logs.

The small test run is a workflow check, not a final performance estimate.
Use the production split/repeat settings before reporting CV and test R² in a
manuscript.
