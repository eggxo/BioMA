# BioMA MAR module

The MAR module wraps `MAR.all.R` as one reproducible command. It runs
`MARPIPELINE`, performs the configured direction test (default `random`),
automatically sets `maxsnps` to the number of variant records in the input VCF,
and creates the official MAR power-law fit and habitat-loss curves.

## Configure and run

Copy `workflow.mar.example.ini`, edit the VCF, sample coordinate table, and
output path, then run:

```bash
bin/bioma mar workflow.mar.ini --dry-run
bin/bioma mar workflow.mar.ini
```

`scheme` may be `random`, `eastwest`, `westeast`, `northsouth`, or
`southnorth`. `maxsnps = auto` is recommended; an integer can be supplied to
cap the number of variants. The selected R environment must provide the
`mar` package; put its `Rscript` on `PATH` or configure it explicitly.

Outputs include the MAR pipeline RDA/GDS objects, `MAR_official_fit_params.tsv`,
`MAR_display_curves.tsv`, `MAR_habitat_loss_curves.png/.pdf`, logs, and
`run_manifest.json`. If `scenario_file` is supplied, scenario points are added
to the curve plot and written to `MAR_scenario_points.tsv`.
