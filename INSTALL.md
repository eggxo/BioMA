# Installing BioMA

BioMA source code is distributed under GPL-3.0-only; see `LICENSE`. The
MaxEnt jar and external research packages are separate dependencies with their
own licenses.

The public source repository is available at
<https://github.com/eggxo/BioMA>. The `v0.8.0` tag is the tested release
candidate; use a tagged release rather than an untracked working tree when
reproducing an analysis.

## Supported platform

BioMA 0.8.0 targets Linux, Python 3.8 or newer, R, GDAL command-line tools, and
Java for MaxEnt. The current production tests run on Linux with Python 3.8.15.
Windows can be used to edit configurations and inspect results, but the full
scientific workflow has not been qualified there.

## Conda environment

`environment.yml` records the common Python, R, geospatial, plotting, and Java
dependencies. Create it and install BioMA from the source directory:

```bash
mamba env create -f environment.yml
mamba activate bioma
python -m pip install .
```

For a reproducible Linux installation, use the checked-in `conda-lock.yml`
instead of solving the unpinned environment again. It contains concrete Conda
and pip package URLs and SHA-256 checksums:

```bash
conda-lock install --name bioma conda-lock.yml
conda activate bioma
python -m pip install --no-deps .
```

The lock targets `linux-64`, which is the validated production platform. After
changing `environment.yml`, regenerate it with:

```bash
conda-lock lock -f environment.yml -p linux-64 --lockfile conda-lock.yml
```

The separately maintained `environment-plot.yml` is a lighter plotting-only
environment. It is useful when calculations already run in module-specific
environments; it is not a substitute for the full lock.

Two research R packages are not available as Conda packages. Install the
reviewed versions used by BioMA:

```bash
Rscript -e 'remotes::install_url("https://download.r-forge.r-project.org/src/contrib/gradientForest_0.1-37.tar.gz")'
Rscript -e 'remotes::install_github("meixilin/mar@f2a60772504a52e518d6827a8d22de1ef4dd11d4")'
```

The GF offset command additionally requires `FNN` in the same R environment as
`gradientForest`, `data.table`, and `geosphere`. Install the Conda binary
package into the environment used by GF (the package name is `r-fnn`):

```bash
conda install -n <gf-environment> -c conda-forge r-fnn
```

Do not mix a system R installation with the Conda libraries used by the GF
environment. Confirm that the selected executable sees all four packages:

```bash
bioma doctor workflow.example.ini --rscript "$CONDA_PREFIX/bin/Rscript"
```

For a project that enables GF, use `--strict` after setting the module's
`rscript` to the same environment. A missing `FNN` is reported before any
offset calculation starts.

Then verify the command and run the bundled Python regression tests:

```bash
bioma --version
python -m unittest discover -s tests -v
```

When running directly from a source checkout, `bin/bioma` honors an explicit
`PYTHON` environment variable, otherwise prefers the full `envs/bioma`
interpreter and then falls back to the legacy `envs/bioma-plot` checkout
environment.  A plotting-only environment should therefore be used explicitly
for plotting commands, not as the calculation environment.

Use `bioma doctor` to diagnose an installation without running a workflow:

```bash
bioma doctor
bioma doctor project.ini
bioma doctor project.ini --json > doctor.json
```

The optional project path lets BioMA mark dependencies for enabled modules as
required and discover module-specific `Rscript`, Python, GDAL, Java, and
`maxent.jar` settings. `--strict` treats optional checks as failures, which is
useful in continuous integration. The command is read-only.

## Portable demo fixtures

The source archive includes two small, server-independent fixtures:

* `tests/data/demo/` generates a seven-module project with placeholder rasters,
  three explicit genomic streams, 19 deterministic RONA LD lists, and no
  private paths. Use `--valid-rasters` when testing GDAL and readable GeoTIFFs.
* `tests/data/gf_frequency_50/` contains a fixed 50-site, 175-sample,
  25-population GF frequency fixture and its expected output. It is intended
  for a fast installation/regression check and is not a substitute for
  species-level validation.

The generated demo files are disposable. BioMA reads the fixture inputs and
does not modify them; generated compatibility tables and manifests are written
under the selected results directory.

## MaxEnt license requirement

The `dismo` R package calls the official `maxent.jar`, which is not distributed
with BioMA. Obtain it under its own license and set `maxent_jar` in the niche
configuration. BioMA records the selected path; redistribution rights for the
jar remain separate from the BioMA source license.

## Existing server environments

Module configurations may point to different R or Python executables. This is
useful on an established cluster where GF/RONA, MaxEnt, load/MAR, and WFmoments
already live in tested environments. Leave project-level `rscript` and
`compute_python` blank in that case so each module keeps its own executable.

Run `bin/bioma project project.ini --dry-run` before a long calculation. It
checks the complete input and scenario plan without changing source data.
