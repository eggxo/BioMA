# Installing BioMA

BioMA source code is distributed under GPL-3.0-only; see `LICENSE`. The
MaxEnt jar and external research packages are separate dependencies with their
own licenses.

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
conda-lock install --name bioma --file conda-lock.yml
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

Then verify the command and run the bundled Python regression tests:

```bash
bioma --version
python -m unittest discover -s tests -v
```

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
