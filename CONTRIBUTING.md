# Contributing to BioMA

Run the local contract suite before opening a change:

```bash
python -m compileall -q bioma tests
python -m unittest discover -s tests -v
```

Changes that affect scientific calculations should include a deterministic
fixture, the relevant input-contract test, and a short note in `CHANGELOG.md`.
Do not commit production VCFs, climate rasters, MaxEnt jars, cluster result
directories, credentials, or files containing private server paths.  The
public demo under `tests/data/demo` is the supported portable fixture.

R and MaxEnt end-to-end checks are run in the reviewed Conda environments. A
MaxEnt jar is licensed software and must be supplied by the test operator; it
is never uploaded to this repository.
