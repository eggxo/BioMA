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

Before publishing a source archive, run `python -m unittest
tests.test_release_boundary`, inspect `MANIFEST.in`, and check the archive
contents (`tar -tzf dist/*.tar.gz`). The release archive must contain the
portable demo and lock file, and must not contain `workflow.*.test.ini`, build
outputs, old archives, private paths, credentials, or research data.

R and MaxEnt end-to-end checks are run in the reviewed Conda environments. A
MaxEnt jar is licensed software and must be supplied by the test operator; it
is never uploaded to this repository.
