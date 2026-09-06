# Changelog

## 0.8.3 - 2026-09-06

- Complete the unified Mamba environment with the Conda/Bioconductor
  dependencies required by MAR and its archived `sars` dependency.
- Install Python numerical and raster packages through Conda so a transient
  GitHub failure while fetching WFmoments does not omit unrelated runtimes.
- Document the reproducible source-install order for `extendedForest`,
  `gradientForest`, `sars`, and `mar`, with an idempotent installation helper.
- Validate one full seven-module server environment with 32 doctor checks,
  84 Python tests, and all bundled R scripts.

## 0.8.2 - 2026-09-06

- Normalize Linux launchers, source, configuration, documentation, and demo
  fixtures to LF through `.gitattributes` while preserving binary files.
- Verify in CI that the source archive retains an executable `bin/bioma` and
  can launch directly after extraction on Linux.

## 0.8.1 - 2026-09-06

- Clarify the public quick-start workflow and the two bundled demo fixtures.
- Document and validate the `r-fnn` and `r-r.utils` dependencies required by GF
  offset calculations and the `bioma doctor --strict` environment check.
- Regenerate the Linux Conda lock after adding `r-r.utils` and verify a real GF
  offset calculation against 54,624 current/future climate cells.

## 0.8.0 - 2026-09-03

- Declare the BioMA source distribution under GPL-3.0-only.
- Add the cross-module input inventory and consolidation boundary document.
- Add a canonical project `[inputs]` layer with generated GF/MAR/load/MaxEnt
  compatibility files and an auditable `input_contract.tsv`.
- Add `bioma project` as a single preflighted, resumable entry point for all
  seven BioMA modules.
- Add standard module directories, shared scenario settings, configuration
  snapshots, checksums, recoverable overwrite backups, project state, TSV
  summary, and HTML report.
- Add automatic GF/RONA/load/niche input wiring for the vulnerability module
  with dependency failure blocking.
- Bundle the loadM/loadD VCF calculator, RF tuning script, and future prediction
  script instead of relying on private absolute script paths.
- Make load tuning `grid_size` and `seed` effective rather than record-only,
  and reject unsupported load schemes.
- Isolate Conda-backed R/Python subprocesses from incompatible login-shell
  libraries and geospatial data paths.
- Allow the niche module to use a licensed, user-supplied `maxent.jar` without
  modifying the installed `dismo` package, and record the jar checksum.
- Keep the revised WFmoments display tail separate from raw calculated tables.
- Make vulnerability RONA variables and summary method configurable, with
  species-agnostic automatic group ordering and no default exclusions.
- Add project, load configuration, and command-status regression tests.
- Add `bioma doctor` for read-only, module-aware checks of Python, R, GDAL,
  Java, MaxEnt, Gradient Forest, MAR, WFmoments, and supporting packages, with
  terminal and JSON reports.
- Record content-level SHA-256 fingerprints for canonical and module-referenced
  inputs; recursively fingerprint directories and shapefile sidecars in
  `00_project/input_manifest.json`.
- Include user-supplied load calculator, RF-tuning, and prediction scripts in
  the project input signature so changing implementation bytes cannot reuse a
  stale result.
- Check project and generated-input fingerprints before materializing
  compatibility files, preventing stale or partially rewritten metadata after
  an input change.
- Validate all 19 RONA LD-pruning lists before computation and record their
  locus counts and content hashes in the RONA manifest.
- Add a reproducible Linux `conda-lock.yml`, a separate plotting environment
  specification, a public server-independent demo generator, and per-module
  end-to-end contract smoke tests.
- Add GitHub Actions checks for Python 3.8/3.10/3.12, source-package contents,
  and the public demo fixture.
- Extend standalone module manifests with recursive input fingerprints and an
  aggregate content signature, matching the project-level provenance contract.
- Keep MaxEnt projections valid for a single configured GCM and make doctor
  GDAL hints follow their command names when only a subset is configured.

## 0.7.0

- Package the independently runnable GF, RONA, MAR, loadM/loadD, niche,
  vulnerability, and WFmoments workflows with configuration examples.
