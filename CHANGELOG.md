# Changelog

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

## 0.7.0

- Package the independently runnable GF, RONA, MAR, loadM/loadD, niche,
  vulnerability, and WFmoments workflows with configuration examples.
