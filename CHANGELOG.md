# Changelog

## 1.0.0

### Fixes

- Pinned the upstream binary Pfam presence/absence encoding with a regression
  test. Upstream's `np.unique(list(pfams), return_counts=True)` cannot produce
  hit counts because it runs on an already deduplicated set; substituting real
  counts changes predictions.
- `docs/results/` was silently excluded from the repository by an unanchored
  `results` pattern in `.gitignore`, so the published artifacts were missing from
  the release. Local run directories are now root-anchored.
- `.gitignore` no longer ignores the `benchmarks/` records, so the published
  runtime table is verifiable by `scripts/build_runtime_table.py --check`.
- Removed the redundant 26 MB `models/models.zip`; `Archiv.zip` is the pinned
  Zenodo release and is checksum-verified.
- `summary.json` now separates samples excluded for a missing reference label
  from those excluded for an unmappable one. Gram previously reported 51
  excluded without distinguishing that 47 simply had no usable label.
- `bacdive-report` raises a clear error for missing trait columns, non-numeric
  quality fields, and non-finite confidences instead of an uncaught
  `ValueError`/`KeyError`.
- Replaced the GNU-only `sed -i` in the prodigal rule with
  `scripts/clean_proteins.py`; BSD and macOS `sed` require an argument to `-i`.
- Lint and format are clean under the pinned `ruff`, so CI runs as written.

### Added

- Fixed sample identity and metadata joins across `.fa`, `.fna`, `.faa`, TSV,
  and prediction filenames, retaining accession version numbers. The metadata
  join now matches 79 of 79 reference records, previously 0.
- Added strict input, model-checksum, class-label, probability, and output
  completeness validation; failed predictions now return failure.
- Added atomic JSON/CSV outputs, optional file logging, and no automatic shared
  timestamp log files.
- Added genome, protein, and annotation workflow entry points, portable sample
  sheets, explicit model/metadata dependencies, and selected-trait support.
- Added pinned model-compatible environments, packaging, and cached batch inference.
- Rebuilt the 120-genome case study with 79 matched reference records, explicit
  evaluation exclusions, provenance, and an offline HTML report.
- Added a runnable upstream example, regression checks, CI, scope documentation,
  attribution, citations, and the verified upstream GPL-3.0-or-later license.
- Added `docs/VALIDATION.md` recording exactly what was executed and measured,
  and `scripts/build_runtime_table.py` to regenerate the runtime table.
