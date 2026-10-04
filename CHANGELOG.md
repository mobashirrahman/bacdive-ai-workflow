# Changelog

## 1.0.0

- Fixed sample identity and metadata joins across `.fa`, `.fna`, `.faa`, TSV,
  and prediction filenames, retaining accession version numbers.
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
