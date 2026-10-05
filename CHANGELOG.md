# Changelog

## 1.2.0

### Added

- Benchmark workflow (`benchmark/`) answering how far BacDive-AI predictions
  can be trusted on unseen genomes, incomplete genomes, and under a different
  annotation route. Full report in `docs/BENCHMARK.md`, generated tables and
  figures in `docs/benchmark/`.
- Leakage-controlled accuracy on 150 seen and 174 unseen genomes (32
  genus-unseen): unseen matches seen within bootstrap intervals on all six
  eligible traits; 35 disagreements listed for follow-up.
- Annotation drift against the published training features at a matched
  E-value threshold: the training-era InterProScan 5.63-95.0 reproduces the
  published Pfam sets (2/150 genomes change a prediction); 5.74-105.0 gives a
  median Jaccard of 0.86 and changes a prediction in 15/150 genomes.
- Completeness degradation by seeded fragment loss with real re-annotation
  shortcut validation (318/320 agree). Reported as the share of full-genome
  positive calls lost: at 50% completeness 97% of spore-forming, 84% of
  motile and 63% of thermophile calls turn negative, while negatives almost
  never turn positive.
- Annotation robustness across E-value, gene caller, and Pfam release
  (class agreement at least 0.94 except motility/spore formation at strict
  E-values).
- Redesigned benchmark figures (colour-blind-checked palette, direct labels,
  byte-identical across reruns) and a README that leads with the findings.
- `benchmark/.test/stage.sh` stages the offline fixture with timestamps in
  dependency order, so the test run never invokes the real InterProScan.
- Frozen selection (`benchmark/config/genomes.tsv`, `label_mapping.yaml`,
  `training_lookup.json`) with calibrated BacDive mappings (all six traits
  pass the 95%-on-500 gate), offline `.test` config, and
  `scripts/build_benchmark_tables.py --check` in CI.

## 1.1.0

### Changed

- Pinned the BacDive-AI v2 release (`10.5281/zenodo.15075932`). The eight
  model files are byte-identical to v1.1, so predictions are unchanged; v2
  additionally publishes the training data. `model_manifest.json` records the
  v2 archive, its nested `models/models.zip`, the training CSV sizes and
  hashes, and the compatible v1.1 `Archiv.zip` entry.
- `python -m bacdive_workflow.models` accepts either archive (identified by
  SHA-256, with the nested v2 zip opened without extracting it to disk) and
  gains `--training-data DIR` to extract and verify the two training CSVs.
  The default archive stays `Archiv.zip` when present, else `v2.zip`;
  `--download` fetches `v2.zip`. `Archiv.zip` remains in git for offline CI.
- Updated README, USAGE, VALIDATION, PROVENANCE, and the report footer to the
  v2 DOI; RESULTS notes the historical run used byte-identical models.

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
