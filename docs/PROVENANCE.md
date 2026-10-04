# Data and model provenance

## Published models and example

The eight serialized models and `examples/1120941.3.faa.tsv` originate from
BacDive-AI v2, released by Julia Koblitz and Lorenz Christian
Reimer: https://doi.org/10.5281/zenodo.15075932. The upstream repository is
https://github.com/LeibnizDSMZ/bacdive-AI. The Zenodo record declares
`gpl-3.0-or-later`; this workflow uses the same license. The original Git history
and bundled upstream archive are preserved for attribution. `Archiv.zip` is the
v1.1 archive (`https://doi.org/10.5281/zenodo.13757323`) kept for offline CI and
quick start; its eight model files are byte-identical to v2, which additionally
publishes the training data. `python -m bacdive_workflow.models --download`
fetches `v2.zip`, the pinned release.

The methods paper is Koblitz et al., *Communications Biology* 8, 897 (2025):
https://doi.org/10.1038/s42003-025-08313-3. This project implements workflow
integration and result reporting using those pretrained models. It does not
train new models or claim the upstream authors' model performance as its own.

The example is the upstream Pfam annotation of *Actinomyces dentalis* DSM 19115,
derived from BV-BRC genome 1120941.3. Reference labels come from the upstream
example/strain description. The model manifest pins the upstream archive and
each extracted model by SHA-256. Pickles are accepted only when they match the
known release; arbitrary downloaded pickle files are not supported.

## Existing project case study

The 120 prediction snapshots originate from the project run on 26 May 2025.
The log evidence identifies InterProScan 5.74-105.0 and Pfam 37.3. Snapshot
prediction values and rounded confidence percentages are preserved. Sample IDs
are normalized to fix the metadata join. Positive-class probabilities in the
rebuilt CSV are reconstructed from the saved label/confidence values and hence
retain their original rounding; they are not fresh full-precision inference.

The public sample table is a reduced export of the supplied
`data/merged_genomes_data.csv`. Reference metadata is a reduced export of
`data/bacdive_curated_matches.csv`. Absolute server paths and unrelated source
columns are excluded. Reference values and taxonomic names are retained as
supplied. Their original retrieval dates, individual BacDive strain identifiers,
and training-set overlap were not recorded in these source tables.

Consequently, comparisons are descriptive agreement with the project's curated
labels. They are not an independent accuracy benchmark or experimental
verification. Missing, variable, and combined Gram labels are excluded from
binary comparison. Oxygen comparisons include only explicit pure
aerobe/aerobic/anaerobe labels; ambiguous, facultative, microaerophilic, and
combined labels are excluded. Each metric records its denominator and label
mapping in `docs/results/summary.json`, and separates samples excluded for a
missing label from those excluded for an unmappable one.

Known quirks in the supplied source values, carried through unchanged:

- `sgbs_user_genome` is the literal `NA` for the 41 records not sourced from
  BacDive, rather than an empty field. They are matched by `sample_id`, which is
  present for every genome, so this does not affect the join.
- Oxygen tolerance mixes casing (`anaerobe/microaerophile` and
  `Anaerobe/Microaerophile`). Labels are lowercased before mapping, so both are
  excluded consistently.
- One label reads `Facultitative anaerobe`, a typo for facultative. It is
  excluded either way because only pure labels are scored.
- An `NA` label is cleaned to empty at read time and is therefore reported as a
  missing label rather than an ambiguous one.

Genome completeness and contamination are carried over from the supplied
tables. Their estimation method is not established here. Quality fields are
shown as context, not used to infer validated probability calibration.

The runtime table summarizes 120 saved benchmark records per stage. Times are
per-job durations from the original server. Summed job seconds are not elapsed
wall time because jobs can overlap. They do not measure the current release.
`scripts/build_runtime_table.py` regenerates it from the tracked `benchmarks/`
records, so CI can verify the published table.

## Feature encoding

Models consume binary Pfam presence/absence. Upstream's feature line reads
`dict(zip(*np.unique(list(pfams), return_counts=True)))`, which appears to build
hit counts but cannot: `pfams` is already a set, so every count is 1. This was
verified against the bundled `models/predict.py` and is pinned by a regression
test. Substituting real hit counts changes predictions, for example the bundled
example's motility confidence from 99.85% to 99.35%.

## Reproduction

`python scripts/reproduce_case_study.py --check` checks that the published CSV,
JSON summary, and HTML report can be rebuilt byte-for-byte from the portable
snapshots. It does not rerun InterProScan or inference. The real model example
and fresh local annotation checks are documented separately in
[VALIDATION.md](VALIDATION.md).
