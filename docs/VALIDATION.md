# Validation record

What was actually executed to check this workflow, and what the results were.
Every figure below was measured on the models and annotations in this repository.

## Environment

| Component | Version |
| --- | --- |
| Python | 3.12.11 |
| numpy | 1.26.4 |
| scipy | 1.12.0 |
| scikit-learn | 1.4.0 |
| Models | BacDive-AI v2, `10.5281/zenodo.15075932` (models byte-identical to v1.1 `10.5281/zenodo.13757323`) |

Library versions are pinned in `pyproject.toml`. Predicted probabilities can
shift slightly between scikit-learn versions, so confidences are not expected to
be bit-identical across environments. Predicted classes are.

## Feature encoding

Upstream `models/predict.py` builds its feature vector with

```python
pfams = set(pfams)
lst = dict(zip(*np.unique(list(pfams), return_counts=True)))
X = [lst.get(k) if k in lst else 0 for k in categories]
```

This reads like a hit count but is not one: `pfams` is already a set, so the list
has no duplicates and every value in `lst` is `1`. The effective upstream
encoding is **presence/absence**. This was verified directly against the bundled
`models/predict.py`.

Using true hit counts instead is a silent behaviour change, not a cleanup. On
the bundled example it moves flagellated-motility confidence from the published
99.85% to 99.35%. `bacdive_workflow/prediction.py` and
`tests/test_core.py::test_features_are_presence_flags_not_hit_counts` pin the
presence/absence encoding so this cannot be "corrected" by mistake.

## Published example

`examples/1120941.3.faa.tsv` is the upstream Pfam annotation of *Actinomyces
dentalis* DSM 19115 (BV-BRC 1120941.3), SHA-256
`89a66a9df4cab7b269ff751006eb87a27917eba7a53c2f08c7f0bca33f8dcc6e`.

```bash
python predict.py all examples/1120941.3.faa.tsv
python scripts/check_example.py results/predictions/1120941.3.json
```

Result: all eight predictions and confidences match the published example.

| Trait | Prediction | Confidence |
| --- | --- | ---: |
| Acidophilic | false | 98.5 |
| Gram-positive | true | 98.93 |
| Spore-forming | false | 96.47 |
| Aerobic | false | 93.91 |
| Anaerobic | true | 63.54 |
| Thermophilic | false | 99.77 |
| Psychrophilic | false | 99.96 |
| Flagellated motility | false | 99.85 |

## Case study reproduction

The 120 saved annotations in `results/interpro_results/` were re-predicted with
the pinned models and compared against `results/predictions/`:

| Measure | Result |
| --- | ---: |
| Genomes | 120 |
| Trait predictions compared | 960 |
| Predicted classes changed | 0 |
| Confidences matching exactly | 956 |
| Confidences differing | 4 |
| Largest confidence difference | 0.01 percentage points |

The four differences are last-digit rounding in the saved run and do not change
any predicted class.

## Model integrity

`bacdive_workflow/model_manifest.json` pins SHA-256 for `v2.zip` (and its nested
`models/models.zip` and training-data members), for each of the eight extracted
models, and for the compatible v1.1 `Archiv.zip`. `python -m bacdive_workflow.models` refuses an
archive or model whose checksum does not match, so a substituted pickle cannot
be loaded silently. `scikit-learn`'s `InconsistentVersionWarning` is promoted to
an error during unpickling.

## Automated checks

Run by CI on every push:

| Check | Command |
| --- | --- |
| Lint | `ruff check bacdive_workflow tests scripts predict.py` |
| Format | `ruff format --check bacdive_workflow tests scripts predict.py` |
| Unit tests | `pytest -q` (25 tests) |
| Artifact reproduction | `python scripts/reproduce_case_study.py --check` |
| Real inference | `snakemake --cores 2` then `scripts/check_example.py` |

`reproduce_case_study.py --check` rebuilds `docs/results/predictions.csv`,
`summary.json`, and `report.html` from the portable snapshots and compares them
byte-for-byte. It validates snapshot-to-report reproducibility. It does not
re-run InterProScan or inference; the table above covers that separately.

## What is not validated

- **No independent accuracy benchmark.** The reference labels come from BacDive,
  the same database the models were trained on. The 69/69 Gram agreement in
  [RESULTS.md](RESULTS.md) is expected from that overlap and is not evidence of
  generalization to new taxa.
- **No runtime or performance guarantee.** `docs/results/runtime.csv` reports
  per-job durations measured on the original server in May 2025. Nothing in CI
  measures runtime.
- **No calibration claim.** Probabilities are model outputs, not measured
  accuracies, and were never calibrated against held-out data here.
- **Genome quality is unverified.** Completeness and contamination are carried
  over from the supplied tables; their estimation method is not established.