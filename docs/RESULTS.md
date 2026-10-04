# Results of the 120-genome case study

These results reproduce the saved project run of 26 May 2025, with corrected
sample identity and metadata joining. They use the published BacDive-AI v1.1
models. Full provenance and interpretation limits are in [PROVENANCE.md](PROVENANCE.md).

## Outcomes

- 120 genomes, with all eight predictions present for every genome.
- 79 reference metadata records matched; the previous aggregation matched zero.
- Re-running the pinned models over the saved annotations reproduces all 960
  trait predictions: no class changes, and 956 of 960 confidence percentages
  match exactly. The 4 that differ are off by 0.01 percentage points, which is
  rounding in the saved run. See [VALIDATION.md](VALIDATION.md).
- 69 genomes have unambiguous binary Gram labels. All 69 predictions agree with
  the supplied labels (35 positive and 34 negative).
- Oxygen tolerance is scored per trait against the explicit single-lifestyle
  reference labels, not pooled:

  | Trait | Eligible records | Agreement | Breakdown |
  | --- | ---: | ---: | --- |
  | Aerobic | 29 | 28/29 (96.6%) | 2 of 3 `aerobe`/`aerobic` predicted aerobic; 26 of 26 `anaerobe` predicted not aerobic |
  | Anaerobic | 29 | 25/29 (86.2%) | 22 of 26 `anaerobe` predicted anaerobic; 3 of 3 `aerobe`/`aerobic` predicted not anaerobic |

  Each row scores both reference groups against one trait, so the denominator is
  29 rather than 3 or 26. Per-group accuracy is in the right-hand column. The
  pooled figures must not be read as "28 of 29 aerobes".
- Completeness and contamination are available for 79 genomes. Median supplied
  completeness is 66.02%; median contamination is 0.83%.

## Why samples are excluded from a comparison

`summary.json` reports the two reasons separately, because combining them
overstates how much data was merely ambiguous:

| Comparison | Eligible | No usable label | Ambiguous label |
| --- | ---: | ---: | ---: |
| Gram | 69 | 47 | 4 |
| Aerobic | 29 | 47 | 44 |
| Anaerobic | 29 | 47 | 44 |

The 47 genomes without a usable label are the 41 sourced outside BacDive, which
have no reference record at all, plus 6 whose source label is `NA`. An empty or
`NA` label is cleaned to empty during reading, so the two cases cannot be told
apart downstream. The 4 ambiguous Gram labels are `variable` and
`negative,positive`.

Agreement is descriptive, and the Gram figure is additionally circular: the
reference labels come from BacDive, the same curated database the models were
trained on, so near-perfect agreement is expected and is not evidence of
generalization. The reference table does not establish strain-level identity,
experimental verification, or separation from training data. These figures must
not be presented as independent model accuracy. Model probabilities are not
measured accuracy, and incomplete genomes can limit the available feature
inventory.

## Trait distribution

| Trait | Positive | Negative |
| --- | ---: | ---: |
| Acidophilic | 0 | 120 |
| Gram-positive | 52 | 68 |
| Spore-forming | 0 | 120 |
| Aerobic | 13 | 107 |
| Anaerobic | 50 | 70 |
| Thermophilic | 0 | 120 |
| Psychrophilic | 0 | 120 |
| Flagellated motility | 7 | 113 |

## Artifacts

- [Predictions and joined metadata](results/predictions.csv)
- [Machine-readable summary and reference comparisons](results/summary.json)
- [Standalone interactive HTML report](results/report.html) (download and open locally)
- [Saved runtime summaries](results/runtime.csv)

The HTML report can be searched by sample, taxon, or result. It is self-contained
and works offline. The saved prediction snapshots have rounded probabilities;
the reconstructed positive-class probabilities retain that rounding.

## Reproduce the artifacts

```bash
python scripts/reproduce_case_study.py
python scripts/reproduce_case_study.py --check
python scripts/build_runtime_table.py --check
```

The first rebuilds the CSV, HTML, and JSON from the saved prediction snapshots;
the checks fail if the published files are stale. This is artifact reproduction
using supplied predictions, not fresh inference. Re-running the models over the
saved annotations is covered by [VALIDATION.md](VALIDATION.md), and the runnable
upstream annotation example by `scripts/check_example.py`.

## Original runtime

| Stage | Jobs | Mean seconds per job | Median seconds per job |
| --- | ---: | ---: | ---: |
| Prodigal | 120 | 18.85 | 14.74 |
| InterProScan | 120 | 176.41 | 166.36 |
| Prediction | 120 | 2.13 | 1.69 |

The source logs identify InterProScan 5.74-105.0 and Pfam 37.3. These are measured
job durations on the original server, not portable performance guarantees or
the wall time of the release workflow. Reusing annotations avoids the largest
recorded stage when those annotations already exist.
