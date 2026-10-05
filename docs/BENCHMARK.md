# Benchmark: how far can BacDive-AI predictions be trusted?

## Question

> How far can BacDive-AI predictions be trusted on genomes the models never
> saw, on incomplete genomes, and under a different annotation route?

The models were trained on complete-ish isolate genomes annotated with
InterProScan and evaluated with random 5-fold cross-validation only. Downstream
use (e.g. metaTraits on 1.1 million MAGs at >50% completeness with
eggNOG-mapper Pfam annotation) leaves that evaluation behind in three ways.
This benchmark measures each one separately with a reproducible Snakemake
workflow (`benchmark/`). Tables below are generated from
`docs/benchmark/*.tsv` by `scripts/build_benchmark_tables.py` and are never
typed by hand.

## Design

- **A. Leakage-controlled accuracy.** 150 seen genomes (training strains, an
  upper bound — training-set performance, not accuracy) against 174 unseen
  genomes with no strain, genome, or species overlap with training, plus a
  32-genome genus-unseen subset. Reference labels come from calibrated BacDive
  mappings, never from BacDive-AI outputs.
- **B. Annotation drift.** Seen-set Pfam sets from our annotation (two
  InterProScan versions) versus the published training feature sets.
- **C. Completeness degradation.** Seeded fragment loss (8 levels, 10
  replicates, contamination at 100% and 70%), with a real re-annotation
  shortcut validation on 20 genomes.
- **D. Annotation robustness.** E-value thresholds, gene caller (Prodigal
  single/meta, PGAP), and Pfam release on a 100-genome subset.

Primary configuration everywhere: InterProScan 5.74-105.0 / Pfam 37.3,
Prodigal single, E-value 1e-20. Full run wall time was about 15 hours on
16 cores with at most 4 concurrent InterProScan jobs.

## Genome selection and leakage control

12,412 training strain IDs and 24,426 candidate records were harvested from
the BacDive API (100-ID batches, cached raw JSON, harvest 2026-10-04). Seen
genomes are training strains with a GCA/GCF accession and published feature
rows. Unseen genomes share no BacDive ID, no accession (INSDC or BV-BRC), and
no species name with training; 32 additionally share no genus. Assemblies are
Complete Genome or Chromosome (17 recorded ≤50-contig top-ups supply rare
positives). Sampling is seeded (20261004), proportional by phylum with a
floor of 3, one genome per species, across 38 phyla; 307 of 324 genomes are
type strains. Frozen in `benchmark/config/genomes.tsv`; exclusion sets in
`benchmark/config/training_lookup.json`, enforced by `tests/test_selection.py`.

## Label calibration and the gate

BacDive field values were cross-tabulated against the 44,747 published
training labels; a value maps to a class only at ≥95% purity over ≥30
strains, and conflicting votes exclude the strain. The mapping reproduces
98.3–99.9% of training labels over thousands of strains per trait, so all
six traits (gram-positive, motile2+, anaerobic, aerobic, thermophile,
spore-forming) enter accuracy analysis. Acidophile and psychrophile have no
training rows and are stability-only. Full cross-tabs in
`docs/benchmark/label_mapping_calibration.csv`.

## A. Leakage-controlled accuracy

Seen numbers are training-set performance (upper bound), not accuracy.

<!-- benchmark:accuracy -->
| Trait | Set | n | Pos | Neg | Sens | Spec | Bal.acc [95% CI] | MCC | AUROC |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| aerobic | seen | 107 | 36 | 71 | 0.944 | 0.944 | 0.944 [0.895, 0.986] | 0.877 | 0.984 |
| aerobic | unseen | 108 | 20 | 88 | 0.900 | 0.966 | 0.933 [0.858, 0.989] | 0.850 | 0.991 |
| aerobic | genus_unseen | 22 | 1 | 21 | — | — | — — | — | — |
| anaerobic | seen | 107 | 46 | 61 | 0.978 | 1.000 | 0.989 [0.967, 1.000] | 0.981 | 1.000 |
| anaerobic | unseen | 87 | 13 | 74 | 1.000 | 0.973 | 0.986 [0.966, 1.000] | 0.918 | 1.000 |
| anaerobic | genus_unseen | 13 | 3 | 10 | — | — | — — | — | — |
| gram-positive | seen | 77 | 26 | 51 | 0.923 | 0.961 | 0.942 [0.884, 0.990] | 0.884 | 0.981 |
| gram-positive | unseen | 98 | 39 | 59 | 0.923 | 0.966 | 0.945 [0.894, 0.983] | 0.893 | 0.984 |
| gram-positive | genus_unseen | 14 | 4 | 10 | — | — | — — | — | — |
| motile2+ | seen | 46 | 5 | 41 | — | — | — — | — | — |
| motile2+ | unseen | 66 | 15 | 51 | 1.000 | 0.882 | 0.941 [0.892, 0.980] | 0.794 | 0.995 |
| motile2+ | genus_unseen | 9 | 2 | 7 | — | — | — — | — | — |
| spore-forming | seen | 28 | 3 | 25 | — | — | — — | — | — |
| spore-forming | unseen | 58 | 21 | 37 | 1.000 | 0.973 | 0.986 [0.959, 1.000] | 0.964 | 1.000 |
| spore-forming | genus_unseen | 8 | 0 | 8 | — | — | — — | — | — |
| thermophile | seen | 146 | 27 | 119 | 0.889 | 1.000 | 0.944 [0.870, 1.000] | 0.931 | 0.999 |
| thermophile | unseen | 169 | 30 | 139 | 0.967 | 0.993 | 0.980 [0.943, 1.000] | 0.959 | 1.000 |
| thermophile | genus_unseen | 31 | 9 | 22 | — | — | — — | — | — |
<!-- /benchmark:accuracy -->

Unseen accuracy matches seen accuracy within bootstrap intervals on every
eligible trait — there is no leakage cliff. The 35 disagreements
(`docs/benchmark/disagreements.tsv`: 11 aerobic, 9 gram-positive, 6
motile2+, 5 thermophile, 3 anaerobic, 1 spore-forming) are listed with
genome, taxon, label, and probability for follow-up. Genus-unseen groups
(9–31 genomes, fewer than 10 per class) report counts only, by design.

## C. Completeness degradation

Flip rate relative to the full-genome prediction (unseen set, no
contamination); direction splits positive→negative and negative→positive.
Losses are overwhelmingly positive→negative: losing genes removes positive
signals rather than inventing them.

<!-- benchmark:flip -->
| Trait | L100 | L70 | L50 | L30 | Direction at L50 |
| --- | ---: | ---: | ---: | ---: | --- |
| acidophile | 0.000 | 0.006 | 0.006 | 0.006 | +→− 0.006, −→+ 0.000 |
| aerobic | 0.000 | 0.027 | 0.078 | 0.219 | +→− 0.075, −→+ 0.002 |
| anaerobic | 0.000 | 0.006 | 0.020 | 0.134 | +→− 0.018, −→+ 0.001 |
| gram-positive | 0.000 | 0.021 | 0.036 | 0.078 | +→− 0.022, −→+ 0.014 |
| motile2+ | 0.000 | 0.194 | 0.339 | 0.394 | +→− 0.339, −→+ 0.000 |
| psychrophile | 0.000 | 0.000 | 0.000 | 0.000 | +→− 0.000, −→+ 0.000 |
| spore-forming | 0.000 | 0.140 | 0.202 | 0.207 | +→− 0.202, −→+ 0.000 |
| thermophile | 0.000 | 0.039 | 0.116 | 0.175 | +→− 0.116, −→+ 0.000 |
<!-- /benchmark:flip -->

Per-trait lowest completeness with flip rate below 5% (bootstrap interval
over genomes):

<!-- benchmark:thresholds -->
| Trait | Lowest completeness, flip rate <5% | Bootstrap interval over genomes |
| --- | ---: | ---: |
| acidophile | 30% | [30%, 30%] |
| aerobic | 60% | [60%, 70%] |
| anaerobic | 50% | [40%, 50%] |
| gram-positive | 40% | [30%, 50%] |
| motile2+ | 100% | [90%, 100%] |
| psychrophile | 30% | [30%, 30%] |
| spore-forming | 90% | [80%, 90%] |
| thermophile | 80% | [70%, 80%] |
<!-- /benchmark:thresholds -->

Motility and spore formation are the fragile traits (5% thresholds at 100%
and 90%); acidophile and psychrophile predictions barely move at any
completeness. Contamination at 5–10% foreign bases shifts flip rates by
roughly 0.5–2 points (e.g. gram-positive at L70: 0.021 clean, 0.025 at 5%,
0.042 at 10%). The cheap simulation is trustworthy: shortcut validation
against real re-annotation agrees on 318/320 comparisons with a mean
probability difference of 0.008.

## B. Annotation drift against the published features

The published features file has an undocumented schema; characterization on
20 seen genomes (`docs/benchmark/features_file_characterization.md`) finds
column 3 consistent with (but not proven to be) hit counts and column 4
consistent with (but not proven to be) best E-values. Drift therefore
compares unfiltered presence sets.

<!-- benchmark:drift -->
| InterProScan | Genomes | Median Jaccard | Min | Genomes with prediction changes |
| --- | ---: | ---: | ---: | ---: |
| 5.63-95.0 | 150 | 0.686 | 0.499 | 26 |
| 5.74-105.0 | 150 | 0.604 | 0.433 | 34 |
<!-- /benchmark:drift -->

Drift is substantial: median Jaccard 0.60 against our 5.74-105.0 annotation
(0.69 against the training-era 5.63-95.0), and 34 of 150 genomes change at
least one prediction when the published set is fed to the models instead of
ours. Annotation route matters as much as the models.

## D. Annotation robustness

Class agreement with the primary configuration on 100 genomes (93 with PGAP
proteins), plus mean absolute probability difference and accuracy on the
labelled subset:

<!-- benchmark:robustness -->
| Trait | Caller meta | Caller PGAP | Pfam 5.63 | E-value 1e-10 | E-value 1e-30 |
| --- | ---: | ---: | ---: | ---: | ---: |
| acidophile | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| aerobic | 1.000 | 1.000 | 1.000 | 0.990 | 1.000 |
| anaerobic | 1.000 | 1.000 | 0.990 | 1.000 | 0.970 |
| gram-positive | 0.990 | 0.989 | 0.990 | 0.970 | 0.980 |
| motile2+ | 0.990 | 0.989 | 0.940 | 0.920 | 0.820 |
| psychrophile | 1.000 | 1.000 | 1.000 | 1.000 | 0.990 |
| spore-forming | 1.000 | 0.989 | 1.000 | 0.990 | 0.890 |
| thermophile | 1.000 | 1.000 | 0.990 | 0.990 | 0.960 |
<!-- /benchmark:robustness -->

Gene caller and Pfam release barely move classes (agreement ≥0.94
everywhere). E-value thresholds move probabilities more than classes, except
for motility (0.82 agreement at 1e-30) and spore formation (0.89) — the same
two traits degradation flagged as fragile.

## Limitations

- Unseen-set labels come from BacDive, the same curation source as training.
  The test controls for strain, genome, and species leakage, not for shared
  curation practice.
- Simulated completeness is fraction of bases retained, not a CheckM
  estimate, and random fragment loss is kinder than real binning, which loses
  genomic islands and plasmids preferentially.
- Acidophile, psychrophile, and any gate-failing trait have stability results
  only. Seen-set numbers are training-set performance.
- The published features file has an undocumented schema; Phase 6
  conclusions depend on the characterization above.

## Reproduce

```bash
snakemake -s benchmark/Snakefile --directory <data_root>/run \
  --workflow-profile benchmark/profiles/local --sdm conda --cores 16 \
  --resources ips=4
python scripts/build_benchmark_tables.py --check
```

`benchmark/.test/` runs the predict, degrade, metrics, and report rules
offline in minutes with the real models. Per-table SHA-256 records are in
`docs/benchmark/MANIFEST.tsv`; the workflow graph is
`docs/benchmark/rulegraph.svg`.
