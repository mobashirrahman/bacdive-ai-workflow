# Published features file: characterization (20 seen genomes)

The v2 `training_data_features.csv` header reads
`ID_strains,sequence_acc,pfam,pfam_description`, but the values are
`strain ID, Pfam accession, integer, E-value` (e.g. `17385,PF00550,14,3.7e-14`).
The schema is undocumented, so before using the file for annotation drift we
checked the two guesses on 20 seen genomes against our own annotation
(InterProScan 5.74-105.0, Prodigal single, E-value 1e-20).

## Column 3: number of hits per Pfam?

Plausible but unconfirmed. Across the 20 genomes, column 3 is 1 in 61% of
rows, 2-5 in 33%, and above 5 in 9% — the shape hit counts would have. Where
a Pfam is found by both routes, the published integer equals our per-Pfam
protein count in about two thirds of cases (e.g. 850/1252, 1025/1462). The
remainder is expected: different assembly versions, different gene calls, and
a different Pfam release all shift counts. Consistent with hit counts, not
proof of them.

## Column 4: best E-value per Pfam?

Plausible but unconfirmed. Values are tiny floats spanning many orders of
magnitude. 73% of rows are at or below 1e-20, 23% sit between 1e-20 and
1e-10, and 10% are above 1e-10 — so the table was not filtered at the 1e-20
threshold the prediction script uses. Applying a 1e-20 cutoff to column 4
moves the published presence set markedly closer to ours (median Jaccard
0.60 unfiltered vs 0.84 filtered over the 20 genomes), which is how an
E-value column should behave. The residual gap reflects genuinely different
annotation routes, not just thresholding.

## Consequence for drift

Analysis B applies the same 1e-20 threshold to both sides: our annotation
through `parse_pfams`, and the published set by keeping Pfams whose fourth
column is at or below 1e-20. Jaccard similarity, gained and lost Pfams, and
prediction changes are all computed on those matched sets. The comparison
against the published sets as shipped (unfiltered) is kept in separate
`_unfiltered` columns of `drift.tsv`, because it mixes annotation drift with
the threshold difference.

The full run supports the E-value reading: with the training-era
InterProScan 5.63-95.0 the matched sets are identical for the median seen
genome (see the drift table in `docs/BENCHMARK.md`).
