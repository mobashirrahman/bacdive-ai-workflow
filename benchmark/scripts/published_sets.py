"""Extract published Pfam presence sets for seen genomes from the training data.

The features CSV header is mislabelled (values are strain, Pfam, count,
E-value); columns are read positionally. Each set is stored unfiltered as
`pfam<TAB>best E-value` so the drift rule can apply the same threshold to both
sides. A missing or unparseable fourth column is stored as `inf`.
"""

import csv
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

genomes = snakemake.params.genomes
# Map genome -> training strain rows via BacDive ID or accession.
# Fixture and real runs share this file; Phase 3 documents the linkage.
index = {}
with open(snakemake.params.labels, newline="", encoding="utf-8") as handle:
    reader = csv.DictReader(handle)
    cols = reader.fieldnames
    for row in reader:
        index.setdefault(row[cols[0]], []).append(row)
        if row[cols[1]].strip():
            index.setdefault(row[cols[1]].strip(), []).append(row)
            index.setdefault(row[cols[1]].strip().split(".")[0], []).append(row)

features: dict = {}
with open(snakemake.params.features, newline="", encoding="utf-8") as handle:
    reader = csv.reader(handle)
    next(reader)  # mislabelled header; columns are positional
    for row in reader:
        if len(row) < 2:
            continue
        try:
            evalue = float(row[3])
        except (IndexError, ValueError):
            evalue = float("inf")
        strain = features.setdefault(row[0], {})
        strain[row[1]] = min(evalue, strain.get(row[1], float("inf")))

meta_rows = []
for genome in genomes:
    candidates = index.get(genome, [])
    match_type = "missing"
    pfams: dict = {}
    if candidates:
        match_type = "exact"
        for row in candidates:
            for pfam, evalue in features.get(row[cols[0]], {}).items():
                pfams[pfam] = min(evalue, pfams.get(pfam, float("inf")))
    out = Path(f"results/published/{genome}.pfams")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("".join(f"{p}\t{pfams[p]!r}\n" for p in sorted(pfams)))
    meta_rows.append({"genome": genome, "match_type": match_type, "n_published": len(pfams)})

with open(snakemake.output.meta, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=["genome", "match_type", "n_published"], delimiter="\t"
    )
    writer.writeheader()
    writer.writerows(meta_rows)
