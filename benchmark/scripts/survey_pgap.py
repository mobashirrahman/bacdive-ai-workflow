"""PGAP availability survey: which robustness genomes have PGAP proteins."""

import csv
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

rows = []
for genome in snakemake.params.genomes:
    staged = Path(f"results/genomes/{genome}.pgap.faa")
    available = staged.is_file() and staged.stat().st_size > 0
    rows.append({"genome": genome, "pgap_available": available})

with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(handle, fieldnames=["genome", "pgap_available"], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows)
