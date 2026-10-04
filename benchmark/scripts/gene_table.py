"""Parse Prodigal FASTA headers into a gene coordinate table."""

import re
import sys

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.common import WorkflowError

HEADER = re.compile(r"^>(\S+)\s+#\s+(\d+)\s+#\s+(\d+)\s+#\s+(-?1)\s+#")
rows = []
with open(snakemake.input.proteins, encoding="utf-8") as handle:
    for line in handle:
        if not line.startswith(">"):
            continue
        match = HEADER.match(line)
        if not match:
            raise WorkflowError(f"Cannot parse Prodigal header: {line.strip()}")
        protein, start, end, strand = match.groups()
        contig = protein.rsplit("_", 1)[0]
        rows.append((protein, contig, start, end, strand))
if not rows:
    raise WorkflowError(f"No genes parsed from {snakemake.input.proteins}.")
with open(snakemake.output.table, "w", encoding="utf-8") as handle:
    handle.write("protein_id\tcontig\tstart\tend\tstrand\n")
    for row in rows:
        handle.write("\t".join(row) + "\n")
