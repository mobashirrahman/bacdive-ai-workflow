"""Snakemake adapter; aggregation itself is also available as a CLI."""

import contextlib
import sys

sys.path.insert(0, snakemake.params.root)
from bacdive_workflow.aggregation import aggregate, save_csv

with open(snakemake.log[0], "w") as log, contextlib.redirect_stdout(log):
    rows = aggregate(
        snakemake.input.predictions,
        snakemake.params.traits,
        snakemake.input.samples,
        snakemake.input.metadata[0] if snakemake.input.metadata else None,
    )
    save_csv(rows, snakemake.output.csv)
    print(
        f"Saved {len(rows)} complete samples; {sum(r['metadata_matched'] for r in rows)} reference matches."
    )
