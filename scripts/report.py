"""Snakemake adapter for the standalone report builder."""

import contextlib
import sys

sys.path.insert(0, snakemake.params.root)
from bacdive_workflow.report import generate_report

with open(snakemake.log[0], "w") as log, contextlib.redirect_stdout(log):
    generate_report(snakemake.input.csv, snakemake.output.html, snakemake.output.summary)
    print("Wrote report and summary.")
