#!/usr/bin/env python3
"""Summarize Snakemake benchmark records into a per-stage runtime table."""

import argparse
import csv
import io
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bacdive_workflow.common import WorkflowError, atomic_text  # noqa: E402

STAGES = ("prodigal", "interproscan", "predict_traits")
FIELDS = ("stage", "jobs", "mean_seconds", "median_seconds", "maximum_seconds", "total_job_seconds")


def collect(benchmarks):
    rows = []
    for stage in STAGES:
        directory = Path(benchmarks) / stage
        if not directory.is_dir():
            raise WorkflowError(f"Missing benchmark directory: {directory}")
        seconds = []
        for record in sorted(directory.glob("*.txt")):
            lines = [
                line for line in record.read_text(encoding="utf-8").splitlines() if line.strip()
            ]
            if len(lines) < 2:
                raise WorkflowError(f"{record}: expected a header and one measurement row.")
            fields = lines[0].split("\t")
            if "s" not in fields:
                raise WorkflowError(f"{record}: no seconds column in header.")
            value = lines[1].split("\t")[fields.index("s")]
            try:
                seconds.append(float(value))
            except ValueError as error:
                raise WorkflowError(f"{record}: invalid seconds value {value!r}.") from error
        if not seconds:
            raise WorkflowError(f"{directory}: no benchmark records found.")
        rows.append(
            {
                "stage": stage,
                "jobs": len(seconds),
                "mean_seconds": round(statistics.fmean(seconds), 2),
                "median_seconds": round(statistics.median(seconds), 2),
                "maximum_seconds": round(max(seconds), 2),
                "total_job_seconds": round(sum(seconds), 2),
            }
        )
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmarks", type=Path, default=ROOT / "benchmarks")
    parser.add_argument("--output", type=Path, default=ROOT / "docs/results/runtime.csv")
    parser.add_argument("--check", action="store_true", help="Fail if the published table is stale")
    args = parser.parse_args(argv)
    try:
        rows = collect(args.benchmarks)
    except (WorkflowError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(FIELDS))
    writer.writeheader()
    writer.writerows(rows)
    if args.check:
        # Compare bytes: csv writes CRLF, which read_text() would normalize away.
        expected = buffer.getvalue().encode("utf-8")
        if not args.output.is_file() or args.output.read_bytes() != expected:
            raise SystemExit(f"Published runtime table is out of date: {args.output}")
        print(f"Runtime table matches {len(rows)} benchmark stages.")
        return 0
    atomic_text(args.output, buffer.getvalue())
    print(f"Wrote runtime table for {len(rows)} stages to {args.output}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
