"""Rebuild the published descriptive case study from portable saved predictions."""

import argparse
import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bacdive_workflow.aggregation import aggregate, save_csv
from bacdive_workflow.report import generate_report


def reproduce(destination):
    case = ROOT / "examples/case-study"
    rows = aggregate(
        sorted((case / "predictions").glob("*.json")),
        samples=case / "samples.tsv",
        metadata=case / "metadata.tsv",
    )
    save_csv(rows, destination / "predictions.csv")
    generate_report(
        destination / "predictions.csv",
        destination / "report.html",
        destination / "summary.json",
        "120-genome case study",
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=ROOT / "docs/results")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        with tempfile.TemporaryDirectory() as directory:
            temporary = Path(directory)
            reproduce(temporary)
            for name in ("predictions.csv", "summary.json", "report.html"):
                if (temporary / name).read_bytes() != (ROOT / "docs/results" / name).read_bytes():
                    raise SystemExit(f"Published artifact is out of date: {name}")
        print("Published CSV, summary, and report reproduce exactly.")
    else:
        args.output_dir.mkdir(parents=True, exist_ok=True)
        reproduce(args.output_dir)
        summary = json.loads((args.output_dir / "summary.json").read_text())
        print(
            f"Reproduced {summary['samples']} samples and {summary['metadata_matches']} reference matches."
        )


if __name__ == "__main__":
    main()
