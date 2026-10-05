"""Regenerate the Markdown tables embedded in docs/BENCHMARK.md.

Tables are generated from docs/benchmark/*.tsv between marker comments and
are never typed by hand. --check fails if the committed file is stale.
"""

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs" / "benchmark"
DOCUMENT = ROOT / "docs" / "BENCHMARK.md"


def read_tsv(name):
    with open(DOCS / name, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fmt(value, digits=3):
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


def accuracy_table():
    lines = [
        "| Trait | Set | n | Pos | Neg | Sens | Spec | Bal.acc [95% CI] | MCC | AUROC |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in read_tsv("accuracy.tsv"):
        if row["phylum"]:
            continue
        ci = (
            f"[{fmt(row['balanced_accuracy_low'])}, {fmt(row['balanced_accuracy_high'])}]"
            if row["balanced_accuracy"]
            else "—"
        )
        lines.append(
            f"| {row['trait']} | {row['group']} | {row['n']} | {row['positives']} | "
            f"{row['negatives']} | {fmt(row['sensitivity'])} | {fmt(row['specificity'])} | "
            f"{fmt(row['balanced_accuracy'])} {ci} | {fmt(row['mcc'])} | {fmt(row['auroc'])} |"
        )
    return "\n".join(lines) + "\n"


def thresholds_table():
    lines = [
        "| Trait | Lowest completeness, flip rate <5% | Bootstrap interval over genomes |",
        "| --- | ---: | ---: |",
    ]
    for row in read_tsv("degrade_thresholds.tsv"):
        lines.append(
            f"| {row['trait']} | {row['lowest_level_flip_below_5pct']}% | "
            f"[{row['bootstrap_low']}%, {row['bootstrap_high']}%] |"
        )
    return "\n".join(lines) + "\n"


def flip_table():
    lines = [
        "| Trait | L100 | L70 | L50 | L30 | Direction at L50 |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ]
    summary = read_tsv("degrade_summary.tsv")
    for trait in sorted({r["trait"] for r in summary}):
        cells = {}
        direction = ""
        for row in summary:
            if row["trait"] != trait or row["set"] != "unseen" or row["contamination"] != "0":
                continue
            cells[row["level"]] = fmt(row["flip_rate"])
            if row["level"] == "50":
                cells[row["level"]] = fmt(row["flip_rate"])
                direction = f"+→− {fmt(row['pos_to_neg'])}, −→+ {fmt(row['neg_to_pos'])}"
        lines.append(
            f"| {trait} | {cells.get('100', '—')} | {cells.get('70', '—')} | "
            f"{cells.get('50', '—')} | {cells.get('30', '—')} | {direction} |"
        )
    return "\n".join(lines) + "\n"


def drift_table():
    lines = [
        "| InterProScan | Genomes | Median Jaccard | Min | Genomes with prediction changes |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    by_version = {}
    for row in read_tsv("drift.tsv"):
        by_version.setdefault(row["ips_version"], []).append(row)
    for version in sorted(by_version):
        group = by_version[version]
        jaccards = sorted(float(r["jaccard"]) for r in group if r["jaccard"])
        changed = sum(1 for r in group if int(r["prediction_changes"]) > 0)
        lines.append(
            f"| {version} | {len(group)} | {statistics.median(jaccards):.3f} | "
            f"{min(jaccards):.3f} | {changed} |"
        )
    return "\n".join(lines) + "\n"


def robustness_table():
    lines = [
        "| Trait | Caller meta | Caller PGAP | Pfam 5.63 | E-value 1e-10 | E-value 1e-30 |",
        "| --- | ---: | ---: | ---: | ---: | ---: |",
    ]
    summary = read_tsv("robustness.tsv")
    for trait in sorted({r["trait"] for r in summary}):
        cells = {r["variant"]: fmt(r["agreement"]) for r in summary if r["trait"] == trait}
        lines.append(
            f"| {trait} | {cells.get('caller-meta', '—')} | {cells.get('caller-pgap', '—')} | "
            f"{cells.get('pfam-5.63', '—')} | {cells.get('evalue-1e-10', '—')} | "
            f"{cells.get('evalue-1e-30', '—')} |"
        )
    return "\n".join(lines) + "\n"


TABLES = {
    "accuracy": accuracy_table,
    "thresholds": thresholds_table,
    "flip": flip_table,
    "drift": drift_table,
    "robustness": robustness_table,
}


def regenerate(check=False):
    text = DOCUMENT.read_text()
    changed = []
    for name, builder in TABLES.items():
        pattern = re.compile(
            rf"<!-- benchmark:{name} -->\n(.*?)<!-- /benchmark:{name} -->", re.DOTALL
        )
        replacement = f"<!-- benchmark:{name} -->\n{builder()}<!-- /benchmark:{name} -->"
        updated, count = pattern.subn(replacement, text)
        if not count:
            raise SystemExit(f"Marker pair benchmark:{name} not found in {DOCUMENT}.")
        if updated != text:
            changed.append(name)
        text = updated
    if check:
        if changed:
            raise SystemExit(
                f"Stale benchmark tables: {', '.join(changed)}. Regenerate without --check."
            )
        print("BENCHMARK.md tables are current.")
    else:
        DOCUMENT.write_text(text)
        print(f"Regenerated tables: {', '.join(TABLES)}.")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Regenerate BENCHMARK.md tables.")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args(argv)
    regenerate(check=args.check)


if __name__ == "__main__":
    sys.exit(main())
