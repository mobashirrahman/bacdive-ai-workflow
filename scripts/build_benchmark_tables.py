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
        "| Trait | Positive genomes | Lowest completeness keeping >95% of positives "
        "| Bootstrap interval | Lowest completeness, all-genome flip rate <5% |",
        "| --- | ---: | ---: | ---: | ---: |",
    ]
    for row in read_tsv("degrade_thresholds.tsv"):
        if row["lowest_level_positive_loss_below_5pct"]:
            kept = f"{row['lowest_level_positive_loss_below_5pct']}%"
            interval = f"[{row['positive_bootstrap_low']}%, {row['positive_bootstrap_high']}%]"
        else:
            kept = interval = "—"
        lines.append(
            f"| {row['trait']} | {row['full_positive_genomes']} of {row['genomes']} | {kept} | "
            f"{interval} | {row['lowest_level_flip_below_5pct']}% |"
        )
    return "\n".join(lines) + "\n"


def flip_table():
    levels = ("90", "70", "50", "30")
    lines = [
        "| Trait | Positive genomes | "
        + " | ".join(f"Lost at {level}%" for level in levels)
        + " | Negatives turned positive at 50% |",
        "| --- | ---: | " + " | ".join("---:" for _ in levels) + " | ---: |",
    ]
    summary = [
        r
        for r in read_tsv("degrade_summary.tsv")
        if r["set"] == "unseen" and r["contamination"] == "0"
    ]
    for trait in sorted({r["trait"] for r in summary}):
        by_level = {r["level"]: r for r in summary if r["trait"] == trait}
        positives = by_level["100"]["full_positive_genomes"]
        genomes = by_level["100"]["genomes"]
        lost = " | ".join(fmt(by_level[level]["positive_loss_rate"], 2) for level in levels)
        lines.append(
            f"| {trait} | {positives} of {genomes} | {lost} | "
            f"{fmt(by_level['50']['negative_gain_rate'])} |"
        )
    return "\n".join(lines) + "\n"


def drift_table():
    lines = [
        "| InterProScan | Genomes | Median Jaccard | Min | Genomes with a changed prediction "
        "| Median Jaccard, published unfiltered | Changed, published unfiltered |",
        "| --- | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    by_version = {}
    for row in read_tsv("drift.tsv"):
        by_version.setdefault(row["ips_version"], []).append(row)
    for version in sorted(by_version):
        group = by_version[version]
        jaccards = [float(r["jaccard"]) for r in group if r["jaccard"]]
        unfiltered = [float(r["jaccard_unfiltered"]) for r in group if r["jaccard_unfiltered"]]
        changed = sum(1 for r in group if int(r["prediction_changes"]) > 0)
        changed_unfiltered = sum(1 for r in group if int(r["prediction_changes_unfiltered"]) > 0)
        lines.append(
            f"| {version} | {len(group)} | {statistics.median(jaccards):.3f} | "
            f"{min(jaccards):.3f} | {changed} | {statistics.median(unfiltered):.3f} | "
            f"{changed_unfiltered} |"
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
