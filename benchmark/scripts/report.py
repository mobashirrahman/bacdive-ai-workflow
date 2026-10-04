"""Report adapter: copy final tables and draw the four benchmark figures."""

import csv
import shutil
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench import figures


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


by_name = {Path(p).name: p for p in snakemake.input.tables}
accuracy = read_tsv(by_name["accuracy.tsv"])
for row in accuracy:
    for key in ("balanced_accuracy", "balanced_accuracy_low", "balanced_accuracy_high"):
        row[key] = as_float(row.get(key))
    row["n"] = int(row["n"] or 0)
degrade = read_tsv(by_name["degrade_summary.tsv"])
flip = [r for r in degrade if r["contamination"] == "0"]
for row in flip:
    row["level"] = int(row["level"])
    row["flip_rate"] = as_float(row["flip_rate"]) or 0.0
drift = read_tsv(by_name["drift.tsv"])
for row in drift:
    row["jaccard"] = as_float(row.get("jaccard"))
robust = read_tsv(by_name["robustness.tsv"])
for row in robust:
    row["agreement"] = as_float(row.get("agreement"))

figures_dir = Path(snakemake.output.figures[0]).parent
figures_dir.mkdir(parents=True, exist_ok=True)
figures.accuracy_figure(accuracy, figures_dir / "accuracy.png", figures_dir / "accuracy.svg")
figures.flip_figure(flip, figures_dir / "flip_rate.png", figures_dir / "flip_rate.svg")
figures.drift_figure(drift, figures_dir / "drift.png", figures_dir / "drift.svg")
figures.robustness_heatmap(robust, figures_dir / "robustness.png", figures_dir / "robustness.svg")

docs_tables = list(snakemake.output.docs_tables)
docs_figures = list(snakemake.output.docs_figures)
for source, dest in zip(sorted(snakemake.input.tables), sorted(docs_tables)):
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
for source, dest in zip(sorted(snakemake.output.figures), sorted(docs_figures)):
    Path(dest).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, dest)
