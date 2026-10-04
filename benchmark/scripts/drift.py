"""Analysis B adapter: Jaccard drift and prediction changes per genome and version."""

import csv
import pickle
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench import drift as drift_lib
from bacdive_workflow.common import TRAITS
from bacdive_workflow.prediction import parse_pfams, predict_one

evalue = float(snakemake.params.evalue)
model_dir = Path("results/models")
bundles = {}
for trait in TRAITS:
    with open(model_dir / f"{trait}_data.p", "rb") as handle:
        bundles[trait] = pickle.load(handle)


def predict_classes(pfams):
    return {trait: predict_one(bundles[trait], pfams)["prediction"] for trait in TRAITS}


def pfam_set(path):
    pfams, _ = parse_pfams(path, evalue)
    return pfams


def published(path):
    text = Path(path).read_text().strip()
    return set(text.split()) if text else set()


match_of = {}
with open("results/published/accession_match.tsv", newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        match_of[row["genome"]] = (row["match_type"], int(row["n_published"]))

rows = []
gain_counter: Counter = Counter()
loss_counter: Counter = Counter()
for path in snakemake.input.primary + snakemake.input.legacy:
    parts = Path(path).parts
    ips = parts[parts.index("interpro") + 1]
    genome = Path(path).stem
    ours = pfam_set(path)
    theirs = published(f"results/published/{genome}.pfams")
    match_type, _ = match_of.get(genome, ("missing", 0))
    ours_classes = predict_classes(ours) if ours else {}
    their_classes = predict_classes(theirs) if theirs else {}
    changes = sum(1 for t in TRAITS if ours_classes.get(t) != their_classes.get(t))
    rows.append(
        drift_lib.drift_row(genome, ips, ours, theirs, changes, match_type),
    )
    model_pfams = set()
    for bundle in bundles.values():
        model_pfams |= set(bundle["categories"])
    for pfam in set(ours) - set(theirs):
        if pfam in model_pfams:
            gain_counter[(ips, pfam)] += 1
    for pfam in set(theirs) - set(ours):
        if pfam in model_pfams:
            loss_counter[(ips, pfam)] += 1

with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "genome",
            "ips_version",
            "jaccard",
            "n_ours",
            "n_published",
            "n_gained",
            "n_lost",
            "prediction_changes",
            "accession_match",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)

with open(snakemake.output.pfams, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle, fieldnames=["ips_version", "pfam", "direction", "genomes"], delimiter="\t"
    )
    writer.writeheader()
    for (ips, pfam), count in gain_counter.most_common(50):
        writer.writerow({"ips_version": ips, "pfam": pfam, "direction": "gained", "genomes": count})
    for (ips, pfam), count in loss_counter.most_common(50):
        writer.writerow({"ips_version": ips, "pfam": pfam, "direction": "lost", "genomes": count})
