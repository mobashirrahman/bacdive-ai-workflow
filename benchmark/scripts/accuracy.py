"""Analysis A adapter: accuracy rows plus the disagreement list."""

import csv
import sys
from pathlib import Path

import yaml

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench import metrics

mapping = yaml.safe_load(open(snakemake.input.mapping, encoding="utf-8"))
eligible = set(mapping.get("accuracy_traits", []))

genomes = {}
with open(snakemake.input.genomes, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        labels = {}
        for column, value in row.items():
            if column.startswith("label_") and value.strip() not in ("", "NA"):
                labels[column[len("label_") :]] = int(value)
        genomes[row["genome_id"]] = {
            "set": row["set"],
            "genus_unseen": row["genus_unseen"].strip().lower() == "true",
            "phylum": row.get("phylum", ""),
            "species": row.get("species", ""),
            "labels": labels,
        }

predictions: dict = {}
scores: dict = {}
classes: dict = {}
with open(snakemake.input.predictions, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        predictions.setdefault(row["genome"], {})[row["trait"]] = int(row["prediction"])
        scores.setdefault(row["genome"], {})[row["trait"]] = float(row["positive_probability"])

traits = sorted(
    {t for labels in (g["labels"] for g in genomes.values()) for t in labels} | eligible
)
records = []
for genome, info in genomes.items():
    records.append(
        {
            "genome": genome,
            "set": info["set"],
            "genus_unseen": info["genus_unseen"],
            "phylum": info["phylum"],
            "labels": info["labels"],
            "predictions": predictions.get(genome, {}),
            "scores": scores.get(genome, {}),
        }
    )

groups = {
    "seen": lambda r: r["set"] == "seen",
    "unseen": lambda r: r["set"] == "unseen",
    "genus_unseen": lambda r: r["set"] == "unseen" and r["genus_unseen"],
}
rows = metrics.accuracy_rows(
    records, traits, groups, snakemake.params.bootstrap, snakemake.params.seed
)
metrics.write_tsv(
    Path(snakemake.output.accuracy),
    rows,
    [
        "trait",
        "group",
        "phylum",
        "n",
        "positives",
        "negatives",
        "sensitivity",
        "specificity",
        "balanced_accuracy",
        "mcc",
        "auroc",
        "balanced_accuracy_low",
        "balanced_accuracy_high",
    ],
)

disagreements = []
with open(snakemake.input.predictions, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        label = genomes.get(row["genome"], {}).get("labels", {}).get(row["trait"])
        if label is not None and int(row["prediction"]) != label:
            info = genomes[row["genome"]]
            disagreements.append(
                {
                    "genome": row["genome"],
                    "taxon": info["species"],
                    "trait": row["trait"],
                    "label": label,
                    "positive_probability": row["positive_probability"],
                }
            )
metrics.write_tsv(
    Path(snakemake.output.disagreements),
    disagreements,
    ["genome", "taxon", "trait", "label", "positive_probability"],
)
