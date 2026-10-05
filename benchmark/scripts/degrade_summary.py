"""Degrade summaries: flip rates, accuracy, probability shift, 5% thresholds."""

import csv
import random
import sys
from pathlib import Path

import yaml

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench.metrics import (
    MIN_CLASS_MEMBERS,
    flip_rates,
    lowest_level_below,
    write_tsv,
)

BOOTSTRAPS = 200

full: dict = {}
with open("results/tables/predictions_full.tsv", newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        full[(row["genome"], row["trait"])] = (
            int(row["prediction"]),
            float(row["positive_probability"]),
        )

genomes = {}
with open(snakemake.input.genomes, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        labels = {
            c[len("label_") :]: int(v)
            for c, v in row.items()
            if c.startswith("label_") and v.strip() not in ("", "NA")
        }
        genomes[row["genome_id"]] = {"set": row["set"], "labels": labels}

mapping = yaml.safe_load(open(snakemake.input.mapping, encoding="utf-8"))
eligible = set(mapping.get("accuracy_traits", []))

# Per (trait, level, contamination) and genome: draws, full-positive draws,
# positives lost, negatives gained, labelled draws, correct draws, summed shift.
cells: dict = {}
for table in snakemake.input.tables:
    with open(table, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["prediction"] == "":
                continue
            trait, genome = row["trait"], row["genome"]
            full_pred, full_prob = full[(genome, trait)]
            pred = int(row["prediction"])
            cell = cells.setdefault(
                (trait, int(row["level"]), int(row["contamination"])), {}
            ).setdefault(genome, [0, 0, 0, 0, 0, 0, 0.0])
            cell[0] += 1
            cell[1] += full_pred
            cell[2] += full_pred == 1 and pred == 0
            cell[3] += full_pred == 0 and pred == 1
            label = genomes[genome]["labels"].get(trait)
            if trait in eligible and label is not None:
                cell[4] += 1
                cell[5] += pred == label
            cell[6] += abs(float(row["positive_probability"]) - full_prob)


def pooled(per_genome, members):
    """Sum the per-genome counters over `members` (repeats count repeatedly)."""
    totals = [0, 0, 0, 0, 0, 0, 0.0]
    for genome in members:
        for index, value in enumerate(per_genome.get(genome, ())):
            totals[index] += value
    return totals


summary = []
for (trait, level, contamination), per_genome in sorted(cells.items()):
    for group in sorted({genomes[g]["set"] for g in per_genome}):
        members = [g for g in per_genome if genomes[g]["set"] == group]
        draws, positives, lost, gained, labelled, correct, shift = pooled(per_genome, members)
        summary.append(
            {
                "trait": trait,
                "level": level,
                "contamination": contamination,
                "set": group,
                "genomes": len(members),
                **flip_rates(draws, positives, lost, gained),
                "full_positive_genomes": sum(1 for g in members if per_genome[g][1]),
                "accuracy": correct / labelled if labelled else "",
                "mean_prob_shift": shift / draws if draws else "",
            }
        )

write_tsv(
    Path(snakemake.output.summary),
    summary,
    [
        "trait",
        "level",
        "contamination",
        "set",
        "genomes",
        "draws",
        "flip_rate",
        "pos_to_neg",
        "neg_to_pos",
        "full_positive_genomes",
        "full_positive_draws",
        "positive_loss_rate",
        "full_negative_draws",
        "negative_gain_rate",
        "accuracy",
        "mean_prob_shift",
    ],
)


# Per-trait lowest completeness with a rate below 5% (no contamination, both
# sets pooled), with a bootstrap interval over genomes. Two rates are reported:
# flips over all genomes, and losses among genomes predicted positive on the
# full genome. The second is left empty below MIN_CLASS_MEMBERS positives.
def threshold(trait, members, rate):
    by_level = {}
    for (cell_trait, level, contamination), per_genome in cells.items():
        if cell_trait == trait and contamination == 0:
            totals = pooled(per_genome, members)
            by_level[level] = flip_rates(*totals[:4])[rate]
    return lowest_level_below(by_level)


def interval(trait, members, rate, rng):
    boots = []
    for _ in range(BOOTSTRAPS):
        pick = [rng.choice(members) for _ in members]
        value = threshold(trait, pick, rate)
        if value != "":
            boots.append(value)
    boots.sort()
    if not boots:
        return "", ""
    return boots[int(0.025 * len(boots))], boots[min(len(boots) - 1, int(0.975 * len(boots)))]


thresholds = []
for trait in sorted({key[0] for key in cells}):
    everyone = sorted({g for key, per in cells.items() if key[0] == trait for g in per})
    positive = [g for g in everyone if full[(g, trait)][0] == 1]
    rng = random.Random(snakemake.params.seed)
    low, high = interval(trait, everyone, "flip_rate", rng)
    row = {
        "trait": trait,
        "genomes": len(everyone),
        "lowest_level_flip_below_5pct": threshold(trait, everyone, "flip_rate"),
        "bootstrap_low": low,
        "bootstrap_high": high,
        "full_positive_genomes": len(positive),
        "lowest_level_positive_loss_below_5pct": "",
        "positive_bootstrap_low": "",
        "positive_bootstrap_high": "",
    }
    if len(positive) >= MIN_CLASS_MEMBERS:
        low, high = interval(trait, positive, "positive_loss_rate", rng)
        row["lowest_level_positive_loss_below_5pct"] = threshold(
            trait, positive, "positive_loss_rate"
        )
        row["positive_bootstrap_low"] = low
        row["positive_bootstrap_high"] = high
    thresholds.append(row)

write_tsv(Path(snakemake.output.thresholds), thresholds, list(thresholds[0]))
