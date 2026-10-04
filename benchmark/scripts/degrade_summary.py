"""Degrade summaries: flip rates, accuracy, probability shift, 5% thresholds."""

import csv
import random
import sys

import yaml

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

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

draws = []
for table in snakemake.input.tables:
    with open(table, newline="", encoding="utf-8") as handle:
        draws.extend(csv.DictReader(handle, delimiter="\t"))

summary = []
grouped: dict = {}
for row in draws:
    grouped.setdefault(
        (row["trait"], row["level"], row["contamination"], genomes[row["genome"]]["set"]), []
    ).append(row)

for (trait, level, contamination, group), rows in sorted(
    grouped.items(), key=lambda kv: str(kv[0])
):
    total = len(rows)
    flips = pos_to_neg = neg_to_pos = 0
    correct = eligible_total = 0
    shifts = []
    for row in rows:
        if row["prediction"] == "":
            continue
        full_pred, full_prob = full[(row["genome"], trait)]
        pred = int(row["prediction"])
        if pred != full_pred:
            flips += 1
            if full_pred == 1:
                pos_to_neg += 1
            else:
                neg_to_pos += 1
        shifts.append(abs(float(row["positive_probability"]) - full_prob))
        label = genomes[row["genome"]]["labels"].get(trait)
        if trait in eligible and label is not None:
            eligible_total += 1
            correct += pred == label
    summary.append(
        {
            "trait": trait,
            "level": level,
            "contamination": contamination,
            "set": group,
            "draws": total,
            "flip_rate": flips / total if total else "",
            "pos_to_neg": pos_to_neg / total if total else "",
            "neg_to_pos": neg_to_pos / total if total else "",
            "accuracy": correct / eligible_total if eligible_total else "",
            "mean_prob_shift": sum(shifts) / len(shifts) if shifts else "",
        }
    )

with open(snakemake.output.summary, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "trait",
            "level",
            "contamination",
            "set",
            "draws",
            "flip_rate",
            "pos_to_neg",
            "neg_to_pos",
            "accuracy",
            "mean_prob_shift",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(summary)

# Per-trait lowest completeness with flip rate below 5% (c=0), with a
# bootstrap interval over genomes.
thresholds = []
for trait in sorted({r["trait"] for r in draws}):
    levels = sorted({int(r["level"]) for r in draws if r["contamination"] == "0"}, reverse=True)
    chosen = ""
    for level in levels:
        subset = [
            r
            for r in draws
            if r["trait"] == trait
            and r["level"] == str(level)
            and r["contamination"] == "0"
            and r["prediction"] != ""
        ]
        if not subset:
            continue
        rate = sum(int(r["prediction"]) != full[(r["genome"], trait)][0] for r in subset) / len(
            subset
        )
        if rate < 0.05:
            chosen = level
        else:
            break
    genomes_for_trait = sorted({r["genome"] for r in draws if r["trait"] == trait})
    rng = random.Random(snakemake.params.seed)
    boots = []
    for _ in range(200):
        pick = [rng.choice(genomes_for_trait) for _ in genomes_for_trait]
        best = ""
        for level in levels:
            subset = [
                r
                for r in draws
                if r["trait"] == trait
                and r["level"] == str(level)
                and r["contamination"] == "0"
                and r["prediction"] != ""
                and r["genome"] in pick
            ]
            if not subset:
                continue
            rate = sum(int(r["prediction"]) != full[(r["genome"], trait)][0] for r in subset) / len(
                subset
            )
            if rate < 0.05:
                best = level
            else:
                break
        if best != "":
            boots.append(best)
    boots.sort()
    thresholds.append(
        {
            "trait": trait,
            "lowest_level_flip_below_5pct": chosen,
            "bootstrap_low": boots[int(0.025 * len(boots))] if boots else "",
            "bootstrap_high": boots[min(len(boots) - 1, int(0.975 * len(boots)))] if boots else "",
        }
    )

with open(snakemake.output.thresholds, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["trait", "lowest_level_flip_below_5pct", "bootstrap_low", "bootstrap_high"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(thresholds)
