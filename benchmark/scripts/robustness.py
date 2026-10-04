"""Analysis D adapter: E-value, gene-caller, and Pfam-release variants."""

import csv
import pickle
import sys
from pathlib import Path

import yaml

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.common import TRAITS, WorkflowError
from bacdive_workflow.prediction import parse_pfams, predict_one

params = snakemake.params
model_dir = Path(params.model_dir)
bundles = {}
for trait in TRAITS:
    with open(model_dir / f"{trait}_data.p", "rb") as handle:
        bundles[trait] = pickle.load(handle)


def variant_prediction(path, evalue):
    pfams, _ = parse_pfams(path, evalue)
    return {trait: predict_one(bundles[trait], pfams) for trait in TRAITS}


by_path = {}
for key in ("single", "meta", "legacy", "pgap"):
    for path in snakemake.input[key]:
        by_path[(key, Path(path).stem)] = path

genomes = {}
with open(snakemake.input.genomes, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        genomes[row["genome_id"]] = {
            c[len("label_") :]: int(v)
            for c, v in row.items()
            if c.startswith("label_") and v.strip() not in ("", "NA")
        }

mapping_name = ".test/label_mapping.yaml" if params.test_mode else "config/label_mapping.yaml"
mapping = yaml.safe_load(open(f"{params.repo}/benchmark/{mapping_name}", encoding="utf-8"))
eligible = set(mapping.get("accuracy_traits", []))

variants: dict = {}
for genome in params.genomes:
    primary = by_path[("single", genome)]
    base = variant_prediction(primary, float(params.primary_evalue))
    variants[(genome, "primary")] = base
    for threshold in params.evalues:
        if float(threshold) == float(params.primary_evalue):
            continue
        variants[(genome, f"evalue-{threshold:g}")] = variant_prediction(primary, float(threshold))
    for key, name in (("meta", "caller-meta"), ("legacy", "pfam-5.63"), ("pgap", "caller-pgap")):
        path = by_path.get((key, genome))
        if path is None:
            continue
        try:
            variants[(genome, name)] = variant_prediction(path, float(params.primary_evalue))
        except WorkflowError:
            continue

rows = []
for (genome, variant), results in sorted(variants.items(), key=lambda kv: str(kv[0])):
    base = variants[(genome, "primary")]
    for trait in TRAITS:
        agreement = int(results[trait]["prediction"]) == int(base[trait]["prediction"])
        difference = abs(
            results[trait]["positive_probability"] - base[trait]["positive_probability"]
        )
        label = genomes.get(genome, {}).get(trait)
        accurate = ""
        if trait in eligible and label is not None:
            accurate = int(results[trait]["prediction"]) == label
        rows.append(
            {
                "genome": genome,
                "variant": variant,
                "trait": trait,
                "agreement": agreement,
                "probability_difference": difference,
                "accuracy": accurate,
            }
        )

summary = []
for trait in TRAITS:
    for variant in sorted({v for _, v in variants}):
        subset = [r for r in rows if r["trait"] == trait and r["variant"] == variant]
        if not subset:
            continue
        labelled = [r for r in subset if r["accuracy"] != ""]
        summary.append(
            {
                "trait": trait,
                "variant": variant,
                "genomes": len(subset),
                "agreement": sum(r["agreement"] for r in subset) / len(subset),
                "mean_probability_difference": sum(r["probability_difference"] for r in subset)
                / len(subset),
                "accuracy": (sum(r["accuracy"] for r in labelled) / len(labelled))
                if labelled
                else "",
            }
        )

with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "trait",
            "variant",
            "genomes",
            "agreement",
            "mean_probability_difference",
            "accuracy",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(summary)
