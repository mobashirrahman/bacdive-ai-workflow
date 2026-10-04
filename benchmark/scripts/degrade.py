"""Analysis C adapter: fragment-loss simulation with predictions per draw.

Level 100 with zero contamination keeps every gene, so it must reproduce the
full-genome prediction; the adapter asserts that.
"""

import csv
import json
import pickle
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench import fragments as frag
from bacdive_workflow.common import TRAITS, WorkflowError
from bacdive_workflow.prediction import predict_one

params = snakemake.params
genome = params.genome
evalue = float(params.evalue)


def read_genes(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return [
            {
                "protein_id": r["protein_id"],
                "contig": r["contig"],
                "start": int(r["start"]),
                "end": int(r["end"]),
            }
            for r in csv.DictReader(handle, delimiter="\t")
        ]


def contig_lengths(path):
    lengths = {}
    name = None
    count = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = count
                name = line[1:].split()[0]
                count = 0
            else:
                count += len(line.strip())
    if name is not None:
        lengths[name] = count
    return lengths


def protein_pfams(path):
    mapping: dict = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.reader(handle, delimiter="\t"):
            if not row or not any(c.strip() for c in row):
                continue
            if len(row) < 11 or row[3].strip().lower() != "pfam":
                continue
            try:
                score = float(row[8])
            except ValueError:
                continue
            if score <= evalue:
                mapping.setdefault(row[0], set()).add(row[4].strip())
    return mapping


genes = read_genes(snakemake.input.genes)
lengths = contig_lengths(snakemake.input.fna)
hits = protein_pfams(snakemake.input.tsv)
gene_pfams = {g["protein_id"]: hits.get(g["protein_id"], set()) for g in genes}

genome_rows = {}
with open(snakemake.input.genomes, newline="", encoding="utf-8") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        genome_rows[row["genome_id"]] = row
phylum = genome_rows[genome]["phylum"]
donors = sorted(g for g, r in genome_rows.items() if r["phylum"] != phylum and g != genome)
donor = donors[frag.draw_seed(params.seed, genome, "donor", 0, 0) % len(donors)] if donors else None
by_id = dict(zip(params.genome_ids, zip(snakemake.input.all_genes, snakemake.input.all_tsvs)))
donor_genes = read_genes(by_id[donor][0]) if donor else []
donor_hits = protein_pfams(by_id[donor][1]) if donor else {}
donor_pfams = {g["protein_id"]: donor_hits.get(g["protein_id"], set()) for g in donor_genes}
for g in donor_genes:
    g["length"] = g["end"] - g["start"]

model_dir = Path(params.model_dir)
bundles = {}
for trait in TRAITS:
    with open(model_dir / f"{trait}_data.p", "rb") as handle:
        bundles[trait] = pickle.load(handle)

full = json.loads(Path(snakemake.input.full).read_text())
full_classes = {trait: full["predictions"][label]["prediction"] for trait, label in TRAITS.items()}

levels = list(params.levels)
replicates = int(params.replicates)
contamination_levels = set(params.contamination_levels)
fp = {
    "fragment_median_bp": params.fragment_median_bp,
    "fragment_sigma": params.fragment_sigma,
    "min_fragment_bp": params.min_fragment_bp,
}

rows = []
for level in levels:
    contaminations = list(params.contamination) if level in contamination_levels else [0]
    for contamination in contaminations:
        for replicate in range(replicates):
            seed = frag.draw_seed(params.seed, genome, level, contamination, replicate)
            kept_ids, retained_bp, foreign_bp, foreign, _ = frag.simulate(
                genes,
                lengths,
                level,
                contamination,
                donor_genes if contamination else None,
                fp,
                seed,
            )
            pfams = set()
            for pid in kept_ids:
                pfams |= gene_pfams.get(pid, set())
            for gene in foreign:
                pfams |= donor_pfams.get(gene["protein_id"], set())
            if (
                level >= 100
                and not contamination
                and set(kept_ids) != {g["protein_id"] for g in genes}
            ):
                raise WorkflowError(f"{genome}: level-100 draw lost genes.")
            for trait in TRAITS:
                try:
                    result = predict_one(bundles[trait], pfams)
                    prediction = int(result["prediction"])
                    probability = result["positive_probability"]
                except WorkflowError:
                    prediction, probability = "", ""
                rows.append(
                    {
                        "genome": genome,
                        "level": level,
                        "contamination": contamination,
                        "replicate": replicate,
                        "seed": seed,
                        "retained_bp": retained_bp,
                        "kept_genes": len(kept_ids),
                        "n_pfams": len(pfams),
                        "trait": trait,
                        "prediction": prediction,
                        "positive_probability": probability,
                    }
                )
            if level >= 100 and not contamination:
                for trait in TRAITS:
                    drawn = predict_one(bundles[trait], pfams)["prediction"]
                    if bool(drawn) != bool(full_classes[trait]):
                        raise WorkflowError(
                            f"{genome}: level-100 shortcut disagrees with full prediction."
                        )

with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "genome",
            "level",
            "contamination",
            "replicate",
            "seed",
            "retained_bp",
            "kept_genes",
            "n_pfams",
            "trait",
            "prediction",
            "positive_probability",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
