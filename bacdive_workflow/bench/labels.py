"""Label mapping helpers shared by selection and accuracy analysis."""

import csv
from pathlib import Path

import yaml

ACCURACY_TRAITS = [
    "gram-positive",
    "motile2+",
    "anaerobic",
    "aerobic",
    "thermophile",
    "spore-forming",
]

STABILITY_ONLY = ["acidophile", "psychrophile"]


def load_mapping(path):
    return yaml.safe_load(Path(path).read_text())


def apply_mapping(mapping, record):
    """Apply per-trait field mappings to one BacDive record; conflicts exclude."""
    labels = {}
    for trait, spec in mapping.get("traits", {}).items():
        votes = set()
        for field in spec.get("fields", []):
            raw_values = record.get(field, [])
            if not isinstance(raw_values, list):
                raw_values = [raw_values]
            for raw in raw_values:
                decision = spec.get("values", {}).get(str(raw).strip(), "exclude")
                if decision != "exclude":
                    votes.add(decision)
        labels[trait] = next(iter(votes)) if len(votes) == 1 else None
    return labels


def read_genomes_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
