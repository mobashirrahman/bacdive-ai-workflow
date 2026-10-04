"""Collect per-genome prediction JSON into a long-format table."""

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.common import TRAITS

LABEL_TO_TRAIT = {label: trait for trait, label in TRAITS.items()}

rows = []
for path in sorted(snakemake.input.predictions):
    payload = json.loads(Path(path).read_text())
    genome = payload["sample_id"]
    for label, entry in payload["predictions"].items():
        rows.append(
            {
                "genome": genome,
                "trait": LABEL_TO_TRAIT[label],
                "prediction": int(entry["prediction"]),
                "positive_probability": entry["positive_probability"],
                "matched_features": entry["matched_features"],
            }
        )
Path(snakemake.output.table).parent.mkdir(parents=True, exist_ok=True)
with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["genome", "trait", "prediction", "positive_probability", "matched_features"],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
