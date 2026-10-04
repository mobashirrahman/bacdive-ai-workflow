"""Calibrate the BacDive-to-trait label mapping from training data.

A field value maps to a class only with >=95% purity over >=30 training
strains; other values map to exclude. A trait enters accuracy analysis only if
the mapping reproduces >=95% of training labels over >=500 strains. No
hand-tuning: thresholds are fixed, outcomes are reported whatever they are.

Command:
    python -m bacdive_workflow.bench.calibrate
        --training-labels <csv> --cache-dir <dir>
        --mapping-out benchmark/config/label_mapping.yaml
        --calibration-csv docs/benchmark/label_mapping_calibration.csv
"""

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

import yaml

from bacdive_workflow.bench.labels import (
    STABILITY_ONLY,
    TRAIT_FIELDS,
    extract_votes,
)
from bacdive_workflow.bench.select import (
    MIN_PURITY,
    MIN_REPRODUCIBILITY,
    MIN_STRAINS_FOR_ACCURACY,
    MIN_SUPPORT,
)
from bacdive_workflow.common import WorkflowError

CALIBRATION_TRAITS = [
    "gram-positive",
    "motile2+",
    "anaerobic",
    "aerobic",
    "thermophile",
    "spore-forming",
]


def load_training_labels(path):
    labels: dict = {}
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            labels.setdefault(row["ID_strains"], {})[row["trait"]] = int(row["ability"])
    return labels


def load_cached_records(cache_dir):
    records = {}
    for batch in sorted(Path(cache_dir, "train").glob("batch-*.json")):
        records.update(json.loads(batch.read_text()).get("results", {}))
    return records


def calibrate(training_labels, records):
    """Cross-tabulate votes vs labels; return (mapping, calibration_rows, scores)."""
    calibration_rows = []
    scores = {}
    traits = {}
    for trait in CALIBRATION_TRAITS:
        by_value: dict = {}
        strain_votes: dict = {}
        skipped_prediction = 0
        for strain, trait_labels in training_labels.items():
            if trait not in trait_labels:
                continue
            record = records.get(str(strain))
            if record is None:
                continue
            votes, skipped = extract_votes(record, trait)
            skipped_prediction += skipped
            strain_votes[strain] = votes
            for vote in votes:
                by_value.setdefault(vote, []).append(trait_labels[trait])
        values = {}
        for value, classes in sorted(by_value.items()):
            support = len(classes)
            positives = sum(classes)
            purity = max(positives, support - positives) / support
            if support < MIN_SUPPORT:
                decision = "exclude"
            elif positives / support >= MIN_PURITY:
                decision = 1
            elif positives / support <= 1 - MIN_PURITY:
                decision = 0
            else:
                decision = "exclude"
            if decision != "exclude":
                values[value] = decision
            calibration_rows.append(
                {
                    "trait": trait,
                    "value": value,
                    "n_total": support,
                    "n_label0": support - positives,
                    "n_label1": positives,
                    "purity_majority": round(purity, 4),
                    "decision": decision,
                }
            )
        covered = reproduced = 0
        for strain, votes in strain_votes.items():
            decisions = {values[vote] for vote in votes if vote in values}
            if len(decisions) != 1:
                continue
            covered += 1
            if next(iter(decisions)) == training_labels[strain][trait]:
                reproduced += 1
        fraction = reproduced / covered if covered else 0.0
        passes = covered >= MIN_STRAINS_FOR_ACCURACY and fraction >= MIN_REPRODUCIBILITY
        scores[trait] = {
            "covered": covered,
            "reproduced": reproduced,
            "fraction": round(fraction, 4),
            "passes_gate": passes,
            "skipped_prediction_refs": skipped_prediction,
        }
        traits[trait] = {
            "fields": ["|".join(path) for path in TRAIT_FIELDS[trait]],
            "values": values,
        }
    return traits, calibration_rows, scores


def accuracy_lists(scores):
    eligible = [t for t in CALIBRATION_TRAITS if scores[t]["passes_gate"]]
    stability = [t for t in CALIBRATION_TRAITS if not scores[t]["passes_gate"]] + STABILITY_ONLY
    return eligible, stability


def main(argv=None):
    parser = argparse.ArgumentParser(description="Calibrate the BacDive label mapping.")
    parser.add_argument("--training-labels", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--mapping-out", type=Path, required=True)
    parser.add_argument("--calibration-csv", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        training_labels = load_training_labels(args.training_labels)
        records = load_cached_records(args.cache_dir)
        if not records:
            raise WorkflowError(f"No cached training records in {args.cache_dir}.")
        traits, rows, scores = calibrate(training_labels, records)
        eligible, stability = accuracy_lists(scores)
        harvest_file = Path(args.cache_dir) / "harvest.json"
        harvest_date = (
            json.loads(harvest_file.read_text())["harvest_date"] if harvest_file.is_file() else ""
        )
        mapping = {
            "harvest_date": harvest_date,
            "thresholds": {"min_support": MIN_SUPPORT, "min_purity": MIN_PURITY},
            "accuracy_traits": eligible,
            "stability_only": stability,
            "per_trait": scores,
            "traits": traits,
        }
        args.mapping_out.parent.mkdir(parents=True, exist_ok=True)
        args.mapping_out.write_text(yaml.safe_dump(mapping, sort_keys=False))
        args.calibration_csv.parent.mkdir(parents=True, exist_ok=True)
        counter = Counter((r["trait"], r["decision"]) for r in rows)
        with open(args.calibration_csv, "w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "trait",
                    "value",
                    "n_total",
                    "n_label0",
                    "n_label1",
                    "purity_majority",
                    "decision",
                ],
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Calibrated {len(rows)} value rows {dict(counter)}.")
        print(f"Accuracy-eligible: {eligible}; stability-only: {stability}.")
        return 0
    except (WorkflowError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
