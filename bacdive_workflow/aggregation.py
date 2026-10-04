"""Validate complete prediction files and join metadata by canonical sample ID."""

import argparse
import csv
import io
import json
import math
from pathlib import Path

from .common import TRAITS, WorkflowError, atomic_text, index_metadata, normalize_id

METADATA_FIELDS = (
    "taxon",
    "gram_stain",
    "oxygen_tolerance",
    "completeness",
    "contamination",
    "data_source",
)


def aggregate(paths, traits=None, samples=None, metadata=None):
    traits = list(TRAITS) if traits is None else list(traits)
    if (
        not traits
        or len(set(traits)) != len(traits)
        or any(trait not in TRAITS for trait in traits)
    ):
        raise WorkflowError("Aggregation requires unique supported traits.")
    sample_rows = index_metadata(samples) if samples else {}
    reference_rows = index_metadata(metadata) if metadata else {}
    rows, seen = [], set()
    expected_labels = {TRAITS[trait] for trait in traits}
    for path in paths:
        try:
            data = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError) as error:
            raise WorkflowError(f"Cannot read prediction JSON {path}: {error}") from error
        if not isinstance(data, dict) or not isinstance(data.get("predictions"), dict):
            raise WorkflowError(f"{path}: missing predictions object.")
        sample_id = data.get("sample_id") or normalize_id(Path(path).name)
        if not isinstance(sample_id, str) or not sample_id or sample_id in seen:
            raise WorkflowError(f"{path}: missing or duplicate sample ID {sample_id!r}.")
        if sample_rows and sample_id not in sample_rows:
            raise WorkflowError(f"{path}: sample {sample_id} is absent from the sample sheet.")
        predictions = data["predictions"]
        if set(predictions) != expected_labels:
            raise WorkflowError(f"{path}: expected exactly {sorted(expected_labels)}.")
        row = {"sample_id": sample_id, "metadata_matched": sample_id in reference_rows}
        for field in METADATA_FIELDS:
            row[field] = reference_rows.get(sample_id, {}).get(field) or sample_rows.get(
                sample_id, {}
            ).get(field, "")
        for trait in traits:
            label = TRAITS[trait]
            prediction = predictions[label]
            if not isinstance(prediction, dict) or type(prediction.get("prediction")) is not bool:
                raise WorkflowError(f"{path}: {label} prediction must be a JSON boolean.")
            confidence = prediction.get("confidence")
            if (
                type(confidence) not in (int, float)
                or not math.isfinite(confidence)
                or not 50 <= confidence <= 100
            ):
                raise WorkflowError(f"{path}: {label} confidence must be between 50 and 100.")
            positive = prediction.get("positive_probability")
            if positive is None:
                positive = confidence / 100 if prediction["prediction"] else 1 - confidence / 100
            if (
                type(positive) not in (int, float)
                or not math.isfinite(positive)
                or not 0 <= positive <= 1
            ):
                raise WorkflowError(f"{path}: invalid {label} positive probability.")
            probability_of_label = positive if prediction["prediction"] else 1 - positive
            if abs(probability_of_label * 100 - confidence) > 0.011:
                raise WorkflowError(f"{path}: inconsistent {label} confidence and probability.")
            row[f"{trait}_prediction"] = prediction["prediction"]
            row[f"{trait}_confidence"] = confidence
            row[f"{trait}_positive_probability"] = round(positive, 8)
        rows.append(row)
        seen.add(sample_id)
    if not rows:
        raise WorkflowError("No prediction files were supplied.")
    if sample_rows and seen != set(sample_rows):
        raise WorkflowError(
            "Predictions are missing for samples: " + ", ".join(sorted(set(sample_rows) - seen))
        )
    return sorted(rows, key=lambda row: row["sample_id"])


def save_csv(rows, path):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    atomic_text(path, buffer.getvalue())


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Combine complete BacDive-AI prediction JSON files."
    )
    parser.add_argument("predictions", nargs="+", type=Path)
    parser.add_argument("--samples", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--traits", nargs="+", choices=list(TRAITS), default=list(TRAITS))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(argv)
    try:
        rows = aggregate(args.predictions, args.traits, args.samples, args.metadata)
        save_csv(rows, args.output)
        print(
            f"Saved {len(rows)} complete samples; {sum(row['metadata_matched'] for row in rows)} reference matches."
        )
        return 0
    except (WorkflowError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
