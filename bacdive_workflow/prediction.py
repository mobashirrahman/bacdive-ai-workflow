"""Strict Pfam parsing and inference with the upstream presence/absence encoding."""

import csv
import importlib.metadata
import logging
import math
import pickle
import re
import warnings
from pathlib import Path

from . import __version__
from .common import DEFAULT_EVALUE, TRAITS, WorkflowError, normalize_id, sha256_file

LOGGER = logging.getLogger("bacdive_workflow")


def parse_pfams(path, evalue=DEFAULT_EVALUE):
    if not math.isfinite(evalue) or evalue <= 0:
        raise WorkflowError("The E-value threshold must be finite and positive.")
    pfams = set()
    counts = {"rows": 0, "pfam_rows": 0, "retained_rows": 0, "other_analysis_rows": 0}
    try:
        with Path(path).open(encoding="utf-8-sig", newline="") as handle:
            for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), start=1):
                if not row or not any(cell.strip() for cell in row):
                    continue
                counts["rows"] += 1
                if len(row) < 11:
                    raise WorkflowError(
                        f"{path}:{line_number}: expected at least 11 InterProScan TSV columns."
                    )
                if row[3].strip().lower() != "pfam":
                    counts["other_analysis_rows"] += 1
                    continue
                counts["pfam_rows"] += 1
                accession = row[4].strip()
                if not re.fullmatch(r"PF\d{5}", accession):
                    raise WorkflowError(
                        f"{path}:{line_number}: invalid Pfam accession {accession!r}."
                    )
                try:
                    score = float(row[8])
                except ValueError as error:
                    raise WorkflowError(
                        f"{path}:{line_number}: invalid Pfam E-value {row[8]!r}."
                    ) from error
                if not math.isfinite(score) or score < 0:
                    raise WorkflowError(
                        f"{path}:{line_number}: E-value must be finite and nonnegative."
                    )
                if score <= evalue:
                    pfams.add(accession)
                    counts["retained_rows"] += 1
    except OSError as error:
        raise WorkflowError(f"Cannot read annotations {path}: {error}") from error
    if not pfams:
        raise WorkflowError(f"{path}: no Pfam matches passed the E-value threshold {evalue:g}.")
    counts["unique_pfams"] = len(pfams)
    LOGGER.info("Retained %d unique Pfams from %d rows", len(pfams), counts["rows"])
    return pfams, counts


def predict_one(bundle, pfams):
    import numpy as np

    if not isinstance(bundle, dict) or "model" not in bundle or "categories" not in bundle:
        raise WorkflowError("Model bundle must contain model and categories.")
    model = bundle["model"]
    categories = list(bundle["categories"])
    if not categories or len(set(categories)) != len(categories):
        raise WorkflowError("Model feature categories are empty or duplicated.")
    # Upstream BacDive-AI encodes each Pfam as a presence flag. Its feature line
    # reads `dict(zip(*np.unique(list(pfams), return_counts=True)))`, which looks
    # like hit counts but is a no-op: `pfams` is already a set, so the list has
    # no duplicates and every count is 1. Verified against models/predict.py.
    # Do not "fix" this into real counts; that changes predictions.
    vector = [int(category in pfams) for category in categories]
    matched = sum(vector)
    if not matched:
        raise WorkflowError("No input Pfams overlap this model's feature categories.")
    classes = list(model.classes_)
    if len(classes) != 2 or set(classes) != {0, 1}:
        raise WorkflowError(f"Expected binary classes 0/1, received {classes!r}.")
    probabilities = np.asarray(model.predict_proba([vector]), dtype=float)
    if probabilities.shape != (1, 2) or not np.isfinite(probabilities).all():
        raise WorkflowError("Model returned invalid probability dimensions or values.")
    if (
        (probabilities < 0).any()
        or (probabilities > 1).any()
        or not np.isclose(probabilities.sum(), 1)
    ):
        raise WorkflowError("Model probabilities must be in [0, 1] and sum to one.")
    index = int(probabilities[0].argmax())
    return {
        "prediction": bool(classes[index]),
        "confidence": round(float(probabilities[0, index]) * 100, 2),
        "positive_probability": round(float(probabilities[0, classes.index(1)]), 8),
        "matched_features": matched,
        "model_features": len(categories),
    }


def predict_file(
    path, selected, model_dir, evalue=DEFAULT_EVALUE, sample_id=None, model_cache=None
):
    from sklearn.exceptions import InconsistentVersionWarning

    from .models import manifest

    if (
        not selected
        or len(set(selected)) != len(selected)
        or any(trait not in TRAITS for trait in selected)
    ):
        raise WorkflowError("Select a nonempty list of unique supported traits.")
    model_paths = {trait: Path(model_dir) / f"{trait}_data.p" for trait in selected}
    missing = [str(p) for p in model_paths.values() if not p.is_file()]
    if missing:
        raise WorkflowError("Required model files are missing: " + ", ".join(missing))
    pfams, counts = parse_pfams(path, evalue)
    output = {
        "schema_version": 2,
        "sample_id": sample_id or normalize_id(path),
        "genome_file": Path(path).name,
        "predictions": {},
        "provenance": {
            "workflow_version": __version__,
            "annotation_sha256": sha256_file(path),
            "evalue_threshold": evalue,
            "feature_encoding": "binary Pfam presence/absence (upstream BacDive-AI encoding)",
            "annotation_counts": counts,
            "software": {
                name: importlib.metadata.version(name)
                for name in ("numpy", "scipy", "scikit-learn")
            },
            "model_sha256": {},
        },
    }
    for trait, model_path in model_paths.items():
        try:
            signature = (
                str(model_path.resolve()),
                model_path.stat().st_size,
                model_path.stat().st_mtime_ns,
            )
            cached = model_cache.get(signature) if model_cache is not None else None
            if cached is None:
                checksum = sha256_file(model_path)
                if checksum != manifest()["models"][model_path.name]["sha256"]:
                    raise WorkflowError(
                        "Model checksum does not match the pinned upstream release."
                    )
                with warnings.catch_warnings():
                    warnings.simplefilter("error", InconsistentVersionWarning)
                    with model_path.open("rb") as handle:
                        bundle = pickle.load(handle)
                if model_cache is not None:
                    model_cache[signature] = (bundle, checksum)
            else:
                bundle, checksum = cached
            output["predictions"][TRAITS[trait]] = predict_one(bundle, pfams)
            output["provenance"]["model_sha256"][trait] = checksum
            LOGGER.info("Predicted %s", TRAITS[trait])
        except Exception as error:
            raise WorkflowError(
                f"Cannot predict {trait} with {model_path.name}: {error}"
            ) from error
    return output
