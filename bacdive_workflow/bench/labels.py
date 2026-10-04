"""BacDive record vote extraction and label-mapping application.

One vote is a normalized BacDive field value for a trait. A strain is labelled
only when all its votes agree; conflicting votes exclude it for that trait.
Entries citing prediction-sourced references are dropped, as is the whole
``Genome-based predictions`` section (BacDive-AI outputs, never labels).
"""

PREDICTION_TITLE_MARKERS = (
    "predicting bacterial phenotypic traits",
    "genome-based prediction",
    "bacdive-ai",
)

# Trait -> record paths (section, subsection, key) tabulated for calibration.
TRAIT_FIELDS = {
    "gram-positive": [
        ("Morphology", "cell morphology", "gram stain"),
    ],
    "motile2+": [
        ("Morphology", "cell morphology", "motility"),
        ("Morphology", "cell morphology", "flagellum arrangement"),
    ],
    "anaerobic": [
        ("Physiology and metabolism", "oxygen tolerance", "oxygen tolerance"),
    ],
    "aerobic": [
        ("Physiology and metabolism", "oxygen tolerance", "oxygen tolerance"),
    ],
    "thermophile": [
        ("Culture and growth conditions", "culture temp", "temperature"),
    ],
    "spore-forming": [
        ("Physiology and metabolism", "spore formation", "spore formation"),
    ],
}

ACCURACY_TRAITS = [
    "gram-positive",
    "motile2+",
    "anaerobic",
    "aerobic",
    "thermophile",
    "spore-forming",
]

STABILITY_ONLY = ["acidophile", "psychrophile"]


def _entries(value):
    if value is None:
        return []
    return value if isinstance(value, list) else [value]


def prediction_refs(record):
    """@ref ids whose Reference entry looks like a genome-based prediction."""
    flagged = set()
    for ref in record.get("Reference", []) or []:
        title = f"{ref.get('title', '')} {ref.get('authors', '')}".lower()
        if any(marker in title for marker in PREDICTION_TITLE_MARKERS):
            flagged.add(ref.get("@id"))
    return flagged


def _normalize(text):
    return " ".join(str(text).strip().lower().split())


def extract_votes(record, trait):
    """Normalized BacDive field values for one trait, sans prediction refs."""
    flagged = prediction_refs(record)
    votes = []
    skipped_prediction = 0
    for section, subsection, key in TRAIT_FIELDS.get(trait, []):
        container = record.get(section, {})
        if not isinstance(container, dict):
            continue
        for entry in _entries(container.get(subsection)):
            if not isinstance(entry, dict):
                continue
            if entry.get("@ref") in flagged:
                skipped_prediction += 1
                continue
            if trait == "thermophile":
                vote = _temperature_vote(entry)
            else:
                raw = entry.get(key)
                if raw is None or _normalize(raw) in ("", "none"):
                    continue
                vote = _normalize(raw)
            if vote is not None:
                votes.append(vote)
    return votes, skipped_prediction


def _temperature_vote(entry):
    """Vote string for a culture-temp entry, or None when unusable."""
    kind = _normalize(entry.get("type", ""))
    temperature = entry.get("temperature")
    if temperature is None or _normalize(temperature) in ("", "none"):
        return None
    if kind == "growth":
        growth = _normalize(entry.get("growth", ""))
        if growth != "positive":
            return None
        return f"growth-positive:{_normalize(temperature)}"
    if kind in ("optimum", "minimum", "maximum"):
        return f"{kind}:{_normalize(temperature)}"
    return None


def apply_mapping(mapping, record):
    """Apply per-trait value mappings; conflicting votes exclude (None)."""
    labels = {}
    for trait, spec in mapping.get("traits", {}).items():
        votes, _ = extract_votes(record, trait)
        decisions = {spec["values"][vote] for vote in votes if vote in spec["values"]}
        labels[trait] = next(iter(decisions)) if len(decisions) == 1 else None
    return labels


def load_mapping(path):
    from pathlib import Path

    import yaml

    return yaml.safe_load(Path(path).read_text())


def read_genomes_tsv(path):
    import csv

    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))
