"""Shared identity, serialization, and input validation helpers."""

import csv
import hashlib
import json
import math
import os
import re
import tempfile
from pathlib import Path


class WorkflowError(ValueError):
    """An actionable error in workflow inputs or model outputs."""


TRAITS = {
    "acidophile": "Acidophilic",
    "gram-positive": "Gram-positive",
    "spore-forming": "Spore-forming",
    "aerobic": "Aerobic",
    "anaerobic": "Anaerobic",
    "thermophile": "Thermophilic",
    "psychrophile": "Psychrophilic",
    "motile2+": "Flagellated motility",
}
DEFAULT_EVALUE = 1e-20
MISSING = {"", "na", "n/a", "nan", "none", "null", "-"}
ALIASES = {
    "gram stain": "gram_stain",
    "oxygen tolerance": "oxygen_tolerance",
    "sgbs_Completeness": "completeness",
    "sgbs_Contamination": "contamination",
}


def normalize_id(value):
    """Remove known data suffixes, preserving accession/version dots."""
    name = Path(str(value)).name
    suffixes = (".gz", ".json", ".tsv", ".fasta", ".faa", ".fna", ".fa")
    while True:
        suffix = next((s for s in suffixes if name.lower().endswith(s)), None)
        if suffix is None:
            return name
        name = name[: -len(suffix)]


def validate_sample_id(value):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", value):
        raise WorkflowError(f"Invalid sample_id {value!r}; use letters, digits, dots, '_' or '-'.")
    return value


def sha256_file(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_text(path, content):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent, prefix=f".{path.name}.", delete=False
        ) as handle:
            temporary = Path(handle.name)
            handle.write(content)
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def write_json(path, payload):
    atomic_text(path, json.dumps(payload, indent=2, allow_nan=False) + "\n")


def clean_value(value):
    value = str(value or "").strip()
    return "" if value.lower() in MISSING else value


def read_table(path):
    path = Path(path)
    try:
        with path.open(encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t" if path.suffix == ".tsv" else ",")
            if not reader.fieldnames or len(set(reader.fieldnames)) != len(reader.fieldnames):
                raise WorkflowError(f"{path}: missing or duplicate column names.")
            rows = []
            for line, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise WorkflowError(f"{path}:{line}: row does not match the table header.")
                normalized = {}
                for key, value in row.items():
                    canonical = ALIASES.get(key.strip(), key.strip())
                    if canonical in normalized:
                        raise WorkflowError(f"{path}: duplicate aliases for {canonical}.")
                    normalized[canonical] = clean_value(value)
                rows.append(normalized)
    except OSError as error:
        raise WorkflowError(f"Cannot read table {path}: {error}") from error
    if not rows:
        raise WorkflowError(f"{path}: table has no samples.")
    return rows


def index_metadata(path):
    indexed = {}
    for row in read_table(path):
        sample_id = row.get("sample_id") or normalize_id(row.get("genome_address", ""))
        validate_sample_id(sample_id)
        if sample_id in indexed:
            raise WorkflowError(f"{path}: duplicate sample_id {sample_id!r}.")
        for field in ("completeness", "contamination"):
            if row.get(field):
                try:
                    value = float(row[field])
                except ValueError as error:
                    raise WorkflowError(f"{sample_id}: invalid {field}.") from error
                if (
                    not math.isfinite(value)
                    or value < 0
                    or (field == "completeness" and value > 100)
                ):
                    raise WorkflowError(f"{sample_id}: invalid {field} {value}.")
        indexed[sample_id] = row
    return indexed


def load_samples(path):
    indexed = index_metadata(path)
    for sample_id, row in indexed.items():
        fields = [
            field for field in ("genome_path", "protein_path", "annotation_path") if row.get(field)
        ]
        if len(fields) != 1:
            raise WorkflowError(
                f"{sample_id}: supply exactly one genome_path, protein_path or annotation_path."
            )
        field = fields[0]
        candidate = Path(row[field])
        if not candidate.is_absolute():
            candidate = Path(path).resolve().parent / candidate
        candidate = candidate.resolve()
        if not candidate.is_file() or candidate.stat().st_size == 0:
            raise WorkflowError(f"{sample_id}: input file is missing or empty: {candidate}.")
        row[field] = str(candidate)
        row["input_type"] = {
            "genome_path": "genome",
            "protein_path": "protein",
            "annotation_path": "annotation",
        }[field]
    return indexed
