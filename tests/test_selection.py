"""Frozen selection checks: schema shape and unseen-set leakage control."""

import csv
import json
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GENOMES = REPO / "benchmark" / "config" / "genomes.tsv"
LOOKUP = REPO / "benchmark" / "config" / "training_lookup.json"

REQUIRED_COLUMNS = {
    "genome_id",
    "bacdive_id",
    "set",
    "genus_unseen",
    "phylum",
    "family",
    "genus",
    "species",
    "type_strain",
    "assembly_level",
    "n_contigs",
    "genome_bp",
    "harvest_date",
}


def load_rows():
    with open(GENOMES, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_frozen_selection_shape():
    rows = load_rows()
    assert rows, "genomes.tsv has no genomes"
    assert REQUIRED_COLUMNS <= set(rows[0])
    label_columns = [c for c in rows[0] if c.startswith("label_")]
    assert label_columns, "no label_<trait> columns"
    assert {r["set"] for r in rows} <= {"seen", "unseen"}
    assert any(r["set"] == "seen" for r in rows)
    assert any(r["set"] == "unseen" for r in rows)
    assert len({r["genome_id"] for r in rows}) == len(rows)


def test_unseen_set_has_no_training_leakage():
    """No unseen row shares a BacDive ID, accession, or species with training."""
    lookup = json.loads(LOOKUP.read_text())
    train_ids = set(lookup["strain_ids"])
    train_accessions = set(lookup["accessions"])
    train_species = set(lookup["species"])

    def forms(accession):
        accession = accession.strip()
        variants = {accession}
        if accession.startswith(("GCA_", "GCF_")) and "." in accession:
            variants.add(accession.split(".")[0])
        return variants

    for row in load_rows():
        if row["set"] != "unseen":
            continue
        assert int(row["bacdive_id"]) not in train_ids, row
        assert not (forms(row["genome_id"]) & train_accessions), row
        assert row["species"] not in train_species, row
