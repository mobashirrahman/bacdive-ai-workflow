"""Unit tests for bacdive_workflow.bench (offline, seeded, deterministic)."""

import numpy as np
import pytest

from bacdive_workflow.bench import derive_seed
from bacdive_workflow.bench.drift import drift_row, jaccard
from bacdive_workflow.bench.fragments import cut_fragments, draw_seed, keep_genes
from bacdive_workflow.bench.metrics import accuracy_rows, stratified_bootstrap_ci, summarize
from bacdive_workflow.bench.select import calibrate, passes_gate, stratified_sample
from bacdive_workflow.bench.selection import (
    assess_assembly,
    normalize_accession,
    parse_record,
    proportional_sample,
    repair_rare_traits,
)


def test_derive_seed_is_stable_and_sensitive():
    assert derive_seed(1, "a", 100, 0, 0) == derive_seed(1, "a", 100, 0, 0)
    assert derive_seed(1, "a", 100, 0, 0) != derive_seed(1, "a", 100, 0, 1)
    assert derive_seed(1, "a", 100, 0, 0) != derive_seed(2, "a", 100, 0, 0)
    assert 0 <= draw_seed(1, "a", 100, 0, 0) < 2**32


def test_summarize_leaves_metrics_empty_for_small_classes():
    row = summarize([1] * 5 + [0] * 20, [1] * 5 + [0] * 20, [0.9] * 25)
    assert row["n"] == 25
    assert row["balanced_accuracy"] is None


def test_summarize_balanced_accuracy_and_mcc():
    labels = [1] * 10 + [0] * 10
    predictions = [1] * 10 + [0] * 10
    scores = [0.9] * 10 + [0.1] * 10
    row = summarize(labels, predictions, scores)
    assert row["balanced_accuracy"] == 1.0
    assert row["mcc"] == 1.0
    assert row["auroc"] == 1.0
    low, high = stratified_bootstrap_ci(labels, predictions, scores, 50, seed=7)
    assert low == high == 1.0


def test_accuracy_rows_groups_and_phylum_floor():
    records = []
    for i in range(20):
        records.append(
            {
                "genome": f"s{i}",
                "set": "seen",
                "genus_unseen": False,
                "phylum": "Firmicutes",
                "labels": {"gram-positive": i % 2},
                "predictions": {"gram-positive": i % 2},
                "scores": {"gram-positive": 0.9 if i % 2 else 0.1},
            }
        )
    for i in range(5):
        records.append(
            {
                "genome": f"u{i}",
                "set": "unseen",
                "genus_unseen": True,
                "phylum": "Tiny",
                "labels": {"gram-positive": 1},
                "predictions": {"gram-positive": 1},
                "scores": {"gram-positive": 0.8},
            }
        )
    groups = {
        "seen": lambda r: r["set"] == "seen",
        "unseen": lambda r: r["set"] == "unseen",
        "genus_unseen": lambda r: r["set"] == "unseen" and r["genus_unseen"],
    }
    rows = accuracy_rows(records, ["gram-positive"], groups, bootstrap=50, seed=1)
    by_group = {(r["group"], r["phylum"]): r for r in rows}
    assert by_group[("seen", "")]["balanced_accuracy"] == 1.0
    # "Tiny" phylum has n=5 < 15, so no per-phylum row.
    assert ("phylum", "Tiny") not in by_group
    assert ("phylum", "Firmicutes") in by_group


def test_fragments_cover_contig_and_keep_interior_genes():
    rng = np.random.default_rng(42)
    fragments = cut_fragments(60000, 20000, 1.0, 1000, rng)
    assert fragments[0][0] == 0
    assert fragments[-1][1] == 60000
    genes = [
        {"protein_id": "g1", "contig": "c", "start": 100, "end": 500},
        {
            "protein_id": "g2",
            "contig": "c",
            "start": fragments[0][1] - 10,
            "end": fragments[0][1] + 10,
        },
    ]
    kept = keep_genes(genes, [("c", s, e) for s, e in fragments[:1]])
    assert [g["protein_id"] for g in kept] == ["g1"]


def test_jaccard_and_drift_row():
    assert jaccard({"a"}, {"a"}) == 1.0
    assert jaccard(set(), set()) is None
    row = drift_row("g", "5.74", {"a", "b"}, {"b", "c"}, 1, "exact")
    assert row["jaccard"] == pytest.approx(1 / 3)
    assert row["n_gained"] == row["n_lost"] == 1


def test_calibration_gate():
    values = {f"s{i}": ("positive" if i < 40 else "negative") for i in range(600)}
    labels = {f"s{i}": (1 if i < 40 else 0) for i in range(600)}
    mapping = calibrate(values, labels)
    assert mapping == {"positive": 1, "negative": 0}
    assert passes_gate(mapping, values, labels)
    noisy = dict(labels)
    for i in range(0, 600, 2):
        noisy[f"s{i}"] = 1 - noisy[f"s{i}"]
    assert not passes_gate(calibrate(values, noisy), values, noisy)


def test_stratified_sample_floor_and_species_dedup():
    candidates = [
        {"species": f"sp{i % 4}", "phylum": "A" if i % 2 else "B", "genome_id": f"g{i}"}
        for i in range(40)
    ]
    chosen = stratified_sample(candidates, 6, seed=3)
    assert len({c["species"] for c in chosen}) == len(chosen)
    phyla = {c["phylum"] for c in chosen}
    assert phyla == {"A", "B"}


def test_parse_record_handles_dict_or_list_genomes():
    base = {
        "Name and taxonomic classification": {
            "phylum": "Firmicutes",
            "family": "F",
            "genus": "G",
            "species": "G s",
            "type strain": "yes",
        },
        "Sequence information": {
            "Genome sequences": {"INSDC accession": "GCA_1", "assembly level": "Scaffold"},
        },
    }
    parsed = parse_record(base)
    assert parsed["species"] == "G s"
    assert parsed["type_strain"] is True
    assert parsed["accessions"][0]["insdc"] == "GCA_1"
    assert normalize_accession("GCA_000429225.1") == {"GCA_000429225.1", "GCA_000429225"}
    assert normalize_accession("1120941.3") == {"1120941.3"}


def test_assess_assembly_prefers_complete_and_filters_flags():
    good = {
        "accession": "GCA_1.1",
        "assembly_info": {"assembly_level": "Complete Genome", "assembly_status": "current"},
        "assembly_stats": {"number_of_contigs": "1", "total_sequence_length": "1000"},
    }
    accepted, reason = assess_assembly(good)
    assert accepted["accession"] == "GCA_1.1"
    assert reason == "preferred"
    relaxed = {
        "accession": "GCA_2.1",
        "assembly_info": {"assembly_level": "Scaffold", "assembly_status": "current"},
        "assembly_stats": {"number_of_contigs": "50", "total_sequence_length": "1000"},
    }
    assert assess_assembly(relaxed)[1] == "relaxed-50-contigs"
    bad = {
        "accession": "GCA_3.1",
        "assembly_info": {
            "assembly_level": "Scaffold",
            "assembly_status": "current",
            "anomalous_list": ["chimeric"],
        },
        "assembly_stats": {"number_of_contigs": "10", "total_sequence_length": "1000"},
    }
    assert assess_assembly(bad) == (None, "anomalous")


def test_proportional_sample_floor_and_repair():
    candidates = [
        {
            "species": f"A{i}",
            "phylum": "A",
            "genome_id": f"a{i}",
            "bacdive_id": i,
            "labels": {"t": 1 if i < 2 else 0},
        }
        for i in range(10)
    ] + [
        {
            "species": f"B{i}",
            "phylum": "B",
            "genome_id": f"b{i}",
            "bacdive_id": 100 + i,
            "labels": {"t": 1 if i < 25 else 0},
        }
        for i in range(50)
    ]
    chosen = proportional_sample(candidates, 20, seed=5)
    assert len(chosen) == 20
    assert len({c["species"] for c in chosen}) == 20
    phyla = [c["phylum"] for c in chosen]
    assert phyla.count("A") >= 3 and phyla.count("B") >= 3
    assert phyla.count("B") > phyla.count("A")  # proportional to pool sizes
    repaired, n_relaxed = repair_rare_traits(chosen[:6], candidates, ["t"], seed=5)
    positives = sum(1 for c in repaired if c["labels"]["t"] == 1)
    assert positives >= 20
    assert n_relaxed == 0
