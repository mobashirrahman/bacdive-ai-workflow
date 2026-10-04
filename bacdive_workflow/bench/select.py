"""BacDive field calibration and stratified genome selection (Phase 3)."""

import csv

MIN_SUPPORT = 30
MIN_PURITY = 0.95
MIN_STRAINS_FOR_ACCURACY = 500
MIN_REPRODUCIBILITY = 0.95


def calibrate(values, labels):
    """Map each field value to 0/1/exclude from training cross-tabulation.

    ``values`` maps strain -> raw field value; ``labels`` maps strain -> 0/1.
    A value maps to a class only with >=95% purity over >=30 strains.
    """
    by_value = {}
    for strain, value in values.items():
        if strain not in labels or value is None or str(value).strip() == "":
            continue
        by_value.setdefault(str(value).strip(), []).append(labels[strain])
    mapping = {}
    for value, classes in by_value.items():
        if len(classes) < MIN_SUPPORT:
            mapping[value] = "exclude"
            continue
        positive = sum(classes) / len(classes)
        if positive >= MIN_PURITY:
            mapping[value] = 1
        elif positive <= 1 - MIN_PURITY:
            mapping[value] = 0
        else:
            mapping[value] = "exclude"
    return mapping


def mapping_reproduces(mapping, values, labels):
    """Fraction of training labels the mapping reproduces on covered strains."""
    matched = total = 0
    for strain, value in values.items():
        if strain not in labels or value is None or str(value).strip() == "":
            continue
        decision = mapping.get(str(value).strip(), "exclude")
        if decision == "exclude":
            continue
        total += 1
        if decision == labels[strain]:
            matched += 1
    return (matched / total if total else 0.0, total)


def passes_gate(mapping, values, labels):
    """A trait enters accuracy analysis at >=95% reproduction over >=500 strains."""
    fraction, total = mapping_reproduces(mapping, values, labels)
    return total >= MIN_STRAINS_FOR_ACCURACY and fraction >= MIN_REPRODUCIBILITY


def stratified_sample(candidates, n, seed, strata_key="phylum", min_per_stratum=3):
    """Proportional sample with a floor per stratum; one genome per species."""
    import random

    rng = random.Random(seed)
    by_species = {}
    for candidate in candidates:
        by_species.setdefault(candidate["species"], candidate)
    pool = list(by_species.values())
    strata = {}
    for candidate in pool:
        strata.setdefault(candidate[strata_key], []).append(candidate)
    chosen = []
    for members in strata.values():
        take = min(len(members), min_per_stratum)
        chosen.extend(rng.sample(members, take))
    remaining = [c for c in pool if c not in chosen]
    rng.shuffle(remaining)
    quota = max(0, n - len(chosen))
    return chosen + remaining[:quota]


def write_genomes_tsv(path, rows):
    columns = [
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
    ]
    trait_columns = sorted({k for row in rows for k in row if k.startswith("label_")})
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns + trait_columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
