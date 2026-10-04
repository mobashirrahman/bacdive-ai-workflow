"""Genome selection for the benchmark: seen and unseen sets (Phase 3, run once).

Seen strains come from the published training labels; unseen strains are
leakage-controlled at strain, genome, and species level with reference labels
from the calibrated BacDive mapping. Output is frozen and committed.

Command:
    python -m bacdive_workflow.bench.select
        --training-labels <csv> --training-features <csv>
        --cache-dir <bacdive_cache> --mapping <label_mapping.yaml>
        --seed 20261004 --n-seen 150 --n-unseen 150
        --genomes-out benchmark/config/genomes.tsv
        --lookup-out benchmark/config/training_lookup.json
        --summary-out docs/benchmark/selection_summary.csv
"""

import argparse
import csv
import json
import random
import subprocess
from pathlib import Path

from bacdive_workflow.bench.labels import apply_mapping
from bacdive_workflow.bench.select import write_genomes_tsv
from bacdive_workflow.common import WorkflowError

RARE_TRAITS = ["thermophile", "spore-forming", "motile2+"]
MIN_CLASS_TARGET = 20


def parse_record(record):
    """Taxonomy and genome accessions from one BacDive record."""
    taxonomy = record.get("Name and taxonomic classification", {})
    if not isinstance(taxonomy, dict):
        return None
    genomes = record.get("Sequence information", {}).get("Genome sequences", [])
    if isinstance(genomes, dict):
        genomes = [genomes]
    accessions = []
    for entry in genomes or []:
        if not isinstance(entry, dict):
            continue
        accessions.append(
            {
                "insdc": (entry.get("INSDC accession") or "").strip(),
                "bvbrc": (entry.get("BV-BRC accession") or "").strip(),
                "level": (entry.get("assembly level") or "").strip(),
            }
        )
    type_strain = str(taxonomy.get("type strain", "")).strip().lower() == "yes"
    return {
        "phylum": (taxonomy.get("phylum") or "").strip(),
        "family": (taxonomy.get("family") or "").strip(),
        "genus": (taxonomy.get("genus") or "").strip(),
        "species": (taxonomy.get("species") or "").strip(),
        "type_strain": type_strain,
        "accessions": accessions,
    }


def normalize_accession(accession):
    """Versioned, base, and raw forms for leakage comparison."""
    accession = accession.strip()
    forms = {accession}
    if accession.startswith(("GCA_", "GCF_")) and "." in accession:
        forms.add(accession.split(".")[0])
    return forms


def training_lookup(labels_path, features_path):
    """Strain IDs, accessions, species, genera, and labels from training data."""
    labels: dict = {}
    strain_accessions: dict = {}
    with open(labels_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            strain = int(row["ID_strains"])
            labels.setdefault(strain, {})[row["trait"]] = int(row["ability"])
            strain_accessions.setdefault(strain, set()).add(row["sequence_acc"].strip())
    with open(features_path, newline="", encoding="utf-8") as handle:
        handle.readline()  # mislabelled header; strain ID is the first column
        strains_with_features = {line.split(",", 1)[0] for line in handle if line.strip()}
    strains_with_features = {int(s) for s in strains_with_features if s.strip().isdigit()}
    accessions: set = set()
    for accs in strain_accessions.values():
        for accession in accs:
            accessions |= normalize_accession(accession)
    return {
        "strain_ids": set(labels),
        "accessions": accessions,
        "labels": labels,
        "strain_accessions": strain_accessions,
        "strains_with_features": strains_with_features,
    }


def load_records(cache_dir, subdir):
    records = {}
    for batch in sorted(Path(cache_dir, subdir).glob("batch-*.json")):
        records.update(json.loads(batch.read_text()).get("results", {}))
    return records


def _cache_path(cache_dir, accession):
    return Path(cache_dir) / "datasets" / f"{accession.replace('.', '_')}.json"


def batch_summarize(queries, cache_dir, chunk=100):
    """Resolve many INSDC accessions with batched datasets calls (cached each)."""
    resolved = {}
    todo = []
    for query in dict.fromkeys(queries):
        cached = _cache_path(cache_dir, query)
        if cached.is_file():
            resolved[query] = json.loads(cached.read_text())
        else:
            todo.append(query)
    for start in range(0, len(todo), chunk):
        group = todo[start : start + chunk]
        print(f"datasets summary {start + 1}-{start + len(group)} of {len(todo)}", flush=True)
        try:
            completed = subprocess.run(
                ["datasets", "summary", "genome", "accession", *group],
                check=True,
                capture_output=True,
                text=True,
                timeout=600,
            )
        except (
            subprocess.CalledProcessError,
            FileNotFoundError,
            subprocess.TimeoutExpired,
        ) as error:
            raise WorkflowError(
                f"datasets summary failed for batch starting {group[0]}: {error}"
            ) from error
        payload = json.loads(completed.stdout)
        by_base = {}
        for report in payload.get("reports", []):
            by_base.setdefault(report.get("accession", "").split(".")[0], report)
            current = report.get("current_accession", "")
            if current:
                by_base.setdefault(current.split(".")[0], report)
        for query in group:
            match = by_base.get(query.split(".")[0])
            single = {"reports": [match] if match else [], "total_count": 1 if match else 0}
            cached = _cache_path(cache_dir, query)
            cached.parent.mkdir(parents=True, exist_ok=True)
            cached.write_text(json.dumps(single, indent=2) + "\n")
            resolved[query] = single
    return resolved


def assess_assembly(report):
    """Accept/reject one datasets assembly report with reasons."""
    info = report.get("assembly_info", {})
    stats = report.get("assembly_stats", {})
    status = str(info.get("assembly_status", "")).strip().lower()
    if status and status != "current":
        return None, f"status-{status}"
    if info.get("anomalous_list") or info.get("anomalousList"):
        return None, "anomalous"
    category = str(info.get("refseq_category", "")).strip().lower()
    organism = str(info.get("organism_name", "")).lower()
    if category == "metagenome" or "metagenome" in organism:
        return None, "metagenome"
    level = str(info.get("assembly_level", "")).strip()

    def _as_int(value):
        try:
            return int(value)
        except (TypeError, ValueError):
            return None

    contigs = _as_int(stats.get("number_of_contigs"))
    total = _as_int(stats.get("total_sequence_length"))
    accession = report.get("accession", "")
    accepted = {"accession": accession, "level": level, "contigs": contigs, "bp": total}
    if level in ("Complete Genome", "Chromosome"):
        return accepted, "preferred"
    if contigs is not None and contigs <= 50:
        return accepted, "relaxed-50-contigs"
    return None, f"level-{level}-contigs-{contigs}"


def proportional_sample(candidates, n, seed, strata_key="phylum", min_per_stratum=3):
    """Proportional allocation with a floor per stratum; one genome per species."""
    rng = random.Random(seed)
    by_species = {}
    for candidate in candidates:
        by_species.setdefault(candidate["species"], candidate)
    pool = list(by_species.values())
    strata = {}
    for candidate in pool:
        strata.setdefault(candidate[strata_key] or "unknown", []).append(candidate)
    quotas = {}
    remaining = n
    for name in sorted(strata):
        quotas[name] = min(len(strata[name]), min_per_stratum)
        remaining -= quotas[name]
    if remaining > 0:
        sizes = {name: len(strata[name]) - quotas[name] for name in strata}
        total = sum(sizes.values())
        if total > 0:
            for name in sorted(strata):
                add = min(sizes[name], round(remaining * sizes[name] / total))
                quotas[name] += add
    chosen = []
    for name in sorted(strata):
        members = strata[name][:]
        rng.shuffle(members)
        chosen.extend(members[: quotas[name]])
    if len(chosen) < n:
        rest = [c for c in pool if c not in chosen]
        rng.shuffle(rest)
        chosen.extend(rest[: n - len(chosen)])
    return chosen


def repair_rare_traits(chosen, preferred_pool, traits, seed, relaxed_pool=()):
    """Add pool strains so each trait reaches 20 positives and 20 negatives.

    Preferred genomes first; only remaining deficits are topped up from the
    relaxed (<=50 contigs) pool. Returns (chosen, n_relaxed_added).
    """
    rng = random.Random(seed)
    chosen = list(chosen)
    have = {c["bacdive_id"] for c in chosen}
    rest_preferred = [c for c in preferred_pool if c["bacdive_id"] not in have]
    rest_relaxed = [c for c in relaxed_pool if c["bacdive_id"] not in have]
    rng.shuffle(rest_preferred)
    rng.shuffle(rest_relaxed)
    n_relaxed = 0
    for trait in traits:
        for wanted in (1, 0):
            while sum(1 for c in chosen if c["labels"].get(trait) == wanted) < MIN_CLASS_TARGET:
                donor = next((c for c in rest_preferred if c["labels"].get(trait) == wanted), None)
                source = rest_preferred
                if donor is None:
                    donor = next(
                        (c for c in rest_relaxed if c["labels"].get(trait) == wanted), None
                    )
                    source = rest_relaxed
                if donor is None:
                    break
                source.remove(donor)
                chosen.append(donor)
                if source is rest_relaxed:
                    n_relaxed += 1
    return chosen, n_relaxed


def main(argv=None):
    parser = argparse.ArgumentParser(description="Select and freeze benchmark genomes.")
    parser.add_argument("--training-labels", type=Path, required=True)
    parser.add_argument("--training-features", type=Path, required=True)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--n-seen", type=int, default=150)
    parser.add_argument("--n-unseen", type=int, default=150)
    parser.add_argument("--genomes-out", type=Path, required=True)
    parser.add_argument("--lookup-out", type=Path, required=True)
    parser.add_argument("--summary-out", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        import yaml

        mapping = yaml.safe_load(args.mapping.read_text())
        eligible = list(mapping.get("accuracy_traits", []))
        print(f"eligible traits: {eligible}", flush=True)
        lookup = training_lookup(args.training_labels, args.training_features)
        print(f"training strains: {len(lookup['strain_ids'])}", flush=True)
        train_species = set()
        train_genera = set()
        records = load_records(args.cache_dir, "train")
        records.update(load_records(args.cache_dir, "candidates"))
        for strain in lookup["strain_ids"]:
            record = records.get(str(strain))
            if record is None:
                continue
            parsed = parse_record(record)
            if parsed is None:
                continue
            if parsed["species"]:
                train_species.add(parsed["species"])
            if parsed["genus"]:
                train_genera.add(parsed["genus"])
        stats = {"excluded": {}}

        def note(reason):
            stats["excluded"][reason] = stats["excluded"].get(reason, 0) + 1

        # Seen pool: training strains with GCA/GCF accessions and feature rows.
        seen_pool = []
        for strain in sorted(lookup["strain_ids"]):
            accs = [
                a
                for a in lookup["strain_accessions"].get(strain, set())
                if a.startswith(("GCA_", "GCF_"))
            ]
            if not accs:
                note("seen-no-insdc-accession")
                continue
            if strain not in lookup["strains_with_features"]:
                note("seen-no-feature-rows")
                continue
            record = records.get(str(strain))
            if record is None:
                note("seen-not-harvested")
                continue
            parsed = parse_record(record)
            if parsed is None or not parsed["species"]:
                note("seen-no-taxonomy")
                continue
            seen_pool.append(
                {
                    "bacdive_id": strain,
                    "set": "seen",
                    "genus_unseen": False,
                    "accession_query": sorted(accs)[0],
                    "labels": {t: lookup["labels"][strain].get(t) for t in eligible},
                    **{
                        k: parsed[k]
                        for k in ("phylum", "family", "genus", "species", "type_strain")
                    },
                }
            )
        # Unseen pool: leakage-controlled at strain, genome, and species level.
        unseen_pool = []
        for bacdive_id, record in records.items():
            strain = int(bacdive_id)
            if strain in lookup["strain_ids"]:
                continue
            parsed = parse_record(record)
            if parsed is None or not parsed["species"]:
                note("unseen-no-taxonomy")
                continue
            insdc = [a["insdc"] for a in parsed["accessions"] if a["insdc"]]
            if not insdc:
                note("unseen-no-insdc-accession")
                continue
            if any(form in lookup["accessions"] for a in insdc for form in normalize_accession(a)):
                note("unseen-accession-in-training")
                continue
            bvbrc = [a["bvbrc"] for a in parsed["accessions"] if a["bvbrc"]]
            if any(b in lookup["accessions"] for b in bvbrc):
                note("unseen-accession-in-training")
                continue
            if parsed["species"] in train_species:
                note("unseen-species-in-training")
                continue
            labels = apply_mapping(mapping, record)
            n_labels = sum(1 for t in eligible if labels.get(t) is not None)
            if n_labels < 2:
                note("unseen-fewer-than-2-labels")
                continue
            unseen_pool.append(
                {
                    "bacdive_id": strain,
                    "set": "unseen",
                    "genus_unseen": parsed["genus"] not in train_genera,
                    "accession_query": sorted(insdc)[0],
                    "labels": {t: labels.get(t) for t in eligible},
                    **{
                        k: parsed[k]
                        for k in ("phylum", "family", "genus", "species", "type_strain")
                    },
                }
            )
        stats["seen_pool"] = len(seen_pool)
        stats["unseen_pool"] = len(unseen_pool)
        print(f"pools: seen {len(seen_pool)}, unseen {len(unseen_pool)}", flush=True)
        # Genome filter via NCBI datasets (batched, cached per accession).
        harvest_date = json.loads((args.cache_dir / "harvest.json").read_text())["harvest_date"]
        queries = [c["accession_query"] for c in seen_pool + unseen_pool]
        print(f"resolving {len(queries)} accessions", flush=True)
        resolved = batch_summarize(queries, args.cache_dir)
        for pool in (seen_pool, unseen_pool):
            for candidate in pool:
                payload = resolved.get(candidate["accession_query"], {})
                reports = payload.get("reports", [])
                if not reports:
                    candidate["genome_exclusion"] = "datasets-no-report"
                    continue
                accepted, reason = assess_assembly(reports[0])
                if accepted is None:
                    candidate["genome_exclusion"] = reason
                else:
                    candidate["assembly"] = accepted
                    candidate["assembly"]["relaxation"] = reason
        for pool in (seen_pool, unseen_pool):
            kept = [c for c in pool if "assembly" in c]
            for c in pool:
                if "assembly" not in c:
                    note(f"genome-{c.get('genome_exclusion', 'unknown')}")
            pool[:] = kept
        stats["seen_after_genome_filter"] = len(seen_pool)
        stats["unseen_after_genome_filter"] = len(unseen_pool)
        # Prefer Complete/Chromosome; relax to <=50 contigs only when needed.
        seen_pref = [c for c in seen_pool if c["assembly"]["relaxation"] == "preferred"]
        unseen_pref = [c for c in unseen_pool if c["assembly"]["relaxation"] == "preferred"]
        stats["seen_preferred"] = len(seen_pref)
        stats["unseen_preferred"] = len(unseen_pref)
        stats["pool_rare_counts"] = {
            name: {
                t: {
                    "positive": sum(1 for c in pool if c["labels"].get(t) == 1),
                    "negative": sum(1 for c in pool if c["labels"].get(t) == 0),
                }
                for t in eligible
            }
            for name, pool in (("seen", seen_pool), ("unseen", unseen_pool))
        }
        # Relaxed (<=50 contigs) unseen candidates, kept only for rare-trait top-ups.
        # Built before the pools are narrowed to preferred assemblies.
        unseen_relaxed = [c for c in unseen_pool if c not in unseen_pref]
        stats["unseen_relaxed_available"] = len(unseen_relaxed)
        if len(seen_pref) >= args.n_seen:
            seen_pool, stats["seen_relaxation"] = seen_pref, False
        else:
            stats["seen_relaxation"] = True
        if len(unseen_pref) >= args.n_unseen:
            unseen_pool, stats["unseen_relaxation"] = unseen_pref, False
        else:
            stats["unseen_relaxation"] = True
        print(
            f"relaxation: seen {stats['seen_relaxation']}, unseen {stats['unseen_relaxation']}",
            flush=True,
        )
        if len(seen_pool) < args.n_seen or len(unseen_pool) < args.n_unseen:
            raise WorkflowError(
                f"Pool too small: seen {len(seen_pool)}, unseen {len(unseen_pool)}."
            )
        chosen_seen = proportional_sample(seen_pool, args.n_seen, args.seed)
        chosen_unseen = proportional_sample(unseen_pool, args.n_unseen, args.seed + 1)
        chosen_unseen, n_relaxed_topup = repair_rare_traits(
            chosen_unseen, unseen_pool, eligible, args.seed + 2, unseen_relaxed
        )
        stats["unseen_relaxed_topup"] = n_relaxed_topup
        stats["relaxed_topup_genomes"] = [
            c["assembly"]["accession"] for c in chosen_unseen if c in unseen_relaxed
        ]
        rows = []
        for candidate in chosen_seen + chosen_unseen:
            row = {
                "genome_id": candidate["assembly"]["accession"],
                "bacdive_id": candidate["bacdive_id"],
                "set": candidate["set"],
                "genus_unseen": candidate["genus_unseen"],
                "phylum": candidate["phylum"],
                "family": candidate["family"],
                "genus": candidate["genus"],
                "species": candidate["species"],
                "type_strain": candidate["type_strain"],
                "assembly_level": candidate["assembly"]["level"],
                "n_contigs": candidate["assembly"]["contigs"],
                "genome_bp": candidate["assembly"]["bp"],
                "harvest_date": harvest_date,
            }
            for trait in eligible:
                value = candidate["labels"].get(trait)
                row[f"label_{trait}"] = "" if value is None else value
            rows.append(row)
        args.genomes_out.parent.mkdir(parents=True, exist_ok=True)
        write_genomes_tsv(args.genomes_out, rows)
        args.lookup_out.parent.mkdir(parents=True, exist_ok=True)
        args.lookup_out.write_text(
            json.dumps(
                {
                    "harvest_date": harvest_date,
                    "strain_ids": sorted(lookup["strain_ids"]),
                    "accessions": sorted(lookup["accessions"]),
                    "species": sorted(train_species),
                    "genera": sorted(train_genera),
                },
                indent=2,
            )
            + "\n"
        )
        summary = []
        for candidate in chosen_seen + chosen_unseen:
            entry = {
                "set": candidate["set"],
                "phylum": candidate["phylum"] or "unknown",
            }
            for trait in eligible:
                value = candidate["labels"].get(trait)
                entry[trait] = (
                    "positive" if value == 1 else ("negative" if value == 0 else "unlabeled")
                )
            summary.append(entry)
        from collections import Counter

        counts = Counter((s["set"], s["phylum"], t, s[t]) for s in summary for t in eligible)
        args.summary_out.parent.mkdir(parents=True, exist_ok=True)
        with open(args.summary_out, "w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["set", "phylum", "trait", "class", "count"])
            for key in sorted(counts):
                writer.writerow([*key, counts[key]])
        stats["chosen_seen"] = len(chosen_seen)
        stats["chosen_unseen"] = len(chosen_unseen)
        stats["rare_trait_counts"] = {
            t: {
                "positive": sum(1 for c in chosen_unseen if c["labels"].get(t) == 1),
                "negative": sum(1 for c in chosen_unseen if c["labels"].get(t) == 0),
            }
            for t in eligible
        }
        stats["seed"] = args.seed
        stats["n_seen_target"] = args.n_seen
        stats["n_unseen_target"] = args.n_unseen
        stats["eligible_traits"] = eligible
        stats_path = args.summary_out.with_name("selection_stats.json")
        stats_path.write_text(json.dumps(stats, indent=2) + "\n")
        print(json.dumps(stats, indent=2))
        return 0
    except (WorkflowError, OSError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
