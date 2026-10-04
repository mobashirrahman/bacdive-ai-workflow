"""Shared helpers: repo paths, genome table, selection subsets."""

import sys
from pathlib import Path

import pandas as pd
from snakemake.utils import validate

BENCH = Path(workflow.basedir)
REPO = BENCH.parent.resolve()
sys.path.insert(0, str(REPO))

from bacdive_workflow.common import TRAITS as TRAITS8

CONFIG_SCHEMA = str(BENCH / "schemas" / "config.schema.yaml")
validate(config, CONFIG_SCHEMA)

SEED = config["seed"]
PRIMARY_IPS = config["interproscan"]["primary"]
PRIMARY_GENES = "prodigal_single"
EVALUE = config["evalue"]
BOOTSTRAP = config["bootstrap"]
TEST_MODE = config.get("test_mode", False)
PUBLISH_DOCS = config.get("publish_docs", True)

SUBDIR = ".test" if TEST_MODE else "config"
GENOMES_FILE = config.get("genomes_file", str(BENCH / SUBDIR / "genomes.tsv"))
MAPPING_FILE = str(BENCH / SUBDIR / "label_mapping.yaml")
FIXTURE_SHORTCUT = str(BENCH / ".test" / "shortcut")
TEST_LABELS = str(BENCH / ".test" / "training_labels.csv")
TEST_FEATURES = str(BENCH / ".test" / "training_features.csv")

BOOL_COLUMNS = {"genus_unseen", "type_strain"}


def _as_bool(value):
    text = str(value).strip().lower()
    if text in ("true", "1", "yes"):
        return True
    if text in ("false", "0", "no"):
        return False
    raise ValueError(f"Not a boolean: {value!r}")


def load_genomes():
    table = pd.read_csv(GENOMES_FILE, sep="\t", dtype=str).fillna("")
    for column in BOOL_COLUMNS:
        if column in table.columns:
            table[column] = table[column].map(_as_bool)
    table["bacdive_id"] = table["bacdive_id"].astype(int)
    table["n_contigs"] = table["n_contigs"].astype(int)
    table["genome_bp"] = table["genome_bp"].astype(int)
    validate(table, str(BENCH / "schemas" / "genomes.schema.yaml"))
    return table


GENOMES = load_genomes()
GENOME_IDS = sorted(GENOMES["genome_id"].tolist())
LABEL_COLUMNS = [c for c in GENOMES.columns if c.startswith("label_")]
TRAITS = [c[len("label_") :] for c in LABEL_COLUMNS]

SEEN_IDS = sorted(GENOMES.loc[GENOMES["set"] == "seen", "genome_id"].tolist())
UNSEEN_IDS = sorted(GENOMES.loc[GENOMES["set"] == "unseen", "genome_id"].tolist())


def robustness_genomes():
    import random

    rng = random.Random(SEED)
    half = max(1, config["robustness_subset"] // 2)
    chosen = rng.sample(SEEN_IDS, min(half, len(SEEN_IDS)))
    chosen += rng.sample(
        UNSEEN_IDS, min(config["robustness_subset"] - len(chosen), len(UNSEEN_IDS))
    )
    return sorted(chosen)


def shortcut_genomes():
    import random

    rng = random.Random(SEED + 1)
    pool = GENOME_IDS
    return sorted(rng.sample(pool, min(config["validate_shortcut_subset"], len(pool))))


ROBUSTNESS_IDS = robustness_genomes()
SHORTCUT_IDS = shortcut_genomes()

RESULT_TABLES = [
    "results/tables/predictions_full.tsv",
    "results/tables/accuracy.tsv",
    "results/tables/disagreements.tsv",
    "results/tables/drift.tsv",
    "results/tables/drift_pfams.tsv",
    "results/tables/degrade_summary.tsv",
    "results/tables/degrade_thresholds.tsv",
    "results/tables/shortcut_validation.tsv",
    "results/tables/robustness.tsv",
]

RESULT_FIGURES = [
    "results/figures/accuracy.png",
    "results/figures/accuracy.svg",
    "results/figures/flip_rate.png",
    "results/figures/flip_rate.svg",
    "results/figures/drift.png",
    "results/figures/drift.svg",
    "results/figures/robustness.png",
    "results/figures/robustness.svg",
]

DOCS_DIR = config.get("docs_dir", str(REPO / "docs/benchmark"))
DOC_TABLES = [f"{DOCS_DIR}/{Path(t).name}" for t in RESULT_TABLES]
DOC_FIGURES = [f"{DOCS_DIR}/{Path(f).name}" for f in RESULT_FIGURES]


def model_archive_path():
    name = config.get("model_archive", "Archiv.zip")
    candidate = Path(name)
    if candidate.is_file():
        return str(candidate)
    return str(REPO / name)


def pgap_inputs(wildcards):
    """InterProScan PGAP outputs, only for genomes whose download has proteins."""
    available = pd.read_csv(checkpoints.survey_pgap.get().output[0], sep="\t")
    genomes = [
        row["genome"]
        for _, row in available.iterrows()
        if str(row["pgap_available"]).strip().lower() == "true"
        and row["genome"] in ROBUSTNESS_IDS
    ]
    return expand(
        f"results/interpro/{PRIMARY_IPS}/pgap/{{genome}}.tsv",
        genome=genomes,
    )
