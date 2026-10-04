"""Analysis B: annotation drift against the published training features (seen set)."""


rule published_sets:
    output:
        sets=expand("results/published/{genome}.pfams", genome=SEEN_IDS),
        meta="results/published/accession_match.tsv",
    log:
        "logs/published_sets.log",
    conda:
        "../envs/stats.yaml"
    params:
        repo=str(REPO),
        labels=(
            TEST_LABELS
            if TEST_MODE
            else str(
                Path(config["data_root"])
                / "training_data"
                / "training_data_labels.csv"
            )
        ),
        features=(
            TEST_FEATURES
            if TEST_MODE
            else str(
                Path(config["data_root"])
                / "training_data"
                / "training_data_features.csv"
            )
        ),
        genomes=SEEN_IDS,
    script:
        "../scripts/published_sets.py"


rule drift:
    input:
        sets=expand("results/published/{genome}.pfams", genome=SEEN_IDS),
        meta="results/published/accession_match.tsv",
        primary=expand(
            f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
            genome=SEEN_IDS,
        ),
        legacy=expand(
            "results/interpro/5.63-95.0/prodigal_single/{genome}.tsv",
            genome=SEEN_IDS,
        ),
        ready="results/models/READY",
    output:
        table="results/tables/drift.tsv",
        pfams="results/tables/drift_pfams.tsv",
    log:
        "logs/drift.log",
    benchmark:
        "benchmarks/drift.tsv"
    conda:
        "../envs/prediction.yaml"
    params:
        repo=str(REPO),
        evalue=EVALUE,
    script:
        "../scripts/drift.py"
