"""Analysis C: completeness degradation with contamination and shortcut validation."""


rule degrade:
    input:
        genes="results/genes/prodigal_single/{genome}.tsv",
        fna="results/genomes/{genome}.fna",
        tsv=f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
        full="results/predictions/full/{genome}.json",
        ready="results/models/READY",
        all_genes=expand(
            "results/genes/prodigal_single/{genome}.tsv", genome=GENOME_IDS
        ),
        all_tsvs=expand(
            f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
            genome=GENOME_IDS,
        ),
        genomes=GENOMES_FILE,
    output:
        table="results/degrade/{genome}.tsv",
    log:
        "logs/degrade/{genome}.log",
    benchmark:
        "benchmarks/degrade/{genome}.tsv"
    conda:
        "../envs/prediction.yaml"
    threads: 2
    resources:
        mem_mb=4096,
        runtime=120,
    params:
        repo=str(REPO),
        genome="{genome}",
        genome_ids=GENOME_IDS,
        levels=config["degrade"]["levels"],
        replicates=config["degrade"]["replicates"],
        fragment_median_bp=config["degrade"]["fragment_median_bp"],
        fragment_sigma=config["degrade"]["fragment_sigma"],
        min_fragment_bp=config["degrade"]["min_fragment_bp"],
        contamination=config["degrade"]["contamination"],
        contamination_levels=config["degrade"]["contamination_levels"],
        evalue=EVALUE,
        seed=SEED,
        model_dir="results/models",
    script:
        "../scripts/degrade.py"


rule degrade_summary:
    input:
        tables=expand("results/degrade/{genome}.tsv", genome=GENOME_IDS),
        predictions="results/tables/predictions_full.tsv",
        genomes=GENOMES_FILE,
        mapping=MAPPING_FILE,
    output:
        summary="results/tables/degrade_summary.tsv",
        thresholds="results/tables/degrade_thresholds.tsv",
    log:
        "logs/degrade_summary.log",
    conda:
        "../envs/stats.yaml"
    params:
        repo=str(REPO),
        seed=SEED,
    script:
        "../scripts/degrade_summary.py"


rule validate_shortcut:
    input:
        tables=expand("results/degrade/{genome}.tsv", genome=SHORTCUT_IDS),
        genes=expand("results/genes/prodigal_single/{genome}.tsv", genome=SHORTCUT_IDS),
        primary=expand(
            f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
            genome=SHORTCUT_IDS,
        ),
        genomes=expand("results/genomes/{genome}.fna", genome=SHORTCUT_IDS),
        ready="results/models/READY",
    output:
        table="results/tables/shortcut_validation.tsv",
    log:
        "logs/validate_shortcut.log",
    benchmark:
        "benchmarks/validate_shortcut.tsv"
    conda:
        "../envs/prodigal.yaml"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=240,
        ips=1,
    params:
        repo=str(REPO),
        genomes=SHORTCUT_IDS,
        test_mode=TEST_MODE,
        fixture_dir=FIXTURE_SHORTCUT,
        ips_executable=config["interproscan"][PRIMARY_IPS],
        ips_args=config["interproscan_args"],
        evalue=EVALUE,
        seed=SEED,
        fragment_median_bp=config["degrade"]["fragment_median_bp"],
        fragment_sigma=config["degrade"]["fragment_sigma"],
        min_fragment_bp=config["degrade"]["min_fragment_bp"],
    script:
        "../scripts/validate_shortcut.py"
