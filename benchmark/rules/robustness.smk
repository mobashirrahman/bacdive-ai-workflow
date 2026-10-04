"""Analysis D: annotation robustness across E-value, gene caller, and Pfam release."""


checkpoint survey_pgap:
    input:
        reports=expand("results/genomes/{genome}.json", genome=ROBUSTNESS_IDS),
    output:
        table="results/tables/pgap_available.tsv",
    log:
        "logs/survey_pgap.log",
    conda:
        "../envs/prediction.yaml"
    params:
        repo=str(REPO),
        genomes=ROBUSTNESS_IDS,
    script:
        "../scripts/survey_pgap.py"


rule stage_pgap:
    input:
        faa="results/genomes/{genome}.pgap.faa",
    output:
        staged="results/proteins/pgap/{genome}.faa",
    log:
        "logs/stage_pgap/{genome}.log",
    shell:
        "cp {input.faa:q} {output.staged:q}"


rule robustness:
    input:
        single=expand(
            f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
            genome=ROBUSTNESS_IDS,
        ),
        meta=expand(
            f"results/interpro/{PRIMARY_IPS}/prodigal_meta/{{genome}}.tsv",
            genome=ROBUSTNESS_IDS,
        ),
        legacy=expand(
            "results/interpro/5.63-95.0/prodigal_single/{genome}.tsv",
            genome=ROBUSTNESS_IDS,
        ),
        pgap=pgap_inputs,
        ready="results/models/READY",
        genomes=GENOMES_FILE,
    output:
        table="results/tables/robustness.tsv",
    log:
        "logs/robustness.log",
    benchmark:
        "benchmarks/robustness.tsv"
    conda:
        "../envs/prediction.yaml"
    params:
        repo=str(REPO),
        genomes=ROBUSTNESS_IDS,
        test_mode=TEST_MODE,
        evalues=config["evalues_robustness"],
        primary_evalue=EVALUE,
        seed=SEED,
        model_dir="results/models",
    script:
        "../scripts/robustness.py"
