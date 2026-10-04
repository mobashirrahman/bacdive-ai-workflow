"""Full-genome prediction with the pinned models and collection to long format."""


rule predict_full:
    input:
        tsv=f"results/interpro/{PRIMARY_IPS}/prodigal_single/{{genome}}.tsv",
        ready="results/models/READY",
    output:
        json="results/predictions/full/{genome}.json",
    log:
        "logs/predict/{genome}.log",
    benchmark:
        "benchmarks/predict/{genome}.tsv"
    conda:
        "../envs/prediction.yaml"
    params:
        repo=str(REPO),
        model_dir="results/models",
        evalue=EVALUE,
        traits=list(TRAITS8),
    script:
        "../scripts/predict_full.py"


rule collect_predictions:
    input:
        predictions=expand("results/predictions/full/{genome}.json", genome=GENOME_IDS),
    output:
        table="results/tables/predictions_full.tsv",
    log:
        "logs/collect_predictions.log",
    conda:
        "../envs/stats.yaml"
    params:
        repo=str(REPO),
    script:
        "../scripts/collect_predictions.py"


rule accuracy:
    input:
        predictions="results/tables/predictions_full.tsv",
        genomes=GENOMES_FILE,
        mapping=MAPPING_FILE,
    output:
        accuracy="results/tables/accuracy.tsv",
        disagreements="results/tables/disagreements.tsv",
    log:
        "logs/accuracy.log",
    benchmark:
        "benchmarks/accuracy.tsv"
    conda:
        "../envs/stats.yaml"
    params:
        repo=str(REPO),
        bootstrap=BOOTSTRAP,
        seed=SEED,
    script:
        "../scripts/accuracy.py"
