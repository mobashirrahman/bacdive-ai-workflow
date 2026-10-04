"""Stage the verified models and fetch genome assemblies."""


rule prepare_models:
    output:
        sentinel=touch("results/models/READY"),
    log:
        "logs/prepare_models.log",
    conda:
        "../envs/prediction.yaml"
    params:
        archive=model_archive_path(),
        model_dir="results/models",
        repo=str(REPO),
    script:
        "../scripts/prepare_models.py"


rule download_genome:
    output:
        fna="results/genomes/{genome}.fna",
        faa="results/genomes/{genome}.pgap.faa",
        report="results/genomes/{genome}.json",
    log:
        "logs/download/{genome}.log",
    benchmark:
        "benchmarks/download/{genome}.tsv"
    retries: 3
    conda:
        "../envs/datasets.yaml"
    params:
        repo=str(REPO),
    script:
        "../scripts/download_genome.py"
