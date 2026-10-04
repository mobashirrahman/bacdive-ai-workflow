"""Prodigal gene calls and version-pure InterProScan Pfam annotation."""


rule prodigal:
    input:
        genome="results/genomes/{genome}.fna",
    output:
        proteins="results/proteins/{genes}/{genome}.faa",
    log:
        "logs/prodigal/{genes}_{genome}.log",
    benchmark:
        "benchmarks/prodigal/{genes}_{genome}.tsv"
    wildcard_constraints:
        genes="prodigal_single|prodigal_meta",
    conda:
        "../envs/prodigal.yaml"
    threads: 1
    resources:
        mem_mb=2048,
        runtime=30,
    params:
        mode=lambda w: "single" if w.genes == "prodigal_single" else "meta",
        strip=str(REPO / "scripts/clean_proteins.py"),
    shell:
        "prodigal -i {input.genome:q} -a {output.proteins:q} -p {params.mode:q} -q > {log:q} 2>&1 "
        "&& python {params.strip:q} {output.proteins:q}"


rule gene_table:
    input:
        proteins="results/proteins/{genes}/{genome}.faa",
    output:
        table="results/genes/{genes}/{genome}.tsv",
    log:
        "logs/genes/{genes}_{genome}.log",
    conda:
        "../envs/prediction.yaml"
    params:
        repo=str(REPO),
    script:
        "../scripts/gene_table.py"


rule interproscan:
    input:
        proteins="results/proteins/{genes}/{genome}.faa",
    output:
        tsv="results/interpro/{ips}/{genes}/{genome}.tsv",
    log:
        "logs/interproscan/{ips}_{genes}_{genome}.log",
    benchmark:
        "benchmarks/interproscan/{ips}_{genes}_{genome}.tsv"
    threads: 4
    resources:
        mem_mb=8192,
        runtime=240,
        ips=1,
    params:
        executable=lambda w: config["interproscan"][w.ips],
        extra=config["interproscan_args"],
    shell:
        "{params.executable:q} -i {input.proteins:q} -f tsv -o {output.tsv:q} "
        "-appl Pfam -cpu {threads} {params.extra} > {log:q} 2>&1"
