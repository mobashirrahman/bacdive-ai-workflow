import os
import glob
import pandas as pd
from pathlib import Path

# Load metadata to get genome file paths
metadata_file = "data/merged_genomes_data.csv"
if os.path.exists(metadata_file):
    metadata = pd.read_csv(metadata_file)
    GENOMES = metadata["genome_address"].tolist()
    # Extract just the filename part for outputs
    GENOME_IDS = [os.path.basename(path).replace('.fa', '') for path in GENOMES]
else:
    # Fallback if no metadata file
    GENOMES = []
    GENOME_IDS = []
    print(f"Warning: {metadata_file} not found. No genomes defined.")

# Output directory configuration
INTERPRO_DIR = "results/interpro_results"
RESULTS_DIR = "results"
PROTEIN_DIR = "results/proteins"

# Create output directories
os.makedirs(INTERPRO_DIR, exist_ok=True)
os.makedirs(RESULTS_DIR, exist_ok=True)
os.makedirs(PROTEIN_DIR, exist_ok=True)

# Define the final target rule
rule all:
    input:
        os.path.join(RESULTS_DIR, "bacdive_ai_predictions.csv")
    benchmark:
        os.path.join("benchmarks", "all.txt")
    message: "Finalizing all predictions."


rule prodigal:
    input:
        genome=lambda w: next((p for p, gid in zip(GENOMES, GENOME_IDS) if gid == w.genome_id), None)
    output:
        proteins=os.path.join(PROTEIN_DIR, "{genome_id}.faa")
    log:
        os.path.join("logs", "prodigal", "{genome_id}.log")
    benchmark:
        os.path.join("benchmarks", "prodigal", "{genome_id}.txt")
    conda:
        "envs/prodigal.yaml"
    message: "Running Prodigal for gene prediction on {wildcards.genome_id}."
    shell:
        """
        mkdir -p {PROTEIN_DIR} logs/prodigal
        {{
            prodigal -i {input.genome} -a {output.proteins} -p meta -q;
            sed -i 's/\\*//g' {output.proteins};
        }} > {log} 2>&1
        """ 


# Rule to run InterProScan on a genome file
rule interproscan:
    input:
        proteins=os.path.join(PROTEIN_DIR, "{genome_id}.faa")
    output:
        tsv=os.path.join(INTERPRO_DIR, "{genome_id}.faa.tsv")
    threads: 4
    log:
        os.path.join("logs", "interproscan", "{genome_id}.log")
    benchmark:
        os.path.join("benchmarks", "interproscan", "{genome_id}.txt")
    message: "Running InterProScan (Pfam) on predicted proteins for {wildcards.genome_id}."
    shell:
        """
        mkdir -p logs/interproscan
        interproscan.sh -i {input.proteins} -f tsv -d {INTERPRO_DIR} \
                        -appl Pfam -cpu {threads} > {log} 2>&1
        """

# Rule to predict traits for a single genome
rule predict_traits:
    input:
        tsv = os.path.join(INTERPRO_DIR, "{genome_id}.faa.tsv")
    output:
        json = os.path.join(RESULTS_DIR, "predictions", "{genome_id}.json")
    log:
        os.path.join("logs", "predict", "{genome_id}.log")
    benchmark:
        os.path.join("benchmarks", "predict_traits", "{genome_id}.txt")
    message: "Predicting traits for {wildcards.genome_id} based on InterProScan results."
    shell:
        """
        mkdir -p {RESULTS_DIR}/predictions logs/predict
        python predict.py all {input.tsv} > {output.json} 2> {log}
        """

# Rule to aggregate all prediction results
rule aggregate_results:
    input:
        predictions = expand(os.path.join(RESULTS_DIR, "predictions", "{genome_id}.json"), genome_id=GENOME_IDS)
    output:
        csv = os.path.join(RESULTS_DIR, "bacdive_ai_predictions.csv")
    log:
        os.path.join("logs", "aggregate_results.log")
    benchmark:
        os.path.join("benchmarks", "aggregate_results.txt")
    message: "Aggregating all individual genome predictions into a final CSV."
    script:
        "scripts/aggregate_results.py"


