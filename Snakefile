"""Genome, protein, or annotation inputs converge on validated trait predictions."""

import math
import sys
from pathlib import Path

ROOT = Path(workflow.basedir).resolve()
sys.path.insert(0, str(ROOT))

from bacdive_workflow.common import TRAITS, WorkflowError, load_samples

configfile: "config/config.yaml"

SAMPLES_FILE = config.get("samples", "examples/samples.tsv")
SAMPLES = load_samples(SAMPLES_FILE)
REFERENCE = config.get("metadata")
if REFERENCE and not Path(REFERENCE).is_file():
    raise WorkflowError(f"Reference metadata does not exist: {REFERENCE}")
SELECTED = config.get("traits", list(TRAITS))
if not isinstance(SELECTED, list) or not SELECTED or len(set(SELECTED)) != len(SELECTED) or any(t not in TRAITS for t in SELECTED):
    raise WorkflowError("traits must be a nonempty list of unique supported traits.")
EVALUE = float(config.get("evalue", 1e-20))
if not math.isfinite(EVALUE) or EVALUE <= 0:
    raise WorkflowError("evalue must be finite and positive.")
MODELS = [str(Path(config.get("model_dir", "models")) / f"{trait}_data.p") for trait in SELECTED]
PRODIGAL_MODE = config.get("prodigal_mode", "meta")
if PRODIGAL_MODE not in ("meta", "single"):
    raise WorkflowError("prodigal_mode must be meta or single.")
PYTHON_SOURCES = [str(p) for p in sorted((ROOT / "bacdive_workflow").glob("*.py"))]


def protein_input(wildcards):
    row = SAMPLES[wildcards.sample]
    return row.get("protein_path") or f"results/proteins/{wildcards.sample}.faa"


def annotation_input(wildcards):
    row = SAMPLES[wildcards.sample]
    return row.get("annotation_path") or f"results/interpro/{wildcards.sample}.tsv"


rule all:
    input:
        "results/bacdive_ai_predictions.csv",
        "results/report.html",
        "results/summary.json",


rule prodigal:
    input:
        genome=lambda w: SAMPLES[w.sample]["genome_path"],
    output:
        proteins="results/proteins/{sample}.faa",
    params:
        mode=PRODIGAL_MODE,
        strip=str(ROOT / "scripts/clean_proteins.py"),
    log:
        "logs/prodigal/{sample}.log",
    benchmark:
        "benchmarks/prodigal/{sample}.tsv"
    conda:
        "envs/prodigal.yaml"
    shell:
        # Portable: BSD/macOS sed requires an argument to -i, so filter via Python
        # instead of an in-place GNU-only substitution.
        "prodigal -i {input.genome:q} -a {output.proteins:q} -p {params.mode:q} -q > {log:q} 2>&1 "
        "&& python {params.strip:q} {output.proteins:q}"


rule interproscan:
    input:
        proteins=protein_input,
    output:
        tsv="results/interpro/{sample}.tsv",
    params:
        executable=config.get("interproscan", "interproscan.sh"),
    threads: 4
    resources:
        mem_mb=8192,
        runtime=120,
    log:
        "logs/interproscan/{sample}.log",
    benchmark:
        "benchmarks/interproscan/{sample}.tsv"
    shell:
        "{params.executable:q} -i {input.proteins:q} -f tsv -o {output.tsv:q} -appl Pfam -cpu {threads} > {log:q} 2>&1"


rule predict_traits:
    input:
        tsv=annotation_input,
        models=MODELS,
        sources=PYTHON_SOURCES,
    output:
        json="results/predictions/{sample}.json",
    params:
        script=str(ROOT / "predict.py"),
        model_dir=config.get("model_dir", "models"),
        evalue=EVALUE,
        traits=SELECTED,
    log:
        "logs/predict/{sample}.log",
    benchmark:
        "benchmarks/predict_traits/{sample}.tsv"
    conda:
        "envs/prediction.yaml"
    shell:
        "python {params.script:q} all {input.tsv:q} --model-dir {params.model_dir:q} --sample-id {wildcards.sample:q} --evalue {params.evalue} --output {output.json:q} --traits {params.traits:q} > {log:q} 2>&1"


rule aggregate_results:
    input:
        predictions=expand("results/predictions/{sample}.json", sample=sorted(SAMPLES)),
        samples=SAMPLES_FILE,
        metadata=[REFERENCE] if REFERENCE else [],
        sources=PYTHON_SOURCES,
    output:
        csv="results/bacdive_ai_predictions.csv",
    params:
        root=str(ROOT),
        traits=SELECTED,
    log:
        "logs/aggregate_results.log",
    conda:
        "envs/report.yaml"
    script:
        "scripts/aggregate_results.py"


rule report:
    input:
        csv="results/bacdive_ai_predictions.csv",
        sources=PYTHON_SOURCES,
    output:
        html="results/report.html",
        summary="results/summary.json",
    params:
        root=str(ROOT),
    log:
        "logs/report.log",
    conda:
        "envs/report.yaml"
    script:
        "scripts/report.py"
