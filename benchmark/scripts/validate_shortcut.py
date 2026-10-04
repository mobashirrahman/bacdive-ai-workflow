"""Shortcut validation: shortcut draws versus a real re-annotation.

Real path writes the kept fragments as FASTA, runs Prodigal meta and
InterProScan for real, and predicts. Test path reads fixture predictions.
Whatever the shortcut shows is reported.
"""

import csv
import json
import pickle
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.bench import fragments as frag
from bacdive_workflow.common import TRAITS, WorkflowError
from bacdive_workflow.prediction import parse_pfams, predict_one

params = snakemake.params
evalue = float(params.evalue)
levels = [70, 50]


def read_genes(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return [
            {
                "protein_id": r["protein_id"],
                "contig": r["contig"],
                "start": int(r["start"]),
                "end": int(r["end"]),
            }
            for r in csv.DictReader(handle, delimiter="\t")
        ]


def read_fasta(path):
    records = {}
    name = None
    chunks = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records[name] = "".join(chunks)
    return records


genes_of = {g: read_genes(p) for g, p in zip(params.genomes, snakemake.input.genes)}
seqs_of = {g: read_fasta(p) for g, p in zip(params.genomes, snakemake.input.genomes)}
lengths_of = {g: {c: len(s) for c, s in seqs.items()} for g, seqs in seqs_of.items()}

shortcut: dict = {}
for table in snakemake.input.tables:
    with open(table, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if (
                row["contamination"] == "0"
                and int(row["level"]) in levels
                and row["replicate"] == "0"
            ):
                shortcut[(row["genome"], int(row["level"]), row["trait"])] = row

model_dir = Path("results/models")
bundles = {}
for trait in TRAITS:
    with open(model_dir / f"{trait}_data.p", "rb") as handle:
        bundles[trait] = pickle.load(handle)

fp = {
    "fragment_median_bp": params.fragment_median_bp,
    "fragment_sigma": params.fragment_sigma,
    "min_fragment_bp": params.min_fragment_bp,
}
rows = []
for genome in params.genomes:
    for level in levels:
        seed = frag.draw_seed(params.seed, genome, level, 0, 0)
        kept = frag.kept_fragments(lengths_of[genome], level, 0, fp, seed)
        genes = genes_of[genome]
        kept_ids, _, _, _, _ = frag.simulate(genes, lengths_of[genome], level, 0, None, fp, seed)
        if params.test_mode:
            fixture = Path(params.fixture_dir) / f"{genome}_{level}.json"
            if not fixture.is_file():
                raise WorkflowError(f"Missing shortcut fixture {fixture}.")
            real = {
                t: (int(v["prediction"]), float(v["positive_probability"]))
                for t, v in json.loads(fixture.read_text()).items()
            }
        else:
            with tempfile.TemporaryDirectory() as directory:
                work = Path(directory)
                fasta = work / "kept.fna"
                with open(fasta, "w", encoding="utf-8") as handle:
                    for i, (contig, start, end) in enumerate(kept):
                        piece = seqs_of[genome][contig][start:end]
                        handle.write(f">kept_{i} {contig}:{start}-{end}\n{piece}\n")
                faa = work / "kept.faa"
                subprocess.run(
                    ["prodigal", "-i", str(fasta), "-a", str(faa), "-p", "meta", "-q"],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                # Strip stop characters: InterProScan rejects '*' in sequences.
                subprocess.run(
                    [sys.executable, params.strip, str(faa)],
                    check=True,
                    capture_output=True,
                    text=True,
                )
                if not any(line.startswith(">") for line in faa.read_text().splitlines()):
                    # No genes called on these fragments; nothing to re-annotate.
                    real = {}
                else:
                    tsv = work / "kept.tsv"
                    subprocess.run(
                        [
                            params.ips_executable,
                            "-i",
                            str(faa),
                            "-f",
                            "tsv",
                            "-o",
                            str(tsv),
                            "-appl",
                            "Pfam",
                            "-cpu",
                            str(snakemake.threads),
                            *params.ips_args.split(),
                        ],
                        check=True,
                        capture_output=True,
                        text=True,
                    )
                    # Only hits on kept proteins count; empty Pfam sets are
                    # reported as missing, never invented.
                    try:
                        pfams, _ = parse_pfams(tsv, evalue)
                    except WorkflowError:
                        pfams = set()
                    real = {}
                    for trait in TRAITS:
                        try:
                            result = predict_one(bundles[trait], pfams)
                            real[trait] = (
                                int(result["prediction"]),
                                result["positive_probability"],
                            )
                        except WorkflowError:
                            real[trait] = ("", "")
        for trait in TRAITS:
            short = shortcut.get((genome, level, trait))
            if short is None or short["prediction"] == "":
                continue
            real_pred, real_prob = real.get(trait, ("", ""))
            if real_pred == "":
                continue
            rows.append(
                {
                    "genome": genome,
                    "level": level,
                    "trait": trait,
                    "shortcut_prediction": short["prediction"],
                    "real_prediction": real_pred,
                    "agreement": int(short["prediction"]) == int(real_pred),
                    "probability_difference": abs(
                        float(short["positive_probability"]) - float(real_prob)
                    ),
                }
            )

with open(snakemake.output.table, "w", newline="", encoding="utf-8") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=[
            "genome",
            "level",
            "trait",
            "shortcut_prediction",
            "real_prediction",
            "agreement",
            "probability_difference",
        ],
        delimiter="\t",
    )
    writer.writeheader()
    writer.writerows(rows)
