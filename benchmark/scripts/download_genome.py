"""Fetch one assembly with NCBI datasets; fail loudly on empty results."""

import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.common import WorkflowError, write_json

accession = snakemake.wildcards.genome
fna = Path(snakemake.output.fna)
faa = Path(snakemake.output.faa)

with tempfile.TemporaryDirectory() as directory:
    package = Path(directory) / "package.zip"
    try:
        subprocess.run(
            [
                "datasets",
                "download",
                "genome",
                "accession",
                accession,
                "--include",
                "genome,protein",
                "--filename",
                str(package),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as error:
        raise WorkflowError(f"Cannot download {accession}: {error}") from error
    try:
        names = zipfile.ZipFile(package).namelist()
    except zipfile.BadZipFile as error:
        raise WorkflowError(f"Download package for {accession} is corrupt.") from error
    genomes = [n for n in names if n.endswith((".fna", ".fna.gz"))]
    proteins = [n for n in names if n.endswith((".faa", ".faa.gz"))]
    if not genomes:
        raise WorkflowError(f"No genome FASTA in the download package for {accession}.")

    def md5_ok(bundle):
        try:
            expected = {
                line.split()[1]: line.split()[0]
                for line in (bundle.read("md5sum.txt").decode().splitlines())
            }
        except KeyError:
            return True
        import hashlib

        for name, digest in expected.items():
            if name in ("md5sum.txt",):
                continue
            try:
                data = bundle.read(name)
            except KeyError:
                continue
            if hashlib.md5(data).hexdigest() != digest:
                return False
        return True

    with zipfile.ZipFile(package) as bundle:
        if not md5_ok(bundle):
            raise WorkflowError(f"md5 verification failed for {accession}.")
        with bundle.open(genomes[0]) as source:
            fna.parent.mkdir(parents=True, exist_ok=True)
            fna.write_bytes(source.read())
        if proteins:
            with bundle.open(proteins[0]) as source:
                faa.write_bytes(source.read())
        else:
            faa.write_bytes(b"")
    report = {
        "accession": accession,
        "genome_file": genomes[0],
        "protein_file": proteins[0] if proteins else None,
        "pgap_available": bool(proteins),
    }
if fna.stat().st_size == 0:
    raise WorkflowError(f"Downloaded genome for {accession} is empty.")
write_json(snakemake.output.report, report)
