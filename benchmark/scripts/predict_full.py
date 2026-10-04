"""Full-genome prediction for one genome with the pinned models."""

import sys

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.common import write_json
from bacdive_workflow.prediction import predict_file

result = predict_file(
    snakemake.input.tsv,
    list(snakemake.params.traits),
    snakemake.params.model_dir,
    evalue=float(snakemake.params.evalue),
    sample_id=snakemake.wildcards.genome,
    model_cache={},
)
result["provenance"]["seed"] = None
write_json(snakemake.output.json, result)
