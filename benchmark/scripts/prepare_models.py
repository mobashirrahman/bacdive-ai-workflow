"""Stage verified models for the benchmark run."""

import sys

sys.path.insert(0, snakemake.params.repo)  # noqa: F821

from bacdive_workflow.models import prepare_models

prepare_models(snakemake.params.archive, snakemake.params.model_dir)
