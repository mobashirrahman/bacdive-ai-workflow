"""Efficient annotation-only batch inference; verified models are loaded once."""

import argparse
from pathlib import Path

from .cli import configure_logging
from .common import DEFAULT_EVALUE, TRAITS, WorkflowError, load_samples, write_json
from .prediction import predict_file


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Predict an annotation sample sheet with cached models."
    )
    parser.add_argument("--samples", required=True, type=Path)
    parser.add_argument("--model-dir", type=Path, default=Path("models"))
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--traits", nargs="+", choices=list(TRAITS), default=list(TRAITS))
    parser.add_argument("--evalue", type=float, default=DEFAULT_EVALUE)
    args = parser.parse_args(argv)
    configure_logging()
    try:
        samples = load_samples(args.samples)
        if any(row["input_type"] != "annotation" for row in samples.values()):
            raise WorkflowError(
                "Batch inference accepts annotation_path inputs; use Snakemake for genomes/proteins."
            )
        cache = {}
        for sample_id, row in sorted(samples.items()):
            result = predict_file(
                row["annotation_path"],
                args.traits,
                args.model_dir,
                args.evalue,
                sample_id,
                model_cache=cache,
            )
            write_json(args.output_dir / f"{sample_id}.json", result)
        print(f"Completed {len(samples)} samples and {len(args.traits)} traits per sample.")
        return 0
    except (WorkflowError, OSError, ImportError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
