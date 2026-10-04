"""Command line interface with stderr logging and atomic file outputs."""

import argparse
import json
import logging
import sys
from pathlib import Path

from . import __version__
from .common import DEFAULT_EVALUE, TRAITS, WorkflowError, validate_sample_id, write_json


def configure_logging(debug=False, log_file=None):
    handlers = [logging.StreamHandler(sys.stderr)]
    if log_file:
        Path(log_file).parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file, encoding="utf-8"))
    logging.basicConfig(
        level=logging.DEBUG if debug else logging.INFO,
        format="%(levelname)s | %(message)s",
        handlers=handlers,
        force=True,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description="Predict traits with published BacDive-AI models.")
    parser.add_argument("trait", choices=["all", *TRAITS])
    parser.add_argument("file", type=Path, help="InterProScan TSV containing Pfam annotations")
    parser.add_argument("--model-dir", type=Path, default=Path("models"))
    parser.add_argument("--sample-id")
    parser.add_argument(
        "--traits", nargs="+", choices=list(TRAITS), help="Select a subset when trait is all"
    )
    parser.add_argument("--evalue", type=float, default=DEFAULT_EVALUE)
    parser.add_argument("--output", type=Path, help="Write complete JSON atomically to this file")
    parser.add_argument("--text-output", action="store_true")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--log-file", type=Path)
    parser.add_argument("--version", action="version", version=__version__)
    args = parser.parse_args(argv)
    if args.output and args.text_output:
        parser.error("--output and --text-output cannot be used together")
    if args.traits and args.trait != "all":
        parser.error("--traits requires the positional trait all")
    try:
        configure_logging(args.debug, args.log_file)
        if args.sample_id:
            validate_sample_id(args.sample_id)
        from .prediction import predict_file

        selected = (args.traits or list(TRAITS)) if args.trait == "all" else [args.trait]
        result = predict_file(args.file, selected, args.model_dir, args.evalue, args.sample_id)
        if args.output:
            write_json(args.output, result)
        elif args.text_output:
            for label, prediction in result["predictions"].items():
                print(f"{label}: {prediction['prediction']} ({prediction['confidence']}%)")
        else:
            print(json.dumps(result, indent=2, allow_nan=False))
        return 0
    except (WorkflowError, OSError, ImportError) as error:
        logging.error("%s", error)
        return 1
