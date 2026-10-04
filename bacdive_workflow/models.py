"""Verify and extract only the eight known upstream model files."""

import argparse
import json
import os
import shutil
import tempfile
import urllib.request
import zipfile
from pathlib import Path

from .common import WorkflowError, sha256_file


def manifest():
    return json.loads(Path(__file__).with_name("model_manifest.json").read_text())


def prepare_models(archive, destination):
    release = manifest()
    if sha256_file(archive) != release["archive_sha256"]:
        raise WorkflowError("Model archive SHA-256 does not match the pinned upstream release.")
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        for name, expected in release["models"].items():
            target = destination / name
            if target.is_file() and sha256_file(target) == expected["sha256"]:
                continue
            member = bundle.getinfo("models/" + name)
            if member.file_size != expected["bytes"]:
                raise WorkflowError(f"Unexpected extracted size for {name}.")
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(
                    dir=destination, prefix=".model-", delete=False
                ) as handle:
                    temporary = Path(handle.name)
                    with bundle.open(member) as source:
                        shutil.copyfileobj(source, handle)
                if sha256_file(temporary) != expected["sha256"]:
                    raise WorkflowError(f"Extracted model checksum mismatch: {name}.")
                os.replace(temporary, target)
            finally:
                if temporary is not None and temporary.exists():
                    temporary.unlink()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare verified BacDive-AI v1.1 models.")
    parser.add_argument("--archive", type=Path, default=Path("Archiv.zip"))
    parser.add_argument("--model-dir", type=Path, default=Path("models"))
    parser.add_argument(
        "--download", action="store_true", help="Download the pinned archive if absent"
    )
    args = parser.parse_args(argv)
    try:
        if not args.archive.is_file() and args.download:
            args.archive.parent.mkdir(parents=True, exist_ok=True)
            with urllib.request.urlopen(manifest()["url"], timeout=120) as response:
                with tempfile.NamedTemporaryFile(dir=args.archive.parent, delete=False) as handle:
                    temporary = Path(handle.name)
                    try:
                        shutil.copyfileobj(response, handle)
                        handle.flush()
                        if sha256_file(temporary) != manifest()["archive_sha256"]:
                            raise WorkflowError("Downloaded archive checksum mismatch.")
                        os.replace(temporary, args.archive)
                    finally:
                        if temporary.exists():
                            temporary.unlink()
        prepare_models(args.archive, args.model_dir)
        print(f"Verified all eight models in {args.model_dir}.")
        return 0
    except (WorkflowError, OSError, KeyError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
