"""Verify and extract only the known upstream model files (v2 pinned, v1.1 compatible)."""

import argparse
import hashlib
import io
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


def _identify_archive(archive_sha):
    release = manifest()
    if archive_sha == release["archive_sha256"]:
        return "v2"
    for entry in release.get("compatible_archives", []):
        if archive_sha == entry["archive_sha256"]:
            return entry["release"]
    raise WorkflowError("Model archive SHA-256 does not match the pinned upstream release.")


def prepare_models(archive, destination):
    archive = Path(archive)
    kind = _identify_archive(sha256_file(archive))
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    release = manifest()
    with zipfile.ZipFile(archive) as bundle:
        if kind == "v2":
            try:
                inner_info = bundle.getinfo(release["models_member"])
            except KeyError as error:
                raise WorkflowError(f"Archive is missing {release['models_member']}.") from error
            if inner_info.file_size != release["models_member_bytes"]:
                raise WorkflowError(f"Unexpected size for {release['models_member']}.")
            with bundle.open(inner_info) as handle:
                raw = handle.read()
            if hashlib.sha256(raw).hexdigest() != release["models_member_sha256"]:
                raise WorkflowError(
                    f"Nested archive checksum mismatch: {release['models_member']}."
                )
            with zipfile.ZipFile(io.BytesIO(raw)) as inner:
                for name, expected in release["models"].items():
                    _extract_one(inner, name, expected, destination)
        else:
            for name, expected in release["models"].items():
                _extract_one(bundle, "models/" + name, expected, destination, target_name=name)


def _extract_one(bundle, member_name, expected, destination, target_name=None):
    target = destination / (target_name or Path(member_name).name)
    if target.is_file() and sha256_file(target) == expected["sha256"]:
        return
    try:
        member = bundle.getinfo(member_name)
    except KeyError as error:
        raise WorkflowError(f"Archive is missing model {member_name}.") from error
    if member.file_size != expected["bytes"]:
        raise WorkflowError(f"Unexpected extracted size for {target.name}.")
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=destination, prefix=".model-", delete=False) as handle:
            temporary = Path(handle.name)
            with bundle.open(member) as source:
                shutil.copyfileobj(source, handle)
        if sha256_file(temporary) != expected["sha256"]:
            raise WorkflowError(f"Extracted model checksum mismatch: {target.name}.")
        os.replace(temporary, target)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def extract_training_data(archive, destination):
    """Extract and verify the two training CSVs from a v2 archive."""
    archive = Path(archive)
    if _identify_archive(sha256_file(archive)) != "v2":
        raise WorkflowError("Training data is only available in the v2 archive.")
    release = manifest()
    training = release["training_data"]
    destination = Path(destination)
    destination.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(archive) as bundle:
        try:
            outer_info = bundle.getinfo(training["member"])
        except KeyError as error:
            raise WorkflowError(f"Archive is missing {training['member']}.") from error
        if outer_info.file_size != training["bytes"]:
            raise WorkflowError(f"Unexpected size for {training['member']}.")
        with bundle.open(outer_info) as handle:
            raw = handle.read()
        if hashlib.sha256(raw).hexdigest() != training["sha256"]:
            raise WorkflowError(f"Nested archive checksum mismatch: {training['member']}.")
        with zipfile.ZipFile(io.BytesIO(raw)) as inner:
            for member_name, expected in training["files"].items():
                target = destination / Path(member_name).name
                if target.is_file() and sha256_file(target) == expected["sha256"]:
                    continue
                try:
                    member = inner.getinfo(member_name)
                except KeyError as error:
                    raise WorkflowError(f"Training archive is missing {member_name}.") from error
                if member.file_size != expected["bytes"]:
                    raise WorkflowError(f"Unexpected extracted size for {member_name}.")
                temporary = None
                try:
                    with tempfile.NamedTemporaryFile(
                        dir=destination, prefix=".training-", delete=False
                    ) as handle:
                        temporary = Path(handle.name)
                        with inner.open(member) as source:
                            shutil.copyfileobj(source, handle)
                    if sha256_file(temporary) != expected["sha256"]:
                        raise WorkflowError(f"Extracted training checksum mismatch: {member_name}.")
                    os.replace(temporary, target)
                finally:
                    if temporary is not None and temporary.exists():
                        temporary.unlink()


def _default_archive():
    if Path("Archiv.zip").is_file():
        return Path("Archiv.zip")
    return Path("v2.zip")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Prepare verified BacDive-AI v2 models.")
    parser.add_argument("--archive", type=Path, default=None)
    parser.add_argument("--model-dir", type=Path, default=Path("models"))
    parser.add_argument(
        "--download", action="store_true", help="Download the pinned archive if absent"
    )
    parser.add_argument(
        "--training-data",
        type=Path,
        default=None,
        help="Extract and verify the v2 training CSVs into DIR",
    )
    args = parser.parse_args(argv)
    try:
        if args.archive is None:
            args.archive = _default_archive()
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
        if args.training_data is not None:
            extract_training_data(args.archive, args.training_data)
        print(f"Verified all eight models in {args.model_dir}.")
        return 0
    except (WorkflowError, OSError, KeyError, zipfile.BadZipFile) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
