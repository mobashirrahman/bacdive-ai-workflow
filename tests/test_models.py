"""Phase 1: v2 manifest provenance and nested-archive extraction."""

import hashlib
import io
import zipfile

import pytest

from bacdive_workflow import models
from bacdive_workflow.common import WorkflowError

V11_MODEL_HASHES = {
    "acidophile_data.p": "e49ef48b6ab1077bb159714b954999b2e61b92aea3a0edb02e28ef9a3660a1f1",
    "gram-positive_data.p": "288e5d618ffb275fb789ba4c0c25f9b85b644948473017302d69335e3a0b5af9",
    "spore-forming_data.p": "6d36821fca9f1573fcea2594f0f774708215758e6e2682bf3b5eebe0ba61f4f1",
    "aerobic_data.p": "495dfa48d91983fdc3456b970424f0509a130e31e9bd49e87deab322d8bb6b95",
    "anaerobic_data.p": "159c1abe9c63f2ec2054125f0fd4300da61c388f468dd3fea9b76b0c9b6197cc",
    "thermophile_data.p": "b2c9dab302cbc262de09b9e533e48431940d10ff65d88e6f2f3bc5668bf4bb95",
    "psychrophile_data.p": "6806852fcee3ffa02ccfd19720cf24484ead18bb3524e1f3fa09bf7cd9ab9ce3",
    "motile2+_data.p": "82b198e8b1b7af03201614891291a5d84065445ebcf212027a265cc4d6458bd7",
}


def test_manifest_models_match_v11_bytes():
    release = models.manifest()
    assert release["release"] == "v2"
    assert release["doi"] == "10.5281/zenodo.15075932"
    assert release["archive_sha256"] == (
        "f48d8a59f36f7925db03796b5557ae48b951e9307af38e0fc064cdf571e59cb2"
    )
    assert release["models_member"] == "models/models.zip"
    assert {name: entry["sha256"] for name, entry in release["models"].items()} == (
        V11_MODEL_HASHES
    )
    assert release["compatible_archives"][0]["archive_sha256"] == (
        "2d4b29a262b6813f88a31b838b050d89de6b5aa7ef24955ce04003e0333ebd11"
    )
    assert "training_data" in release


def _make_nested_archive(path, payloads, inner_name="models/models.zip"):
    inner_buffer = io.BytesIO()
    with zipfile.ZipFile(inner_buffer, "w", zipfile.ZIP_STORED) as inner:
        for name, data in payloads.items():
            inner.writestr(name, data)
    inner_bytes = inner_buffer.getvalue()
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as outer:
        outer.writestr(inner_name, inner_bytes)
    return hashlib.sha256(inner_bytes).hexdigest()


def test_nested_zip_extraction_without_disk_roundtrip(tmp_path, monkeypatch):
    payloads = {"acidophile_data.p": b"model-bytes-a"}
    archive = tmp_path / "synthetic-v2.zip"
    inner_sha = _make_nested_archive(archive, payloads)
    release = models.manifest()
    fake = {
        "release": "v2",
        "doi": "x",
        "url": "x",
        "archive_sha256": hashlib.sha256(archive.read_bytes()).hexdigest(),
        "models_member": "models/models.zip",
        "models_member_sha256": inner_sha,
        "models_member_bytes": None,  # filled below
        "models": {
            "acidophile_data.p": {
                "sha256": hashlib.sha256(b"model-bytes-a").hexdigest(),
                "bytes": len(b"model-bytes-a"),
            }
        },
        "training_data": release["training_data"],
        "compatible_archives": [],
    }
    with zipfile.ZipFile(archive) as bundle:
        fake["models_member_bytes"] = bundle.getinfo("models/models.zip").file_size
    monkeypatch.setattr(models, "manifest", lambda: fake)
    destination = tmp_path / "models"
    models.prepare_models(archive, destination)
    assert (destination / "acidophile_data.p").read_bytes() == b"model-bytes-a"


def test_wrong_hash_archive_is_rejected(tmp_path, monkeypatch):
    archive = tmp_path / "bad.zip"
    with zipfile.ZipFile(archive, "w") as bundle:
        bundle.writestr("models/acidophile_data.p", b"nope")
    monkeypatch.setattr(models, "manifest", models.manifest)
    with pytest.raises(WorkflowError, match="SHA-256"):
        models.prepare_models(archive, tmp_path / "out")
