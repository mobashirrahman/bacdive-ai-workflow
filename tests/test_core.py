import json

import numpy as np
import pytest

from bacdive_workflow.aggregation import aggregate
from bacdive_workflow.cli import main
from bacdive_workflow.common import (
    TRAITS,
    WorkflowError,
    index_metadata,
    load_samples,
    normalize_id,
)
from bacdive_workflow.prediction import parse_pfams, predict_one
from bacdive_workflow.report import comparison, render_html, summarize


def annotation(accession="PF00001", evalue="1e-30", analysis="Pfam"):
    return f"protein\tmd5\t100\t{analysis}\t{accession}\tdescription\t1\t100\t{evalue}\tT\t01-01-2024\n"


@pytest.mark.parametrize(
    "value,expected",
    [
        ("/a/GCF_0123.1.fna", "GCF_0123.1"),
        ("AB_Plaque.25.fa", "AB_Plaque.25"),
        ("GCF_0123.1.fna.faa.tsv", "GCF_0123.1"),
        ("sample.fasta.gz", "sample"),
        ("sample.json", "sample"),
        ("version.1", "version.1"),
    ],
)
def test_identity_preserves_accession_versions(value, expected):
    assert normalize_id(value) == expected


def test_pfam_filtering_and_presence_encoding(tmp_path):
    path = tmp_path / "input.tsv"
    path.write_text(
        annotation()
        + annotation()
        + annotation("PF00002", "1e-10")
        + annotation("PS00001", "-", "Prosite")
    )
    pfams, counts = parse_pfams(path)
    assert pfams == {"PF00001"}
    assert counts == {
        "rows": 4,
        "pfam_rows": 3,
        "retained_rows": 2,
        "other_analysis_rows": 1,
        "unique_pfams": 1,
    }


@pytest.mark.parametrize(
    "content",
    [
        "",
        "bad\trow\n",
        annotation(evalue="nan"),
        annotation(evalue="-1"),
        annotation(evalue="text"),
        annotation(accession="unexpected"),
    ],
)
def test_invalid_annotation_is_rejected(tmp_path, content):
    path = tmp_path / "input.tsv"
    path.write_text(content)
    with pytest.raises(WorkflowError):
        parse_pfams(path)


class ReversedClasses:
    classes_ = np.array([1, 0])

    def predict_proba(self, data):
        assert data == [[1, 0]]
        return np.array([[0.9, 0.1]])


def test_argmax_is_mapped_through_class_labels():
    result = predict_one(
        {"model": ReversedClasses(), "categories": ["PF00001", "PF00002"]}, {"PF00001"}
    )
    assert result["prediction"] is True
    assert result["confidence"] == 90.0
    assert result["positive_probability"] == 0.9


class RecordingModel:
    """Records the feature vector so the upstream encoding stays pinned."""

    classes_ = np.array([0, 1])

    def __init__(self):
        self.seen = None

    def predict_proba(self, data):
        self.seen = data
        return np.array([[0.25, 0.75]])


def test_features_are_presence_flags_not_hit_counts():
    # Upstream's `np.unique(list(pfams), return_counts=True)` runs on an already
    # deduplicated set, so its counts are always 1. Real hit counts shift
    # predictions: the published example reports 99.85% for motility, not 99.35%.
    model = RecordingModel()
    bundle = {"model": model, "categories": ["PF00001", "PF00002", "PF00003"]}
    result = predict_one(bundle, {"PF00001", "PF00002"})
    assert model.seen == [[1, 1, 0]]
    assert result["matched_features"] == 2
    assert result["model_features"] == 3


def test_missing_models_fail_without_creating_output(tmp_path, capsys):
    source = tmp_path / "source.tsv"
    source.write_text(annotation())
    output = tmp_path / "output.json"
    assert main(["all", str(source), "--model-dir", str(tmp_path), "--output", str(output)]) == 1
    assert not output.exists()
    assert capsys.readouterr().out == ""


def legacy_result(tmp_path, name="sample"):
    path = tmp_path / f"{name}.json"
    path.write_text(
        json.dumps(
            {
                "predictions": {
                    label: {"prediction": False, "confidence": 90.0} for label in TRAITS.values()
                }
            }
        )
    )
    return path


def test_metadata_join_supports_legacy_extensions_and_columns(tmp_path):
    metadata = tmp_path / "metadata.csv"
    metadata.write_text(
        "genome_address,taxon,gram stain,oxygen tolerance\n/a/sample.fa,Example,positive,anaerobe\n"
    )
    rows = aggregate([legacy_result(tmp_path)], metadata=metadata)
    assert rows[0]["metadata_matched"] is True
    assert rows[0]["gram_stain"] == "positive"
    assert rows[0]["oxygen_tolerance"] == "anaerobe"
    assert rows[0]["taxon"] == "Example"


def test_incomplete_and_corrupt_predictions_fail(tmp_path):
    path = legacy_result(tmp_path)
    data = json.loads(path.read_text())
    data["predictions"].pop("Aerobic")
    path.write_text(json.dumps(data))
    with pytest.raises(WorkflowError, match="expected exactly"):
        aggregate([path])
    path.write_text("not json")
    with pytest.raises(WorkflowError, match="Cannot read prediction"):
        aggregate([path])


def test_duplicate_samples_and_missing_predictions_fail(tmp_path):
    path = legacy_result(tmp_path)
    with pytest.raises(WorkflowError, match="duplicate"):
        aggregate([path, path])
    metadata = tmp_path / "metadata.tsv"
    metadata.write_text("sample_id\ttaxon\nsample\tExample\nmissing\tAnother\n")
    with pytest.raises(WorkflowError, match="missing for samples"):
        aggregate([path], samples=metadata)
    metadata.write_text("sample_id\ttaxon\nsample\tA\nsample\tB\n")
    with pytest.raises(WorkflowError, match="duplicate"):
        index_metadata(metadata)


def test_sample_paths_are_relative_to_sample_sheet(tmp_path, monkeypatch):
    sheet_dir = tmp_path / "directory with spaces"
    sheet_dir.mkdir()
    (sheet_dir / "annotation.tsv").write_text(annotation())
    sheet = sheet_dir / "samples.tsv"
    sheet.write_text("sample_id\tannotation_path\nsample\tannotation.tsv\n")
    monkeypatch.chdir(tmp_path)
    assert load_samples(sheet)["sample"]["annotation_path"] == str(sheet_dir / "annotation.tsv")


@pytest.mark.parametrize(
    "content",
    [
        "sample_id\tannotation_path\n",
        "sample_id\tannotation_path\nx\tmissing.tsv\n",
        "sample_id\tannotation_path\nx\tx\textra\n",
    ],
)
def test_bad_sample_sheets_fail(tmp_path, content):
    sheet = tmp_path / "samples.tsv"
    sheet.write_text(content)
    with pytest.raises(WorkflowError):
        load_samples(sheet)


def test_ambiguous_reference_labels_are_excluded():
    rows = [
        {"sample_id": "a", "gram_stain": "positive", "gram-positive_prediction": True},
        {"sample_id": "b", "gram_stain": "variable", "gram-positive_prediction": False},
        {"sample_id": "c", "gram_stain": "", "gram-positive_prediction": False},
    ]
    result = comparison(rows, "gram-positive", "gram_stain", {"positive": True, "negative": False})
    assert result["eligible_samples"] == 1
    assert result["excluded_samples"] == 2
    assert result["agreement"] == 1


def test_exclusion_reasons_are_reported_separately():
    # A sample with no usable label and a sample with an ambiguous label are both
    # excluded, but for different reasons and must not be conflated.
    rows = [
        {"sample_id": "a", "gram_stain": "positive", "gram-positive_prediction": True},
        {"sample_id": "b", "gram_stain": "variable", "gram-positive_prediction": False},
        {"sample_id": "c", "gram_stain": "", "gram-positive_prediction": False},
        {"sample_id": "d", "gram_stain": "negative,positive", "gram-positive_prediction": False},
    ]
    result = comparison(rows, "gram-positive", "gram_stain", {"positive": True, "negative": False})
    assert result["eligible_samples"] == 1
    assert result["samples_missing_reference_label"] == 1
    assert result["samples_with_unmapped_reference_label"] == 2
    assert (
        result["excluded_samples"]
        == result["samples_missing_reference_label"]
        + result["samples_with_unmapped_reference_label"]
    )


def test_reference_labels_are_case_insensitive():
    # The source table mixes "anaerobe/microaerophile" and "Anaerobe/Microaerophile".
    mapping = {"aerobe": True, "anaerobe": False}
    rows = [
        {"sample_id": "a", "oxygen_tolerance": "Aerobe", "aerobic_prediction": True},
        {"sample_id": "b", "oxygen_tolerance": "anaerobe", "aerobic_prediction": False},
        {
            "sample_id": "c",
            "oxygen_tolerance": "Anaerobe/Microaerophile",
            "aerobic_prediction": True,
        },
    ]
    result = comparison(rows, "aerobic", "oxygen_tolerance", mapping)
    assert result["eligible_samples"] == 2
    assert result["samples_with_unmapped_reference_label"] == 1
    assert result["agreement"] == 1


def test_interproscan_header_row_is_tolerated(tmp_path):
    # Real `interproscan.sh -f tsv` output starts with a commented header line.
    header = (
        "#EVIDENCE_DATABASE\tSEQUENCE_LENGTH\tSEQUENCE_ACC\tANALYSIS\tSIGNATURE_ACCESSION\t"
        "SIGNATURE_DESCRIPTION\tSTART_LOCATION\tSTOP_LOCATION\tSCORE\tSTATUS\tDATE\t"
        "INTERPRO_ANNOTATION\tINTERPRO_DESCRIPTION\tGO_ANNOTATION\tPATHWAYS\n"
    )
    path = tmp_path / "input.tsv"
    path.write_text(header + annotation())
    pfams, counts = parse_pfams(path)
    assert pfams == {"PF00001"}
    # The header is counted as a row but classified as a non-Pfam analysis, so it
    # never reaches the Pfam accession check.
    assert counts["rows"] == 2
    assert counts["pfam_rows"] == 1
    assert counts["retained_rows"] == 1
    assert counts["other_analysis_rows"] == 1


def test_report_rejects_non_numeric_quality_fields(tmp_path):
    rows = aggregate([legacy_result(tmp_path)])
    rows[0]["completeness"] = "not-a-number"
    with pytest.raises(WorkflowError, match="completeness"):
        summarize(rows)


def test_report_escapes_external_metadata(tmp_path):
    rows = aggregate([legacy_result(tmp_path)])
    rows[0]["taxon"] = "<script>alert('x')</script>"
    output = render_html(rows, summarize(rows), "Example <title>")
    assert "&lt;script&gt;alert" in output
    assert "<script>alert" not in output
    assert "Example &lt;title&gt;" in output
