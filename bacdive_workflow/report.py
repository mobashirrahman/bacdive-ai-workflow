"""Portable descriptive HTML reports and conservative reference comparisons."""

import argparse
import csv
import html
import math
import statistics
from pathlib import Path

from .common import TRAITS, WorkflowError, atomic_text, sha256_file, write_json


def as_bool(value):
    if value in (True, "True", "true"):
        return True
    if value in (False, "False", "false"):
        return False
    raise WorkflowError(f"Expected a boolean, received {value!r}.")


def comparison(rows, trait, field, mapping):
    tp = tn = fp = fn = 0
    disagreements = []
    missing_label = 0
    unmapped_label = 0
    for row in rows:
        label = row.get(field, "").strip().lower()
        if not label:
            missing_label += 1
            continue
        if label not in mapping:
            unmapped_label += 1
            continue
        truth = mapping[label]
        prediction = as_bool(row[f"{trait}_prediction"])
        if truth and prediction:
            tp += 1
        elif not truth and not prediction:
            tn += 1
        elif prediction:
            fp += 1
        else:
            fn += 1
        if truth != prediction:
            disagreements.append(row["sample_id"])
    count = tp + tn + fp + fn
    return {
        "eligible_samples": count,
        "excluded_samples": missing_label + unmapped_label,
        # Kept separate so the two reasons for exclusion are not conflated. A
        # source "NA" is cleaned to empty upstream and lands in the first bucket.
        "samples_missing_reference_label": missing_label,
        "samples_with_unmapped_reference_label": unmapped_label,
        "true_positive": tp,
        "true_negative": tn,
        "false_positive": fp,
        "false_negative": fn,
        "agreement": (tp + tn) / count if count else None,
        "precision": tp / (tp + fp) if tp + fp else None,
        "recall": tp / (tp + fn) if tp + fn else None,
        "f1": 2 * tp / (2 * tp + fp + fn) if 2 * tp + fp + fn else None,
        "disagreement_sample_ids": disagreements,
        "reference_mapping": mapping,
        "interpretation": "Descriptive agreement; training overlap and strain identity are unverified.",
    }


def summarize(rows):
    if not rows or len({r["sample_id"] for r in rows}) != len(rows):
        raise WorkflowError("Report input is empty or contains duplicate sample IDs.")
    selected = [trait for trait in TRAITS if f"{trait}_prediction" in rows[0]]
    if not selected:
        raise WorkflowError("Report input has no supported trait columns.")
    output = {
        "samples": len(rows),
        "metadata_matches": sum(as_bool(r["metadata_matched"]) for r in rows),
        "traits": {},
        "reference_comparisons": {},
        "interpretation": "Published model predictions, not experimental phenotype measurements.",
    }
    required = {"metadata_matched", *(f"{t}_prediction" for t in selected)}
    missing_columns = required - set(rows[0])
    if missing_columns:
        raise WorkflowError(
            "Report input is missing columns: " + ", ".join(sorted(missing_columns))
        )
    for trait in selected:
        positives = sum(as_bool(r[f"{trait}_prediction"]) for r in rows)
        try:
            confidences = [float(r[f"{trait}_confidence"]) for r in rows]
        except (KeyError, ValueError) as error:
            raise WorkflowError(f"Invalid or missing confidence for {trait}: {error}") from error
        if any(not math.isfinite(c) or not 50 <= c <= 100 for c in confidences):
            raise WorkflowError(f"Invalid confidence values for {trait}.")
        output["traits"][trait] = {
            "label": TRAITS[trait],
            "positive": positives,
            "negative": len(rows) - positives,
            "mean_label_probability_percent": round(statistics.mean(confidences), 2),
        }
    if "gram-positive" in selected:
        output["reference_comparisons"]["gram-positive"] = comparison(
            rows, "gram-positive", "gram_stain", {"positive": True, "negative": False}
        )
    # Compare only explicit pure oxygen labels. Intermediate/combined labels are excluded.
    pure_oxygen = {"aerobe": True, "aerobic": True, "anaerobe": False, "anaerobic": False}
    for trait in ("aerobic", "anaerobic"):
        if trait in selected:
            mapping = (
                pure_oxygen if trait == "aerobic" else {k: not v for k, v in pure_oxygen.items()}
            )
            output["reference_comparisons"][trait] = comparison(
                rows, trait, "oxygen_tolerance", mapping
            )
    for field in ("completeness", "contamination"):
        values = []
        for row in rows:
            raw = row.get(field, "").strip()
            if not raw:
                continue
            try:
                values.append(float(raw))
            except ValueError as error:
                raise WorkflowError(
                    f"{row.get('sample_id', '?')}: {field} is not a number: {raw!r}."
                ) from error
        output[field] = {
            "available": len(values),
            "median": round(statistics.median(values), 2) if values else None,
        }
    return output


def render_html(rows, summary, title):
    def escape(value):
        return html.escape(str(value))

    trait_rows = []
    for trait, data in summary["traits"].items():
        width = data["positive"] / summary["samples"] * 100
        trait_rows.append(
            f"<tr><th>{escape(data['label'])}</th><td>{data['positive']}</td>"
            f"<td>{data['negative']}</td><td><div class='bar'><span style='width:{width:.2f}%'></span></div></td>"
            f"<td>{data['mean_label_probability_percent']:.2f}%</td></tr>"
        )
    comparison_rows = []
    for trait, data in summary["reference_comparisons"].items():
        agreement = (
            f"{data['agreement'] * 100:.1f}%" if data["agreement"] is not None else "Unavailable"
        )
        comparison_rows.append(
            f"<tr><th>{escape(TRAITS[trait])}</th><td>{data['eligible_samples']}</td>"
            f"<td>{agreement}</td><td>{data['true_positive']}</td><td>{data['true_negative']}</td>"
            f"<td>{data['false_positive']}</td><td>{data['false_negative']}</td></tr>"
        )
    selected = list(summary["traits"])
    headers = [
        "Sample",
        "Taxon",
        "Completeness %",
        "Contamination %",
        *[TRAITS[t] for t in selected],
    ]
    sample_rows = []
    for row in rows:
        cells = [
            escape(row.get(field) or "—")
            for field in ("sample_id", "taxon", "completeness", "contamination")
        ]
        for trait in selected:
            positive = as_bool(row[f"{trait}_prediction"])
            label = "Yes" if positive else "No"
            cells.append(
                f"<span class='{'yes' if positive else 'no'}'>{label}</span> "
                f"<small>{float(row[f'{trait}_confidence']):.2f}%</small>"
            )
        sample_rows.append("<tr>" + "".join(f"<td>{cell}</td>" for cell in cells) + "</tr>")
    return f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>{escape(title)}</title><style>
:root{{font-family:system-ui,sans-serif;color:#172b40;background:#f4f7fa}}body{{max-width:1200px;margin:auto;padding:32px}}
h1{{font-size:2.3rem}}h2{{margin-top:32px}}.cards{{display:flex;gap:16px;flex-wrap:wrap}}.card{{background:white;padding:20px;border-radius:12px;flex:1}}
.card strong{{display:block;font-size:2rem;color:#176b79}}table{{border-collapse:collapse;width:100%;background:white}}td,th{{text-align:left;padding:12px;border-bottom:1px solid #e1e7ed}}
th{{font-size:.9rem}}.scroll{{overflow:auto}}.bar{{height:12px;background:#e3e9ef;width:160px;border-radius:6px}}.bar span{{display:block;height:100%;background:#176b79;border-radius:6px}}
.note{{background:#fff2d8;padding:18px;border-radius:10px;line-height:1.6}}input{{padding:12px;width:min(90%,480px);margin:12px 0;border:1px solid #abbac5;border-radius:6px}}
.yes{{color:#0b6b48;font-weight:600}}.no{{color:#5d6670}}small{{white-space:nowrap}}a{{color:#176b79}}footer{{margin-top:32px;font-size:.9rem;line-height:1.6}}
</style></head><body>
<p>BacDive-AI workflow · reproducible results</p><h1>{escape(title)}</h1>
<div class="cards"><div class="card"><strong>{summary["samples"]}</strong>genomes</div>
<div class="card"><strong>{len(selected)}</strong>traits per genome</div>
<div class="card"><strong>{summary["metadata_matches"]}</strong>reference metadata matches</div></div>
<p class="note">These are predictions from published BacDive-AI models. The percentages show the probability assigned to the predicted label;
they are not measured accuracy. Reference comparisons describe this dataset. Model training overlap and strain-level identity have not been established.</p>
<h2>Trait distribution</h2><div class="scroll"><table><thead><tr><th>Trait</th><th>Positive</th><th>Negative</th><th>Positive share</th><th>Mean label probability</th></tr></thead>
<tbody>{"".join(trait_rows)}</tbody></table></div>
<h2>Reference agreement</h2><p>Gram comparisons use only explicit positive/negative labels. Oxygen comparisons use only explicit aerobe/anaerobe labels.
Missing, variable, combined, facultative, and ambiguous labels are excluded. TP/TN/FP/FN are relative to the named trait.</p>
<div class="scroll"><table><thead><tr><th>Trait</th><th>Eligible</th><th>Agreement</th><th>TP</th><th>TN</th><th>FP</th><th>FN</th></tr></thead><tbody>{"".join(comparison_rows)}</tbody></table></div>
<h2>Explore samples</h2><label for="search">Filter by sample, taxon, or result</label><br><input id="search" type="search" placeholder="Search the table">
<p id="visible" aria-live="polite">{len(rows)} samples shown</p><div class="scroll"><table id="samples"><thead><tr>{"".join(f"<th>{escape(h)}</th>" for h in headers)}</tr></thead><tbody>{"".join(sample_rows)}</tbody></table></div>
<footer>Feature encoding: binary Pfam presence/absence. Model origin: <a href="https://doi.org/10.5281/zenodo.15075932">BacDive-AI v2</a>.
Methods: <a href="https://doi.org/10.1038/s42003-025-08313-3">Koblitz et al., Communications Biology (2025)</a>.<br>
Source CSV SHA-256: <code>{escape(summary.get("source_sha256", ""))}</code>. This report contains no external scripts or tracking.</footer>
<script>document.getElementById('search').addEventListener('input',function(){{const q=this.value.toLowerCase();let n=0;document.querySelectorAll('#samples tbody tr').forEach(r=>{{r.hidden=!r.textContent.toLowerCase().includes(q);if(!r.hidden)n++;}});document.getElementById('visible').textContent=n+' samples shown';}});</script>
</body></html>"""


def generate_report(csv_path, html_path, summary_path, title="Bacterial phenotype predictions"):
    with Path(csv_path).open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    summary = summarize(rows)
    summary["source_sha256"] = sha256_file(csv_path)
    write_json(summary_path, summary)
    atomic_text(html_path, render_html(rows, summary, title))
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description="Build a standalone HTML prediction report.")
    parser.add_argument("csv", type=Path)
    parser.add_argument("--html", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--title", default="Bacterial phenotype predictions")
    args = parser.parse_args(argv)
    try:
        generate_report(args.csv, args.html, args.summary, args.title)
        return 0
    except (WorkflowError, OSError, KeyError, ValueError) as error:
        parser.exit(1, f"Error: {error}\n")


if __name__ == "__main__":
    raise SystemExit(main())
