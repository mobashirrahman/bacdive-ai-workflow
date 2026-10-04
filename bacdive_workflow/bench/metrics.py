"""Leakage-controlled accuracy with stratified bootstrap intervals.

Seen-set numbers are training-set performance (an upper bound), never
accuracy; callers must keep that label on every table and figure.
"""

import csv
import math

from bacdive_workflow.bench import derive_seed

MIN_CLASS_MEMBERS = 10


def _mcc(tp, tn, fp, fn):
    denominator = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    if denominator == 0:
        return None
    return (tp * tn - fp * fn) / denominator


def _auroc(labels, scores):
    pairs = [(score, label) for label, score in zip(labels, scores)]
    positives = sum(labels)
    negatives = len(labels) - positives
    if positives == 0 or negatives == 0:
        return None
    order = sorted(range(len(pairs)), key=lambda i: pairs[i][0])
    ranks = [0.0] * len(pairs)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and pairs[order[j + 1]][0] == pairs[order[i]][0]:
            j += 1
        average = (i + j) / 2 + 1
        for k in range(i, j + 1):
            ranks[order[k]] = average
        i = j + 1
    rank_sum = sum(rank for rank, (_, label) in zip(ranks, pairs) if label == 1)
    return (rank_sum - positives * (positives + 1) / 2) / (positives * negatives)


def confusion(labels, predictions):
    tp = sum(1 for truth, pred in zip(labels, predictions) if truth == 1 and pred == 1)
    tn = sum(1 for truth, pred in zip(labels, predictions) if truth == 0 and pred == 0)
    fp = sum(1 for truth, pred in zip(labels, predictions) if truth == 0 and pred == 1)
    fn = sum(1 for truth, pred in zip(labels, predictions) if truth == 1 and pred == 0)
    return tp, tn, fp, fn


def summarize(labels, predictions, scores):
    """Point estimates for one trait/group; empty metrics when a class is small."""
    tp, tn, fp, fn = confusion(labels, predictions)
    positives = tp + fn
    negatives = tn + fp
    row = {
        "n": len(labels),
        "positives": positives,
        "negatives": negatives,
        "sensitivity": None,
        "specificity": None,
        "balanced_accuracy": None,
        "mcc": None,
        "auroc": None,
    }
    if positives < MIN_CLASS_MEMBERS or negatives < MIN_CLASS_MEMBERS:
        return row
    sensitivity = tp / positives
    specificity = tn / negatives
    row.update(
        {
            "sensitivity": sensitivity,
            "specificity": specificity,
            "balanced_accuracy": (sensitivity + specificity) / 2,
            "mcc": _mcc(tp, tn, fp, fn),
            "auroc": _auroc(labels, scores),
        }
    )
    return row


def stratified_bootstrap_ci(
    labels, predictions, scores, resamples, seed, metric="balanced_accuracy"
):
    """Percentile 95% interval over stratified resamples; None when not computable."""
    import random

    positives = [i for i, truth in enumerate(labels) if truth == 1]
    negatives = [i for i, truth in enumerate(labels) if truth == 0]
    if len(positives) < MIN_CLASS_MEMBERS or len(negatives) < MIN_CLASS_MEMBERS:
        return None, None
    rng = random.Random(seed)
    estimates = []
    for _ in range(resamples):
        draw = [rng.choice(positives) for _ in positives] + [
            rng.choice(negatives) for _ in negatives
        ]
        sub = summarize(
            [labels[i] for i in draw],
            [predictions[i] for i in draw],
            [scores[i] for i in draw],
        )
        if sub[metric] is not None:
            estimates.append(sub[metric])
    if not estimates:
        return None, None
    estimates.sort()
    low = estimates[max(0, int(0.025 * len(estimates)))]
    high = estimates[min(len(estimates) - 1, int(0.975 * len(estimates)))]
    return low, high


def accuracy_rows(records, traits, groups, bootstrap, seed):
    """Yield one row per trait x group plus per-phylum rows with n >= 15."""
    rows = []
    for trait in traits:
        for group_name, predicate in groups.items():
            subset = [r for r in records if predicate(r)]
            labelled = [r for r in subset if r["labels"].get(trait) in (0, 1)]
            if not labelled:
                continue
            labels = [r["labels"][trait] for r in labelled]
            predictions = [r["predictions"][trait] for r in labelled]
            scores = [r["scores"][trait] for r in labelled]
            summary = summarize(labels, predictions, scores)
            low, high = stratified_bootstrap_ci(
                labels,
                predictions,
                scores,
                bootstrap,
                derive_seed(seed, trait, group_name),
            )
            rows.append(
                {
                    "trait": trait,
                    "group": group_name,
                    "phylum": "",
                    **summary,
                    "balanced_accuracy_low": low,
                    "balanced_accuracy_high": high,
                }
            )
        phyla = sorted({r["phylum"] for r in records if r["phylum"]})
        for phylum in phyla:
            subset = [r for r in records if r["phylum"] == phylum]
            labelled = [r for r in subset if r["labels"].get(trait) in (0, 1)]
            if len(labelled) < 15:
                continue
            labels = [r["labels"][trait] for r in labelled]
            predictions = [r["predictions"][trait] for r in labelled]
            scores = [r["scores"][trait] for r in labelled]
            summary = summarize(labels, predictions, scores)
            low, high = stratified_bootstrap_ci(
                labels,
                predictions,
                scores,
                bootstrap,
                derive_seed(seed, trait, "phylum", phylum),
            )
            rows.append(
                {
                    "trait": trait,
                    "group": "phylum",
                    "phylum": phylum,
                    **summary,
                    "balanced_accuracy_low": low,
                    "balanced_accuracy_high": high,
                }
            )
    return rows


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in columns})
