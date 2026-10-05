"""Annotation drift: our Pfam sets versus the published training features."""


def jaccard(first, second):
    if not first and not second:
        return None
    union = set(first) | set(second)
    if not union:
        return None
    return len(set(first) & set(second)) / len(union)


def threshold_set(evalues, threshold):
    """Pfams of a published {pfam: E-value} map that pass the E-value threshold."""
    return {pfam for pfam, evalue in evalues.items() if evalue <= threshold}


def drift_row(
    genome,
    ips_version,
    ours,
    published,
    published_unfiltered,
    prediction_changes,
    prediction_changes_unfiltered,
    accession_match,
):
    """One drift record.

    `ours` and `published` carry the same E-value threshold, so `jaccard`,
    the gained/lost counts and `prediction_changes` measure annotation drift
    alone. The `_unfiltered` fields compare against the published set as
    shipped and therefore also contain the threshold effect.
    """
    ours, published = set(ours), set(published)
    return {
        "genome": genome,
        "ips_version": ips_version,
        "jaccard": jaccard(ours, published),
        "jaccard_unfiltered": jaccard(ours, published_unfiltered),
        "n_ours": len(ours),
        "n_published": len(published),
        "n_published_unfiltered": len(set(published_unfiltered)),
        "n_gained": len(ours - published),
        "n_lost": len(published - ours),
        "prediction_changes": prediction_changes,
        "prediction_changes_unfiltered": prediction_changes_unfiltered,
        "accession_match": accession_match,
    }
