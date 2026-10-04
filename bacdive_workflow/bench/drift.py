"""Annotation drift: our Pfam sets versus the published training features."""


def jaccard(first, second):
    if not first and not second:
        return None
    union = set(first) | set(second)
    if not union:
        return None
    return len(set(first) & set(second)) / len(union)


def drift_row(genome, ips_version, ours, published, prediction_changes, accession_match):
    gained = sorted(set(ours) - set(published))
    lost = sorted(set(published) - set(ours))
    return {
        "genome": genome,
        "ips_version": ips_version,
        "jaccard": jaccard(ours, published),
        "n_ours": len(set(ours)),
        "n_published": len(set(published)),
        "n_gained": len(gained),
        "n_lost": len(lost),
        "prediction_changes": prediction_changes,
        "accession_match": accession_match,
    }
