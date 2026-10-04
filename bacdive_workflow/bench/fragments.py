"""Completeness degradation by seeded fragment loss.

Fragments are log-normal lengths cut per contig; a gene survives only when one
kept fragment covers it entirely. Contamination adds donor fragments from a
different phylum. Level 100 with zero contamination keeps every gene, so it
must reproduce the full-genome prediction exactly.
"""

import math

from bacdive_workflow.bench import derive_seed


def cut_fragments(length, median_bp, sigma, min_bp, rng):
    """Cover a contig of ``length`` bases with log-normal fragment lengths."""
    fragments = []
    position = 0
    while position < length:
        draw = float(rng.lognormal(math.log(median_bp), sigma))
        size = max(int(draw), min_bp)
        size = min(size, length - position)
        fragments.append((position, position + size))
        position += size
    return fragments


def keep_genes(genes, kept):
    """Return genes lying entirely inside one kept (contig, start, end) fragment."""
    survivors = []
    for gene in genes:
        for fragment_contig, start, end in kept:
            if gene["contig"] == fragment_contig and gene["start"] >= start and gene["end"] <= end:
                survivors.append(gene)
                break
    return survivors


def kept_fragments(contig_lengths, level, contamination, params, seed):
    """Fragment coordinates kept for a draw; level 100 keeps everything.

    Level 100 with zero contamination applies no fragmentation: each contig
    is one fragment, so every gene survives and the draw must reproduce the
    full-genome prediction exactly.
    """
    if level >= 100 and not contamination:
        return [(contig, 0, length) for contig, length in contig_lengths.items()]
    import numpy as np

    rng = np.random.default_rng(seed)
    fragments = []
    for contig, length in contig_lengths.items():
        for start, end in cut_fragments(
            length,
            params["fragment_median_bp"],
            params["fragment_sigma"],
            params["min_fragment_bp"],
            rng,
        ):
            fragments.append((contig, start, end))
    if level >= 100 and not contamination:
        return fragments
    order = rng.permutation(len(fragments))
    target = contig_total(contig_lengths) * level / 100
    kept = []
    retained = 0
    for index in order:
        if retained >= target:
            break
        kept.append(fragments[index])
        retained += fragments[index][2] - fragments[index][1]
    return kept


def simulate(genome_genes, contig_lengths, level, contamination, donor_genes, params, seed):
    """Simulate one (level, contamination, replicate) draw.

    Returns ``(kept_gene_ids, retained_bp, foreign_bp, foreign_genes,
    kept_fragments)``. Level 100 with no contamination keeps every gene by
    construction.
    """
    import numpy as np

    rng = np.random.default_rng(seed)
    kept = kept_fragments(contig_lengths, level, contamination, params, seed)
    kept_genes = keep_genes(genome_genes, kept)
    retained_bp = sum(end - start for _, start, end in kept)
    foreign_bp = 0
    foreign_genes = []
    if contamination and donor_genes:
        target_foreign = retained_bp * contamination / 100
        permutation = rng.permutation(len(donor_genes))
        for index in permutation:
            if foreign_bp >= target_foreign:
                break
            foreign_genes.append(donor_genes[index])
            foreign_bp += donor_genes[index].get("length", params["fragment_median_bp"])
    return [g["protein_id"] for g in kept_genes], retained_bp, foreign_bp, foreign_genes, kept


def contig_total(contig_lengths):
    return sum(contig_lengths.values())


def draw_seed(config_seed, genome_id, level, contamination, replicate):
    return derive_seed(config_seed, genome_id, level, contamination, replicate)
