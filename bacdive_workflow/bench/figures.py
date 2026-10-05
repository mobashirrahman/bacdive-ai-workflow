"""Matplotlib figures for the benchmark report (PNG and SVG)."""

PALETTE = ["#0173b2", "#de8f05", "#029e73", "#d55e00", "#cc78bc", "#949494"]


def _setup():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 150,
            "axes.spines.top": False,
            "axes.spines.right": False,
            # Fixed salt: SVG element ids are otherwise random per run.
            "svg.hashsalt": "bacdive-benchmark",
        }
    )
    return plt


def _save(fig, path_png, path_svg):
    """Write both formats without the creation date, so reruns are byte-identical."""
    fig.savefig(path_png, metadata={"Software": None})
    fig.savefig(path_svg, metadata={"Date": None})


def accuracy_figure(rows, path_png, path_svg):
    plt = _setup()
    traits = sorted({r["trait"] for r in rows})
    groups = sorted({r["group"] for r in rows})
    width = max(6, len(traits) * len(groups) * 0.6)
    fig, ax = plt.subplots(figsize=(width, 4))
    positions, values, errors, colors, labels = [], [], [], [], []
    for gi, group in enumerate(groups):
        for ti, trait in enumerate(traits):
            match = [r for r in rows if r["trait"] == trait and r["group"] == group]
            if not match or match[0]["balanced_accuracy"] is None:
                continue
            row = match[0]
            positions.append(ti * (len(groups) + 1) + gi)
            values.append(row["balanced_accuracy"])
            low = row.get("balanced_accuracy_low") or row["balanced_accuracy"]
            high = row.get("balanced_accuracy_high") or row["balanced_accuracy"]
            errors.append([row["balanced_accuracy"] - low, high - row["balanced_accuracy"]])
            colors.append(PALETTE[gi % len(PALETTE)])
            labels.append(f"{trait} (n={row['n']})")
    if values:
        import numpy as np

        ax.errorbar(positions, values, yerr=np.array(errors).T, fmt="o", ecolor="black")
        for x, y, color in zip(positions, values, colors):
            ax.scatter([x], [y], color=color, s=60, zorder=3)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Balanced accuracy (seen = training-set upper bound)")
    ax.set_title("Accuracy by trait and set")
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)


def flip_figure(rows, path_png, path_svg):
    """Share of full-genome positives lost per completeness level, one line per trait.

    `rows` are one set's uncontaminated summary rows. Traits with too few
    full-genome positives have no `positive_loss_rate` and are not drawn.
    """
    plt = _setup()
    fig, ax = plt.subplots(figsize=(7, 4))
    drawn = [r for r in rows if r["positive_loss_rate"] is not None]
    for ti, trait in enumerate(sorted({r["trait"] for r in drawn})):
        series = sorted((r for r in drawn if r["trait"] == trait), key=lambda r: r["level"])
        ax.plot(
            [r["level"] for r in series],
            [r["positive_loss_rate"] for r in series],
            marker="o",
            color=PALETTE[ti % len(PALETTE)],
            label=f"{trait} (n={series[0]['full_positive_genomes']})",
        )
    ax.invert_xaxis()
    ax.set_ylim(0, 1)
    ax.set_xlabel("Completeness level (% bases retained)")
    ax.set_ylabel("Share of full-genome positives called negative")
    ax.set_title("Positive calls lost with completeness (unseen set, n = positive genomes)")
    ax.legend(fontsize="small")
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)


def drift_figure(rows, path_png, path_svg):
    plt = _setup()
    fig, ax = plt.subplots(figsize=(6, 4))
    values = [r["jaccard"] for r in rows if r["jaccard"] is not None]
    ax.hist(values, bins=20, color=PALETTE[0])
    ax.set_xlabel("Jaccard similarity, ours vs published Pfam set (same E-value threshold)")
    ax.set_ylabel(f"Genomes (n={len(values)})")
    ax.set_title("Annotation drift")
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)


def robustness_heatmap(rows, path_png, path_svg):
    plt = _setup()
    traits = sorted({r["trait"] for r in rows})
    variants = sorted({r["variant"] for r in rows})
    import numpy as np

    matrix = np.full((len(traits), len(variants)), float("nan"))
    for r in rows:
        if r.get("agreement") is not None:
            matrix[traits.index(r["trait"]), variants.index(r["variant"])] = r["agreement"]
    fig, ax = plt.subplots(figsize=(max(6, len(variants) * 1.2), max(3, len(traits) * 0.6)))
    im = ax.imshow(matrix, vmin=0, vmax=1, cmap="Blues")
    ax.set_xticks(range(len(variants)), variants, rotation=30, ha="right")
    ax.set_yticks(range(len(traits)), traits)
    ax.set_title("Robustness: class agreement with primary configuration (n in table)")
    fig.colorbar(im, ax=ax, label="Agreement")
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)
