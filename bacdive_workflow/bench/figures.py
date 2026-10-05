"""Matplotlib figures for the benchmark report (PNG and SVG).

Colours are a fixed, colour-blind-checked categorical order; a trait keeps its
colour in every figure. Text stays in ink colours, never in a series colour.
"""

from bacdive_workflow.common import TRAITS

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK_SECONDARY = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
SEQUENTIAL = ["#fcfcfb", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"]

# Fixed colour per trait, in this order, wherever traits are series.
TRAIT_ORDER = ["spore-forming", "motile2+", "thermophile", "aerobic", "gram-positive", "anaerobic"]
TRAIT_COLOR = dict(zip(TRAIT_ORDER, SERIES))

VARIANTS = [
    ("caller-meta", "Prodigal\nmeta mode"),
    ("caller-pgap", "PGAP\nproteins"),
    ("pfam-5.63", "Pfam 35.0\n(InterProScan 5.63)"),
    ("evalue-1e-10", "E-value\n1e-10"),
    ("evalue-1e-30", "E-value\n1e-30"),
]


def _setup():
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "figure.dpi": 200,
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.spines.left": False,
            "axes.edgecolor": AXIS,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlecolor": INK,
            "axes.titlesize": 11,
            "axes.titleweight": "bold",
            "axes.titlelocation": "left",
            "axes.labelsize": 9,
            "xtick.color": MUTED,
            "ytick.color": MUTED,
            "xtick.labelcolor": INK_SECONDARY,
            "ytick.labelcolor": INK_SECONDARY,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "ytick.left": False,
            "grid.color": GRID,
            "grid.linewidth": 0.8,
            "legend.frameon": False,
            "legend.fontsize": 9,
            "text.color": INK,
            # Fixed salt: SVG element ids are otherwise random per run.
            "svg.hashsalt": "bacdive-benchmark",
        }
    )
    return plt


def _save(fig, path_png, path_svg):
    """Write both formats without the creation date, so reruns are byte-identical."""
    fig.savefig(path_png, metadata={"Software": None})
    fig.savefig(path_svg, metadata={"Date": None})


def _subtitle(ax, text):
    ax.text(0, 1.02, text, transform=ax.transAxes, fontsize=9, color=INK_SECONDARY, va="bottom")


def accuracy_figure(rows, path_png, path_svg):
    """Balanced accuracy with bootstrap interval per trait: unseen against seen."""
    plt = _setup()
    groups = [
        ("unseen", "Unseen species (held out)", SERIES[0]),
        ("seen", "Training genomes (upper bound)", SERIES[1]),
    ]
    rows = [r for r in rows if not r.get("phylum") and r["balanced_accuracy"] is not None]
    traits = [t for t in TRAIT_ORDER if any(r["trait"] == t for r in rows)]
    fig, ax = plt.subplots(figsize=(7.5, 0.62 * max(len(traits), 1) + 1.6))
    for gi, (group, label, color) in enumerate(groups):
        offset = 0.16 if gi == 0 else -0.16
        drawn = False
        for ti, trait in enumerate(traits):
            match = [r for r in rows if r["trait"] == trait and r["group"] == group]
            if not match:
                continue
            row = match[0]
            y = len(traits) - 1 - ti + offset
            low = row.get("balanced_accuracy_low") or row["balanced_accuracy"]
            high = row.get("balanced_accuracy_high") or row["balanced_accuracy"]
            ax.plot([low, high], [y, y], color=color, linewidth=2, solid_capstyle="round")
            ax.plot(
                [row["balanced_accuracy"]],
                [y],
                marker="o",
                markersize=7,
                color=color,
                markeredgecolor=SURFACE,
                markeredgewidth=1.5,
                linestyle="none",
                label=None if drawn else label,
            )
            drawn = True
            ax.text(1.004, y, f"n={row['n']}", fontsize=8, color=MUTED, va="center")
    ax.set_yticks(
        [len(traits) - 1 - i for i in range(len(traits))], [TRAITS.get(t, t) for t in traits]
    )
    ax.set_ylim(-0.6, len(traits) - 0.4)
    ax.set_xlim(0.8, 1.0)
    ax.set_xticks([0.80, 0.85, 0.90, 0.95, 1.00])
    ax.set_xlabel("Balanced accuracy, 0.80 to 1.00 (line: 95% bootstrap interval)")
    ax.grid(axis="x")
    ax.set_axisbelow(True)
    ax.set_title("Accuracy holds on species the models never saw", pad=22)
    _subtitle(ax, "Groups with fewer than 10 genomes in a class are not scored")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=2)
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)


def flip_figure(rows, path_png, path_svg):
    """Share of full-genome positives lost per completeness level, one line per trait.

    `rows` are one set's uncontaminated summary rows. Traits with too few
    full-genome positives have no `positive_loss_rate` and are not drawn.
    """
    plt = _setup()
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    drawn = [r for r in rows if r["positive_loss_rate"] is not None]
    present = {r["trait"] for r in drawn}
    traits = [t for t in TRAIT_ORDER if t in present] + sorted(present - set(TRAIT_ORDER))
    ends = []
    for ti, trait in enumerate(traits):
        series = sorted((r for r in drawn if r["trait"] == trait), key=lambda r: -r["level"])
        color = TRAIT_COLOR.get(trait, SERIES[ti % len(SERIES)])
        ax.plot(
            [r["level"] for r in series],
            [100 * r["positive_loss_rate"] for r in series],
            color=color,
            linewidth=2,
            marker="o",
            markersize=5,
            markeredgecolor=SURFACE,
            markeredgewidth=1,
        )
        label = f"{TRAITS.get(trait, trait)} (n={series[0]['full_positive_genomes']})"
        ends.append([100 * series[-1]["positive_loss_rate"], label, color, series[-1]["level"]])
    # Direct labels at the line ends, pushed apart so they never overlap.
    ends.sort(key=lambda end: end[0])
    position = None
    for end in ends:
        position = end[0] if position is None else max(end[0], position + 6.5)
        end.append(position)
    overshoot = max((end[4] for end in ends), default=100) - 100
    for value, label, color, level, y in ends:
        y -= max(overshoot, 0)
        ax.plot([level - 0.6, level - 1.8], [value, y], color=color, linewidth=0.8, clip_on=False)
        ax.text(level - 2.2, y, label, fontsize=8.5, color=INK, va="center", clip_on=False)
    levels = sorted({r["level"] for r in drawn}, reverse=True)
    if levels:
        ax.set_xlim(levels[0] + 2, levels[-1] - 2)
        ax.set_xticks(levels, [f"{level}%" for level in levels])
    ax.set_ylim(0, 102)
    ax.set_yticks(range(0, 101, 25), [f"{v}%" for v in range(0, 101, 25)])
    ax.grid(axis="y")
    ax.set_axisbelow(True)
    ax.set_xlabel("Genome completeness (share of bases retained)")
    ax.set_title("Incomplete genomes lose positive calls", pad=22)
    _subtitle(ax, "Share of full-genome positive calls that turn negative (unseen set)")
    fig.subplots_adjust(left=0.08, right=0.72, top=0.86, bottom=0.13)
    _save(fig, path_png, path_svg)
    plt.close(fig)


def drift_figure(rows, path_png, path_svg):
    """Jaccard similarity to the published Pfam sets, one panel per InterProScan version."""
    plt = _setup()
    versions = sorted({r["ips_version"] for r in rows})
    fig, axes = plt.subplots(
        len(versions), 1, figsize=(7.5, 1.9 * max(len(versions), 1) + 1.1), sharex=True
    )
    axes = [axes] if len(versions) == 1 else list(axes)
    bins = [0.70 + 0.01 * i for i in range(31)]
    for ax, version in zip(axes, versions):
        group = [r for r in rows if r["ips_version"] == version and r["jaccard"] is not None]
        values = [min(max(r["jaccard"], 0.70), 1.0) for r in group]
        changed = sum(1 for r in group if int(r.get("prediction_changes") or 0) > 0)
        ax.hist(values, bins=bins, color=SERIES[0], edgecolor=SURFACE, linewidth=1)
        ax.set_title(f"InterProScan {version}", fontsize=10)
        ax.text(
            0.0,
            0.82,
            f"{changed} of {len(group)} genomes change a prediction",
            transform=ax.transAxes,
            fontsize=9,
            color=INK_SECONDARY,
        )
        ax.set_ylabel("Genomes")
        ax.grid(axis="y")
        ax.set_axisbelow(True)
    axes[-1].set_xlim(0.70, 1.0)
    axes[-1].set_xlabel(
        "Jaccard similarity to the published training Pfam set (both at E-value 1e-20)"
    )
    fig.suptitle(
        "The training-era annotation reproduces the published features",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)


def robustness_heatmap(rows, path_png, path_svg):
    """Share of calls that change against the primary configuration, per trait and variant."""
    plt = _setup()
    import numpy as np
    from matplotlib.colors import LinearSegmentedColormap

    present = {r["variant"] for r in rows}
    variants = [(key, label) for key, label in VARIANTS if key in present]
    known = {key for key, _ in VARIANTS} | {"primary"}
    variants += [(key, key) for key in sorted(present - known)]
    traits = [t for t in TRAITS if any(r["trait"] == t for r in rows)]
    matrix = np.full((len(traits), len(variants)), float("nan"))
    keys = [key for key, _ in variants]
    for r in rows:
        if r.get("agreement") is not None and r["variant"] in keys and r["trait"] in traits:
            matrix[traits.index(r["trait"]), keys.index(r["variant"])] = 100 * (1 - r["agreement"])
    fig, ax = plt.subplots(
        figsize=(max(6.5, 1.35 * len(variants) + 2.2), 0.5 * max(len(traits), 1) + 2.0)
    )
    top = 20
    cmap = LinearSegmentedColormap.from_list("blue", SEQUENTIAL)
    im = ax.imshow(matrix, vmin=0, vmax=top, cmap=cmap, aspect="auto")
    for i in range(len(traits)):
        for j in range(len(variants)):
            if not np.isnan(matrix[i, j]):
                ax.text(
                    j,
                    i,
                    f"{matrix[i, j]:.0f}%",
                    ha="center",
                    va="center",
                    fontsize=9,
                    color=SURFACE if matrix[i, j] > 0.55 * top else INK,
                )
    ax.set_xticks(range(len(variants)), [label for _, label in variants])
    ax.set_yticks(range(len(traits)), [TRAITS.get(t, t) for t in traits])
    ax.tick_params(length=0)
    ax.xaxis.tick_top()
    for spine in ax.spines.values():
        spine.set_visible(False)
    # 2px surface gaps between cells.
    ax.set_xticks([x - 0.5 for x in range(1, len(variants))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(traits))], minor=True)
    ax.grid(which="minor", color=SURFACE, linewidth=2)
    ax.tick_params(which="minor", length=0)
    bar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.03)
    bar.outline.set_visible(False)
    bar.set_ticks([0, 10, 20], labels=["0%", "10%", "20%"])
    bar.set_label("Calls that change")
    fig.suptitle(
        "Which annotation choices change predictions",
        x=0.02,
        ha="left",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    fig.tight_layout()
    _save(fig, path_png, path_svg)
    plt.close(fig)
