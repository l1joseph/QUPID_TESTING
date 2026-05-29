#!/usr/bin/env python3
"""
Publication-quality 4-panel benchmark figure for Qupid.

Reads:  benchmark_agp/benchmark_results.tsv  (real barnacle2 data)
Writes: figures/qupid_benchmark_runtime.{pdf,png}

Panels:
  a  Runtime vs. background pool size (log-log)
  b  Runtime vs. number of matchings  (log-log)
  c  Matching success rate vs. background pool size
  d  Qupid speedup factor over naive rebuild
"""
import string
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Paul Tol colorblind-safe qualitative palette (max 9 categories)
PALETTE = [
    "#332288",
    "#88CCEE",
    "#44AA99",
    "#117733",
    "#999933",
    "#DDCC77",
    "#CC6677",
    "#882255",
    "#AA4499",
]

NATURE_WIDTHS = {"single": 3.50, "1.5col": 5.04, "double": 7.09}
NATURE_MAX_H = 6.69

NATURE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 7,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 5,
    "legend.title_fontsize": 6,
    "figure.titlesize": 8,
    "axes.linewidth": 0.8,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "text.color": "black",
    "xtick.color": "black",
    "ytick.color": "black",
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.minor.width": 0.3,
    "ytick.minor.width": 0.3,
    "lines.linewidth": 0.75,
    "patch.linewidth": 0.5,
    "grid.linewidth": 0.3,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
    "xtick.minor.size": 1.5,
    "ytick.minor.size": 1.5,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.pad": 2,
    "ytick.major.pad": 2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "axes.axisbelow": True,
    "legend.frameon": False,
    "legend.borderpad": 0.3,
    "legend.handlelength": 1.2,
    "legend.handletextpad": 0.4,
    "legend.labelspacing": 0.3,
    "scatter.edgecolors": "white",
    "errorbar.capsize": 2,
    "figure.dpi": 150,
    "savefig.dpi": 450,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.constrained_layout.use": False,
    "image.interpolation": "nearest",
    "image.cmap": "viridis",
    "mathtext.fontset": "dejavusans",
}

NATURE_PANEL = {"fontsize": 8, "fontweight": "bold"}
NATURE_FORMATS = ["pdf", "png"]
NATURE_DPI = 450

METHOD_COLORS = {
    "qupid": "#332288",
    "graph_once": "#44AA99",
    "naive_loop": "#CC6677",
}
METHOD_LABELS = {
    "qupid": "Qupid (cached graph + Hopcroft–Karp)",
    "graph_once": "Cached graph + greedy",
    "naive_loop": "Rebuild graph every iteration",
}
METHOD_MARKERS = {
    "qupid": "o",
    "graph_once": "s",
    "naive_loop": "^",
}
METHOD_ORDER = ["naive_loop", "graph_once", "qupid"]


def _add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.06,
        label,
        fontsize=NATURE_PANEL["fontsize"],
        fontweight=NATURE_PANEL["fontweight"],
        va="top",
        ha="right",
        transform=ax.transAxes,
    )


def main() -> None:
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)

    df = pd.read_csv("benchmark_agp/benchmark_results.tsv", sep="\t")

    bg = df[df["experiment"] == "background_size"].sort_values("n_controls")
    it = df[df["experiment"] == "iterations"].sort_values("n_iterations")

    # Metadata for titles
    k_bg = int(bg[bg["method"] == "qupid"]["n_iterations"].iloc[0])
    n_cases_bg = int(bg[bg["method"] == "qupid"]["n_cases"].iloc[0])
    n_cases_it = int(it[it["method"] == "qupid"]["n_cases"].iloc[0])
    n_bg_it = int(it[it["method"] == "qupid"]["n_controls"].iloc[0])

    w = NATURE_WIDTHS["double"]
    h = min(w * 1.05, NATURE_MAX_H)  # extra height for bottom legend
    fig, axes = plt.subplots(2, 2, figsize=(w, h))
    fig.subplots_adjust(
        left=0.12, right=0.97, top=0.93, bottom=0.16, hspace=0.52, wspace=0.38
    )

    # ── Panel a: Runtime vs background size ──────────────────────────────
    ax = axes[0, 0]
    for method in METHOD_ORDER:
        sub = bg[bg["method"] == method]
        if sub.empty:
            continue
        ax.plot(
            sub["n_controls"],
            sub["elapsed_sec"],
            marker=METHOD_MARKERS[method],
            color=METHOD_COLORS[method],
            label=METHOD_LABELS[method],
            linewidth=1.0,
            markersize=3.5,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Background pool size")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(f"Runtime vs. background size\n({n_cases_bg} cases, k = {k_bg})")
    ax.grid(True, which="both", alpha=0.2, linewidth=0.3)

    # ── Panel b: Runtime vs iterations ───────────────────────────────────
    ax = axes[0, 1]
    for method in METHOD_ORDER:
        sub = it[it["method"] == method]
        if sub.empty:
            continue
        ax.plot(
            sub["n_iterations"],
            sub["elapsed_sec"],
            marker=METHOD_MARKERS[method],
            color=METHOD_COLORS[method],
            label=METHOD_LABELS[method],
            linewidth=1.0,
            markersize=3.5,
        )

    # Annotate speedup at largest shared k
    shared_k = sorted(
        set(it[it["method"] == "qupid"]["n_iterations"])
        & set(it[it["method"] == "naive_loop"]["n_iterations"])
    )
    if shared_k:
        k_ann = shared_k[-1]
        q_t = it[(it["method"] == "qupid") & (it["n_iterations"] == k_ann)][
            "elapsed_sec"
        ].values[0]
        n_t = it[(it["method"] == "naive_loop") & (it["n_iterations"] == k_ann)][
            "elapsed_sec"
        ].values[0]
        speedup = n_t / q_t
        ax.annotate(
            f"{speedup:.0f}× speedup\nat k = {k_ann}",
            xy=(k_ann, q_t),
            xytext=(k_ann * 0.22, q_t * 8),
            fontsize=5.5,
            color="#333333",
            arrowprops=dict(arrowstyle="->", color="gray", lw=0.6),
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(
        f"Runtime vs. number of matchings\n({n_cases_it} cases, {n_bg_it} controls)"
    )
    ax.grid(True, which="both", alpha=0.2, linewidth=0.3)

    # ── Panel c: Matching success rate ───────────────────────────────────
    ax = axes[1, 0]
    for method in METHOD_ORDER:
        sub = bg[bg["method"] == method].copy()
        if sub.empty:
            continue
        sub["success_rate"] = sub["successful_iterations"] / sub["n_iterations"]
        ax.plot(
            sub["n_controls"],
            sub["success_rate"],
            marker=METHOD_MARKERS[method],
            color=METHOD_COLORS[method],
            label=METHOD_LABELS[method],
            linewidth=1.0,
            markersize=3.5,
        )
    ax.set_xscale("log")
    ax.set_xlabel("Background pool size")
    ax.set_ylabel("Fraction of matchings fully solved")
    ax.set_title(f"Matching success rate\n({n_cases_bg} cases, k = {k_bg})")
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, which="both", alpha=0.2, linewidth=0.3)

    # ── Panel d: Speedup factor ───────────────────────────────────────────
    ax = axes[1, 1]
    qupid_t = it[it["method"] == "qupid"].set_index("n_iterations")["elapsed_sec"]
    naive_t = it[it["method"] == "naive_loop"].set_index("n_iterations")["elapsed_sec"]
    shared_k = sorted(set(qupid_t.index) & set(naive_t.index))
    if shared_k:
        speedups = [naive_t[k] / qupid_t[k] for k in shared_k]
        ax.plot(
            shared_k,
            speedups,
            marker="o",
            color=METHOD_COLORS["qupid"],
            linewidth=1.0,
            markersize=3.5,
        )
        ax.fill_between(shared_k, 1, speedups, alpha=0.12, color=METHOD_COLORS["qupid"])
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=0.6, alpha=0.7)
        # label the reference line directly to avoid legend overlap
        ax.text(
            shared_k[0] * 1.1,
            1.0,
            "no speedup",
            va="bottom",
            ha="left",
            fontsize=5,
            color="gray",
        )
    ax.set_xscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Speedup over naive (×)")
    ax.set_title(f"Qupid speedup over naive\n({n_cases_it} cases, {n_bg_it} controls)")
    ax.grid(True, which="both", alpha=0.2, linewidth=0.3)

    # Panel labels
    for i, ax in enumerate(axes.flatten()):
        _add_panel_label(ax, string.ascii_lowercase[i])

    # Shared legend for the three methods, placed below all panels
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=3,
        bbox_to_anchor=(0.5, 0.01),
        fontsize=5.5,
        frameon=False,
        handlelength=1.4,
        handletextpad=0.4,
        columnspacing=1.0,
    )

    # Save
    out = pathlib.Path("figures")
    out.mkdir(exist_ok=True)
    OUTNAME = "qupid_benchmark_runtime"
    for fmt in NATURE_FORMATS:
        fig.savefig(
            out / f"{OUTNAME}.{fmt}",
            dpi=NATURE_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
            format=fmt,
            facecolor="white",
        )
        print(f"Saved: figures/{OUTNAME}.{fmt}")
    plt.close(fig)


if __name__ == "__main__":
    main()
