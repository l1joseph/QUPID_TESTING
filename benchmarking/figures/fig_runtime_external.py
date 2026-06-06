#!/usr/bin/env python3
"""
fig_runtime_external — qupid vs. external matching tools (2-panel, Nature double column).
Panels: a) qupid vs MatchIt on the AGP IBD cohort (real benchmark data),
        b) historical benchmarks (SPSS FUZZY, R Matching, qupid; Patel et al.).
Run from benchmarking/  (reads benchmark_real/agp_external_results.tsv).
"""
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ── Palette ──────────────────────────────────────────────────────────────────
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
PALETTE_TWO = {
    "blue_orange": ["#0077BB", "#EE7733"],
    "blue_red": ["#4477AA", "#CC3311"],
    "teal_wine": ["#44AA99", "#882255"],
    "blue_gold": ["#2E86AB", "#D4A03C"],
}

# ── Nature RC ────────────────────────────────────────────────────────────────
NATURE_WIDTHS = {"single": 3.50, "1.5col": 5.04, "double": 7.09}
NATURE_MAX_H = 6.69
NATURE_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.titlesize": 8,
    "axes.labelsize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 5,
    "legend.title_fontsize": 6,
    "figure.titlesize": 8,
    "axes.linewidth": 0.5,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
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
    "legend.handlelength": 1.5,
    "legend.handletextpad": 0.4,
    "legend.labelspacing": 0.3,
    "scatter.edgecolors": "white",
    "errorbar.capsize": 2,
    "figure.dpi": 150,
    "savefig.dpi": 450,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "figure.constrained_layout.use": True,
    "image.interpolation": "nearest",
    "image.cmap": "viridis",
    "mathtext.fontset": "dejavusans",
}
NATURE_PANEL = {"fontsize": 8, "fontweight": "bold"}
NATURE_FORMATS = ["pdf", "png"]
NATURE_DPI = 450

# ── Series styling (qupid blue kept consistent across both panels) ───────────
STYLE = {
    "qupid": {"color": "#0077BB", "ls": "-", "marker": "o"},
    "MatchIt": {"color": "#EE7733", "ls": "--", "marker": "s"},
    "R Matching": {"color": "#009988", "ls": "--", "marker": "^"},
    "SPSS FUZZY": {"color": "#CC3311", "ls": ":", "marker": "D"},
}

# ── Historical benchmarks (Patel et al.; Wisconsin 16S, sex + age ±4.5 yr) ───
# Legend label intentionally "qupid" (no "(Lucas)" suffix).
HISTORICAL = {
    "SPSS FUZZY": {10: 11.75, 100: 498.6, 1000: 39654.0},
    "R Matching": {10: 0.112, 100: 1.162, 1000: 11.84, 10000: 113.7},
    "qupid": {10: 0.015, 100: 0.069, 1000: 0.68, 10000: 6.48, 100000: 63.91},
}


def _plot_series(ax, k, t, label):
    s = STYLE[label]
    ax.plot(
        k,
        t,
        marker=s["marker"],
        markersize=3.5,
        linewidth=1.0,
        linestyle=s["ls"],
        color=s["color"],
        label=label,
        clip_on=False,
    )


def _add_panel_label(ax, label):
    ax.text(
        -0.16,
        1.05,
        label,
        fontsize=NATURE_PANEL["fontsize"],
        fontweight=NATURE_PANEL["fontweight"],
        va="top",
        ha="right",
        transform=ax.transAxes,
    )


def main():
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)

    w = NATURE_WIDTHS["double"]
    h = min(w * 0.618 * (1 / 2) ** 0.5, NATURE_MAX_H)
    fig, axes = plt.subplots(1, 2, figsize=(w, h))

    # ── Panel a: real AGP benchmark (qupid vs MatchIt) ───────────────────────
    ax = axes[0]
    df = pd.read_csv("benchmark_real/agp_external_results.tsv", sep="\t")
    for tool in ["qupid", "MatchIt"]:
        sub = df[df["tool"] == tool].sort_values("k")
        if not sub.empty:
            _plot_series(ax, sub["k"].values, sub["elapsed_sec"].values, tool)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings ($k$)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("AGP IBD cohort\n(sex + age category + BMI category)")
    ax.legend(loc="upper left")

    # ── Panel b: historical benchmarks ───────────────────────────────────────
    ax = axes[1]
    for label in ["SPSS FUZZY", "R Matching", "qupid"]:
        data = HISTORICAL[label]
        ks = sorted(data)
        _plot_series(ax, ks, [data[k] for k in ks], label)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings ($k$)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("Wisconsin 16S\n(sex + age ±4.5 yr; Patel et al.)")
    ax.legend(loc="upper left")

    for ax, lab in zip(axes, ("a", "b")):
        _add_panel_label(ax, lab)

    out = pathlib.Path("figures")
    out.mkdir(exist_ok=True)
    for fmt in NATURE_FORMATS:
        fig.savefig(
            out / f"fig_runtime_external.{fmt}",
            dpi=NATURE_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
            format=fmt,
            facecolor="white",
        )
        print(f"Saved: figures/fig_runtime_external.{fmt}")
    plt.close(fig)


if __name__ == "__main__":
    main()
