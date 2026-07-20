#!/usr/bin/env python3
"""
fig_runtime_external — qupid vs. external matching tools.

Two-panel figure.

Panel (a): new AGP IBD benchmark. Qupid vs. four widely used host-covariate
matching tools (MatchIt, R Matching, CEM, q2-matchmaker) under categorical
matching (sex + age_cat + bmi_cat). Reads benchmark_real/agp_external_results.tsv.

Panel (b): historical Wisconsin 16S dementia-AD benchmark.
SPSS FUZZY, R Matching, and qupid timings at k ∈ {10, 100, 1000, 10000},
sex + age ±4.5 yr. Hard-coded so the figure is self-contained.

Both panels share log-log axes for direct visual comparison of scale; SPSS
FUZZY hitting ~40,000 s at k = 1,000 vs qupid at < 1 s is the headline
order-of-magnitude gap motivating Qupid.

Run from benchmarking/  (reads benchmark_real/agp_external_results.tsv).
"""
import pathlib

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

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
    "legend.fontsize": 6,
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
NATURE_FORMATS = ["pdf", "png"]
NATURE_DPI = 450

# ── Per-tool styling (Paul Tol colorblind-safe palette) ──────────────────────
STYLE = {
    "qupid": {"color": "#0077BB", "ls": "-", "marker": "o"},
    "MatchIt": {"color": "#EE7733", "ls": "--", "marker": "s"},
    "R Matching": {"color": "#009988", "ls": "-.", "marker": "^"},
    "CEM": {"color": "#CC3311", "ls": ":", "marker": "D"},
    "q2-matchmaker": {"color": "#EE3377", "ls": (0, (3, 1, 1, 1)), "marker": "v"},
    "SPSS FUZZY": {"color": "#CC3311", "ls": ":", "marker": "D"},
}

# ── Historical Wisconsin 16S benchmark ────────────────────────────────────────
#    Dementia-AD cohort, sex + age ±4.5 yr.
HISTORICAL = {
    "SPSS FUZZY": {10: 11.75, 100: 498.6, 1000: 39654.0},
    "R Matching": {10: 0.112, 100: 1.162, 1000: 11.84, 10000: 113.7},
    "qupid": {10: 0.015, 100: 0.069, 1000: 0.68, 10000: 6.48},
}


def _plot_panel(ax, sources, title):
    """Plot one panel; `sources` is a list of (label, k_array, t_array)."""
    for label, k_arr, t_arr in sources:
        if not len(k_arr):
            continue
        s = STYLE[label]
        ax.plot(
            k_arr,
            t_arr,
            marker=s["marker"],
            markersize=3.5,
            linewidth=1.0,
            linestyle=s["ls"],
            color=s["color"],
            label=label,
            clip_on=False,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings ($k$)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(title)
    ax.legend(loc="upper left")


def _add_panel_label(ax, label):
    ax.text(
        -0.16,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=8,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def main():
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)

    # Two-panel layout, Nature double-column width
    w = NATURE_WIDTHS["double"]
    h = w * 0.45
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(w, h))

    # ── Panel a: new AGP IBD benchmark ───────────────────────────────────────
    panel_a_sources = []
    agp_tsv = pathlib.Path("benchmark_real/agp_external_results.tsv")
    if agp_tsv.exists():
        df = pd.read_csv(agp_tsv, sep="\t")
        for tool in ["qupid", "MatchIt", "R Matching", "CEM", "q2-matchmaker"]:
            sub = df[df["tool"] == tool].sort_values("k")
            if sub.empty:
                continue
            panel_a_sources.append((tool, sub["k"].values, sub["elapsed_sec"].values))
    _plot_panel(
        ax_a,
        panel_a_sources,
        "AGP IBD cohort\n(sex + age category + BMI category)",
    )

    # ── Panel b: historical Wisconsin 16S ────────────────────────────────────
    panel_b_sources = []
    for tool in ["qupid", "R Matching", "SPSS FUZZY"]:
        kv = sorted(HISTORICAL[tool].keys())
        tv = [HISTORICAL[tool][k] for k in kv]
        panel_b_sources.append((tool, kv, tv))
    _plot_panel(
        ax_b,
        panel_b_sources,
        "Wisconsin 16S dementia-AD\n(historical; sex + age ±4.5 yr)",
    )

    _add_panel_label(ax_a, "a")
    _add_panel_label(ax_b, "b")

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
