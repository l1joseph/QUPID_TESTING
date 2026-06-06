#!/usr/bin/env python3
"""
fig_runtime_external — qupid vs. external matching tools on the AGP IBD cohort.
Single-panel, Nature single column. Runtime vs. number of matchings (k) for
qupid, MatchIt, and R Matching, all run on the same AGP IBD cohort
(categorical sex + age_cat + bmi_cat).
Run from benchmarking/  (reads benchmark_real/agp_external_results.tsv).
"""
import pathlib

import matplotlib.pyplot as plt
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

# ── Series styling (qupid blue; one distinct color per external tool) ────────
STYLE = {
    "qupid": {"color": "#0077BB", "ls": "-", "marker": "o"},
    "MatchIt": {"color": "#EE7733", "ls": "--", "marker": "s"},
    "R Matching": {"color": "#009988", "ls": "-.", "marker": "^"},
    "CEM": {"color": "#CC3311", "ls": ":", "marker": "D"},
}


def main():
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)

    w = NATURE_WIDTHS["single"]
    h = min(w * 0.78, NATURE_MAX_H)
    fig, ax = plt.subplots(1, 1, figsize=(w, h))

    df = pd.read_csv("benchmark_real/agp_external_results.tsv", sep="\t")
    for tool in ["qupid", "MatchIt", "R Matching", "CEM"]:
        sub = df[df["tool"] == tool].sort_values("k")
        if sub.empty:
            continue
        s = STYLE[tool]
        ax.plot(
            sub["k"].values,
            sub["elapsed_sec"].values,
            marker=s["marker"],
            markersize=3.5,
            linewidth=1.0,
            linestyle=s["ls"],
            color=s["color"],
            label=tool,
            clip_on=False,
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings ($k$)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("AGP IBD cohort\n(sex + age category + BMI category)")
    ax.legend(loc="upper left")

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
