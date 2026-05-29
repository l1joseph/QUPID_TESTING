#!/usr/bin/env python3
"""
Figure 2 — AGP + HMP2 CCM effect size results (4-panel, Nature double column).
Panels: a) AGP distribution, b) AGP scatter, c) HMP2 distribution, d) HMP2 scatter.
Run from benchmarking/.
"""
import string
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
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

COLOR_PRE = PALETTE_TWO["blue_red"][1]  # warm red  #CC3311
COLOR_POST = PALETTE_TWO["blue_red"][0]  # cool blue #4477AA
COLOR_FIT = PALETTE[0]  # dark blue #332288


def _add_panel_label(ax, label):
    ax.text(
        -0.18,
        1.06,
        label,
        fontsize=NATURE_PANEL["fontsize"],
        fontweight=NATURE_PANEL["fontweight"],
        va="top",
        ha="right",
        transform=ax.transAxes,
    )


def _violin_panel(ax, pre_vals, post_vals, pre_label, post_label, title, ylim=None):
    """Draw violin + strip for two groups."""
    df = pd.DataFrame(
        {
            "R²": np.concatenate([pre_vals, post_vals]),
            "Group": [pre_label] * len(pre_vals) + [post_label] * len(post_vals),
        }
    )
    order = [pre_label, post_label]
    palette = {pre_label: COLOR_PRE, post_label: COLOR_POST}

    sns.violinplot(
        data=df,
        x="Group",
        y="R²",
        order=order,
        palette=palette,
        inner=None,
        linewidth=0.6,
        cut=0,
        ax=ax,
    )
    sns.stripplot(
        data=df,
        x="Group",
        y="R²",
        order=order,
        hue="Group",
        hue_order=order,
        palette=palette,
        size=1.8,
        jitter=True,
        alpha=0.5,
        legend=False,
        ax=ax,
        zorder=3,
    )

    # means as horizontal lines
    for i, (vals, color) in enumerate([(pre_vals, COLOR_PRE), (post_vals, COLOR_POST)]):
        m = vals.mean()
        ax.hlines(
            m,
            i - 0.35,
            i + 0.35,
            colors=color,
            linewidths=1.2,
            linestyles="--",
            zorder=4,
        )
        ax.text(
            i + 0.38, m, f"{m:.2f}%", va="center", ha="left", fontsize=5, color=color
        )

    ax.set_xlabel("")
    ax.set_ylabel("PERMANOVA R² (%)")
    ax.set_title(title, fontsize=7)
    if ylim:
        ax.set_ylim(ylim)

    # significance annotation
    pre_m, post_m = pre_vals.mean(), post_vals.mean()
    _, pval = stats.ttest_ind(pre_vals, post_vals)
    stars = (
        "***"
        if pval < 0.001
        else ("**" if pval < 0.01 else ("*" if pval < 0.05 else "ns"))
    )
    y_top = max(pre_vals.max(), post_vals.max()) * 1.05
    ax.annotate(
        "",
        xy=(1, y_top),
        xytext=(0, y_top),
        arrowprops=dict(arrowstyle="-", color="black", lw=0.6),
    )
    ax.text(0.5, y_top * 1.01, stars, ha="center", va="bottom", fontsize=6)


def _scatter_panel(
    ax,
    pre,
    post,
    label_col,
    highlight_label,
    title,
    pearson_r,
    color_by_direction=False,
):
    """Scatter of pre vs post means with error bars and identity line."""
    x = pre["r_squared_mean"].values
    y = post["r_squared_mean"].values
    xe = pre["r_squared_std"].values
    ye = post["r_squared_std"].values

    if color_by_direction:
        colors = [COLOR_POST if yi > xi else COLOR_PRE for xi, yi in zip(x, y)]
    else:
        colors = [PALETTE[2]] * len(x)

    lims = [min(x.min(), y.min()) * 0.85, max(x.max(), y.max()) * 1.08]
    ax.plot(lims, lims, color="gray", lw=0.6, ls="--", zorder=0, label="y = x")

    for xi, yi, xe_i, ye_i, c in zip(x, y, xe, ye, colors):
        ax.errorbar(
            xi,
            yi,
            xerr=xe_i,
            yerr=ye_i,
            fmt="none",
            ecolor=c,
            elinewidth=0.4,
            capsize=1.5,
            alpha=0.6,
            zorder=1,
        )
    ax.scatter(x, y, c=colors, s=12, zorder=2, linewidths=0.3, edgecolors="white")

    # regression line
    slope, intercept, *_ = stats.linregress(x, y)
    xfit = np.linspace(lims[0], lims[1], 200)
    ax.plot(xfit, slope * xfit + intercept, color=COLOR_FIT, lw=0.9, zorder=3)

    # label highlight point
    cats = pre[label_col].values
    for xi, yi, cat in zip(x, y, cats):
        if cat == highlight_label:
            ax.annotate(
                cat,
                (xi, yi),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=4.5,
                color="black",
            )

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Pre-CCM R² (%)")
    ax.set_ylabel("Post-CCM R² (%)")
    ax.set_title(title, fontsize=7)
    ax.text(
        0.05,
        0.95,
        f"r = {pearson_r:.3f}\nn = {len(x)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=5.5,
    )


def main():
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)

    # ── Load data ─────────────────────────────────────────────────────────
    agp_pre_all = pd.read_csv("../agp_ccm_analysis_controlled/AGP_cc_pre_all.csv")
    agp_post_all = pd.read_csv("../agp_ccm_analysis_controlled/AGP_cc_post_all.csv")
    agp_pre_sum = pd.read_csv("../agp_ccm_analysis_controlled/AGP_pre_bootstrap.csv")
    agp_post_sum = pd.read_csv("../agp_ccm_analysis_controlled/AGP_post_ccm.csv")
    hmp2_pre_sum = pd.read_csv(
        "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_pre_bootstrap.csv"
    )
    hmp2_post_sum = pd.read_csv(
        "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_post_ccm.csv"
    )

    agp_pre_vals = agp_pre_all[agp_pre_all["category"] == "is_case"]["r_squared"].values
    agp_post_vals = agp_post_all[agp_post_all["category"] == "is_case"][
        "r_squared"
    ].values

    # HMP2 IBD — real per-iteration R² values from CCM analysis
    hmp2_pre_all = pd.read_csv(
        "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_cc_pre_all.csv"
    )
    hmp2_post_all = pd.read_csv(
        "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_cc_post_all.csv"
    )
    hmp2_pre_vals = hmp2_pre_all[hmp2_pre_all["category"] == "is_case"][
        "r_squared"
    ].values
    hmp2_post_vals = hmp2_post_all[hmp2_post_all["category"] == "is_case"][
        "r_squared"
    ].values

    # Merge summary for scatter panels
    agp_merged = agp_pre_sum.merge(agp_post_sum, on="category", suffixes=("", "_post"))
    agp_merged = agp_merged.rename(
        columns={"r_squared_mean": "r_squared_mean", "r_squared_std": "r_squared_std"}
    )
    agp_pre_s = agp_merged[["category", "r_squared_mean", "r_squared_std"]].copy()
    agp_post_s = agp_merged[
        ["category", "r_squared_mean_post", "r_squared_std_post"]
    ].rename(
        columns={
            "r_squared_mean_post": "r_squared_mean",
            "r_squared_std_post": "r_squared_std",
        }
    )

    hmp2_merged = hmp2_pre_sum.merge(
        hmp2_post_sum, on="category", suffixes=("", "_post")
    )
    hmp2_pre_s = hmp2_merged[["category", "r_squared_mean", "r_squared_std"]].copy()
    hmp2_post_s = hmp2_merged[
        ["category", "r_squared_mean_post", "r_squared_std_post"]
    ].rename(
        columns={
            "r_squared_mean_post": "r_squared_mean",
            "r_squared_std_post": "r_squared_std",
        }
    )

    # ── Layout ────────────────────────────────────────────────────────────
    w = NATURE_WIDTHS["double"]
    h = min(w * 1.0, NATURE_MAX_H)
    fig, axes = plt.subplots(2, 2, figsize=(w, h))
    fig.subplots_adjust(
        left=0.11, right=0.97, top=0.93, bottom=0.10, hspace=0.55, wspace=0.40
    )

    # ── Panel a: AGP distribution ─────────────────────────────────────────
    _violin_panel(
        axes[0, 0],
        agp_pre_vals,
        agp_post_vals,
        "Pre-CCM\n(bootstrap)",
        "Post-CCM\n(matched)",
        "AGP IBD effect size\n(169 cases, Bray-Curtis)",
    )

    # ── Panel b: AGP scatter ──────────────────────────────────────────────
    _scatter_panel(
        axes[0, 1],
        agp_pre_s,
        agp_post_s,
        "category",
        "ibd",
        "AGP effect sizes across categories\n(Pearson r = 0.982)",
        pearson_r=0.982,
        color_by_direction=False,
    )

    # ── Panel c: HMP2 distribution ────────────────────────────────────────
    _violin_panel(
        axes[1, 0],
        hmp2_pre_vals,
        hmp2_post_vals,
        "Pre-CCM\n(bootstrap)",
        "Post-CCM\n(matched)",
        "HMP2 IBD effect size\n(259 cases, is_case, Unweighted UniFrac)",
    )
    axes[1, 0].text(
        0.5,
        0.02,
        "Mean: 2.08% → 1.80%  |  SD: 0.30% → 0.04% (~8×)",
        transform=axes[1, 0].transAxes,
        ha="center",
        va="bottom",
        fontsize=4.5,
        color="gray",
        style="italic",
    )

    # ── Panel d: HMP2 scatter ─────────────────────────────────────────────
    _scatter_panel(
        axes[1, 1],
        hmp2_pre_s,
        hmp2_post_s,
        "category",
        "diagnosis",
        "HMP2 effect sizes across categories\n(Pearson r = 0.960)",
        pearson_r=0.960,
        color_by_direction=True,
    )

    # ── Panel labels ──────────────────────────────────────────────────────
    for i, ax in enumerate(axes.flatten()):
        _add_panel_label(ax, string.ascii_lowercase[i])

    # ── Save ──────────────────────────────────────────────────────────────
    out = pathlib.Path("figures")
    out.mkdir(exist_ok=True)
    for fmt in NATURE_FORMATS:
        fig.savefig(
            out / f"fig2_agp_hmp2.{fmt}",
            dpi=NATURE_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
            format=fmt,
            facecolor="white",
        )
        print(f"Saved: figures/fig2_agp_hmp2.{fmt}")
    plt.close(fig)


if __name__ == "__main__":
    main()
