#!/usr/bin/env python3
"""
Figure 3 — THDMI CCM confounding exposure (2-panel, Nature 1.5col).
Panels: a) distribution pre vs post, b) scatter across categories.
Run from benchmarking/.
"""
import string
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns

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

COLOR_PRE = PALETTE_TWO["blue_red"][1]  # #CC3311 warm red
COLOR_POST = PALETTE_TWO["blue_red"][0]  # #4477AA cool blue
COLOR_FIT = PALETTE[0]  # #332288


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


def main():
    sns.set_theme(rc=NATURE_RC)
    plt.rcParams.update(NATURE_RC)
    plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)

    # ── Load data ─────────────────────────────────────────────────────────
    pre_all = pd.read_csv(
        "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_cc_pre_all.csv"
    )
    post_all = pd.read_csv(
        "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_cc_post_all.csv"
    )
    pre_sum = pd.read_csv(
        "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_pre_bootstrap.csv"
    )
    post_sum = pd.read_csv(
        "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_post_ccm.csv"
    )

    pre_vals = pre_all[pre_all["category"] == "is_case"]["r_squared"].values
    post_vals = post_all[post_all["category"] == "is_case"]["r_squared"].values

    # Pre-CCM: standard p-value (bootstrap arm has no pair structure)
    # Post-CCM: pair-aware restricted-permutation p-value (primary significance criterion)
    pre_col = "p_value"
    post_col = (
        "p_value_restricted" if "p_value_restricted" in post_all.columns else "p_value"
    )
    pre_sig = (pre_all[pre_all["category"] == "is_case"][pre_col] < 0.05).sum()
    post_sig = (post_all[post_all["category"] == "is_case"][post_col] < 0.05).sum()

    # Summary scatter data
    merged = pre_sum.merge(post_sum, on="category", suffixes=("", "_post"))

    # ── Layout ────────────────────────────────────────────────────────────
    w = NATURE_WIDTHS["1.5col"]
    h = min(w * 0.88, NATURE_MAX_H)
    fig, axes = plt.subplots(1, 2, figsize=(w, h))
    fig.subplots_adjust(left=0.13, right=0.97, top=0.88, bottom=0.14, wspace=0.45)

    # ── Panel a: THDMI distribution ───────────────────────────────────────
    ax = axes[0]
    df_dist = pd.DataFrame(
        {
            "R²": np.concatenate([pre_vals, post_vals]),
            "Group": ["Pre-CCM\n(bootstrap)"] * len(pre_vals)
            + ["Post-CCM\n(matched)"] * len(post_vals),
        }
    )
    order = ["Pre-CCM\n(bootstrap)", "Post-CCM\n(matched)"]
    palette = {"Pre-CCM\n(bootstrap)": COLOR_PRE, "Post-CCM\n(matched)": COLOR_POST}

    sns.violinplot(
        data=df_dist,
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
        data=df_dist,
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
            i + 0.38, m, f"{m:.3f}%", va="center", ha="left", fontsize=5, color=color
        )

    # significance annotation
    _, pval = stats.ttest_ind(pre_vals, post_vals)
    y_top = pre_vals.max() * 1.05
    ax.annotate(
        "",
        xy=(1, y_top),
        xytext=(0, y_top),
        arrowprops=dict(arrowstyle="-", color="black", lw=0.6),
    )
    ax.text(0.5, y_top * 1.01, "***", ha="center", va="bottom", fontsize=6)

    # significance counts (pre: standard p, post: pair-aware restricted p)
    ax.text(
        0,
        -0.14,
        f"{pre_sig}/100\nsig. (p<0.05)",
        transform=ax.get_xaxis_transform(),
        ha="center",
        fontsize=5,
        color=COLOR_PRE,
    )
    ax.text(
        1,
        -0.14,
        f"{post_sig}/100\nsig. (restr. p<0.05)",
        transform=ax.get_xaxis_transform(),
        ha="center",
        fontsize=5,
        color=COLOR_POST,
    )

    ax.set_xlabel("")
    ax.set_ylabel("PERMANOVA R² (%)")
    ax.set_title(
        "THDMI unhealthy effect size\n(563 cases, Unweighted UniFrac)", fontsize=7
    )

    # ── Panel b: THDMI scatter ────────────────────────────────────────────
    ax = axes[1]
    x = merged["r_squared_mean"].values
    y = merged["r_squared_mean_post"].values
    xe = merged["r_squared_std"].values
    ye = merged["r_squared_std_post"].values
    cats = merged["category"].values

    # highlight matched-upon variables
    matched_vars = {"bmi_cat", "sex", "thdmi_cohort", "cosmetics_frequency", "host_age"}
    colors = ["#882255" if c in matched_vars else PALETTE[2] for c in cats]

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

    # regression
    slope, intercept, r_val, p_val, _ = stats.linregress(x, y)
    xfit = np.linspace(lims[0], lims[1], 200)
    ax.plot(xfit, slope * xfit + intercept, color=COLOR_FIT, lw=0.9, zorder=3)

    # label matching variables
    for xi, yi, cat in zip(x, y, cats):
        if cat in matched_vars:
            ax.annotate(
                cat.replace("_", " "),
                (xi, yi),
                xytext=(3, 3),
                textcoords="offset points",
                fontsize=4,
                color="#882255",
            )

    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("Pre-CCM R² (%)")
    ax.set_ylabel("Post-CCM R² (%)")
    r_pearson, _ = stats.pearsonr(x, y)
    ax.set_title("THDMI effect sizes across categories", fontsize=7)
    ax.text(
        0.05,
        0.95,
        f"r = {r_pearson:.3f}\nn = {len(x)}",
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=5.5,
    )

    # legend for matched variable highlight
    from matplotlib.lines import Line2D

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor="#882255",
            markersize=4,
            label="Matched variable",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            markerfacecolor=PALETTE[2],
            markersize=4,
            label="Other covariate",
        ),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=4.5)

    # ── Panel labels ──────────────────────────────────────────────────────
    for i, ax in enumerate(axes):
        _add_panel_label(ax, string.ascii_lowercase[i])

    # ── Save ──────────────────────────────────────────────────────────────
    out = pathlib.Path("figures")
    out.mkdir(exist_ok=True)
    for fmt in NATURE_FORMATS:
        fig.savefig(
            out / f"fig3_thdmi.{fmt}",
            dpi=NATURE_DPI,
            bbox_inches="tight",
            pad_inches=0.02,
            format=fmt,
            facecolor="white",
        )
        print(f"Saved: figures/fig3_thdmi.{fmt}")
    plt.close(fig)


if __name__ == "__main__":
    main()
