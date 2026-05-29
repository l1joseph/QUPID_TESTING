"""Supplementary Figure 1 — k=100 convergence check.

Running mean R² ± 1 SD for the primary case-control (is_case) comparison,
as a function of the number of matching iterations k, for all three cohorts.
Demonstrates that distribution statistics stabilise well before k=100.
"""

import pathlib
import string

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# ── Style (Nature / double-column) ───────────────────────────────────────────
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
NATURE_MAX_H = 9.0

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
    "figure.constrained_layout.use": True,
    "image.interpolation": "nearest",
    "image.cmap": "viridis",
    "mathtext.fontset": "dejavusans",
}
PANEL_CFG = {"fontsize": 8, "fontweight": "bold"}

sns.set_theme(rc=NATURE_RC)
plt.rcParams.update(NATURE_RC)
plt.rcParams["axes.prop_cycle"] = plt.cycler(color=PALETTE)


# ── Data ─────────────────────────────────────────────────────────────────────
BASE = pathlib.Path(__file__).parent.parent  # …/benchmarking/

COHORTS = [
    {
        "name": "AGP",
        "label": "AGP IBD\n(Bray-Curtis, 169 pairs)",
        "pre_csv": BASE / "../agp_ccm_analysis_controlled/AGP_cc_pre_all.csv",
        "post_csv": BASE / "../agp_ccm_analysis_controlled/AGP_cc_post_all.csv",
    },
    {
        "name": "HMP2",
        "label": "HMP2 IBD\n(UniFrac, 259 pairs)",
        "pre_csv": BASE
        / "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_cc_pre_all.csv",
        "post_csv": BASE
        / "../final_notebooks/hmp2_ccm_analysis_controlled/HMP2_cc_post_all.csv",
    },
    {
        "name": "THDMI",
        "label": "THDMI unhealthy\n(UniFrac, 563 pairs)",
        "pre_csv": BASE
        / "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_cc_pre_all.csv",
        "post_csv": BASE
        / "../final_notebooks/thdmi_ccm_analysis_controlled/THDMI_cc_post_all.csv",
    },
]

COL_PRE = PALETTE_TWO["blue_red"][1]  # red  — pre-CCM
COL_POST = PALETTE_TWO["blue_red"][0]  # blue — post-CCM


def _running_stats(values: np.ndarray):
    """Return running mean and SD (ddof=1) vectors for k=1..n."""
    n = len(values)
    rmean = np.array([values[:k].mean() for k in range(1, n + 1)])
    rsd = np.array([values[:k].std(ddof=1) if k > 1 else 0.0 for k in range(1, n + 1)])
    return rmean, rsd


def _add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.14,
        1.06,
        label,
        fontsize=PANEL_CFG["fontsize"],
        fontweight=PANEL_CFG["fontweight"],
        va="top",
        ha="right",
        transform=ax.transAxes,
    )


# ── Figure ───────────────────────────────────────────────────────────────────
w = NATURE_WIDTHS["double"]
h = 2.6
fig, axes = plt.subplots(1, 3, figsize=(w, h))

for ax, cohort, panel_label in zip(axes, COHORTS, string.ascii_lowercase):
    pre_df = pd.read_csv(cohort["pre_csv"])
    post_df = pd.read_csv(cohort["post_csv"])

    pre_vals = (
        pre_df[pre_df["category"] == "is_case"]
        .sort_values("iteration")["r_squared"]
        .values
    )
    post_vals = (
        post_df[post_df["category"] == "is_case"]
        .sort_values("iteration")["r_squared"]
        .values
    )

    n = len(pre_vals)
    ks = np.arange(1, n + 1)

    pre_mean, pre_sd = _running_stats(pre_vals)
    post_mean, post_sd = _running_stats(post_vals)

    # Pre-CCM line + band
    ax.plot(ks, pre_mean, color=COL_PRE, lw=0.9, label="Pre-CCM (bootstrap)")
    ax.fill_between(
        ks, pre_mean - pre_sd, pre_mean + pre_sd, color=COL_PRE, alpha=0.15, linewidth=0
    )

    # Post-CCM line + band
    ax.plot(ks, post_mean, color=COL_POST, lw=0.9, label="Post-CCM (matched)")
    ax.fill_between(
        ks,
        post_mean - post_sd,
        post_mean + post_sd,
        color=COL_POST,
        alpha=0.15,
        linewidth=0,
    )

    ax.set_xlabel("Number of iterations (k)", labelpad=3)
    ax.set_ylabel("Running mean R² (%)", labelpad=3)
    ax.set_title(cohort["label"], pad=4)
    ax.set_xlim(1, n)
    ax.tick_params(axis="both", which="major")

    if ax is axes[0]:
        ax.legend(loc="lower right", fontsize=5)

    _add_panel_label(ax, panel_label)

OUTNAME = "figS1_convergence"
FORMATS = ["pdf", "png"]
DPI = 450
out_dir = pathlib.Path(__file__).parent
out_dir.mkdir(exist_ok=True)

for fmt in FORMATS:
    fig.savefig(
        out_dir / f"{OUTNAME}.{fmt}",
        dpi=DPI,
        bbox_inches="tight",
        pad_inches=0.02,
        format=fmt,
        facecolor="white",
    )
    print(f"Saved: figures/{OUTNAME}.{fmt}")

plt.close(fig)
