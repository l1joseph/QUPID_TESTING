"""
hmp2_covariate_balance.py
==========================

Demonstrate Qupid's `compute_covariate_balance` (Supplementary Figure 4)
on the HMP2 IBD cohort.

Compute standardized mean differences (SMD; Cohen's d) for each matching
covariate before and after 1:1 case-control matching, then render a
Love plot — the standard covariate-balance visualization in causal
inference.

Output:
    figures/figS4_hmp2_covariate_balance.{png,pdf}
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from qupid import (
    compute_covariate_balance,
    match_by_multiple,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

HMP2_METADATA = Path(
    "/Users/leojo/Developer/knight_lab/ms_thesis/microbiome_mechinterp/"
    "prelim_analysis/data/metadata/72996_72996_analysis_mapping.txt"
)
OUT_DIR = Path(__file__).parent / "figures"
OUT_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Load + prep
# ---------------------------------------------------------------------------

print("Loading full HMP2 metadata...")
md = pd.read_table(HMP2_METADATA, index_col=0, low_memory=False)
print(f"  raw shape: {md.shape}")

# Encode categorical matching variables as binary indicators so SMD
# (Cohen's d on a 0/1 variable) is well-defined for each level.
md["is_female"] = (md["sex"] == "female").astype(float)
md["is_male"] = (md["sex"] == "male").astype(float)
md["site_Cincinnati"] = (md["site_name"] == "Cincinnati").astype(float)
md["site_MGH"] = (md["site_name"] == "MGH").astype(float)
md["site_MGH_Pediatrics"] = (md["site_name"] == "MGH Pediatrics").astype(float)

# Restrict to samples with complete matching-variable metadata
matching_cols = ["sex", "race", "site_name", "consent_age"]
md = md.dropna(subset=matching_cols)
print(f"  after dropna on matching cols: {md.shape}")

# Case/control split
md["is_ibd"] = md["diagnosis"].isin({"CD", "UC"}).astype(int)
focus = md[md["diagnosis"].isin({"CD", "UC"})].copy()
background = md[md["diagnosis"] == "nonIBD"].copy()
print(f"  full HMP2: {len(focus)} IBD cases, {len(background)} nonIBD controls")


# ---------------------------------------------------------------------------
# Match (same configuration as the main analysis)
# ---------------------------------------------------------------------------

print("\nRunning Qupid matching (sex + race + site + consent_age ±5 yr)...")
cm_one_to_many = match_by_multiple(
    focus=focus,
    background=background,
    categories=["sex", "race", "site_name", "consent_age"],
    tolerance_map={"consent_age": 5.0},
    on_failure="continue",
)
print(f"  matched {len(cm_one_to_many.cases)} cases with at least one valid control")

# One representative 1:1 matching (consistent with Supp Fig 3's framing of
# the marginal distribution being preserved under Qupid's sampler)
collection = cm_one_to_many.create_matched_pairs(iterations=1, strict=False, seed=42)
matching = collection[0]
print(f"  representative matched pair count: {len(matching.cases)}")


# ---------------------------------------------------------------------------
# Compute SMD pre vs post matching
# ---------------------------------------------------------------------------

# Cohen's d is only defined for numeric/boolean columns, so we evaluate the
# matching variables in their numeric/binary-encoded form.
balance_cols = [
    "is_female",
    "is_male",
    "site_Cincinnati",
    "site_MGH",
    "site_MGH_Pediatrics",
    "consent_age",
]

print("\nComputing covariate balance...")
balance = compute_covariate_balance(
    focus=focus,
    background=background,
    casematch=matching,
    categories=balance_cols,
)
print(balance.to_string(index=False))

# ---------------------------------------------------------------------------
# Love plot
# ---------------------------------------------------------------------------

# Pretty display names
display_names = {
    "is_female": "Sex (female)",
    "is_male": "Sex (male)",
    "site_Cincinnati": "Site: Cincinnati",
    "site_MGH": "Site: MGH",
    "site_MGH_Pediatrics": "Site: MGH Pediatrics",
    "consent_age": "Consent age",
}
balance["display"] = balance["covariate"].map(display_names)

# Order from largest |smd_pre| to smallest
balance = balance.reindex(balance["smd_pre"].abs().sort_values(ascending=True).index)

sns.set_style("whitegrid")
plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

fig, ax = plt.subplots(figsize=(6.5, 4.0))

y_pos = np.arange(len(balance))
PRE_COLOR = "#CC3311"
POST_COLOR = "#0077BB"

ax.scatter(
    np.abs(balance["smd_pre"]),
    y_pos,
    s=70,
    color=PRE_COLOR,
    label="Pre-CCM (unmatched)",
    zorder=3,
    edgecolor="white",
    linewidth=1,
)
ax.scatter(
    np.abs(balance["smd_post"]),
    y_pos,
    s=70,
    color=POST_COLOR,
    label="Post-CCM (matched)",
    zorder=3,
    edgecolor="white",
    linewidth=1,
)

# Connect pre/post with thin lines so the eye tracks the improvement
for i, row in enumerate(balance.itertuples()):
    ax.plot(
        [abs(row.smd_pre), abs(row.smd_post)],
        [i, i],
        color="gray",
        linewidth=0.5,
        alpha=0.5,
        zorder=2,
    )

# Reference line at the conventional well-balanced threshold
ax.axvline(0.1, color="black", linestyle="--", linewidth=0.75, alpha=0.6)
ax.text(
    0.105,
    len(balance) - 0.3,
    "|SMD| = 0.1\n(well-balanced)",
    fontsize=8,
    color="black",
    alpha=0.7,
    va="top",
)

ax.set_yticks(y_pos)
ax.set_yticklabels(balance["display"])
ax.set_xlabel("|Standardized mean difference| (Cohen's d)")
ax.set_xlim(0, max(0.3, balance["smd_pre"].abs().max() * 1.1))
ax.set_title(
    f"HMP2 covariate balance pre- vs. post-CCM\n"
    f"({len(focus)} IBD cases, {len(background)} nonIBD controls; "
    f"{len(matching.cases)} matched pairs)",
    fontsize=10,
)
ax.legend(loc="lower right", frameon=True, fontsize=9)
ax.grid(True, axis="x", alpha=0.3)

fig.tight_layout()
for ext in ("png", "pdf"):
    fig.savefig(
        OUT_DIR / f"figS4_hmp2_covariate_balance.{ext}",
        dpi=200,
        bbox_inches="tight",
    )
plt.close(fig)

print(f"\nFigure written to {OUT_DIR}/figS4_hmp2_covariate_balance.{{png,pdf}}")

# Save the numeric summary too
summary_path = OUT_DIR.parent / "benchmark_real" / "hmp2_covariate_balance_summary.tsv"
summary_path.parent.mkdir(parents=True, exist_ok=True)
balance.to_csv(summary_path, sep="\t", index=False)
print(f"Summary written to {summary_path}")
