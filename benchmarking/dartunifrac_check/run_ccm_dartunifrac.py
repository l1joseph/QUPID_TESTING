"""DartUniFrac sanity check for HMP2 and THDMI CCM analyses.

Replicates the published HMP2 + THDMI is_case CCM results, but loads the
distance matrix from a DartUniFrac TSV instead of computing UniFrac via
QIIME2/scikit-bio. All matching/bootstrap/PERMANOVA logic is copied verbatim
from final_notebooks/{HMP2,THDMI}_effect_size_controlled.ipynb so the
comparison is apples-to-apples.

Run via the qiime2-amplicon-2023.9 env which has qupid + skbio + qiime2.
"""

from __future__ import annotations

import argparse
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.distance import DistanceMatrix, permanova

from qupid import match_by_multiple

warnings.filterwarnings("ignore")

# --------------------------------------------------------------------------- #
#  PERMANOVA / restricted-permutation helpers (verbatim from notebooks)
# --------------------------------------------------------------------------- #


def calculate_effect_size_permanova(dm, metadata_df, category, permutations=999):
    common_ids = list(set(dm.ids) & set(metadata_df.index))
    if len(common_ids) < 10:
        return None
    dm_filtered = dm.filter(common_ids)
    md_filtered = metadata_df.loc[common_ids, [category]].dropna()
    common_ids = list(set(dm_filtered.ids) & set(md_filtered.index))
    if len(common_ids) < 10:
        return None
    dm_filtered = dm_filtered.filter(common_ids)
    md_filtered = md_filtered.loc[common_ids]
    grouping = md_filtered[category]
    n_groups = grouping.nunique()
    n_samples = len(common_ids)
    if n_groups < 2:
        return None
    try:
        results = permanova(dm_filtered, grouping, permutations=permutations)
        pseudo_f = results["test statistic"]
        df_between = n_groups - 1
        df_within = n_samples - n_groups
        r_squared = (pseudo_f * df_between) / (pseudo_f * df_between + df_within) * 100
        return {
            "r_squared": r_squared,
            "pseudo_f": pseudo_f,
            "p_value": results["p-value"],
            "n_samples": n_samples,
            "n_groups": n_groups,
        }
    except Exception:
        return None


def _pseudo_f_from_dm(dm_array, labels):
    n = len(labels)
    unique_groups = np.unique(labels)
    k = len(unique_groups)
    if k < 2 or n <= k:
        return np.nan
    total_ss = np.sum(dm_array**2) / (2 * n)
    within_ss = 0.0
    for g in unique_groups:
        idx = np.where(labels == g)[0]
        n_g = len(idx)
        if n_g >= 2:
            sub = dm_array[np.ix_(idx, idx)]
            within_ss += np.sum(sub**2) / (2 * n_g)
    between_ss = total_ss - within_ss
    df_between = k - 1
    df_within = n - k
    if df_within <= 0 or within_ss <= 0:
        return np.nan
    return (between_ss / df_between) / (within_ss / df_within)


def restricted_permanova_pvalue(dm, pair_map, is_case_series, permutations=999, seed=0):
    rng = np.random.default_rng(seed)
    valid_pairs = [
        (case, ctrl)
        for case, ctrl in pair_map.items()
        if (
            case in dm.ids
            and ctrl in dm.ids
            and case in is_case_series.index
            and ctrl in is_case_series.index
        )
    ]
    if len(valid_pairs) < 5:
        return np.nan
    sample_ids: list[str] = []
    labels: list[float] = []
    for case, ctrl in valid_pairs:
        sample_ids.extend([case, ctrl])
        labels.extend([float(is_case_series[case]), float(is_case_series[ctrl])])
    dm_sub = dm.filter(sample_ids)
    dm_arr = np.array(dm_sub.data)
    obs_labels = np.array(labels)
    obs_f = _pseudo_f_from_dm(dm_arr, obs_labels)
    if np.isnan(obs_f):
        return np.nan
    n_pairs = len(valid_pairs)
    n_exceed = 0
    for _ in range(permutations):
        perm_labels = obs_labels.copy()
        flip = rng.random(n_pairs) < 0.5
        for p_idx in np.where(flip)[0]:
            i = p_idx * 2
            perm_labels[i], perm_labels[i + 1] = perm_labels[i + 1], perm_labels[i]
        null_f = _pseudo_f_from_dm(dm_arr, perm_labels)
        if not np.isnan(null_f) and null_f >= obs_f:
            n_exceed += 1
    return (n_exceed + 1) / (permutations + 1)


# --------------------------------------------------------------------------- #
#  DM loader
# --------------------------------------------------------------------------- #


def load_dartunifrac_dm(tsv_path: str) -> DistanceMatrix:
    """Read a square TSV (first column = sample IDs, first row = headers)."""
    # Force the index column to be parsed as string — sample IDs like
    # ``10317.000114575`` get coerced to float64 otherwise, which then breaks
    # ``reindex(columns=ids)`` because the column labels remain strings.
    df = pd.read_csv(tsv_path, sep="\t", index_col=0, dtype={0: str})
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    ids = df.index.tolist()
    df = df.reindex(index=ids, columns=ids)
    arr = df.to_numpy(dtype=float)
    # DartUniFrac TSV may have small floating-point asymmetry; symmetrize.
    arr = (arr + arr.T) / 2.0
    np.fill_diagonal(arr, 0.0)
    return DistanceMatrix(arr, ids=ids)


# --------------------------------------------------------------------------- #
#  Main per-cohort runner
# --------------------------------------------------------------------------- #


def run_cohort(
    cohort: str,
    dm_path: str,
    metadata_df: pd.DataFrame,
    matching_categories: list[str],
    numeric_categories: list[str],
    tolerance_map: dict,
    case_column: str,
    case_values: list,
    control_value,
    out_dir: Path,
    n_iterations: int = 100,
    n_bootstrap: int = 100,
    n_permutations: int = 999,
    seed: int = 42,
):
    print(f"\n{'=' * 70}\n  {cohort}  (DartUniFrac unweighted)\n{'=' * 70}")
    np.random.seed(seed)

    # Case / control split (mirrors notebook prepare_case_control_groups)
    if not isinstance(case_values, list):
        case_values = [case_values]
    all_values = case_values + [control_value]
    relevant = metadata_df[metadata_df[case_column].isin(all_values)].copy()
    relevant["is_case"] = relevant[case_column].isin(case_values).astype(int)

    # Load DartUniFrac DM
    dm_full = load_dartunifrac_dm(dm_path)
    print(f"DartUniFrac DM: {dm_full.shape[0]} samples")

    # Restrict to common samples
    common = list(set(relevant.index) & set(dm_full.ids))
    relevant = relevant.loc[relevant.index.isin(common)]
    dm_full = dm_full.filter(common)
    focus = relevant[relevant["is_case"] == 1]
    background = relevant[relevant["is_case"] == 0]
    print(f"After filtering to DM: {len(focus)} cases, {len(background)} controls")

    # Run CCM
    all_cats = matching_categories + [
        c for c in numeric_categories if c not in matching_categories
    ]
    focus_clean = focus.dropna(subset=all_cats)
    background_clean = background.dropna(subset=all_cats)
    print(
        f"After dropping NaN in matching covariates: "
        f"{len(focus_clean)} cases, {len(background_clean)} controls"
    )
    cm = match_by_multiple(
        focus=focus_clean,
        background=background_clean,
        categories=all_cats,
        tolerance_map=tolerance_map,
        on_failure="warn",
    )
    matched_pairs = cm.create_matched_pairs(
        iterations=n_iterations, strict=False, seed=seed
    )

    cm_single = matched_pairs[0]
    n_cases = len(list(cm_single.cases))
    n_controls = len(list(cm_single.controls))
    print(f"Post-CCM per iter: {n_cases} cases + {n_controls} controls")

    # ---- Pre-CCM bootstrap on is_case ----
    case_ids = [s for s in relevant[relevant["is_case"] == 1].index if s in dm_full.ids]
    control_ids = [
        s for s in relevant[relevant["is_case"] == 0].index if s in dm_full.ids
    ]
    print(
        f"Available for bootstrap: {len(case_ids)} cases, {len(control_ids)} controls"
    )

    rng = np.random.RandomState(seed + 1)  # notebook uses seed+1 for cc_pre bootstrap
    pre_r2: list[float] = []
    pre_p: list[float] = []
    for i in range(n_bootstrap):
        sampled_cases = rng.choice(case_ids, size=n_cases, replace=False)
        sampled_controls = rng.choice(control_ids, size=n_controls, replace=False)
        sampled = list(sampled_cases) + list(sampled_controls)
        dm_sub = dm_full.filter(sampled)
        md_sub = relevant.loc[sampled]
        eff = calculate_effect_size_permanova(
            dm_sub, md_sub, "is_case", permutations=n_permutations
        )
        if eff is not None:
            pre_r2.append(eff["r_squared"])
            pre_p.append(eff["p_value"])
        if (i + 1) % 20 == 0:
            print(f"  bootstrap {i + 1}/{n_bootstrap}")

    # ---- Post-CCM matched is_case ----
    post_r2: list[float] = []
    post_p: list[float] = []
    post_p_restricted: list[float] = []
    for i in range(n_iterations):
        cm_i = matched_pairs[i]
        cases_i = list(cm_i.cases)
        ctrls_i = list(cm_i.controls)
        matched = [s for s in cases_i + ctrls_i if s in dm_full.ids]
        if len(matched) < 20:
            continue
        dm_matched = dm_full.filter(matched)
        md_matched = relevant.loc[relevant.index.isin(matched)]
        eff = calculate_effect_size_permanova(
            dm_matched, md_matched, "is_case", permutations=n_permutations
        )
        if eff is None:
            continue
        post_r2.append(eff["r_squared"])
        post_p.append(eff["p_value"])
        pair_map = {
            case: list(ctrl_set)[0] for case, ctrl_set in cm_i.case_control_map.items()
        }
        p_rest = restricted_permanova_pvalue(
            dm_matched,
            pair_map,
            relevant["is_case"],
            permutations=n_permutations,
            seed=i,
        )
        post_p_restricted.append(p_rest)
        if (i + 1) % 20 == 0:
            print(f"  post-CCM iter {i + 1}/{n_iterations}")

    pre_r2_arr = np.array(pre_r2)
    post_r2_arr = np.array(post_r2)
    pre_p_arr = np.array(pre_p)
    post_p_arr = np.array(post_p)
    post_p_restricted_arr = np.array(post_p_restricted)

    t_stat, p_val = stats.ttest_ind(pre_r2_arr, post_r2_arr)

    n_sig_pre = int(np.sum(pre_p_arr < 0.05))
    n_sig_post_pool = int(np.sum(post_p_arr < 0.05))
    n_sig_post_restricted = int(
        np.sum(post_p_restricted_arr[~np.isnan(post_p_restricted_arr)] < 0.05)
    )

    print(
        f"\n{cohort} RESULTS:\n"
        f"  Pre-CCM (bootstrap, n={len(pre_r2_arr)}):  "
        f"R² = {pre_r2_arr.mean():.4f} ± {pre_r2_arr.std():.4f} %\n"
        f"  Post-CCM (matched, n={len(post_r2_arr)}):  "
        f"R² = {post_r2_arr.mean():.4f} ± {post_r2_arr.std():.4f} %\n"
        f"  t-test: t={t_stat:.3f}, p={p_val:.4e}\n"
        f"  Significant iterations:\n"
        f"    pre  (pooled PERMANOVA p < 0.05):   {n_sig_pre}/{len(pre_p_arr)}\n"
        f"    post (pooled PERMANOVA p < 0.05):   {n_sig_post_pool}/{len(post_p_arr)}\n"
        f"    post (restricted PERMANOVA p<0.05): {n_sig_post_restricted}"
        f"/{int((~np.isnan(post_p_restricted_arr)).sum())}"
    )

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"r_squared": pre_r2_arr, "p_value": pre_p_arr}).to_csv(
        out_dir / f"{cohort}_pre_ccm.csv", index=False
    )
    pd.DataFrame(
        {
            "r_squared": post_r2_arr,
            "p_value": post_p_arr,
            "p_value_restricted": post_p_restricted_arr,
        }
    ).to_csv(out_dir / f"{cohort}_post_ccm.csv", index=False)
    pd.DataFrame(
        [
            {
                "cohort": cohort,
                "pre_r2_mean": pre_r2_arr.mean(),
                "pre_r2_std": pre_r2_arr.std(),
                "post_r2_mean": post_r2_arr.mean(),
                "post_r2_std": post_r2_arr.std(),
                "t_stat": t_stat,
                "p_value": p_val,
                "n_sig_pre_pooled": n_sig_pre,
                "n_sig_post_pooled": n_sig_post_pool,
                "n_sig_post_restricted": n_sig_post_restricted,
            }
        ]
    ).to_csv(out_dir / f"{cohort}_summary.csv", index=False)


# --------------------------------------------------------------------------- #
#  Cohort configs
# --------------------------------------------------------------------------- #


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cohort", choices=["hmp2", "thdmi", "both"])
    args = parser.parse_args()

    out_dir = Path(__file__).parent
    if args.cohort in ("hmp2", "both"):
        hmp2_md_path = (
            "/Users/leojo/Developer/knight_lab/ms_thesis/microbiome_mechinterp/"
            "prelim_analysis/data/metadata/72996_72996_analysis_mapping.txt"
        )
        # Metadata is QIIME2-style: first row may have a #q2:types row.
        md = pd.read_csv(hmp2_md_path, sep="\t", index_col=0, low_memory=False)
        if str(md.index[0]).startswith("#"):
            md = md.iloc[1:]
        run_cohort(
            cohort="HMP2",
            dm_path=str(out_dir / "hmp2_dartunifrac.tsv"),
            metadata_df=md,
            matching_categories=["sex", "race", "site_name"],
            numeric_categories=["consent_age"],
            tolerance_map={"consent_age": 5},
            case_column="diagnosis",
            case_values=["UC", "CD"],
            control_value="nonIBD",
            out_dir=out_dir,
        )

    if args.cohort in ("thdmi", "both"):
        thdmi_md_path = (
            "/Users/leojo/Developer/knight_lab/QUPID_TESTING/thdmi/data/"
            "thdmi_metadata_filtered_samples_clean_variables.tsv"
        )
        md = pd.read_csv(thdmi_md_path, sep="\t", index_col=0, low_memory=False)
        if str(md.index[0]).startswith("#"):
            md = md.iloc[1:]
        # Boolify unhealthy column (matches notebook)
        if md["unhealthy"].dtype == object:
            md["unhealthy"] = md["unhealthy"].map(
                {"True": True, "False": False, True: True, False: False}
            )
        # Numeric host_age
        md["host_age"] = pd.to_numeric(md["host_age"], errors="coerce")
        run_cohort(
            cohort="THDMI",
            dm_path=str(out_dir / "thdmi_dartunifrac.tsv"),
            metadata_df=md,
            matching_categories=[
                "sex",
                "thdmi_cohort",
                "bmi_cat",
                "cosmetics_frequency",
            ],
            numeric_categories=["host_age"],
            tolerance_map={"host_age": 5},
            case_column="unhealthy",
            case_values=[True],
            control_value=False,
            out_dir=out_dir,
        )


if __name__ == "__main__":
    main()
