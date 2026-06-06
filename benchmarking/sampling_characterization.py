"""
sampling_characterization.py
============================

Empirical characterization of Qupid's randomized Hopcroft-Karp sampler
relative to uniform sampling over the space of maximum-cardinality matchings.

On a small bipartite case-control graph where all maximum matchings can be
enumerated exhaustively, we:

1. Enumerate every distinct maximum-cardinality 1:1 matching by backtracking.
2. Run Qupid's create_matched_pairs() N = 100,000 times with random seeds,
   recording the matching produced on each call.
3. Compare the empirical frequency of each matching under Qupid's sampler to
   the uniform expectation (1/M, where M is the number of maximum matchings).
4. Compute total-variation distance and a chi-square goodness-of-fit test
   between the empirical and uniform distributions.
5. Also compute the marginal distribution of PERMANOVA pseudo-F across (a) all
   matchings under uniform weight and (b) all matchings under Qupid's empirical
   weight, on a synthetic distance matrix. If the two pseudo-F distributions
   are similar, this empirically validates that Qupid's effect-size
   distribution is a faithful sample of the over-all matching distribution
   even though the per-matching probabilities are not exactly uniform.

Produces Supplementary Figure 3.
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

# Qupid imports
from qupid.casematch import CaseMatchOneToMany

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

CASES = ["c0", "c1", "c2", "c3", "c4"]
CONTROLS = [f"t{i}" for i in range(8)]
N_RUNS = 100_000
SEED = 42

OUT = Path(__file__).parent / "figures"
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Test bipartite graph (realistic partial connectivity)
# ---------------------------------------------------------------------------
#
# c0 has 4 valid controls; c1 has 3; c2 has 4; c3 has 3; c4 has 4.
# Designed so that the number of maximum matchings is non-trivial but tractable
# to enumerate exhaustively.
#
CASE_CONTROL_MAP: dict[str, set[str]] = {
    "c0": {"t0", "t1", "t2", "t3"},
    "c1": {"t1", "t3", "t4"},
    "c2": {"t2", "t4", "t5", "t6"},
    "c3": {"t0", "t5", "t7"},
    "c4": {"t3", "t4", "t6", "t7"},
}


# ---------------------------------------------------------------------------
# Step 1: enumerate all maximum-cardinality matchings
# ---------------------------------------------------------------------------


def enumerate_max_matchings(
    ccm: dict[str, set[str]],
) -> list[frozenset[tuple[str, str]]]:
    """Backtracking enumeration of all maximum-cardinality 1:1 matchings."""
    cases = sorted(ccm.keys())
    n_cases = len(cases)

    best_size = [0]
    best_matchings: set[frozenset[tuple[str, str]]] = set()

    def recurse(idx: int, used: set[str], current: dict[str, str]) -> None:
        # Upper-bound pruning: even if we match every remaining case, we
        # cannot exceed best_size if we already trail by too much.
        remaining = n_cases - idx
        if len(current) + remaining < best_size[0]:
            return

        if idx == n_cases:
            if len(current) > best_size[0]:
                best_size[0] = len(current)
                best_matchings.clear()
            if len(current) == best_size[0]:
                best_matchings.add(frozenset(current.items()))
            return

        case = cases[idx]
        valid = ccm[case] - used
        # Option 1: match this case to each valid control
        for ctrl in sorted(valid):
            current[case] = ctrl
            recurse(idx + 1, used | {ctrl}, current)
            del current[case]
        # Option 2: leave this case unmatched
        recurse(idx + 1, used, current)

    recurse(0, set(), {})
    return sorted(best_matchings, key=lambda m: sorted(m))


# ---------------------------------------------------------------------------
# Step 2: Qupid empirical sampling
# ---------------------------------------------------------------------------


def qupid_sample(
    ccm: dict[str, set[str]], n_runs: int, seed: int
) -> Counter[frozenset[tuple[str, str]]]:
    """Run Qupid's create_matched_pairs n_runs times and count each matching."""
    cm = CaseMatchOneToMany(ccm)
    counts: Counter[frozenset[tuple[str, str]]] = Counter()

    # Run in batches to amortize the joblib overhead
    batch_size = 1000
    rng = np.random.default_rng(seed)
    for batch_start in range(0, n_runs, batch_size):
        n_this_batch = min(batch_size, n_runs - batch_start)
        batch_seed = int(rng.integers(0, 2**31 - 1))
        coll = cm.create_matched_pairs(
            iterations=n_this_batch, strict=True, seed=batch_seed
        )
        # NOTE: create_matched_pairs dedupes via set() so frequencies are
        # not preserved across iterations within one call. We instead run
        # iterations=1 calls so each is independently sampled.
        # Adjust the loop:
        pass

    # Proper sampling: one independent call per run with a fresh seed.
    rng = np.random.default_rng(seed)
    for _ in range(n_runs):
        run_seed = int(rng.integers(0, 2**31 - 1))
        coll = cm.create_matched_pairs(iterations=1, strict=True, seed=run_seed)
        matching = next(iter(coll))
        # Extract (case, control) pairs from CaseMatchOneToOne
        pairs = frozenset(
            (case, next(iter(ctrls)))
            for case, ctrls in matching.case_control_map.items()
        )
        counts[pairs] += 1
    return counts


# ---------------------------------------------------------------------------
# Step 3: PERMANOVA pseudo-F under each matching
# ---------------------------------------------------------------------------


def compute_pseudo_f_per_matching(
    matchings: list[frozenset[tuple[str, str]]], dm: np.ndarray, sample_ids: list[str]
) -> np.ndarray:
    """Compute PERMANOVA pseudo-F for each matching on a fixed distance matrix."""
    from qupid.stats import _fast_permanova

    id_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
    f_values = np.empty(len(matchings))
    for i, matching in enumerate(matchings):
        sample_ids_used = []
        is_case = []
        for case, ctrl in matching:
            sample_ids_used.append(case)
            is_case.append(True)
        for case, ctrl in matching:
            sample_ids_used.append(ctrl)
            is_case.append(False)
        idx = [id_to_idx[s] for s in sample_ids_used]
        sub_dm = dm[np.ix_(idx, idx)].copy()
        is_case_arr = np.array(is_case)
        f, _ = _fast_permanova(sub_dm, is_case_arr, 0, np.random.default_rng())
        f_values[i] = f
    return f_values


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def main() -> None:
    print("Step 1: enumerating all maximum-cardinality matchings...")
    all_matchings = enumerate_max_matchings(CASE_CONTROL_MAP)
    M = len(all_matchings)
    n_cases = len(CASES)
    max_size = len(next(iter(all_matchings)))
    print(f"  found M = {M} matchings, each of size {max_size} (n_cases = {n_cases})")
    assert max_size == n_cases, "graph does not admit perfect matching"

    print(f"\nStep 2: running Qupid N = {N_RUNS:,} times...")
    counts = qupid_sample(CASE_CONTROL_MAP, N_RUNS, SEED)
    print(f"  Qupid produced {len(counts)} unique matchings out of {M} possible")

    # Build empirical and uniform distributions in a fixed order
    uniform_prob = 1.0 / M
    empirical_prob = np.array([counts[m] / N_RUNS for m in all_matchings])
    uniform_dist = np.full(M, uniform_prob)

    # Total-variation distance
    tv_distance = 0.5 * np.abs(empirical_prob - uniform_dist).sum()
    print(f"\n  Total-variation distance: {tv_distance:.4f}")

    # Chi-square goodness-of-fit
    observed = np.array([counts[m] for m in all_matchings])
    expected = np.full(M, N_RUNS / M)
    chi2_stat, chi2_p = stats.chisquare(observed, expected)
    print(f"  Chi-square statistic: {chi2_stat:.2f}, p-value: {chi2_p:.4f}")

    # KL divergence (Qupid || uniform)
    eps = 1e-12
    kl = np.sum(empirical_prob * np.log((empirical_prob + eps) / uniform_prob))
    print(f"  KL divergence (Qupid || uniform): {kl:.4f}")

    # Step 3: pseudo-F distribution under uniform vs. Qupid sampling
    print("\nStep 3: computing pseudo-F under uniform vs. Qupid sampling...")
    sample_ids = CASES + CONTROLS
    n_total = len(sample_ids)
    rng = np.random.default_rng(SEED)
    raw = rng.random((n_total, n_total))
    dm = np.triu(raw, 1) + np.triu(raw, 1).T
    f_per_matching = compute_pseudo_f_per_matching(all_matchings, dm, sample_ids)

    # Uniform: each matching contributes 1/M to the pseudo-F distribution.
    # Qupid: each matching contributes count/N to the pseudo-F distribution.
    # The two distributions are weighted averages of the same f_per_matching
    # vector, so comparing them is a fair test of whether the sampler's
    # non-uniformity is consequential for downstream pseudo-F summaries.
    uniform_mean = float(f_per_matching.mean())
    uniform_std = float(f_per_matching.std())
    qupid_weights = np.array([counts[m] / N_RUNS for m in all_matchings])
    qupid_mean = float(np.average(f_per_matching, weights=qupid_weights))
    qupid_var = float(
        np.average((f_per_matching - qupid_mean) ** 2, weights=qupid_weights)
    )
    qupid_std = qupid_var**0.5
    print(f"  Uniform pseudo-F: mean = {uniform_mean:.4f}, SD = {uniform_std:.4f}")
    print(f"  Qupid pseudo-F:   mean = {qupid_mean:.4f}, SD = {qupid_std:.4f}")
    mean_relative_diff = abs(qupid_mean - uniform_mean) / uniform_mean
    print(f"  Relative difference in mean pseudo-F: {mean_relative_diff:.2%}")

    # ── Plot ──────────────────────────────────────────────────────────────
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # Panel a: empirical frequency of each matching, sorted, vs. uniform
    ax = axes[0]
    sorted_idx = np.argsort(empirical_prob)[::-1]
    rank = np.arange(1, M + 1)
    ax.bar(
        rank,
        empirical_prob[sorted_idx],
        width=1.0,
        color="#4477AA",
        edgecolor="none",
        label="Qupid empirical",
    )
    ax.axhline(
        uniform_prob,
        color="#CC3311",
        linewidth=1.5,
        linestyle="--",
        label=f"Uniform (1/M = {uniform_prob:.4f})",
    )
    ax.set_xlabel("Matching rank (sorted by Qupid frequency)")
    ax.set_ylabel("Sampling probability")
    ax.set_title(
        f"a. Per-matching sampling distribution\n"
        f"M = {M} matchings, N = {N_RUNS:,} runs,\n"
        f"TV = {tv_distance:.3f}, KL = {kl:.3f}",
        fontsize=9,
    )
    ax.legend(frameon=True, fontsize=8)
    ax.grid(True, alpha=0.3)

    # Panel b: pseudo-F distribution under uniform vs. Qupid weights
    ax = axes[1]
    bins = np.linspace(f_per_matching.min() - 0.1, f_per_matching.max() + 0.1, 30)
    # Uniform: weight each matching equally
    ax.hist(
        f_per_matching,
        bins=bins,
        weights=np.full(M, 1.0 / M),
        color="#CC3311",
        alpha=0.5,
        label=f"Uniform (mean = {uniform_mean:.3f}, SD = {uniform_std:.3f})",
    )
    # Qupid: weight by empirical frequency
    ax.hist(
        f_per_matching,
        bins=bins,
        weights=qupid_weights,
        color="#4477AA",
        alpha=0.5,
        label=f"Qupid (mean = {qupid_mean:.3f}, SD = {qupid_std:.3f})",
    )
    ax.set_xlabel("PERMANOVA pseudo-F")
    ax.set_ylabel("Density")
    ax.set_title(
        f"b. pseudo-F distribution\n"
        f"|mean(Qupid) − mean(uniform)| / mean(uniform) = {mean_relative_diff:.2%}",
        fontsize=9,
    )
    ax.legend(frameon=True, fontsize=8)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Qupid's randomized Hopcroft-Karp sampler vs. uniform sampling",
        fontsize=10,
        y=1.02,
    )
    fig.tight_layout()

    for ext in ("png", "pdf"):
        fig.savefig(
            OUT / f"figS3_sampling_characterization.{ext}", dpi=200, bbox_inches="tight"
        )
    plt.close(fig)
    print(f"\nFigure written to {OUT}/figS3_sampling_characterization.{{png,pdf}}")

    # Also write a one-line summary TSV for the manuscript
    summary_path = OUT.parent / "benchmark_real/sampling_characterization_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(
        [
            {
                "M_matchings": M,
                "N_runs": N_RUNS,
                "unique_observed": len(counts),
                "total_variation_distance": tv_distance,
                "chi2_statistic": chi2_stat,
                "chi2_pvalue": chi2_p,
                "kl_divergence": kl,
                "uniform_mean_pseudoF": uniform_mean,
                "uniform_sd_pseudoF": uniform_std,
                "qupid_mean_pseudoF": qupid_mean,
                "qupid_sd_pseudoF": qupid_std,
                "mean_pseudoF_rel_diff": mean_relative_diff,
            }
        ]
    ).to_csv(summary_path, sep="\t", index=False)
    print(f"Summary TSV written to {summary_path}")


if __name__ == "__main__":
    main()
