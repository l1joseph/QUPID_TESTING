"""
sampling_characterization.py
============================

Empirical characterization of Qupid's randomized Hopcroft-Karp sampler
relative to uniform sampling over the space of maximum-cardinality matchings.

We test multiple graph sizes (spanning more than an order of magnitude in the
number of distinct maximum matchings M) and ask whether the sampler's bias
(measured by total-variation distance against uniform) and the downstream
pseudo-F shift remain bounded as the matching space grows.

For each test graph, we:
1. Enumerate every distinct maximum-cardinality 1:1 matching by backtracking.
2. Run Qupid's create_matched_pairs() N times with independent random seeds,
   recording the matching produced on each call.
3. Compare the empirical frequency of each matching to the uniform expectation
   (1/M). Report total-variation distance and KL divergence.
4. Compute PERMANOVA pseudo-F for each matching on a synthetic distance matrix
   and report the marginal distribution under uniform vs. Qupid weights.

The four test graphs span M = 165 to M ≈ 10⁵ matchings, covering more than
two orders of magnitude in matching-space size.

Produces Supplementary Figure 3.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats

from qupid.casematch import CaseMatchOneToMany
from qupid.stats import _fast_permanova


# ---------------------------------------------------------------------------
# Test graphs (designed for tractable enumeration but spanning >2 orders of M)
# ---------------------------------------------------------------------------


def _graph_small() -> dict[str, set[str]]:
    """5 cases × 8 controls; M = 165."""
    return {
        "c0": {"t0", "t1", "t2", "t3"},
        "c1": {"t1", "t3", "t4"},
        "c2": {"t2", "t4", "t5", "t6"},
        "c3": {"t0", "t5", "t7"},
        "c4": {"t3", "t4", "t6", "t7"},
    }


def _graph_medium() -> dict[str, set[str]]:
    """7 cases × 12 controls, denser connectivity."""
    return {
        "c0": {"t0", "t1", "t2", "t3", "t4"},
        "c1": {"t1", "t2", "t5", "t6"},
        "c2": {"t0", "t3", "t6", "t7", "t8"},
        "c3": {"t2", "t4", "t7", "t9"},
        "c4": {"t1", "t5", "t8", "t10"},
        "c5": {"t3", "t6", "t9", "t10", "t11"},
        "c6": {"t0", "t4", "t7", "t11"},
    }


def _graph_large() -> dict[str, set[str]]:
    """9 cases × 16 controls, denser still."""
    return {
        "c0": {"t0", "t1", "t2", "t3", "t4"},
        "c1": {"t1", "t2", "t5", "t6", "t7"},
        "c2": {"t0", "t3", "t6", "t8", "t9"},
        "c3": {"t2", "t4", "t7", "t10", "t11"},
        "c4": {"t1", "t5", "t8", "t12"},
        "c5": {"t3", "t9", "t10", "t13"},
        "c6": {"t6", "t11", "t12", "t14"},
        "c7": {"t4", "t7", "t13", "t14", "t15"},
        "c8": {"t0", "t5", "t8", "t15"},
    }


TEST_GRAPHS = {
    "small (5×8)": _graph_small(),
    "medium (7×12)": _graph_medium(),
    "large (9×16)": _graph_large(),
}


# ---------------------------------------------------------------------------
# Step 1: enumerate all maximum-cardinality matchings
# ---------------------------------------------------------------------------


def enumerate_max_matchings(
    ccm: dict[str, set[str]],
) -> list[frozenset[tuple[str, str]]]:
    """Backtracking enumeration with branch-and-bound pruning."""
    cases = sorted(ccm.keys())
    n_cases = len(cases)
    best_size = [0]
    best: set[frozenset[tuple[str, str]]] = set()

    def recurse(idx: int, used: set[str], current: dict[str, str]) -> None:
        remaining = n_cases - idx
        if len(current) + remaining < best_size[0]:
            return
        if idx == n_cases:
            sz = len(current)
            if sz > best_size[0]:
                best_size[0] = sz
                best.clear()
            if sz == best_size[0]:
                best.add(frozenset(current.items()))
            return
        case = cases[idx]
        for ctrl in sorted(ccm[case] - used):
            current[case] = ctrl
            recurse(idx + 1, used | {ctrl}, current)
            del current[case]
        recurse(idx + 1, used, current)

    recurse(0, set(), {})
    return sorted(best, key=lambda m: sorted(m))


# ---------------------------------------------------------------------------
# Step 2: Qupid sampling
# ---------------------------------------------------------------------------


def qupid_sample(
    ccm: dict[str, set[str]], n_runs: int, seed: int
) -> Counter[frozenset[tuple[str, str]]]:
    cm = CaseMatchOneToMany(ccm)
    counts: Counter[frozenset[tuple[str, str]]] = Counter()
    rng = np.random.default_rng(seed)
    for _ in range(n_runs):
        run_seed = int(rng.integers(0, 2**31 - 1))
        coll = cm.create_matched_pairs(iterations=1, strict=True, seed=run_seed)
        matching = next(iter(coll))
        pairs = frozenset(
            (case, next(iter(ctrls)))
            for case, ctrls in matching.case_control_map.items()
        )
        counts[pairs] += 1
    return counts


# ---------------------------------------------------------------------------
# Step 3: PERMANOVA pseudo-F per matching
# ---------------------------------------------------------------------------


def compute_pseudo_f_per_matching(
    matchings: list[frozenset[tuple[str, str]]],
    dm: np.ndarray,
    sample_ids: list[str],
) -> np.ndarray:
    id_to_idx = {sid: i for i, sid in enumerate(sample_ids)}
    f_values = np.empty(len(matchings))
    for i, matching in enumerate(matchings):
        cases_used = [case for case, _ in matching]
        ctrls_used = [ctrl for _, ctrl in matching]
        ordered = cases_used + ctrls_used
        idx = [id_to_idx[s] for s in ordered]
        sub_dm = dm[np.ix_(idx, idx)].copy()
        is_case = np.array([True] * len(cases_used) + [False] * len(ctrls_used))
        f, _ = _fast_permanova(sub_dm, is_case, 0, np.random.default_rng())
        f_values[i] = f
    return f_values


# ---------------------------------------------------------------------------
# Per-graph analysis
# ---------------------------------------------------------------------------


@dataclass
class GraphResult:
    label: str
    n_cases: int
    n_controls: int
    M: int
    n_runs: int
    counts: Counter
    matchings: list[frozenset[tuple[str, str]]]
    f_per_matching: np.ndarray
    tv_distance: float
    tv_noise_floor: float
    kl_divergence: float
    chi2_stat: float
    chi2_p: float
    uniform_mean_f: float
    uniform_sd_f: float
    qupid_mean_f: float
    qupid_sd_f: float
    mean_f_rel_diff: float
    ks_pseudo_f: float


def analyze_graph(label: str, ccm: dict[str, set[str]], n_runs: int) -> GraphResult:
    print(f"\n=== {label} ===")
    print("  enumerating max matchings...")
    matchings = enumerate_max_matchings(ccm)
    M = len(matchings)
    n_cases = len(ccm)
    n_controls = len({c for ctrls in ccm.values() for c in ctrls})
    max_size = len(next(iter(matchings)))
    print(f"  M = {M}; matching size = {max_size}/{n_cases}")
    assert max_size == n_cases, f"{label}: graph does not admit perfect matching"

    print(f"  running Qupid N = {n_runs:,} times...")
    counts = qupid_sample(ccm, n_runs, seed=42)
    print(f"  Qupid reached {len(counts)}/{M} unique matchings")

    uniform_p = 1.0 / M
    empirical_p = np.array([counts[m] / n_runs for m in matchings])
    tv = 0.5 * float(np.abs(empirical_p - uniform_p).sum())
    # Sampling-noise floor on TV for a true uniform distribution with M bins
    # and N samples: TV ≈ sqrt((M-1)/N)/2 (multinomial standard error scaled).
    tv_noise_floor = 0.5 * float(np.sqrt((M - 1) / n_runs))
    observed = np.array([counts[m] for m in matchings])
    expected = np.full(M, n_runs / M)
    chi2_stat, chi2_p = stats.chisquare(observed, expected)
    eps = 1e-12
    kl = float(np.sum(empirical_p * np.log((empirical_p + eps) / uniform_p)))
    print(
        f"  TV = {tv:.4f} (noise floor {tv_noise_floor:.4f}); "
        f"KL = {kl:.4f}; chi^2 = {chi2_stat:.1f}"
    )

    print("  computing pseudo-F under uniform vs. Qupid weights...")
    sample_ids = sorted(ccm.keys()) + sorted(
        {c for ctrls in ccm.values() for c in ctrls}
    )
    n_total = len(sample_ids)
    rng = np.random.default_rng(42)
    raw = rng.random((n_total, n_total))
    dm = np.triu(raw, 1) + np.triu(raw, 1).T
    f_per_matching = compute_pseudo_f_per_matching(matchings, dm, sample_ids)

    qupid_w = np.array([counts[m] / n_runs for m in matchings])
    uniform_mean = float(f_per_matching.mean())
    uniform_sd = float(f_per_matching.std())
    qupid_mean = float(np.average(f_per_matching, weights=qupid_w))
    qupid_sd = float(
        np.sqrt(np.average((f_per_matching - qupid_mean) ** 2, weights=qupid_w))
    )
    rel_diff = (
        abs(qupid_mean - uniform_mean) / abs(uniform_mean) if uniform_mean else 0.0
    )
    order = np.argsort(f_per_matching)
    uniform_cdf = np.cumsum(np.full(M, 1.0 / M)[order])
    qupid_cdf = np.cumsum(qupid_w[order])
    ks_pseudo_f = float(np.max(np.abs(uniform_cdf - qupid_cdf)))
    print(f"  uniform pseudo-F: mean = {uniform_mean:.3f}, SD = {uniform_sd:.3f}")
    print(f"  Qupid pseudo-F:   mean = {qupid_mean:.3f}, SD = {qupid_sd:.3f}")
    print(
        f"  mean shift = {abs(qupid_mean - uniform_mean):.4f} absolute "
        f"({rel_diff:.2%} relative); KS = {ks_pseudo_f:.4f}"
    )

    return GraphResult(
        label=label,
        n_cases=n_cases,
        n_controls=n_controls,
        M=M,
        n_runs=n_runs,
        counts=counts,
        matchings=matchings,
        f_per_matching=f_per_matching,
        tv_distance=tv,
        tv_noise_floor=tv_noise_floor,
        kl_divergence=kl,
        chi2_stat=chi2_stat,
        chi2_p=chi2_p,
        uniform_mean_f=uniform_mean,
        uniform_sd_f=uniform_sd,
        qupid_mean_f=qupid_mean,
        qupid_sd_f=qupid_sd,
        mean_f_rel_diff=rel_diff,
        ks_pseudo_f=ks_pseudo_f,
    )


# ---------------------------------------------------------------------------
# Plot + summary
# ---------------------------------------------------------------------------


def main() -> None:
    OUT = Path(__file__).parent / "figures"
    OUT.mkdir(exist_ok=True)
    SUMMARY = (
        Path(__file__).parent
        / "benchmark_real"
        / "sampling_characterization_summary.tsv"
    )
    SUMMARY.parent.mkdir(parents=True, exist_ok=True)

    # Use enough runs so the expected hits per matching are comparable across
    # scales (~600 hits per matching for the small graph). This keeps the
    # sampling-noise floor on TV roughly proportional rather than letting
    # finite-N noise dominate the larger graphs.
    runs_by_label = {
        "small (5×8)": 100_000,
        "medium (7×12)": 1_000_000,
        "large (9×16)": 3_000_000,
    }

    results: list[GraphResult] = []
    for label, ccm in TEST_GRAPHS.items():
        n_runs = runs_by_label[label]
        results.append(analyze_graph(label, ccm, n_runs))

    # --- Save summary ---
    pd.DataFrame(
        [
            {
                "label": r.label,
                "n_cases": r.n_cases,
                "n_controls": r.n_controls,
                "M": r.M,
                "N_runs": r.n_runs,
                "unique_observed": len(r.counts),
                "coverage_fraction": len(r.counts) / r.M,
                "TV_distance": r.tv_distance,
                "TV_noise_floor": r.tv_noise_floor,
                "KS_pseudoF": r.ks_pseudo_f,
                "KL_divergence": r.kl_divergence,
                "chi2_statistic": r.chi2_stat,
                "chi2_p": r.chi2_p,
                "uniform_mean_pseudoF": r.uniform_mean_f,
                "uniform_sd_pseudoF": r.uniform_sd_f,
                "qupid_mean_pseudoF": r.qupid_mean_f,
                "qupid_sd_pseudoF": r.qupid_sd_f,
                "mean_pseudoF_rel_diff": r.mean_f_rel_diff,
            }
            for r in results
        ]
    ).to_csv(SUMMARY, sep="\t", index=False)
    print(f"\nSummary written to {SUMMARY}")

    # --- Plot ---
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.5))

    colors = ["#332288", "#117733", "#882255"]

    # Row 1: per-matching frequency distributions
    for ax, r, color in zip(axes[0], results, colors):
        sorted_idx = np.argsort([r.counts[m] / r.n_runs for m in r.matchings])[::-1]
        rank = np.arange(1, r.M + 1)
        sorted_probs = np.array(
            [r.counts[r.matchings[i]] / r.n_runs for i in sorted_idx]
        )
        uniform = 1.0 / r.M
        ax.bar(rank, sorted_probs, width=1.0, color=color, edgecolor="none")
        ax.axhline(uniform, color="#CC3311", linewidth=1.5, linestyle="--")
        ax.set_xlabel("Matching rank")
        ax.set_ylabel("Sampling probability")
        ax.set_title(
            f"{r.label}\n"
            f"M = {r.M:,}; TV = {r.tv_distance:.3f} (noise floor {r.tv_noise_floor:.3f})",
            fontsize=9,
        )
        ax.grid(True, alpha=0.3)

    # Row 2: pseudo-F distributions under uniform vs. Qupid weights
    for ax, r, color in zip(axes[1], results, colors):
        qupid_w = np.array([r.counts[m] / r.n_runs for m in r.matchings])
        bins = np.linspace(
            r.f_per_matching.min() - 0.1, r.f_per_matching.max() + 0.1, 30
        )
        ax.hist(
            r.f_per_matching,
            bins=bins,
            weights=np.full(r.M, 1.0 / r.M),
            color="#CC3311",
            alpha=0.5,
            label=f"Uniform (μ={r.uniform_mean_f:.2f}, σ={r.uniform_sd_f:.2f})",
        )
        ax.hist(
            r.f_per_matching,
            bins=bins,
            weights=qupid_w,
            color=color,
            alpha=0.5,
            label=f"Qupid (μ={r.qupid_mean_f:.2f}, σ={r.qupid_sd_f:.2f})",
        )
        ax.set_xlabel("PERMANOVA pseudo-F")
        ax.set_ylabel("Density")
        ax.set_title(
            f"|Δmean| = {abs(r.qupid_mean_f - r.uniform_mean_f):.3f}; "
            f"KS = {r.ks_pseudo_f:.3f}",
            fontsize=9,
        )
        ax.legend(frameon=True, fontsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Qupid's randomized Hopcroft-Karp sampler vs. uniform sampling across graph scales",
        fontsize=11,
        y=1.01,
    )
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(
            OUT / f"figS3_sampling_characterization.{ext}",
            dpi=200,
            bbox_inches="tight",
        )
    plt.close(fig)
    print(f"Figure written to {OUT}/figS3_sampling_characterization.{{png,pdf}}")


if __name__ == "__main__":
    main()
