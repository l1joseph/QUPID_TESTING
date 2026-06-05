"""
qupid_benchmark.py
==================

Runtime benchmarks for Qupid multi-matching, designed to support the
"Qupid scales to large background pools" section of the manuscript.

Three experiments:
  1. Wall-clock vs. background pool size (cases fixed, controls vary)
  2. Wall-clock vs. number of iterations (graph fixed, k varies)
  3. Wall-clock vs. number of matching covariates (graph density varies)

In each experiment we compare three approaches:
  - qupid       : match_by_multiple + create_matched_pairs(iterations=k)
  - naive_loop  : for each iteration, for each case, randomly pick a valid
                  control without replacement (no graph caching)
  - graph_once  : build the valid-control map once with vanilla pandas, then
                  for each iteration randomly assign one control per case
                  without replacement (caches the valid-control sets but does
                  not solve the bipartite max-cardinality problem properly).
                  This isolates the cost of the Hopcroft-Karp step.

The naive_loop and graph_once baselines are *not* equivalent to Qupid; they
do not guarantee a maximum-cardinality matching and may fail to assign some
cases. They exist purely to demonstrate the algorithmic value of building
the bipartite graph once and applying randomized Hopcroft-Karp on top.

Usage
-----
    python qupid_benchmark.py --metadata /path/to/metadata.tsv \
        --case-control-column ibd \
        --case-value 'Diagnosed by a medical professional' \
        --control-value 'I do not have this condition' \
        --discrete-cats sex age_cat bmi_cat \
        --output-dir benchmark_results

Requirements: qupid, numpy, pandas, matplotlib, seaborn.
"""

from __future__ import annotations

import argparse
import json
import os
import time
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")


# ---------------------------------------------------------------------------
# Baselines for comparison against Qupid
# ---------------------------------------------------------------------------


def _build_valid_control_map(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None = None,
) -> dict[str, list[str]]:
    """Build the valid-control map once, mirroring what Qupid caches.

    For each case, returns the list of control IDs satisfying all matching
    constraints. This is intentionally written in plain pandas/numpy so it can
    be benchmarked outside of Qupid's internals.
    """
    numeric_tolerances = numeric_tolerances or {}

    valid: dict[str, list[str]] = {}
    bg = background.copy()
    for case_id, case_row in focus.iterrows():
        mask = pd.Series(True, index=bg.index)
        for cat in discrete_cats:
            mask &= bg[cat] == case_row[cat]
        for cat, tol in numeric_tolerances.items():
            mask &= (bg[cat] - case_row[cat]).abs() <= tol
        valid[case_id] = bg.index[mask].tolist()
    return valid


def naive_loop_match(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None,
    n_iterations: int,
    rng: np.random.Generator,
) -> int:
    """Recompute the valid-control map every iteration (worst case)."""
    successful = 0
    for _ in range(n_iterations):
        valid = _build_valid_control_map(
            focus, background, discrete_cats, numeric_tolerances
        )
        used: set[str] = set()
        ok = True
        for case_id in focus.index:
            candidates = [c for c in valid[case_id] if c not in used]
            if not candidates:
                ok = False
                break
            chosen = rng.choice(candidates)
            used.add(chosen)
        if ok:
            successful += 1
    return successful


def graph_once_match(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None,
    n_iterations: int,
    rng: np.random.Generator,
) -> tuple[int, dict[str, list[str]]]:
    """Build the valid-control map once, then random-greedy assign per iter.

    Note: this is greedy, not max-cardinality. It will fail to assign some
    cases that Qupid (Hopcroft-Karp) would successfully assign.
    """
    valid = _build_valid_control_map(
        focus, background, discrete_cats, numeric_tolerances
    )
    successful = 0
    for _ in range(n_iterations):
        used: set[str] = set()
        ok = True
        for case_id in focus.index:
            candidates = [c for c in valid[case_id] if c not in used]
            if not candidates:
                ok = False
                break
            chosen = rng.choice(candidates)
            used.add(chosen)
        if ok:
            successful += 1
    return successful, valid


def qupid_match(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None,
    n_iterations: int,
) -> int:
    """Run Qupid's match_by_multiple + create_matched_pairs."""
    from qupid import match_by_multiple

    cm = match_by_multiple(
        focus=focus,
        background=background,
        categories=discrete_cats + list((numeric_tolerances or {}).keys()),
        tolerance_map=numeric_tolerances or {},
        on_failure="continue",
    )
    pairs = cm.create_matched_pairs(iterations=n_iterations, strict=False)
    return len(pairs.case_matches)


# ---------------------------------------------------------------------------
# Experiment runner
# ---------------------------------------------------------------------------


@dataclass
class TimingResult:
    method: str
    n_cases: int
    n_controls: int
    n_iterations: int
    n_covariates: int
    elapsed_sec: float
    successful_iterations: int


def time_block(fn, *args, **kwargs):
    """Median of 3 wall-clock measurements."""
    times = []
    result = None
    for _ in range(3):
        t0 = time.perf_counter()
        result = fn(*args, **kwargs)
        times.append(time.perf_counter() - t0)
    return float(np.median(times)), result


def experiment_vs_background_size(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None,
    background_sizes: list[int],
    n_iterations: int,
    seed: int,
) -> list[TimingResult]:
    rng = np.random.default_rng(seed)
    n_covs = len(discrete_cats) + len(numeric_tolerances or {})
    results: list[TimingResult] = []

    for n_bg in background_sizes:
        if n_bg > len(background):
            continue
        bg_sub = background.sample(n=n_bg, random_state=seed)

        # qupid
        elapsed, n_ok = time_block(
            qupid_match,
            focus,
            bg_sub,
            discrete_cats,
            numeric_tolerances,
            n_iterations,
        )
        results.append(
            TimingResult(
                "qupid",
                len(focus),
                n_bg,
                n_iterations,
                n_covs,
                elapsed,
                n_ok,
            )
        )

        # graph_once
        elapsed, (n_ok, _) = time_block(
            graph_once_match,
            focus,
            bg_sub,
            discrete_cats,
            numeric_tolerances,
            n_iterations,
            rng,
        )
        results.append(
            TimingResult(
                "graph_once",
                len(focus),
                n_bg,
                n_iterations,
                n_covs,
                elapsed,
                n_ok,
            )
        )

        # naive_loop (skip on huge backgrounds; quadratic in k)
        if n_bg <= 5000:
            elapsed, n_ok = time_block(
                naive_loop_match,
                focus,
                bg_sub,
                discrete_cats,
                numeric_tolerances,
                n_iterations,
                rng,
            )
            results.append(
                TimingResult(
                    "naive_loop",
                    len(focus),
                    n_bg,
                    n_iterations,
                    n_covs,
                    elapsed,
                    n_ok,
                )
            )

        print(f"  n_bg={n_bg}: done")

    return results


def experiment_vs_iterations(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    discrete_cats: list[str],
    numeric_tolerances: dict[str, float] | None,
    iteration_counts: list[int],
    seed: int,
) -> list[TimingResult]:
    rng = np.random.default_rng(seed)
    n_covs = len(discrete_cats) + len(numeric_tolerances or {})
    results: list[TimingResult] = []

    for k in iteration_counts:
        elapsed, n_ok = time_block(
            qupid_match,
            focus,
            background,
            discrete_cats,
            numeric_tolerances,
            k,
        )
        results.append(
            TimingResult(
                "qupid",
                len(focus),
                len(background),
                k,
                n_covs,
                elapsed,
                n_ok,
            )
        )

        elapsed, (n_ok, _) = time_block(
            graph_once_match,
            focus,
            background,
            discrete_cats,
            numeric_tolerances,
            k,
            rng,
        )
        results.append(
            TimingResult(
                "graph_once",
                len(focus),
                len(background),
                k,
                n_covs,
                elapsed,
                n_ok,
            )
        )

        if k <= 200:
            elapsed, n_ok = time_block(
                naive_loop_match,
                focus,
                background,
                discrete_cats,
                numeric_tolerances,
                k,
                rng,
            )
            results.append(
                TimingResult(
                    "naive_loop",
                    len(focus),
                    len(background),
                    k,
                    n_covs,
                    elapsed,
                    n_ok,
                )
            )

        print(f"  k={k}: done")

    return results


def experiment_vs_n_covariates(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    candidate_cats: list[str],
    n_iterations: int,
    seed: int,
) -> list[TimingResult]:
    """Add covariates one at a time; graph density falls as more constraints apply."""
    rng = np.random.default_rng(seed)
    results: list[TimingResult] = []

    for k in range(1, len(candidate_cats) + 1):
        cats = candidate_cats[:k]

        elapsed, n_ok = time_block(
            qupid_match,
            focus,
            background,
            cats,
            None,
            n_iterations,
        )
        results.append(
            TimingResult(
                "qupid",
                len(focus),
                len(background),
                n_iterations,
                k,
                elapsed,
                n_ok,
            )
        )

        elapsed, (n_ok, _) = time_block(
            graph_once_match,
            focus,
            background,
            cats,
            None,
            n_iterations,
            rng,
        )
        results.append(
            TimingResult(
                "graph_once",
                len(focus),
                len(background),
                n_iterations,
                k,
                elapsed,
                n_ok,
            )
        )

        elapsed, n_ok = time_block(
            naive_loop_match,
            focus,
            background,
            cats,
            None,
            n_iterations,
            rng,
        )
        results.append(
            TimingResult(
                "naive_loop",
                len(focus),
                len(background),
                n_iterations,
                k,
                elapsed,
                n_ok,
            )
        )

        print(f"  n_cov={k}: done")

    return results


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def make_plots(df: pd.DataFrame, output_dir: Path) -> None:
    import matplotlib.pyplot as plt
    import seaborn as sns

    sns.set_style("whitegrid")
    palette = {"qupid": "#1f77b4", "graph_once": "#ff7f0e", "naive_loop": "#d62728"}

    # Background size
    bg_df = df[df["experiment"] == "background_size"]
    if not bg_df.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for method in ["naive_loop", "graph_once", "qupid"]:
            sub = bg_df[bg_df["method"] == method].sort_values("n_controls")
            if sub.empty:
                continue
            ax.plot(
                sub["n_controls"],
                sub["elapsed_sec"],
                marker="o",
                label=method,
                color=palette[method],
                linewidth=2,
            )
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Background pool size (n_controls)")
        ax.set_ylabel("Wall-clock time (s)")
        ax.set_title(
            "Runtime vs. background pool size\n(cases fixed, k = 100 iterations)"
        )
        ax.legend(frameon=True)
        fig.tight_layout()
        fig.savefig(output_dir / "fig_runtime_vs_background.png", dpi=200)
        fig.savefig(output_dir / "fig_runtime_vs_background.pdf")
        plt.close(fig)

    # Iterations
    it_df = df[df["experiment"] == "iterations"]
    if not it_df.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for method in ["naive_loop", "graph_once", "qupid"]:
            sub = it_df[it_df["method"] == method].sort_values("n_iterations")
            if sub.empty:
                continue
            ax.plot(
                sub["n_iterations"],
                sub["elapsed_sec"],
                marker="o",
                label=method,
                color=palette[method],
                linewidth=2,
            )
        ax.set_xlabel("Number of matchings (k)")
        ax.set_ylabel("Wall-clock time (s)")
        ax.set_title("Runtime vs. number of matchings\n(graph fixed)")
        ax.legend(frameon=True)
        fig.tight_layout()
        fig.savefig(output_dir / "fig_runtime_vs_iterations.png", dpi=200)
        fig.savefig(output_dir / "fig_runtime_vs_iterations.pdf")
        plt.close(fig)

    # Covariates
    cov_df = df[df["experiment"] == "n_covariates"]
    if not cov_df.empty:
        fig, ax = plt.subplots(figsize=(6.5, 4.5))
        for method in ["naive_loop", "graph_once", "qupid"]:
            sub = cov_df[cov_df["method"] == method].sort_values("n_covariates")
            if sub.empty:
                continue
            ax.plot(
                sub["n_covariates"],
                sub["elapsed_sec"],
                marker="o",
                label=method,
                color=palette[method],
                linewidth=2,
            )
        ax.set_xlabel("Number of matching covariates")
        ax.set_ylabel("Wall-clock time (s)")
        ax.set_title("Runtime vs. matching covariates\n(k = 100 iterations)")
        ax.legend(frameon=True)
        fig.tight_layout()
        fig.savefig(output_dir / "fig_runtime_vs_covariates.png", dpi=200)
        fig.savefig(output_dir / "fig_runtime_vs_covariates.pdf")
        plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--metadata", required=True, type=Path)
    p.add_argument("--case-control-column", required=True)
    p.add_argument("--case-value", required=True)
    p.add_argument("--control-value", required=True)
    p.add_argument("--discrete-cats", nargs="+", default=[])
    p.add_argument(
        "--numeric-cat",
        action="append",
        default=[],
        help="Numeric matching: column:tolerance, e.g. age_years:5",
    )
    p.add_argument(
        "--background-sizes",
        nargs="+",
        type=int,
        default=[100, 250, 500, 1000, 2500, 5000, 10000],
    )
    p.add_argument(
        "--iteration-counts", nargs="+", type=int, default=[10, 50, 100, 250, 500, 1000]
    )
    p.add_argument("--default-iterations", type=int, default=100)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--output-dir", required=True, type=Path)
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    md = pd.read_csv(args.metadata, sep="\t", index_col=0)
    if md.index[0].startswith("#"):
        md = pd.read_csv(args.metadata, sep="\t", index_col=0, skiprows=1)

    focus = md[md[args.case_control_column] == args.case_value].dropna(
        subset=args.discrete_cats
    )
    background = md[md[args.case_control_column] == args.control_value].dropna(
        subset=args.discrete_cats
    )

    numeric_tol: dict[str, float] = {}
    for spec in args.numeric_cat:
        col, tol = spec.split(":")
        numeric_tol[col] = float(tol)
        focus = focus.dropna(subset=[col])
        background = background.dropna(subset=[col])

    print(f"Cases: {len(focus)}")
    print(f"Available controls: {len(background)}")
    print(f"Discrete categories: {args.discrete_cats}")
    print(f"Numeric tolerances: {numeric_tol}")
    print()

    all_results: list[tuple[str, TimingResult]] = []

    print("Experiment 1: vs. background pool size")
    res = experiment_vs_background_size(
        focus,
        background,
        args.discrete_cats,
        numeric_tol,
        args.background_sizes,
        args.default_iterations,
        args.seed,
    )
    all_results.extend(("background_size", r) for r in res)

    print("Experiment 2: vs. number of iterations")
    res = experiment_vs_iterations(
        focus,
        background,
        args.discrete_cats,
        numeric_tol,
        args.iteration_counts,
        args.seed,
    )
    all_results.extend(("iterations", r) for r in res)

    if len(args.discrete_cats) >= 2:
        print("Experiment 3: vs. number of covariates")
        res = experiment_vs_n_covariates(
            focus,
            background,
            args.discrete_cats,
            args.default_iterations,
            args.seed,
        )
        all_results.extend(("n_covariates", r) for r in res)

    rows = []
    for exp, r in all_results:
        rows.append({"experiment": exp, **r.__dict__})
    df = pd.DataFrame(rows)
    df.to_csv(args.output_dir / "benchmark_results.tsv", sep="\t", index=False)
    print(f"\nResults written to {args.output_dir / 'benchmark_results.tsv'}")

    make_plots(df, args.output_dir)
    print(f"Figures written to {args.output_dir}")


if __name__ == "__main__":
    main()
