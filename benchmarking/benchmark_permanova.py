"""
benchmark_permanova.py
======================

Benchmark the vectorized _fast_permanova implementation against scikit-bio's
permanova across a range of sample sizes and permutation counts.

Produces a two-panel figure:
  - Left:  wall-clock time vs. n_samples for old (skbio) and new (_fast_permanova)
           at 999 permutations
  - Right: speedup factor (skbio / fast) vs. n_samples

Usage (from benchmarking/):
    python benchmark_permanova.py [--output-dir figures]

Requirements: qupid (fork), skbio, numpy, matplotlib, seaborn.
"""

from __future__ import annotations

import argparse
import timeit
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from skbio import DistanceMatrix
from skbio.stats.distance import permanova as skbio_permanova

# Import after qupid is on the path
import sys

sys.path.insert(0, str(Path(__file__).parent.parent / "qupid"))
from qupid.stats import _fast_permanova


N_REPEAT = 10  # timeit repetitions per (n, perms) cell
PERM_COUNTS = [99, 499, 999]
N_SAMPLES = [20, 50, 90, 150, 300, 500]
RNG_SEED = 42

PALETTE = {
    "skbio": "#d62728",
    "fast": "#1f77b4",
}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def make_random_dm(n: int, seed: int) -> tuple[np.ndarray, np.ndarray, DistanceMatrix]:
    rng = np.random.default_rng(seed)
    raw = rng.random((n, n))
    dm_arr = np.triu(raw, 1) + np.triu(raw, 1).T
    ids = [f"s{i}" for i in range(n)]
    is_case = np.array([True] * (n // 2) + [False] * (n - n // 2))
    grouping = pd.Series(["case"] * (n // 2) + ["control"] * (n - n // 2), index=ids)
    dm = DistanceMatrix(dm_arr, ids=ids)
    return dm_arr, is_case, dm, grouping


def time_fn(fn, n_repeat: int = N_REPEAT) -> float:
    """Return median wall-clock time in seconds over n_repeat calls."""
    times = timeit.repeat(fn, number=1, repeat=n_repeat)
    return float(np.median(times))


def run_benchmarks() -> pd.DataFrame:
    rows = []
    for n in N_SAMPLES:
        dm_arr, is_case, dm, grouping = make_random_dm(n, RNG_SEED)
        dm_arr_copy = dm_arr.copy()
        print(f"n={n}", end="", flush=True)

        for perms in PERM_COUNTS:
            rng = np.random.default_rng(RNG_SEED)

            t_fast = time_fn(
                lambda: _fast_permanova(
                    dm_arr_copy.copy(), is_case.copy(), perms, np.random.default_rng()
                )
            )
            t_skbio = time_fn(lambda: skbio_permanova(dm, grouping, permutations=perms))

            rows.append(
                {
                    "n_samples": n,
                    "permutations": perms,
                    "method": "fast",
                    "elapsed_sec": t_fast,
                }
            )
            rows.append(
                {
                    "n_samples": n,
                    "permutations": perms,
                    "method": "skbio",
                    "elapsed_sec": t_skbio,
                }
            )
            speedup = t_skbio / t_fast
            print(
                f"  perms={perms}: fast={t_fast:.4f}s  skbio={t_skbio:.4f}s  speedup={speedup:.1f}x"
            )

        print()

    return pd.DataFrame(rows)


def plot(df: pd.DataFrame, out: Path) -> None:
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.3))

    # ── Panel a: wall-clock time at 999 permutations ──────────────────────
    ax = axes[0]
    sub = df[df["permutations"] == 999]
    for method, label, color in [
        ("skbio", "scikit-bio (current)", PALETTE["skbio"]),
        ("fast", "Vectorized (this work)", PALETTE["fast"]),
    ]:
        s = sub[sub["method"] == method].sort_values("n_samples")
        ax.plot(
            s["n_samples"],
            s["elapsed_sec"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=label,
            color=color,
        )
    ax.set_xlabel("Samples per matching (n)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("a. PERMANOVA runtime\n(999 permutations)", fontsize=10)
    ax.legend(frameon=True, fontsize=9)
    ax.grid(True, alpha=0.3)

    # ── Panel b: speedup factor across permutation counts ─────────────────
    ax = axes[1]
    linestyles = {99: ":", 499: "--", 999: "-"}
    for perms in PERM_COUNTS:
        sub_p = df[df["permutations"] == perms]
        fast = sub_p[sub_p["method"] == "fast"].set_index("n_samples")["elapsed_sec"]
        slow = sub_p[sub_p["method"] == "skbio"].set_index("n_samples")["elapsed_sec"]
        shared_n = sorted(set(fast.index) & set(slow.index))
        speedups = [slow[n] / fast[n] for n in shared_n]
        ax.plot(
            shared_n,
            speedups,
            marker="o",
            linewidth=2,
            markersize=6,
            linestyle=linestyles[perms],
            color=PALETTE["fast"],
            label=f"{perms} permutations",
        )
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xlabel("Samples per matching (n)")
    ax.set_ylabel("Speedup factor (skbio / vectorized)")
    ax.set_title("b. Speedup across permutation counts", fontsize=10)
    ax.legend(frameon=True, fontsize=9)
    ax.grid(True, alpha=0.3)

    fig.suptitle(
        "Vectorized PERMANOVA vs. scikit-bio (median of 10 runs)",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out / "fig_permanova_speedup.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "fig_permanova_speedup.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"Figures written to {out}/fig_permanova_speedup.{{png,pdf}}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--output-dir", default="figures", type=Path, help="Output directory"
    )
    p.add_argument(
        "--results-out",
        default=None,
        type=Path,
        help="Save raw timing data to this TSV path",
    )
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    print("Running PERMANOVA speedup benchmark...")
    print(f"  n_samples: {N_SAMPLES}")
    print(f"  permutations: {PERM_COUNTS}")
    print(f"  repeats per cell: {N_REPEAT}")
    print()

    df = run_benchmarks()

    if args.results_out:
        args.results_out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.results_out, sep="\t", index=False)
        print(f"Timing data written to {args.results_out}")

    plot(df, args.output_dir)


if __name__ == "__main__":
    main()
