"""
Generate manuscript figures from real qupid_benchmark.py output.

Reads benchmark_results.tsv produced on barnacle2 and writes:
  - fig_runtime_combined.{png,pdf}   -- 3-panel overview figure
  - fig_speedup_vs_iterations.{png,pdf}  -- qupid speedup over naive

Usage (from benchmarking/):
    python generate_benchmark_figures.py [--results benchmark_agp/benchmark_results.tsv]
                                         [--output-dir benchmark_agp]
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PALETTE = {
    "qupid": "#1f77b4",
    "graph_once": "#ff7f0e",
    "naive_loop": "#d62728",
}
LABELS = {
    "qupid": "Qupid (cached graph + Hopcroft-Karp)",
    "graph_once": "Cached graph + greedy assignment",
    "naive_loop": "Rebuild graph every iteration",
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


def plot_combined(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))

    # ── Panel a: runtime vs background size ──────────────────────────────
    ax = axes[0]
    bg = df[df["experiment"] == "background_size"].sort_values("n_controls")
    for method in ["naive_loop", "graph_once", "qupid"]:
        sub = bg[bg["method"] == method]
        if sub.empty:
            continue
        ax.plot(
            sub["n_controls"],
            sub["elapsed_sec"],
            marker="o",
            color=PALETTE[method],
            label=LABELS[method],
            linewidth=2,
            markersize=5,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Background pool size (n_controls)")
    ax.set_ylabel("Wall-clock time (s)")
    k_bg = int(bg[bg["method"] == "qupid"]["n_iterations"].iloc[0])
    n_cases = int(bg[bg["method"] == "qupid"]["n_cases"].iloc[0])
    ax.set_title(
        f"a. Runtime vs. background size\n({n_cases} cases, k = {k_bg} matchings)",
        fontsize=10,
    )
    ax.legend(loc="upper left", frameon=True, fontsize=8.5)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel b: runtime vs iterations ───────────────────────────────────
    ax = axes[1]
    it = df[df["experiment"] == "iterations"].sort_values("n_iterations")
    for method in ["naive_loop", "graph_once", "qupid"]:
        sub = it[it["method"] == method]
        if sub.empty:
            continue
        ax.plot(
            sub["n_iterations"],
            sub["elapsed_sec"],
            marker="o",
            color=PALETTE[method],
            label=LABELS[method],
            linewidth=2,
            markersize=5,
        )

    # annotate speedup at the largest shared k
    shared_k = sorted(
        set(it[it["method"] == "qupid"]["n_iterations"])
        & set(it[it["method"] == "naive_loop"]["n_iterations"])
    )
    if shared_k:
        k_ann = shared_k[-1]
        q_t = it[(it["method"] == "qupid") & (it["n_iterations"] == k_ann)][
            "elapsed_sec"
        ].values[0]
        n_t = it[(it["method"] == "naive_loop") & (it["n_iterations"] == k_ann)][
            "elapsed_sec"
        ].values[0]
        speedup = n_t / q_t
        ax.annotate(
            f"{speedup:.0f}× speedup\nat k = {k_ann}",
            xy=(k_ann, q_t),
            xytext=(k_ann * 0.3, q_t * 5),
            fontsize=9,
            arrowprops=dict(arrowstyle="->", color="gray", lw=1),
        )

    n_bg_it = int(it[it["method"] == "qupid"]["n_controls"].iloc[0])
    n_cases_it = int(it[it["method"] == "qupid"]["n_cases"].iloc[0])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(
        f"b. Runtime vs. number of matchings\n({n_cases_it} cases, {n_bg_it} controls)",
        fontsize=10,
    )
    ax.legend(loc="upper left", frameon=True, fontsize=8.5)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel c: matching correctness vs background size ─────────────────
    ax = axes[2]
    for method in ["qupid", "graph_once", "naive_loop"]:
        sub = bg[bg["method"] == method].copy()
        if sub.empty:
            continue
        sub["success_rate"] = sub["successful_iterations"] / sub["n_iterations"]
        ax.plot(
            sub["n_controls"],
            sub["success_rate"],
            marker="o" if method == "qupid" else "s",
            color=PALETTE[method],
            label=LABELS[method],
            linewidth=2,
            markersize=6,
        )

    ax.set_xscale("log")
    ax.set_xlabel("Background pool size (n_controls)")
    ax.set_ylabel("Fraction of iterations fully matched")
    ax.set_title(
        f"c. Matching success rate\n({n_cases} cases, k = {k_bg})",
        fontsize=10,
    )
    ax.legend(loc="lower right", frameon=True, fontsize=8.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, which="both", alpha=0.3)

    fig.suptitle(
        "Qupid runtime characteristics (AGP IBD dataset, real benchmarks)",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    fig.savefig(out / "fig_runtime_combined.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "fig_runtime_combined.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  fig_runtime_combined.{{png,pdf}}")


def plot_speedup(df: pd.DataFrame, out: Path) -> None:
    it = df[df["experiment"] == "iterations"].sort_values("n_iterations")
    qupid_t = it[it["method"] == "qupid"].set_index("n_iterations")["elapsed_sec"]
    naive_t = it[it["method"] == "naive_loop"].set_index("n_iterations")["elapsed_sec"]

    shared_k = sorted(set(qupid_t.index) & set(naive_t.index))
    if not shared_k:
        print("  (no shared k values for speedup plot — skipped)")
        return

    speedups = [naive_t[k] / qupid_t[k] for k in shared_k]

    fig, ax = plt.subplots(figsize=(6.0, 4.2))
    ax.plot(
        shared_k,
        speedups,
        marker="o",
        color=PALETTE["qupid"],
        linewidth=2,
        markersize=6,
    )
    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.6)

    n_cases = int(it[it["method"] == "qupid"]["n_cases"].iloc[0])
    n_bg = int(it[it["method"] == "qupid"]["n_controls"].iloc[0])
    ax.set_xscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Speedup factor over naive")
    ax.set_title(
        f"Qupid speedup over naive multi-matching\n({n_cases} cases, {n_bg} controls)",
        fontsize=10,
    )
    ax.grid(True, which="both", alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "fig_speedup_vs_iterations.png", dpi=200, bbox_inches="tight")
    fig.savefig(out / "fig_speedup_vs_iterations.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  fig_speedup_vs_iterations.{{png,pdf}}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument(
        "--results",
        default="benchmark_agp/benchmark_results.tsv",
        type=Path,
        help="Path to benchmark_results.tsv",
    )
    p.add_argument(
        "--output-dir",
        default=None,
        type=Path,
        help="Output directory (defaults to same dir as --results)",
    )
    args = p.parse_args()

    out = args.output_dir or args.results.parent
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.results, sep="\t")
    print(f"Loaded {len(df)} rows from {args.results}")
    print(f"Writing figures to {out}/")

    plot_combined(df, out)
    plot_speedup(df, out)
    print("Done.")


if __name__ == "__main__":
    main()
