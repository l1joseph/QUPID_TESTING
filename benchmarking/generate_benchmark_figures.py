"""
Generate manuscript figures from real qupid_benchmark.py output.

Reads benchmark_results.tsv produced by qupid_benchmark.py (single dataset)
or benchmark_results_real.tsv (multi-dataset, with a leading `dataset` column)
and writes:
  - fig_runtime_combined.{png,pdf}      3-panel overview figure
  - fig_speedup_vs_iterations.{png,pdf} qupid speedup over naive

When multiple datasets are present, each panel shows one line per
(method, dataset) combination using the same color per method and different
linestyles per dataset.

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
DATASET_LINESTYLES = ["solid", "dashed", "dotted", "dashdot"]
MARKERS = {"qupid": "o", "graph_once": "s", "naive_loop": "^"}

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _iter_series(df: pd.DataFrame, experiment: str):
    """Yield (method, dataset, subset_df) tuples for an experiment."""
    sub = df[df["experiment"] == experiment]
    datasets = sorted(sub["dataset"].unique()) if "dataset" in sub.columns else [None]
    for method in ["naive_loop", "graph_once", "qupid"]:
        for i, ds in enumerate(datasets):
            if ds is not None:
                rows = sub[(sub["method"] == method) & (sub["dataset"] == ds)]
            else:
                rows = sub[sub["method"] == method]
            if rows.empty:
                continue
            label = LABELS[method]
            if ds is not None and len(datasets) > 1:
                label = f"{LABELS[method]} ({ds})"
            yield method, ds, i, rows.sort_values(
                "n_controls"
                if experiment == "background_size"
                else "n_iterations" if experiment == "iterations" else "n_covariates"
            ), label


def plot_combined(df: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))

    # ── Panel a: runtime vs background size ──────────────────────────────
    ax = axes[0]
    for method, ds, i, rows, label in _iter_series(df, "background_size"):
        ax.plot(
            rows["n_controls"],
            rows["elapsed_sec"],
            marker=MARKERS.get(method, "o"),
            color=PALETTE[method],
            linestyle=DATASET_LINESTYLES[i % len(DATASET_LINESTYLES)],
            label=label,
            linewidth=2,
            markersize=5,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Background pool size (n_controls)")
    ax.set_ylabel("Wall-clock time (s)")
    bg = df[df["experiment"] == "background_size"]
    k_bg = (
        int(bg[bg["method"] == "qupid"]["n_iterations"].iloc[0])
        if not bg.empty
        else 100
    )
    n_cases = int(bg[bg["method"] == "qupid"]["n_cases"].iloc[0]) if not bg.empty else 0
    ax.set_title(
        f"a. Runtime vs. background size\n({n_cases} cases, k = {k_bg} matchings)",
        fontsize=10,
    )
    ax.legend(loc="upper left", frameon=True, fontsize=7.5)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel b: runtime vs iterations ───────────────────────────────────
    ax = axes[1]
    for method, ds, i, rows, label in _iter_series(df, "iterations"):
        ax.plot(
            rows["n_iterations"],
            rows["elapsed_sec"],
            marker=MARKERS.get(method, "o"),
            color=PALETTE[method],
            linestyle=DATASET_LINESTYLES[i % len(DATASET_LINESTYLES)],
            label=label,
            linewidth=2,
            markersize=5,
        )

    it = df[df["experiment"] == "iterations"]
    # Annotate speedup at the largest shared k for the first dataset
    ds0 = sorted(it["dataset"].unique())[0] if "dataset" in it.columns else None
    it0 = it[it["dataset"] == ds0] if ds0 is not None else it
    shared_k = sorted(
        set(it0[it0["method"] == "qupid"]["n_iterations"])
        & set(it0[it0["method"] == "naive_loop"]["n_iterations"])
    )
    if shared_k:
        k_ann = shared_k[-1]
        q_t = it0[(it0["method"] == "qupid") & (it0["n_iterations"] == k_ann)][
            "elapsed_sec"
        ].values[0]
        n_t = it0[(it0["method"] == "naive_loop") & (it0["n_iterations"] == k_ann)][
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

    n_bg_it = (
        int(it0[it0["method"] == "qupid"]["n_controls"].iloc[0]) if not it0.empty else 0
    )
    n_cases_it = (
        int(it0[it0["method"] == "qupid"]["n_cases"].iloc[0]) if not it0.empty else 0
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(
        f"b. Runtime vs. number of matchings\n({n_cases_it} cases, {n_bg_it} controls)",
        fontsize=10,
    )
    ax.legend(loc="upper left", frameon=True, fontsize=7.5)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel c: matching success rate ────────────────────────────────────
    ax = axes[2]
    for method, ds, i, rows, label in _iter_series(df, "background_size"):
        rows = rows.copy()
        rows["success_rate"] = rows["successful_iterations"] / rows["n_iterations"]
        ax.plot(
            rows["n_controls"],
            rows["success_rate"],
            marker=MARKERS.get(method, "o"),
            color=PALETTE[method],
            linestyle=DATASET_LINESTYLES[i % len(DATASET_LINESTYLES)],
            label=label,
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
    ax.legend(loc="lower right", frameon=True, fontsize=7.5)
    ax.set_ylim(-0.05, 1.05)
    ax.grid(True, which="both", alpha=0.3)

    datasets = sorted(df["dataset"].unique()) if "dataset" in df.columns else []
    ds_str = " + ".join(datasets) if datasets else "real benchmarks"
    fig.suptitle(
        f"Qupid runtime characteristics ({ds_str})",
        fontsize=11,
        y=1.02,
    )
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(out / f"fig_runtime_combined.{fmt}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("  fig_runtime_combined.{png,pdf}")


def plot_speedup(df: pd.DataFrame, out: Path) -> None:
    it = df[df["experiment"] == "iterations"]
    datasets = sorted(it["dataset"].unique()) if "dataset" in it.columns else [None]

    fig, ax = plt.subplots(figsize=(6.0, 4.2))

    for i, ds in enumerate(datasets):
        it_ds = it[it["dataset"] == ds] if ds is not None else it
        qupid_t = it_ds[it_ds["method"] == "qupid"].set_index("n_iterations")[
            "elapsed_sec"
        ]
        naive_t = it_ds[it_ds["method"] == "naive_loop"].set_index("n_iterations")[
            "elapsed_sec"
        ]

        shared_k = sorted(set(qupid_t.index) & set(naive_t.index))
        if not shared_k:
            continue

        speedups = [naive_t[k] / qupid_t[k] for k in shared_k]
        label = ds if ds is not None else "qupid"
        ax.plot(
            shared_k,
            speedups,
            marker="o",
            color=PALETTE["qupid"],
            linestyle=DATASET_LINESTYLES[i % len(DATASET_LINESTYLES)],
            linewidth=2,
            markersize=6,
            label=label,
        )

    ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.6)
    ax.set_xscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Speedup factor over naive")
    ax.set_title("Qupid speedup over naive multi-matching", fontsize=10)
    if len(datasets) > 1:
        ax.legend(frameon=True, fontsize=9)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(
            out / f"fig_speedup_vs_iterations.{fmt}", dpi=200, bbox_inches="tight"
        )
    plt.close(fig)
    print("  fig_speedup_vs_iterations.{png,pdf}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--results",
        default="benchmark_agp/benchmark_results.tsv",
        type=Path,
    )
    p.add_argument("--output-dir", default=None, type=Path)
    args = p.parse_args()

    out = args.output_dir or args.results.parent
    out.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.results, sep="\t")
    if "dataset" not in df.columns:
        df.insert(0, "dataset", args.results.parent.name)
    print(f"Loaded {len(df)} rows from {args.results}")
    print(f"Datasets: {sorted(df['dataset'].unique())}")
    print(f"Writing figures to {out}/")

    plot_combined(df, out)
    plot_speedup(df, out)
    print("Done.")


if __name__ == "__main__":
    main()
