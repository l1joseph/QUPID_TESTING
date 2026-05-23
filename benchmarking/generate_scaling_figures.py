"""
Generate illustrative scaling figures for the Qupid manuscript.

These curves are not real measurements; they reflect the expected asymptotic
behavior of each method:

  - naive_loop  : O(k * n_cases * n_controls) - rebuilds the valid-control
                  graph every iteration.
  - qupid       : O(n_cases * n_controls + k * E * sqrt(V)) - cache + run
                  randomized Hopcroft-Karp on the cached graph each iteration.
  - greedy      : same caching as qupid but uses random-greedy assignment
                  instead of Hopcroft-Karp; faster but does not guarantee a
                  maximum-cardinality matching, so it fails as constraints
                  tighten.

Replace with real benchmark data from qupid_benchmark.py before submitting.
"""

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).parent / "benchmark_agp"
OUT.mkdir(exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 10,
        "axes.linewidth": 1.0,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)

PALETTE = {
    "qupid": "#1f77b4",
    "greedy": "#ff7f0e",
    "naive_loop": "#d62728",
}

LABELS = {
    "qupid": "Qupid (cached graph + Hopcroft-Karp)",
    "greedy": "Cached graph + greedy assignment",
    "naive_loop": "Rebuild graph every iteration",
}


def t_graph_build(n_cases, n_controls):
    return 4e-7 * n_cases * n_controls + 5e-3


def t_hopcroft_karp_per_iter(n_cases, avg_valid):
    E = n_cases * avg_valid
    V = 2 * n_cases
    return 5e-8 * E * np.sqrt(V) + 1e-3


def time_naive_loop(n_cases, n_controls, k, valid_frac):
    avg_valid = max(n_controls * valid_frac, 1)
    return k * (
        t_graph_build(n_cases, n_controls)
        + t_hopcroft_karp_per_iter(n_cases, avg_valid)
    )


def time_qupid(n_cases, n_controls, k, valid_frac):
    avg_valid = max(n_controls * valid_frac, 1)
    return t_graph_build(n_cases, n_controls) + k * t_hopcroft_karp_per_iter(
        n_cases, avg_valid
    )


def greedy_success_rate(n_cases, avg_valid):
    if avg_valid <= 1:
        return 0.0
    ratio = avg_valid / n_cases
    return float(np.clip(1 - np.exp(-1.5 * ratio), 0.0, 1.0))


# ---------------------------------------------------------------------------
# Combined figure
# ---------------------------------------------------------------------------

bg_sizes = np.array([100, 250, 500, 1000, 2500, 5000, 10000, 25000, 50000])
iter_counts = np.array([10, 25, 50, 100, 250, 500, 1000, 2500])
n_cases_default = 200
valid_frac = 0.05

fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.3))

# Panel a
ax = axes[0]
for method, fn in [("naive_loop", time_naive_loop), ("qupid", time_qupid)]:
    times = [fn(n_cases_default, n, 100, valid_frac) for n in bg_sizes]
    ax.plot(
        bg_sizes,
        times,
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
ax.set_title(
    "a. Runtime vs. background size\n(200 cases, k = 100 matchings)", fontsize=10
)
ax.legend(loc="upper left", frameon=True, fontsize=8.5)
ax.grid(True, which="both", alpha=0.3)

# Panel b
ax = axes[1]
for method, fn in [("naive_loop", time_naive_loop), ("qupid", time_qupid)]:
    times = [fn(n_cases_default, 5000, kk, valid_frac) for kk in iter_counts]
    ax.plot(
        iter_counts,
        times,
        marker="o",
        color=PALETTE[method],
        label=LABELS[method],
        linewidth=2,
        markersize=5,
    )

naive_at_100 = time_naive_loop(n_cases_default, 5000, 100, valid_frac)
qupid_at_100 = time_qupid(n_cases_default, 5000, 100, valid_frac)
speedup_100 = naive_at_100 / qupid_at_100
ax.annotate(
    f"{speedup_100:.0f}x speedup\nat k = 100",
    xy=(100, qupid_at_100),
    xytext=(15, 60),
    fontsize=9,
    arrowprops=dict(arrowstyle="->", color="gray", lw=1),
)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of matchings (k)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(
    "b. Runtime vs. number of matchings\n(200 cases, 5000 controls)", fontsize=10
)
ax.legend(loc="upper left", frameon=True, fontsize=8.5)
ax.grid(True, which="both", alpha=0.3)

# Panel c
ax = axes[2]
n_cases_panel_c = 200
avg_valid_range = np.array([10, 25, 50, 100, 200, 400, 800, 1600])
success_qupid = [1.0] * len(avg_valid_range)
success_greedy = [greedy_success_rate(n_cases_panel_c, v) for v in avg_valid_range]

ax.plot(
    avg_valid_range,
    success_qupid,
    marker="o",
    color=PALETTE["qupid"],
    label="Qupid (Hopcroft-Karp)",
    linewidth=2,
    markersize=6,
)
ax.plot(
    avg_valid_range,
    success_greedy,
    marker="s",
    color=PALETTE["greedy"],
    label="Random-greedy assignment",
    linewidth=2,
    markersize=6,
)
ax.axvline(n_cases_panel_c, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax.text(
    n_cases_panel_c * 1.15,
    0.05,
    "n_cases",
    color="gray",
    fontsize=8,
    rotation=90,
    va="bottom",
)

ax.set_xscale("log")
ax.set_xlabel("Avg. valid controls per case")
ax.set_ylabel("Fraction of iterations all cases matched")
ax.set_title("c. Matching correctness\n(200 cases, illustrative)", fontsize=10)
ax.legend(loc="lower right", frameon=True, fontsize=8.5)
ax.set_ylim(-0.05, 1.05)
ax.grid(True, which="both", alpha=0.3)

fig.suptitle(
    "Qupid runtime characteristics (illustrative; replace with real benchmarks)",
    fontsize=11,
    y=1.02,
)
fig.tight_layout()
fig.savefig(OUT / "fig_runtime_combined.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "fig_runtime_combined.pdf", bbox_inches="tight")
plt.close(fig)


# ---------------------------------------------------------------------------
# Standalone panels
# ---------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6.0, 4.2))
for method, fn in [("naive_loop", time_naive_loop), ("qupid", time_qupid)]:
    times = [fn(n_cases_default, n, 100, valid_frac) for n in bg_sizes]
    ax.plot(
        bg_sizes,
        times,
        marker="o",
        color=PALETTE[method],
        label=LABELS[method],
        linewidth=2,
        markersize=6,
    )
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Background pool size (n_controls)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(
    "Runtime vs. background pool size\n" "(200 cases, k = 100, illustrative)",
    fontsize=10,
)
ax.legend(loc="upper left", frameon=True, fontsize=9)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_runtime_vs_background.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "fig_runtime_vs_background.pdf", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(6.0, 4.2))
for method, fn in [("naive_loop", time_naive_loop), ("qupid", time_qupid)]:
    times = [fn(n_cases_default, 5000, kk, valid_frac) for kk in iter_counts]
    ax.plot(
        iter_counts,
        times,
        marker="o",
        color=PALETTE[method],
        label=LABELS[method],
        linewidth=2,
        markersize=6,
    )
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Number of matchings (k)")
ax.set_ylabel("Wall-clock time (s)")
ax.set_title(
    "Runtime vs. number of matchings\n" "(200 cases, 5000 controls, illustrative)",
    fontsize=10,
)
ax.legend(loc="upper left", frameon=True, fontsize=9)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_runtime_vs_iterations.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "fig_runtime_vs_iterations.pdf", bbox_inches="tight")
plt.close(fig)

speedup_qupid = [
    time_naive_loop(n_cases_default, 5000, kk, valid_frac)
    / time_qupid(n_cases_default, 5000, kk, valid_frac)
    for kk in iter_counts
]
fig, ax = plt.subplots(figsize=(6.0, 4.2))
ax.plot(
    iter_counts,
    speedup_qupid,
    marker="o",
    color=PALETTE["qupid"],
    linewidth=2,
    markersize=6,
)
ax.axhline(1.0, color="gray", linestyle="--", linewidth=1, alpha=0.6)
ax.set_xscale("log")
ax.set_xlabel("Number of matchings (k)")
ax.set_ylabel("Speedup factor over naive")
ax.set_title(
    "Qupid speedup over naive multi-matching\n"
    "(200 cases, 5000 controls, illustrative)",
    fontsize=10,
)
ax.grid(True, which="both", alpha=0.3)
fig.tight_layout()
fig.savefig(OUT / "fig_speedup_vs_iterations.png", dpi=200, bbox_inches="tight")
fig.savefig(OUT / "fig_speedup_vs_iterations.pdf", bbox_inches="tight")
plt.close(fig)

print("Figures written to:", OUT)
for f in sorted(OUT.glob("*.png")):
    print(f" - {f.name}")

print("\nKey numbers (illustrative):")
print(f"  Speedup at k=100 (200 cases, 5000 controls): {speedup_100:.0f}x")
print(f"  Qupid wall-clock at k=100: {qupid_at_100:.2f} s")
print(f"  Naive wall-clock at k=100: {naive_at_100:.2f} s")
print(
    f"  AGP (169 cases, 4672 controls), k=100:  "
    f"{time_qupid(169, 4672, 100, 0.05):.2f} s"
)
print(
    f"  HMP2 (259 cases, 413 controls), k=100: "
    f"{time_qupid(259, 413, 100, 0.10):.2f} s"
)
print(
    f"  THDMI (563 cases, 918 controls), k=100: "
    f"{time_qupid(563, 918, 100, 0.05):.2f} s"
)
