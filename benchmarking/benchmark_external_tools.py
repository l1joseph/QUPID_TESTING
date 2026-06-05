"""
benchmark_external_tools.py
============================

Compare qupid against modern case-control matching tools on the AGP IBD cohort.

Tools benchmarked:
  - qupid    : match_by_multiple + create_matched_pairs(k)  [Python]
  - MatchIt  : matchit(..., method='nearest', caliper=...) looped k times  [R/subprocess]
  - miMatch  : if importable; otherwise skipped with a note  [Python]

Additionally loads Lucas Patel's pre-computed timing data (LNP_01) for historical
context (SPSS FUZZY and R Matching on Wisconsin 16S; different dataset, different
hardware — presented on a separate panel).

AGP metadata lives on barnacle2. Set AGP_METADATA env variable or pass --metadata:
  AGP_METADATA=/projects/tmi-public-results/22Oct2025/human-gut/WGS/10317/metadata-by-status/All_good.tsv

Usage:
    python benchmark_external_tools.py \\
        --metadata /path/to/agp_metadata.tsv \\
        --output-dir figures

    # Dry-run with synthetic data (for testing):
    python benchmark_external_tools.py --dry-run --output-dir figures

Outputs:
    fig_runtime_external.{png,pdf}
    benchmark_external_results.tsv
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import tempfile
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

AGP_METADATA_DEFAULT = (
    "/projects/tmi-public-results/22Oct2025/human-gut/WGS/10317/"
    "metadata-by-status/All_good.tsv"
)
AGP_CASE_COL = "ibd"
AGP_CASE_VAL = "Diagnosed by a medical professional (doctor, physician assistant)"
AGP_CTRL_VAL = "I do not have this condition"
AGP_DISCRETE_CATS = ["sex"]
AGP_NUMERIC_TOLS = {"age_years": 10.0}

K_VALUES = [1, 5, 10, 25, 50, 100]
N_REPEAT = 3  # median of N runs

PALETTE = {
    "qupid": "#1f77b4",
    "MatchIt": "#2ca02c",
    "miMatch": "#9467bd",
    "SPSS FUZZY": "#d62728",
    "R Matching": "#ff7f0e",
}
LINESTYLES = {
    "qupid": "-",
    "MatchIt": "--",
    "miMatch": "-.",
    "SPSS FUZZY": ":",
    "R Matching": "--",
}

# Lucas's pre-computed results (LNP_01_Runtime_Analysis.ipynb)
# Dataset: Wisconsin 16S dementia-AD, sex + marsage ±4.5 yr, ~75 cases
LUCAS_DATA = {
    "SPSS FUZZY": {1: None, 10: 11.75, 100: 498.6, 1000: 39654.0},
    "R Matching": {10: 0.112, 100: 1.162, 1000: 11.84, 10000: 113.7},
    "qupid (Lucas)": {10: 0.015, 100: 0.069, 1000: 0.68, 10000: 6.48, 100000: 63.91},
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


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_agp(metadata_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    md = pd.read_table(metadata_path, index_col=0, low_memory=False)
    if str(md.index[0]).startswith("#"):
        md = pd.read_table(metadata_path, index_col=0, skiprows=1, low_memory=False)

    focus = md[md[AGP_CASE_COL] == AGP_CASE_VAL].copy()
    background = md[md[AGP_CASE_COL] == AGP_CTRL_VAL].copy()

    keep_cols = AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())
    focus = focus.dropna(subset=keep_cols)
    background = background.dropna(subset=keep_cols)
    return focus, background


def make_synthetic(
    n_cases: int = 50, n_controls: int = 500
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Small synthetic cohort for --dry-run testing."""
    rng = np.random.default_rng(42)
    sexes = ["male", "female"]
    focus = pd.DataFrame(
        {
            "sex": rng.choice(sexes, n_cases),
            "age_years": rng.integers(20, 70, n_cases).astype(float),
        },
        index=[f"case_{i}" for i in range(n_cases)],
    )
    background = pd.DataFrame(
        {
            "sex": rng.choice(sexes, n_controls),
            "age_years": rng.integers(20, 70, n_controls).astype(float),
        },
        index=[f"ctrl_{i}" for i in range(n_controls)],
    )
    return focus, background


# ---------------------------------------------------------------------------
# qupid timing
# ---------------------------------------------------------------------------


def time_qupid(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k: int,
    n_repeat: int = N_REPEAT,
) -> float:
    from qupid import match_by_multiple

    cats = AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())
    times = []
    for _ in range(n_repeat):
        t0 = time.perf_counter()
        cm = match_by_multiple(
            focus=focus,
            background=background,
            categories=cats,
            tolerance_map=AGP_NUMERIC_TOLS,
            on_failure="continue",
        )
        cm.create_matched_pairs(iterations=k, strict=False)
        times.append(time.perf_counter() - t0)
    return float(np.median(times))


# ---------------------------------------------------------------------------
# MatchIt timing (via R subprocess)
# ---------------------------------------------------------------------------

MATCHIT_R_TEMPLATE = """\
suppressPackageStartupMessages({{
    library(MatchIt)
    library(data.table)
}})

focus    <- fread("{focus_path}", data.table=FALSE)
background <- fread("{background_path}", data.table=FALSE)
focus$group    <- 1L
background$group <- 0L
df <- rbind(focus, background)

t_total <- 0
for (i in seq_len({k})) {{
    t0 <- proc.time()["elapsed"]
    m <- matchit(
        group ~ sex + age_years,
        data       = df,
        method     = "nearest",
        distance   = "mahalanobis",
        exact      = ~ sex,
        caliper    = c(age_years = {age_tol}),
        ratio      = 1,
        replace    = FALSE
    )
    t_total <- t_total + (proc.time()["elapsed"] - t0)
}}
cat(t_total, "\\n")
"""


def time_matchit(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k: int,
    n_repeat: int = N_REPEAT,
) -> float | None:
    if not shutil.which("Rscript"):
        return None

    times = []
    with tempfile.TemporaryDirectory() as tmp:
        fp = Path(tmp) / "focus.csv"
        bp = Path(tmp) / "background.csv"
        focus[AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())].to_csv(fp, index=False)
        background[AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())].to_csv(
            bp, index=False
        )

        script = MATCHIT_R_TEMPLATE.format(
            focus_path=fp,
            background_path=bp,
            k=k,
            age_tol=AGP_NUMERIC_TOLS["age_years"],
        )
        rscript_path = Path(tmp) / "matchit_timing.R"
        rscript_path.write_text(script)

        for _ in range(n_repeat):
            result = subprocess.run(
                ["Rscript", str(rscript_path)],
                capture_output=True,
                text=True,
                timeout=300,
            )
            if result.returncode != 0:
                return None
            try:
                elapsed = float(result.stdout.strip().split()[-1])
                times.append(elapsed)
            except (ValueError, IndexError):
                return None

    return float(np.median(times))


# ---------------------------------------------------------------------------
# miMatch timing (Python, if available)
# ---------------------------------------------------------------------------


def time_mimatch(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k: int,
    n_repeat: int = N_REPEAT,
) -> float | None:
    try:
        import mimatch  # noqa: F401
    except ImportError:
        return None

    # miMatch API varies by version — attempt a generic call
    try:
        from mimatch import match as mm_match

        cats = AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())
        times = []
        for _ in range(n_repeat):
            t0 = time.perf_counter()
            for _ in range(k):
                mm_match(
                    focus=focus,
                    background=background,
                    categories=AGP_DISCRETE_CATS,
                    tolerance_map=AGP_NUMERIC_TOLS,
                )
            times.append(time.perf_counter() - t0)
        return float(np.median(times))
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Results assembly
# ---------------------------------------------------------------------------


def run_benchmarks(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k_values: list[int] = K_VALUES,
) -> pd.DataFrame:
    rows = []

    for k in k_values:
        print(f"  k={k}", end="", flush=True)

        t = time_qupid(focus, background, k)
        rows.append({"tool": "qupid", "k": k, "elapsed_sec": t, "dataset": "AGP"})
        print(f"  qupid={t:.3f}s", end="", flush=True)

        t = time_matchit(focus, background, k)
        if t is not None:
            rows.append({"tool": "MatchIt", "k": k, "elapsed_sec": t, "dataset": "AGP"})
            print(f"  MatchIt={t:.3f}s", end="", flush=True)
        else:
            print("  MatchIt=N/A", end="", flush=True)

        t = time_mimatch(focus, background, k)
        if t is not None:
            rows.append({"tool": "miMatch", "k": k, "elapsed_sec": t, "dataset": "AGP"})
            print(f"  miMatch={t:.3f}s", end="", flush=True)
        else:
            print("  miMatch=N/A", end="", flush=True)

        print()

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------


def plot(agp_df: pd.DataFrame, out: Path) -> None:
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(12.0, 4.5))

    # ── Panel a: new benchmarks on AGP ────────────────────────────────────
    ax = axes[0]
    for tool in agp_df["tool"].unique():
        sub = agp_df[agp_df["tool"] == tool].sort_values("k")
        color = PALETTE.get(tool, "#333333")
        ls = LINESTYLES.get(tool, "-")
        ax.plot(
            sub["k"],
            sub["elapsed_sec"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=tool,
            color=color,
            linestyle=ls,
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(
        "a. Runtime vs. k — modern tools\n(AGP IBD cohort, sex + age ±10 yr)",
        fontsize=10,
    )
    ax.legend(frameon=True, fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    # ── Panel b: Lucas's historical benchmarks ────────────────────────────
    ax = axes[1]
    lucas_palette = {
        "SPSS FUZZY": PALETTE["SPSS FUZZY"],
        "R Matching": PALETTE["R Matching"],
        "qupid (Lucas)": PALETTE["qupid"],
    }
    lucas_ls = {"SPSS FUZZY": ":", "R Matching": "--", "qupid (Lucas)": "-"}
    for tool, data in LUCAS_DATA.items():
        ks = sorted(k for k, v in data.items() if v is not None)
        ts = [data[k] for k in ks]
        ax.plot(
            ks,
            ts,
            marker="o",
            linewidth=2,
            markersize=6,
            label=tool,
            color=lucas_palette.get(tool, "#333333"),
            linestyle=lucas_ls.get(tool, "-"),
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title(
        "b. Runtime vs. k — historical benchmarks\n"
        "(Wisconsin 16S, sex + age ±4.5 yr; Patel et al.)",
        fontsize=10,
    )
    ax.legend(frameon=True, fontsize=9)
    ax.grid(True, which="both", alpha=0.3)

    fig.tight_layout()
    for fmt in ("png", "pdf"):
        fig.savefig(out / f"fig_runtime_external.{fmt}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Figures written to {out}/fig_runtime_external.{{png,pdf}}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--metadata",
        default=os.environ.get("AGP_METADATA", AGP_METADATA_DEFAULT),
        type=Path,
        help="Path to AGP metadata TSV (default: $AGP_METADATA env var or barnacle2 path)",
    )
    p.add_argument("--output-dir", default="figures", type=Path)
    p.add_argument(
        "--results-out",
        default=None,
        type=Path,
        help="Save raw timing data to this TSV",
    )
    p.add_argument(
        "--dry-run",
        action="store_true",
        help="Use synthetic data instead of AGP (for testing)",
    )
    p.add_argument(
        "--k-values",
        nargs="+",
        type=int,
        default=K_VALUES,
    )
    args = p.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    if args.dry_run:
        print("Dry-run mode: using synthetic data (50 cases, 500 controls)")
        focus, background = make_synthetic()
    else:
        if not args.metadata.exists():
            raise FileNotFoundError(
                f"AGP metadata not found: {args.metadata}\n"
                "Run on barnacle2 or set AGP_METADATA env variable."
            )
        focus, background = load_agp(args.metadata)

    print(f"Cases: {len(focus)}  Controls: {len(background)}")
    print(f"k values: {args.k_values}")
    print()

    df = run_benchmarks(focus, background, args.k_values)

    if args.results_out:
        args.results_out.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(args.results_out, sep="\t", index=False)
        print(f"Results written to {args.results_out}")

    plot(df, args.output_dir)


if __name__ == "__main__":
    main()
