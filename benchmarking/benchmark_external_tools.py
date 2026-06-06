"""
benchmark_external_tools.py
============================

Compare qupid against modern case-control matching tools on the AGP IBD cohort.

Tools benchmarked (all on the same AGP IBD cohort, categorical matching):
  - qupid       : match_by_multiple + create_matched_pairs(k)  [Python]
  - MatchIt     : matchit(..., method='nearest', exact=...) looped k times  [R/subprocess]
  - R Matching  : Matching::Match(..., exact=...) looped k times  [R/subprocess]
  - CEM         : cem::cem(...) looped k times; set CEM_RSCRIPT  [R/subprocess]
  - miMatch     : if importable; otherwise skipped with a note  [Python]

AGP metadata lives on barnacle. Set AGP_METADATA env variable or pass --metadata:
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
# The public tmi-public-results AGP metadata redacts numeric age_years
# ("not provided"/"not collected"), so matching follows the published
# qupid_agp analysis: categorical sex + age_cat + bmi_cat, no numeric tolerance.
AGP_DISCRETE_CATS = ["sex", "age_cat", "bmi_cat"]
AGP_NUMERIC_TOLS: dict[str, float] = {}
# AGP missingness sentinels to drop from categorical matching columns
AGP_JUNK_VALUES = {
    "not provided",
    "not collected",
    "not applicable",
    "unspecified",
    "nan",
    "",
}

K_VALUES = [1, 5, 10, 25, 50, 100]
N_REPEAT = 3  # median of N runs

PALETTE = {
    "qupid": "#0077BB",
    "MatchIt": "#EE7733",
    "R Matching": "#009988",
    "CEM": "#CC3311",
    "miMatch": "#9467bd",
}
LINESTYLES = {
    "qupid": "-",
    "MatchIt": "--",
    "R Matching": "-.",
    "CEM": ":",
    "miMatch": (0, (1, 1)),
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

    # Drop AGP missingness sentinels (e.g. "not provided") from categorical cols
    def _drop_junk(df: pd.DataFrame) -> pd.DataFrame:
        for col in AGP_DISCRETE_CATS:
            mask = df[col].astype(str).str.strip().str.lower().isin(AGP_JUNK_VALUES)
            df = df[~mask]
        return df

    focus = _drop_junk(focus)
    background = _drop_junk(background)

    # Coerce any numeric matching columns to float (continuous covariates)
    for col in AGP_NUMERIC_TOLS:
        focus[col] = pd.to_numeric(focus[col], errors="coerce")
        background[col] = pd.to_numeric(background[col], errors="coerce")
        focus = focus.dropna(subset=[col])
        background = background.dropna(subset=[col])

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

# MatchIt analogue of qupid's categorical matching: 1:1 nearest matching
# constrained to exact strata on all categorical covariates (sex, age_cat,
# bmi_cat), looped k times. A logistic propensity score is fit for the
# nearest-neighbour ordering within each exact stratum.
MATCHIT_R_TEMPLATE = """\
suppressPackageStartupMessages({{
    library(MatchIt)
    library(data.table)
}})

focus      <- fread("{focus_path}", data.table=FALSE, colClasses="character")
background <- fread("{background_path}", data.table=FALSE, colClasses="character")
focus$group      <- 1L
background$group  <- 0L
df <- rbind(focus, background)

cats <- strsplit("{cats}", ",")[[1]]
for (c in cats) df[[c]] <- as.factor(df[[c]])
form <- as.formula(paste("group ~", paste(cats, collapse=" + ")))

t_total <- 0
for (i in seq_len({k})) {{
    t0 <- proc.time()["elapsed"]
    m <- matchit(
        form,
        data     = df,
        method   = "nearest",
        distance = "glm",
        exact    = cats,
        ratio    = 1,
        replace  = FALSE
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
            cats=",".join(AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())),
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
#
# NOTE: miMatch (Ye et al., Gut Microbes 2024; github.com/dfwlab/miMatch) is
# intentionally NOT benchmarked on this runtime-vs-k panel, and the import
# below is expected to fail (-> N/A) on any standard install. It is a
# fundamentally different tool and a fair side-by-side comparison is not
# possible here:
#   1. It matches on *microbial metabolic background* (PCs of inferred
#      metabolic pathway profiles), not host covariates (sex/age_cat/bmi_cat),
#      so it would require a HUMAnN-style pathway table the AGP cohort here
#      does not carry.
#   2. It produces a *single* deterministic matched cohort, so there is no
#      "k matchings" axis to place it on against qupid/MatchIt.
#   3. It ships as a config.ini-driven `miMatch.py` CLI script, not an
#      importable module with the `match()` API assumed below.
# The manuscript frames miMatch accordingly as an indirect, microbiome-
# intrinsic approach distinct from direct host-covariate matching. This stub
# is kept only so the tool degrades gracefully to N/A rather than erroring.


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
# R Matching (Sekhon `Matching` package) timing — via R subprocess
# ---------------------------------------------------------------------------
#
# Replicates Patel et al.'s perform_matching.R: Matching::Match looped k times
# with reshuffling. The AGP covariates are all categorical, so every covariate
# is integer-encoded and matched exactly (exact=TRUE, no caliper).
RMATCHING_R_TEMPLATE = """\
suppressPackageStartupMessages(library(Matching))

focus      <- read.csv("{focus_path}", colClasses="character")
background <- read.csv("{background_path}", colClasses="character")
focus$group      <- 1L
background$group  <- 0L
df <- rbind(focus, background)

cats <- strsplit("{cats}", ",")[[1]]
Xmat <- sapply(cats, function(cc) as.integer(as.factor(df[[cc]])))
Xmat <- matrix(as.numeric(Xmat), nrow=nrow(df))
Tr <- df$group

t_total <- 0
for (i in seq_len({k})) {{
    ord <- sample(nrow(df))
    t0 <- proc.time()["elapsed"]
    m <- Match(Tr=Tr[ord], X=Xmat[ord, , drop=FALSE], M=1,
               replace=FALSE, ties=FALSE, exact=rep(TRUE, length(cats)))
    t_total <- t_total + (proc.time()["elapsed"] - t0)
}}
cat(t_total, "\\n")
"""


def time_r_matching(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k: int,
    n_repeat: int = N_REPEAT,
) -> float | None:
    if not shutil.which("Rscript"):
        return None

    cats = AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())
    times = []
    with tempfile.TemporaryDirectory() as tmp:
        fp = Path(tmp) / "focus.csv"
        bp = Path(tmp) / "background.csv"
        focus[cats].to_csv(fp, index=False)
        background[cats].to_csv(bp, index=False)

        script = RMATCHING_R_TEMPLATE.format(
            focus_path=fp,
            background_path=bp,
            k=k,
            cats=",".join(cats),
        )
        rscript_path = Path(tmp) / "rmatching_timing.R"
        rscript_path.write_text(script)

        for _ in range(n_repeat):
            result = subprocess.run(
                ["Rscript", str(rscript_path)],
                capture_output=True,
                text=True,
                timeout=600,
            )
            if result.returncode != 0:
                return None
            try:
                times.append(float(result.stdout.strip().split()[-1]))
            except (ValueError, IndexError):
                return None

    return float(np.median(times))


# ---------------------------------------------------------------------------
# CEM (Coarsened Exact Matching) timing — via R subprocess
# ---------------------------------------------------------------------------
#
# Uses the original `cem` package (Iacus, King & Porro), not MatchIt's
# method="cem" — the latter errors on all-categorical covariate sets in
# MatchIt 4.7.x ("cutpoints must be ... for each numeric variable"). The cem
# package only ships a conda build for R<=4.3, so it lives in its own `r-cem`
# env; point CEM_RSCRIPT at that env's Rscript (falls back to PATH "Rscript").
# Factor covariates are matched exactly (no coarsening needed).
CEM_RSCRIPT = os.environ.get("CEM_RSCRIPT", "Rscript")

CEM_R_TEMPLATE = """\
suppressPackageStartupMessages(library(cem))

focus      <- read.csv("{focus_path}", colClasses="character")
background <- read.csv("{background_path}", colClasses="character")
focus$group      <- 1L
background$group  <- 0L
df <- rbind(focus, background)

cats <- strsplit("{cats}", ",")[[1]]
for (c in cats) df[[c]] <- as.factor(df[[c]])

t_total <- 0
for (i in seq_len({k})) {{
    t0 <- proc.time()["elapsed"]
    m <- cem(treatment="group", data=df, drop=NULL, keep.all=FALSE, verbose=0)
    t_total <- t_total + (proc.time()["elapsed"] - t0)
}}
cat(t_total, "\\n")
"""


def time_cem(
    focus: pd.DataFrame,
    background: pd.DataFrame,
    k: int,
    n_repeat: int = N_REPEAT,
) -> float | None:
    if not shutil.which(CEM_RSCRIPT):
        return None

    cats = AGP_DISCRETE_CATS + list(AGP_NUMERIC_TOLS.keys())
    times = []
    with tempfile.TemporaryDirectory() as tmp:
        fp = Path(tmp) / "focus.csv"
        bp = Path(tmp) / "background.csv"
        focus[cats].to_csv(fp, index=False)
        background[cats].to_csv(bp, index=False)

        script = CEM_R_TEMPLATE.format(
            focus_path=fp,
            background_path=bp,
            k=k,
            cats=",".join(cats),
        )
        rscript_path = Path(tmp) / "cem_timing.R"
        rscript_path.write_text(script)

        for _ in range(n_repeat):
            result = subprocess.run(
                [CEM_RSCRIPT, str(rscript_path)],
                capture_output=True,
                text=True,
                timeout=600,
            )
            if result.returncode != 0:
                return None
            try:
                times.append(float(result.stdout.strip().split()[-1]))
            except (ValueError, IndexError):
                return None

    return float(np.median(times))


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

        t = time_r_matching(focus, background, k)
        if t is not None:
            rows.append(
                {"tool": "R Matching", "k": k, "elapsed_sec": t, "dataset": "AGP"}
            )
            print(f"  RMatching={t:.3f}s", end="", flush=True)
        else:
            print("  RMatching=N/A", end="", flush=True)

        t = time_cem(focus, background, k)
        if t is not None:
            rows.append({"tool": "CEM", "k": k, "elapsed_sec": t, "dataset": "AGP"})
            print(f"  CEM={t:.3f}s", end="", flush=True)
        else:
            print("  CEM=N/A", end="", flush=True)

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
    """Quick single-panel runtime figure (qupid vs external tools on AGP).

    The publication figure is produced by figures/fig_runtime_external.py; this
    is a fast sanity-check plot written alongside the benchmark run.
    """
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(1, 1, figsize=(6.0, 4.5))

    for tool in agp_df["tool"].unique():
        sub = agp_df[agp_df["tool"] == tool].sort_values("k")
        ax.plot(
            sub["k"],
            sub["elapsed_sec"],
            marker="o",
            linewidth=2,
            markersize=6,
            label=tool,
            color=PALETTE.get(tool, "#333333"),
            linestyle=LINESTYLES.get(tool, "-"),
        )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Number of matchings (k)")
    ax.set_ylabel("Wall-clock time (s)")
    ax.set_title("Runtime vs. k — modern tools\n(AGP IBD cohort, sex + age_cat + bmi_cat)")
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
