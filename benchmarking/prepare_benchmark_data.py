"""
prepare_benchmark_data.py
=========================

Preprocess HMP2 and THDMI metadata into benchmark-ready TSVs for use with
qupid_benchmark.py.  Each output TSV has a binary case/control column so the
benchmark script can be invoked identically for all datasets.

Outputs (relative to benchmarking/):
  benchmark_real/hmp2_benchmark.tsv
  benchmark_real/thdmi_benchmark.tsv

Usage:
    python prepare_benchmark_data.py
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

BASE = Path(__file__).parent.parent  # QUPID_TESTING/
OUT = Path(__file__).parent / "benchmark_real"
OUT.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# HMP2
# ---------------------------------------------------------------------------

HMP2_PATH = BASE / "hmp2_ccm_analysis/HMP2_matched_metadata.tsv"

hmp2 = pd.read_table(HMP2_PATH, index_col=0, low_memory=False)

# Create binary case/control column: CD or UC → 'IBD'; nonIBD → 'nonIBD'
hmp2["is_ibd"] = hmp2["diagnosis"].apply(
    lambda x: "IBD" if x in ("CD", "UC") else "nonIBD"
)

# Drop rows missing the matching covariates we'll use
hmp2 = hmp2.dropna(subset=["sex", "consent_age"])

out_hmp2 = OUT / "hmp2_benchmark.tsv"
hmp2.to_csv(out_hmp2, sep="\t", index=True)
cases = (hmp2["is_ibd"] == "IBD").sum()
ctrls = (hmp2["is_ibd"] == "nonIBD").sum()
print(f"HMP2: {cases} cases, {ctrls} controls → {out_hmp2}")
print(f"  Matching cols: sex (discrete), consent_age (numeric ±5 yr)")


# ---------------------------------------------------------------------------
# THDMI
# ---------------------------------------------------------------------------

THDMI_PATH = BASE / "thdmi/data/thdmi_metadata_filtered_samples_clean_variables.tsv"

thdmi = pd.read_table(THDMI_PATH, index_col=0, low_memory=False)

# Binary column: True → 'unhealthy'; False → 'healthy'
thdmi["health_status"] = thdmi["unhealthy"].map({True: "unhealthy", False: "healthy"})

thdmi = thdmi.dropna(subset=["sex", "host_age", "bmi_cat", "thdmi_cohort"])

out_thdmi = OUT / "thdmi_benchmark.tsv"
thdmi.to_csv(out_thdmi, sep="\t", index=True)
cases = (thdmi["health_status"] == "unhealthy").sum()
ctrls = (thdmi["health_status"] == "healthy").sum()
print(f"THDMI: {cases} cases, {ctrls} controls → {out_thdmi}")
print(
    f"  Matching cols: sex, bmi_cat, thdmi_cohort (discrete); host_age (numeric ±5 yr)"
)
