#!/usr/bin/env bash
# run_real_benchmarks.sh
# ======================
# Run qupid_benchmark.py on HMP2 and THDMI with real metadata.
# Must be run from benchmarking/ with the qiime2-metagenome-2024.10 conda env active.
#
# Usage:
#   conda activate qiime2-metagenome-2024.10
#   cd benchmarking/
#   bash run_real_benchmarks.sh
#
# Outputs:
#   benchmark_real/hmp2/benchmark_results.tsv
#   benchmark_real/thdmi/benchmark_results.tsv

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/benchmark_real"

# Prepare preprocessed TSVs if not already done
if [ ! -f "$DATA_DIR/hmp2_benchmark.tsv" ] || [ ! -f "$DATA_DIR/thdmi_benchmark.tsv" ]; then
    echo "Preparing benchmark data..."
    python "$SCRIPT_DIR/prepare_benchmark_data.py"
fi

# ---------------------------------------------------------------------------
# HMP2  (259 cases, 259 controls)
# ---------------------------------------------------------------------------
echo ""
echo "=== HMP2 benchmark ==="
python "$SCRIPT_DIR/qupid_benchmark.py" \
    --metadata "$DATA_DIR/hmp2_benchmark.tsv" \
    --case-control-column is_ibd \
    --case-value IBD \
    --control-value nonIBD \
    --discrete-cats sex \
    --numeric-cat consent_age:5 \
    --background-sizes 50 100 150 200 259 \
    --iteration-counts 10 25 50 100 250 500 1000 \
    --default-iterations 100 \
    --seed 42 \
    --output-dir "$DATA_DIR/hmp2"

# ---------------------------------------------------------------------------
# THDMI  (1058 cases, 918 controls)
# ---------------------------------------------------------------------------
echo ""
echo "=== THDMI benchmark ==="
python "$SCRIPT_DIR/qupid_benchmark.py" \
    --metadata "$DATA_DIR/thdmi_benchmark.tsv" \
    --case-control-column health_status \
    --case-value unhealthy \
    --control-value healthy \
    --discrete-cats sex bmi_cat thdmi_cohort \
    --numeric-cat host_age:5 \
    --background-sizes 100 250 500 918 \
    --iteration-counts 10 25 50 100 250 500 1000 \
    --default-iterations 100 \
    --seed 42 \
    --output-dir "$DATA_DIR/thdmi"

# ---------------------------------------------------------------------------
# Merge results and generate combined figures
# ---------------------------------------------------------------------------
echo ""
echo "=== Merging results and generating figures ==="
python - <<'EOF'
import pandas as pd
from pathlib import Path

data_dir = Path("benchmark_real")
dfs = []
for dataset, path in [("HMP2", data_dir / "hmp2/benchmark_results.tsv"),
                       ("THDMI", data_dir / "thdmi/benchmark_results.tsv")]:
    if path.exists():
        df = pd.read_csv(path, sep="\t")
        df.insert(0, "dataset", dataset)
        dfs.append(df)

if dfs:
    combined = pd.concat(dfs, ignore_index=True)
    out = data_dir / "benchmark_results_real.tsv"
    combined.to_csv(out, sep="\t", index=False)
    print(f"Combined {len(combined)} rows → {out}")
EOF

python "$SCRIPT_DIR/generate_benchmark_figures.py" \
    --results "$DATA_DIR/benchmark_results_real.tsv" \
    --output-dir "$DATA_DIR"

echo ""
echo "Done. Figures written to $DATA_DIR/"
