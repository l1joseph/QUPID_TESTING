#!/bin/bash
#SBATCH --partition=short
#SBATCH --job-name=qupid_bench
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=4:00:00
#SBATCH --output=/home/l1joseph/QUPID_TESTING/benchmarking/logs/qupid_benchmark_%j.log
#SBATCH --error=/home/l1joseph/QUPID_TESTING/benchmarking/logs/qupid_benchmark_%j.err

set -euo pipefail

cd /home/l1joseph/QUPID_TESTING/benchmarking

source ~/miniforge3/etc/profile.d/conda.sh
conda activate qiime2-shotgun-2024.2

python -u qupid_benchmark.py \
    --metadata /projects/tmi-public-results/22Oct2025/human-gut/WGS/10317/metadata-by-status/All_good.tsv \
    --case-control-column ibd \
    --case-value "Diagnosed by a medical professional (doctor, physician assistant)" \
    --control-value "I do not have this condition" \
    --discrete-cats sex age_cat bmi_cat \
    --output-dir benchmark_agp
