# Barnacle2 Re-run Instructions (AGP + THDMI)

These steps re-run the AGP analysis with the peer-review fixes applied:
- **Major 3**: restricted-permutation PERMANOVA (pair-aware p-values for is_case)
- **Major 4**: Benjamini–Hochberg FDR across categories within each iteration

## Prerequisites

- Conda environment `qiime2-shotgun-2024.2` with qupid installed
- AGP data at `/projects/tmi-public-results/22Oct2025/human-gut/WGS/10317/`
- GG2 phylogeny at `/databases/gg/2024.09/2024.09.phylogeny.id.nwk.qza`
- `statsmodels` installed (`pip install statsmodels` if not present in env)

## Steps

```bash
# 1. Connect
ssh l1joseph@barnacle2.ucsd.edu

# 2. Pull the latest notebooks
cd ~/path/to/QUPID_TESTING
git pull origin main

# 3. Activate environment
conda activate qiime2-shotgun-2024.2

# 4. Check statsmodels is available
python -c "from statsmodels.stats.multitest import multipletests; print('statsmodels OK')"
# If not: pip install statsmodels

# 5. Run AGP notebook
cd final_notebooks
jupyter nbconvert --to notebook --execute \
    --ExecutePreprocessor.timeout=7200 \
    AGP_effect_size_controlled.ipynb \
    --output AGP_effect_size_controlled.ipynb

# 6. Copy results back to local machine (run on local):
scp l1joseph@barnacle2.ucsd.edu:~/path/to/QUPID_TESTING/final_notebooks/agp_ccm_analysis_controlled/AGP_cc_pre_all.csv \
    /Users/leojo/Developer/knight_lab/QUPID_TESTING/agp_ccm_analysis_controlled/
scp l1joseph@barnacle2.ucsd.edu:~/path/to/QUPID_TESTING/final_notebooks/agp_ccm_analysis_controlled/AGP_cc_post_all.csv \
    /Users/leojo/Developer/knight_lab/QUPID_TESTING/agp_ccm_analysis_controlled/
scp l1joseph@barnacle2.ucsd.edu:~/path/to/QUPID_TESTING/final_notebooks/agp_ccm_analysis_controlled/AGP_post_ccm.csv \
    /Users/leojo/Developer/knight_lab/QUPID_TESTING/agp_ccm_analysis_controlled/
scp l1joseph@barnacle2.ucsd.edu:~/path/to/QUPID_TESTING/final_notebooks/agp_ccm_analysis_controlled/AGP_pre_bootstrap.csv \
    /Users/leojo/Developer/knight_lab/QUPID_TESTING/agp_ccm_analysis_controlled/
```

## What the updated AGP notebook produces

The notebook now writes these additional columns in the exported CSVs:
- `p_value` — original pooled PERMANOVA p-value (unrestricted permutations)
- `p_value_restricted` — pair-aware p-value (within-pair restricted permutations, is_case only)
- `q_value` — Benjamini–Hochberg FDR q-value (corrected across all categories within each iteration)

The `*_cc_pre_all.csv` and `*_cc_post_all.csv` files are needed for Figure 2 panel a (AGP violin).

## After copying results back

Re-run the figure scripts from `benchmarking/`:
```bash
python figures/fig_agp_hmp2.py
python figures/fig_thdmi.py
```

Then update the manuscript significance counts to match the new FDR-corrected numbers.
