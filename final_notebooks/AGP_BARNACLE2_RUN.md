# AGP Re-run at N=999 permutations (barnacle2)

The AGP notebook has been updated to use N=999 permutations (both standard PERMANOVA
and within-pair restricted permutation).  Because the AGP data lives only on barnacle2,
you must run it there and copy the output CSVs back.

## Run on barnacle2

```bash
ssh l1joseph@barnacle2.ucsd.edu

# Activate the same environment used previously
conda activate qiime2-shotgun-2024.2

# Navigate to the repo (adjust path as needed)
cd ~/Developer/knight_lab/QUPID_TESTING/final_notebooks

# Execute the notebook in-place (~30–60 min with 999 permutations)
jupyter nbconvert --to notebook --execute AGP_effect_size_controlled.ipynb \
    --output AGP_effect_size_controlled.ipynb \
    --ExecutePreprocessor.timeout=7200
```

All file paths are read from environment variables with barnacle2 defaults:

| Variable | Default (barnacle2) |
|---|---|
| `AGP_BASE_PATH` | `/projects/tmi-public-results/22Oct2025/human-gut/WGS/10317` |
| `AGP_FEATURE_TABLE` | `$AGP_BASE_PATH/raw.minfeat.mindepth.biom.qza` |
| `AGP_METADATA` | `$AGP_BASE_PATH/metadata-by-status/All_good.tsv` |
| `GG2_PHYLOGENY` | `/databases/gg/2024.09/2024.09.phylogeny.id.nwk.qza` |

No environment variables need to be set if the paths above are correct.

## Copy output CSVs back to your laptop

The notebook writes its output to `final_notebooks/agp_ccm_analysis_controlled/`
(relative to where nbconvert runs).  Copy the CSVs to the project root's
`agp_ccm_analysis_controlled/` directory so that `benchmarking/figures/fig2_agp_hmp2.py`
can find them:

```bash
# Run this on your local machine (after the barnacle2 job finishes)
scp -r l1joseph@barnacle2.ucsd.edu:~/Developer/knight_lab/QUPID_TESTING/final_notebooks/agp_ccm_analysis_controlled/ \
    /Users/leojo/Developer/knight_lab/QUPID_TESTING/agp_ccm_analysis_controlled/
```

Confirm these four files are present locally before re-rendering fig2:

```
agp_ccm_analysis_controlled/AGP_cc_pre_all.csv
agp_ccm_analysis_controlled/AGP_cc_post_all.csv
agp_ccm_analysis_controlled/AGP_pre_bootstrap.csv
agp_ccm_analysis_controlled/AGP_post_ccm.csv
```

## After copying the CSVs

Re-render Figure 2 from the `benchmarking/` directory:

```bash
conda activate qiime2-shotgun-2024.2
cd benchmarking
python figures/fig2_agp_hmp2.py
```

Then update the AGP numbers in `manuscript.md` to match the regenerated CSV values
(mean R², SD, p-value from t-test, and per-iteration significance counts).
