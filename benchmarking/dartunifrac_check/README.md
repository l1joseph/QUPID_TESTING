# DartUniFrac reproduction check (HMP2 + THDMI)

Verification that the HMP2 and THDMI cohorts' published PERMANOVA R²
distributions reproduce when the unweighted-UniFrac distance computation is
swapped from scikit-bio's reference implementation to DartUniFrac (Zhao et al.
2026; https://doi.org/10.64898/2026.03.01.708916), a fast sketching-based
algorithm for the same metric.

## Headline results

| cohort | metric | manuscript | DartUniFrac | drift |
|---|---|---|---|---|
| HMP2 | t-statistic (pre vs post R²) | 9.43 | 12.07 | +28% (more significant) |
| HMP2 | p-value | < 1 × 10⁻¹⁵ | 1.67 × 10⁻²⁵ | reproduces |
| HMP2 | SD collapse (pre / post) | 0.30 → 0.04 | 0.272 → 0.036 | ~7.5×, matches "~8-fold" |
| THDMI | relative R² change pre→post | −35% | −38% | reproduces |
| THDMI | # significant iterations (pre / post) | 97 / 67 | 97 / 65 | reproduces |

Both cohorts reproduce robustly. This matters because the AGP cohort did *not*
reproduce when switched from Bray-Curtis to unweighted UniFrac (p = 0.694 vs.
0.006; covariate balance worsened post-CCM). HMP2/THDMI reproducing here
confirms the AGP non-reproducibility is AGP-specific (likely
abundance-vs-phylogenetic confounder structure in that cohort), not a
metric, tree, sketching, or pipeline bug. See `manuscript/manuscript.md`
Methods §"Beta diversity" for the per-cohort metric-selection rationale.

## How to regenerate the distance matrices

Install DartUniFrac on macOS via Homebrew (cargo install fails without
`--features macos-accelerate`):

```bash
brew tap jianshu93/DartUniFrac
brew install DartUniFrac
```

Unweighted UniFrac with DartMinHash sketching:

```bash
dartunifrac -t <gg2_tree.nwk> -b <feature_table.biom> -m dmh -s 2048 -o <out.tsv>
```

Sketch params `-m dmh -s 2048` follow the DartUniFrac README example.

## Files

| file | tracked? | content |
|---|---|---|
| `README.md` | yes | this file |
| `run_ccm_dartunifrac.py` | yes | adapted CCM analysis (DM-source swap from final_notebooks; preserves matching/bootstrap/PERMANOVA/restricted-permutation logic) |
| `{HMP2,THDMI}_pre_ccm.csv` | yes | per-iteration pre-CCM R² (100 rows) |
| `{HMP2,THDMI}_post_ccm.csv` | yes | per-iteration post-CCM R² (100 rows) |
| `{HMP2,THDMI}_summary.csv` | yes | mean/SD/t/p summary stats |
| `{hmp2,thdmi}_dartunifrac.tsv` | **gitignored** | DartUniFrac distance matrices (~26 MB + 77 MB; regenerable from BIOM + tree) |
| `{hmp2,thdmi}_table.biom` | **gitignored** | extracted from QZAs (`*.biom` globally ignored) |

## Source notebooks

- `final_notebooks/HMP2_effect_size_controlled.ipynb`
- `final_notebooks/THDMI_effect_size_controlled.ipynb`

`run_ccm_dartunifrac.py` copies their matching + bootstrap + PERMANOVA logic
verbatim, only swapping the distance-matrix source. Compare results against the
notebook outputs as the regression check.
