# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains validation analyses for Qupid, a Python library for generating and statistically evaluating multiple case-control matchings of microbiome data. The goal is to demonstrate Qupid's utility for publication.

**Qupid library**: Located in `qupid/` subdirectory (has its own CLAUDE.md with detailed API documentation).

**Manuscript**: `manuscript` is a Word document (.docx) containing the paper draft.

## Repository Structure

- `qupid/` - The Qupid library source code
- `final_notebooks/` - Cleaned analysis notebooks ready for publication
  - `AGP_effect_size_controlled.ipynb` - American Gut Project analysis
  - `HMP2_effect_size_controlled.ipynb` - HMP2 (Human Microbiome Project 2) analysis
- `agp_ccm_analysis_controlled/` - AGP analysis outputs (CSVs, PNGs)
- `hmp2_ccm_analysis_controlled_v2/` - HMP2 analysis outputs (CSVs, PNGs)
- Scratch notebooks at root level (older iterations)

## Environment Setup

Use conda environment `qiime2-shotgun-2024.2` (preferred) or `qiime2-metagenome-2024.10`:

```bash
conda activate qiime2-shotgun-2024.2  # preferred, qupid pre-installed

# OR
conda activate qiime2-metagenome-2024.10
pip install -e qupid/  # must install qupid manually in this env
```

## Documentation

Use **context7** MCP server to pull QIIME2 documentation when needed.

## Build & Development Commands

```bash
# Install qupid in development mode
cd qupid && pip install -e .[dev]

# Run all tests
cd qupid && make test_all

# Run standalone tests only (no QIIME 2 dependency)
cd qupid && make standalone_test

# Run QIIME 2 plugin tests (requires QIIME 2 environment)
cd qupid && make q2_test

# Run a single test file
cd qupid && pytest qupid/tests/test_casematch.py -v

# Linting
cd qupid && make stylecheck_all
```

## Key Qupid Concepts

- **Focus**: Case samples to be matched
- **Background**: Pool of potential control samples
- **Tolerance**: For continuous metadata, acceptable difference between case/control values (e.g., `age_years+-10`)
- **Discrete matching**: Exact value match required

## Qupid Workflow

1. Match each case to all valid controls → `CaseMatchOneToMany`
2. Generate multiple one-to-one matchings → `CaseMatchCollection`
3. Evaluate statistical differences across all matchings

## Quick API Reference

```python
import qupid
from qupid.stats import bulk_univariate_test, bulk_permanova

# Generate multiple case-control matchings
matches = qupid.shuffle(
    focus=focus_df,
    background=background_df,
    categories=["sex", "age_years"],
    tolerance_map={"age_years": 10},
    iterations=100
)

# Statistical evaluation
test_results = bulk_univariate_test(casematches=results, values=sample_values, test="t")
```

## Validation Analysis Pattern

The notebooks follow this pattern for comparing effect sizes:
1. Load metadata and feature tables (QIIME 2 artifacts or TSVs)
2. Define focus (cases) and background (controls) populations
3. Run Qupid matching with demographic confounders (sex, age, BMI, etc.)
4. Compare effect sizes pre- and post-matching using bootstrap analysis
5. Generate visualizations (distribution plots, scatter comparisons)

## Output Files Convention

Analysis outputs use this naming pattern:
- `{DATASET}_pre_bootstrap.csv` - Effect sizes before matching
- `{DATASET}_post_ccm.csv` - Effect sizes after case-control matching
- `{DATASET}_comparison_controlled.csv` - Pre/post comparison summary
- `{DATASET}_effect_size_controlled.png` - Effect size visualization
- `{DATASET}_scatter_comparison.png` - Pre vs post scatter plot
