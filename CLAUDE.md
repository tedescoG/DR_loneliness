# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This R-based causal inference project analyzes COVID-19's impact on loneliness using survey data from Italy. The analysis implements doubly robust estimation with propensity score weighting and outcome modeling to estimate the Average Treatment Effect on the Treated (ATT) of increased remote social contacts on severe loneliness outcomes.

## Development Commands

This is an R research project with no formal build system. Scripts are executed directly in R/RStudio:

```r
# Each script sets its own working directory
source("data_manipulation.R")  # Prepare data
source("IPTW.R")               # Estimate propensity scores
source("analysis.R")           # Run doubly robust estimation

# Or for utility functions only
source("utils.R")
```

## Data Structure

- `data/` - All datasets
  - `pooled_ita_w.dta` - Raw survey data (Stata format)
  - `data_clean.rds` - Processed dataset for analysis
  - `inc_dec.rds`, `inc_mix.rds` - Datasets with propensity score weights
- `results/` - Analysis outputs
  - `weighting/` - Propensity score diagnostics and balance tables
  - `outcome/` - Treatment effect estimates and bootstrap results

## High-Level Architecture

### Core Analysis Pipeline

1. **Data Preparation** (`data_manipulation.R`)
   - Imports raw survey data and constructs 20 analysis variables
   - Creates multinomial treatment variable (`remote_contact`: increase/mix/decrease)
   - Defines binary outcome (`severe_loneliness`)
   - Filters to ages 50+, removes students, complete cases only
   - Output: `data_clean.rds` (~1,700 observations)

2. **Utility Functions** (`utils.R` - 996 lines)
   Core statistical infrastructure:
   - **Balance metrics**: `std.diff()`, `es.mean()`, `compute_KS()`, `ks.max()`
   - **Hyperparameter tuning**: `tune.gbm()` - parallel grid search for optimal GBM parameters
   - **Doubly robust estimators**:
     - `aipw_att()` - Augmented IPW (Moodie 2018)
     - `drs_att()` - Doubly Robust Standardization (g-computation)
   - **Bootstrap inference**: `DR_att()` - main function with parallel processing
     - Supports stratified/unstratified bootstrap
     - Returns estimates, SEs, CIs, p-values, diagnostic plots

3. **Propensity Score Weighting** (`IPTW.R`)
   - Uses gradient boosting (GBM) via `twang` package
   - Hyperparameter grid search: 128 combinations (n.trees, depth, shrinkage, bag.fraction, minobs)
   - Optimizes covariate balance (not prediction accuracy)
   - Two binary comparisons: increase vs decrease, increase vs mix
   - Achieves excellent balance: max ASD < 0.05, all p-values > 0.5

4. **Outcome Modeling** (`analysis.R`)
   - Runs 8 models: 4 outcome specifications × 2 comparisons
   - Model complexity ranges from no covariates to fully saturated
   - 500 bootstrap replications per model
   - Both AIPW and DRS estimators for robustness
   - Key finding: 3-5% reduction in severe loneliness from increasing remote contacts

### Key Statistical Methods

- **Causal Estimand**: Average Treatment Effect on Treated (ATT)
- **Identification**: Doubly robust estimation combining propensity scores with outcome regression
- **Propensity Scores**: Gradient Boosting Machine optimizing balance metrics
- **Inference**: Stratified nonparametric bootstrap with parallel processing
- **Balance Assessment**: Absolute standardized differences (ASD) and Kolmogorov-Smirnov statistics

### Treatment and Outcome Variables

- **Treatment**: Change in non-physical contact patterns during COVID
  - Binary: increase (1) vs decrease/unchanged (0)
  - Multinomial: increase vs mix vs decrease
- **Outcome**: Severe loneliness (felt lonely "very often" during COVID)
- **Covariates**: 18 baseline variables including demographics, health, socioeconomic status, pre-COVID mental health

## Dependencies

Primary R packages (loaded in `utils.R`):
- `twang` - Propensity score estimation with GBM
- `gbm` - Gradient boosting implementation
- `survey` - Survey-weighted regression
- `marginaleffects` - Average marginal effects and predictions
- `tidyverse` - Data manipulation and visualization
- `haven` - Stata file import
- `parallel` - Parallel processing for bootstrap

## Common Workflow

1. **Initial data preparation**: Run `data_manipulation.R` once to create cleaned dataset
2. **Propensity score estimation**: Run `IPTW.R` to estimate weights and assess balance
3. **Causal estimation**: Run `analysis.R` for final doubly robust ATT estimates
4. **Results inspection**: Check `results/outcome/bootstrap_DR_summary.csv` for main findings

## Important Notes

- Scripts use `rstudioapi::getActiveDocumentContext()$path` to set working directory
- Random seed fixed at 123 for reproducibility
- Parallel processing auto-detects available cores
- Some scripts contain hardcoded absolute paths that may need updating
- Bootstrap analysis is computationally intensive (~30-60 minutes with full 500 replications)