# Supplementary Code

R code for heterogeneous treatment effect (HTE) estimation using three causal machine learning methods:

- **R-Lasso**: R-learner with Lasso regularization (glmnet)
- **R-Boost**: R-learner with gradient boosting (xgboost)
- **UpliftRF**: Uplift Random Forest via transformed outcome (randomForest)

Model discrimination is evaluated using a custom Qini curve implementation, C-for-benefit, and Harrell's C-index.

## Requirements

R >= 4.2.1 with the following packages:

```
glmnet, xgboost, randomForest, grf, dplyr, mice, ggplot2, survival, readxl
```

Optional: `EpiForsk` (for C-for-benefit)

Install all at once:

```r
install.packages(c("glmnet", "xgboost", "randomForest", "grf",
                    "dplyr", "mice", "ggplot2", "survival", "readxl"))
```

## Usage

1. Open `run_analysis.R` and set `data_file` and `external_data_file` to your data paths.
2. Run:

```bash
Rscript run_analysis.R
```

Results (CSV tables and PNG figures) are saved to the `results/` directory.

## File Structure

```
run_analysis.R        Main script (configuration and execution)
R/
  utils.R             Package loading and helper functions
  data_prep.R         Data loading, MICE imputation, train/test split
  r_learner.R         R-learner framework (R-Lasso and R-Boost)
  uplift_rf.R         UpliftRF (transformed outcome random forest)
  qini.R              Custom Qini curve, bootstrap CI, AIPW variant
  evaluate.R          Calibration, C-for-benefit, survival concordance
  plots.R             Qini curves, ITE distributions, KM plots
```

## Key Parameters

Configurable in `run_analysis.R`:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `horizon` | 1.0 | Time horizon in years |
| `age_min` | 65 | Minimum age filter |
| `rboost_nrounds` | 500 | XGBoost iterations |
| `upliftrf_ntree` | 1000 | Number of trees for UpliftRF |
| `use_calibration` | TRUE | Post-hoc calibration |
| `qini_bootstraps` | 200 | Bootstrap replications for CI |
| `seed` | 2026 | Random seed |
