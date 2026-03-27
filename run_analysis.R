###############################################################################
# run_analysis.R - Main analysis script
#
# Supplementary code for: [Paper Title]
#
# This script implements heterogeneous treatment effect (HTE) estimation
# using three causal machine learning methods:
#   1. UpliftRF - Uplift Random Forest (transformed outcome approach)
#   2. R-Lasso  - R-learner with Lasso regularization
#   3. R-Boost  - R-learner with XGBoost
#
# The Qini curve implementation is custom (not from an external package).
#
# Required R packages:
#   glmnet, xgboost, randomForest, grf, EpiForsk, tools4uplift,
#   readxl, dplyr, mice, ggplot2, survival
#
# Usage:
#   1. Set the data_file and external_data_file paths below
#   2. Set the output directory
#   3. Run: source("run_analysis.R")
#
# All results (CSV tables and PNG figures) are saved to the output directory.
###############################################################################

# --- Source all modules ---
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || !nzchar(script_dir)) {
  script_dir <- getwd()
}
r_dir <- file.path(script_dir, "R")
source(file.path(r_dir, "utils.R"))
source(file.path(r_dir, "data_prep.R"))
source(file.path(r_dir, "r_learner.R"))
source(file.path(r_dir, "uplift_rf.R"))
source(file.path(r_dir, "qini.R"))
source(file.path(r_dir, "evaluate.R"))
source(file.path(r_dir, "plots.R"))

# --- Load packages ---
load_required_packages()

# ============================================================================
# CONFIGURATION - Modify these paths for your environment
# ============================================================================
cfg <- build_config(
  # Data files (modify paths as needed)
  data_file          = "",
  external_data_file = "",
  output_dir         = file.path(script_dir, "results"),

  # Analysis parameters
  horizon  = 1.0,       # Time horizon (years) for event indicator
  age_min  = 65,        # Minimum age filter
  seed     = 2026,      # Random seed for reproducibility

  # Feature variables
  feature_vars = c(
    "age", "male", "current_smoke", "paroxysmal_af", "hypertension",
    "diabetes_mellitus", "ckd", "heart_failure", "cad", "pad",
    "cardiomyopathy", "liver_dysfunction", "prior_is_se"
  ),

  # Outcome / treatment variables
  y_time_var  = "composite_out_fu",
  y_event_var = "composite_out",
  w_var       = "intervention",
  split_col   = "cluster_site",

  # R-learner nuisance estimation
  rlearner_cf_folds    = 5L,
  rlearner_e_trim      = 0.02,
  rlearner_min_residual = 1e-3,

  # R-Boost (XGBoost) hyperparameters
  rboost_nrounds          = 500L,
  rboost_eta              = 0.05,
  rboost_max_depth        = 3L,
  rboost_subsample        = 0.8,
  rboost_colsample        = 0.8,
  rboost_min_child_weight = 5,

  # UpliftRF hyperparameters
  upliftrf_ntree    = 1000L,
  upliftrf_nodesize = 5L,
  upliftrf_mtry     = NA_integer_,   # default: sqrt(p)

  # Post-hoc calibration
  use_calibration  = TRUE,
  calib_bins       = 10L,
  calib_min_group  = 20L,

  # Evaluation
  qini_bootstraps      = 200L,
  qini_bootstrap_seed  = 2026L,
  enable_qini_aipw     = TRUE,
  enable_qini_bootstrap = TRUE,
  cfb_bootstraps       = 200L
)

# ============================================================================
# STEP 1: Data Preparation
# ============================================================================
cat("\n================================================================\n")
cat("STEP 1: Data Preparation\n")
cat("================================================================\n")

state <- load_and_prepare_data(cfg)

# ============================================================================
# STEP 2: Model Training
# ============================================================================
cat("\n================================================================\n")
cat("STEP 2: Model Training\n")
cat("================================================================\n")

# Train R-learner models (R-Lasso and R-Boost)
r_results <- train_r_learner(state)

# Train UpliftRF
uf_results <- train_uplift_rf(state)

# Save diagnostics
write_output_csv(r_results$diagnostics, "r_learner_diagnostics.csv",
                 cfg, row.names = FALSE)
write_output_csv(uf_results$diagnostics, "upliftrf_diagnostics.csv",
                 cfg, row.names = FALSE)

# ============================================================================
# STEP 3: Model Evaluation
# ============================================================================
cat("\n================================================================\n")
cat("STEP 3: Model Evaluation\n")
cat("================================================================\n")

results <- evaluate_all_models(state, r_results, uf_results)

# Save results
write_output_csv(results$summary, "all_models_summary.csv",
                 cfg, row.names = FALSE)
write_output_csv(results$predictions, "all_models_predictions.csv",
                 cfg, row.names = FALSE)

if (!is.null(results$bootstrap_draws) && nrow(results$bootstrap_draws) > 0) {
  write_output_csv(results$bootstrap_draws,
                   "all_models_qini_bootstrap_draws.csv",
                   cfg, row.names = FALSE)
}

# Save calibration parameters
if (length(results$calibration_params) > 0) {
  calib_df <- data.frame(
    method = names(results$calibration_params),
    a = sapply(results$calibration_params, `[[`, "a"),
    b = sapply(results$calibration_params, `[[`, "b"),
    message = sapply(results$calibration_params, `[[`, "message"),
    stringsAsFactors = FALSE
  )
  write_output_csv(calib_df, "calibration_params.csv", cfg, row.names = FALSE)
}

# Ranked summary
ranked <- results$summary %>%
  dplyr::group_by(dataset) %>%
  dplyr::arrange(dplyr::desc(qini), .by_group = TRUE) %>%
  dplyr::mutate(qini_rank = dplyr::row_number()) %>%
  dplyr::ungroup()
write_output_csv(ranked, "all_models_qini_ranked.csv", cfg, row.names = FALSE)

# Save Qini curve values for each model/dataset
for (key in names(results$qini_curves)) {
  qc <- results$qini_curves[[key]]
  if (!is.null(qc)) {
    fname <- paste0("qini_curve_", tolower(gsub("[^A-Za-z0-9]", "_", key)),
                    "_values.csv")
    write_output_csv(qc, fname, cfg, row.names = FALSE)
  }
}

# ============================================================================
# STEP 4: Generate Plots
# ============================================================================
cat("\n================================================================\n")
cat("STEP 4: Generating Plots\n")
cat("================================================================\n")

generate_all_plots(results, state, cfg)

# ============================================================================
# Summary
# ============================================================================
cat("\n================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("================================================================\n")
cat("\nModel comparison (sorted by Qini):\n")
print(ranked[, c("method", "dataset", "qini", "qini_ci_lower", "qini_ci_upper",
                  "qini_bootstrap_p", "qini_aipw", "c_for_benefit",
                  "harrell_c", "qini_rank")])
cat("\nAll results saved to:", cfg$output_dir, "\n")
