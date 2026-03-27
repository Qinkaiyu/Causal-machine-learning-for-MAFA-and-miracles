###############################################################################
# utils.R - Utility functions and configuration
#
# Helper functions for package management, configuration, and file I/O.
###############################################################################

options(repos = c(CRAN = "https://cran.r-project.org/"))

install_if_missing <- function(packages) {
  missing_pkgs <- packages[!(packages %in% installed.packages()[, "Package"])]
  if (length(missing_pkgs) > 0) {
    cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
    install.packages(missing_pkgs, dependencies = TRUE, quiet = TRUE)
  }
}

load_required_packages <- function() {
  required_pkgs <- c(
    "readxl", "dplyr", "mice", "ggplot2", "survival",
    "glmnet", "xgboost", "randomForest", "grf"
  )
  optional_pkgs <- c("EpiForsk")
  install_if_missing(required_pkgs)

  suppressPackageStartupMessages({
    library(readxl)
    library(dplyr)
    library(mice)
    library(ggplot2)
    library(survival)
    library(glmnet)
    library(xgboost)
    library(randomForest)
    library(grf)
  })

  if (all(optional_pkgs %in% installed.packages()[, "Package"])) {
    suppressPackageStartupMessages(library(EpiForsk))
  } else {
    cat("Optional package EpiForsk not installed. C-for-benefit will be skipped.\n")
  }
}

build_config <- function(
    data_file,
    external_data_file = NA_character_,
    output_dir = "results",
    horizon = 1.0,
    age_min = 65,
    seed = 2026,
    feature_vars = c(
      "age", "male", "current_smoke", "paroxysmal_af", "hypertension",
      "diabetes_mellitus", "ckd", "heart_failure", "cad", "pad",
      "cardiomyopathy", "liver_dysfunction", "prior_is_se"
    ),
    y_time_var = "composite_out_fu",
    y_event_var = "composite_out",
    w_var = "intervention",
    split_col = "cluster_site",
    # R-learner parameters
    rlearner_cf_folds = 5L,
    rlearner_e_trim = 0.02,
    rlearner_min_residual = 1e-3,
    # R-Boost parameters
    rboost_nrounds = 200L,
    rboost_eta = 0.05,
    rboost_max_depth = 3L,
    rboost_subsample = 0.8,
    rboost_colsample = 0.8,
    rboost_min_child_weight = 5,
    # UpliftRF parameters
    upliftrf_ntree = 1000L,
    upliftrf_nodesize = 5L,
    upliftrf_mtry = NA_integer_,
    # Calibration
    use_calibration = TRUE,
    calib_bins = 10L,
    calib_min_group = 20L,
    # Evaluation
    qini_bootstraps = 200L,
    qini_bootstrap_seed = 2026L,
    enable_qini_aipw = TRUE,
    enable_qini_bootstrap = TRUE,
    cfb_bootstraps = 200L
) {
  set.seed(seed)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  list(
    data_file = data_file,
    external_data_file = external_data_file,
    output_dir = output_dir,
    horizon = horizon,
    age_min = age_min,
    apply_age_filter = is.finite(age_min),
    seed = seed,
    feature_vars = feature_vars,
    y_time_var = y_time_var,
    y_event_var = y_event_var,
    w_var = w_var,
    split_col = split_col,
    rlearner_cf_folds = rlearner_cf_folds,
    rlearner_e_trim = rlearner_e_trim,
    rlearner_min_residual = rlearner_min_residual,
    rboost_nrounds = rboost_nrounds,
    rboost_eta = rboost_eta,
    rboost_max_depth = rboost_max_depth,
    rboost_subsample = rboost_subsample,
    rboost_colsample = rboost_colsample,
    rboost_min_child_weight = rboost_min_child_weight,
    upliftrf_ntree = upliftrf_ntree,
    upliftrf_nodesize = upliftrf_nodesize,
    upliftrf_mtry = upliftrf_mtry,
    use_calibration = use_calibration,
    calib_bins = calib_bins,
    calib_min_group = calib_min_group,
    qini_bootstraps = qini_bootstraps,
    qini_bootstrap_seed = qini_bootstrap_seed,
    enable_qini_aipw = enable_qini_aipw,
    enable_qini_bootstrap = enable_qini_bootstrap,
    cfb_bootstraps = cfb_bootstraps
  )
}

output_path <- function(cfg, filename) {
  file.path(cfg$output_dir, filename)
}

write_output_csv <- function(df, filename, cfg, ...) {
  path <- output_path(cfg, filename)
  utils::write.csv(df, path, ...)
  cat("Saved:", path, "\n")
}

save_output_plot <- function(plot, filename, cfg, ...) {
  path <- output_path(cfg, filename)
  ggplot2::ggsave(filename = path, plot = plot, ...)
  cat("Saved:", path, "\n")
}

coerce_binary01 <- function(x, name) {
  ux <- sort(unique(x[!is.na(x)]))
  if (!all(ux %in% c(0, 1))) {
    stop(sprintf("%s must be coded as 0/1. Found values: %s", name, paste(ux, collapse = ", ")))
  }
  as.numeric(x)
}

log_section <- function(title) {
  cat("\n== ", title, " ==\n", sep = "")
}
