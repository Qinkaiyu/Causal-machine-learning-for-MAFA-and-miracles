###############################################################################
# evaluate.R - Model evaluation and comparison
#
# Evaluation metrics for heterogeneous treatment effect models:
#   - Qini curve (custom implementation, see qini.R)
#   - C-for-benefit (van Klaveren et al., via EpiForsk package)
#   - Harrell's C-index for treatment benefit
#   - Survival benefit concordance (treated/control arms)
#
# Also includes post-hoc calibration for tau estimates.
###############################################################################

# --- Post-hoc calibration ---

fit_posthoc_calibration <- function(tau_raw, time, event, w, horizon_val,
                                     n_bins = 10L, min_group = 20L) {
  y_surv <- as.numeric(!(event == 1 & time <= horizon_val))
  ok <- is.finite(tau_raw) & !is.na(y_surv) & !is.na(w)
  tau_raw <- as.numeric(tau_raw[ok])
  y_surv  <- as.numeric(y_surv[ok])
  w       <- as.numeric(w[ok])

  if (length(tau_raw) < 100 || length(unique(w)) < 2) {
    return(list(a = 0, b = 1, bin_table = NULL, message = "insufficient data"))
  }

  breaks <- unique(as.numeric(quantile(tau_raw, probs = seq(0, 1,
                    length.out = n_bins + 1), na.rm = TRUE, type = 8)))
  if (length(breaks) < 3) {
    return(list(a = 0, b = 1, bin_table = NULL,
                message = "insufficient tau variation"))
  }

  bins <- cut(tau_raw, breaks = breaks, include.lowest = TRUE,
              ordered_result = TRUE)
  rows <- lapply(levels(bins), function(bin_name) {
    idx <- which(bins == bin_name)
    if (length(idx) < (2 * min_group)) return(NULL)
    n_t <- sum(w[idx] == 1)
    n_c <- sum(w[idx] == 0)
    if (n_t < min_group || n_c < min_group) return(NULL)
    data.frame(
      bin = bin_name, n = length(idx),
      n_treated = n_t, n_control = n_c,
      mean_tau_raw = mean(tau_raw[idx]),
      obs_survival_gain = mean(y_surv[idx][w[idx] == 1]) -
        mean(y_surv[idx][w[idx] == 0]),
      stringsAsFactors = FALSE
    )
  })
  bin_table <- do.call(rbind, rows)

  if (is.null(bin_table) || nrow(bin_table) < 2) {
    return(list(a = 0, b = 1, bin_table = bin_table,
                message = "too few valid bins"))
  }

  fit <- lm(obs_survival_gain ~ mean_tau_raw, data = bin_table, weights = n)
  co <- coef(fit)
  list(
    a = ifelse(is.finite(co[1]), as.numeric(co[1]), 0),
    b = ifelse(is.finite(co[2]), as.numeric(co[2]), 1),
    bin_table = bin_table, message = "ok"
  )
}

apply_posthoc_calibration <- function(tau_raw, a, b, clip = TRUE) {
  tau_cal <- as.numeric(a + b * tau_raw)
  if (clip) tau_cal <- pmax(pmin(tau_cal, 1), -1)
  tau_cal
}

calibrate_tau_event <- function(tau_event_train_raw, tau_event_all_raw,
                                 time_train, event_train, w_train,
                                 horizon, n_bins = 10L, min_group = 20L) {
  calib_fit <- fit_posthoc_calibration(
    tau_raw = -as.numeric(tau_event_train_raw),
    time = as.numeric(time_train), event = as.numeric(event_train),
    w = as.numeric(w_train), horizon_val = horizon,
    n_bins = n_bins, min_group = min_group
  )
  tau_survival_cal <- apply_posthoc_calibration(
    -as.numeric(tau_event_all_raw), calib_fit$a, calib_fit$b, clip = TRUE
  )
  list(tau_event_all = -tau_survival_cal, calib_fit = calib_fit)
}

# --- C-for-benefit ---

calculate_c_for_benefit <- function(tau_event, X_eval, Y_eval, W_eval,
                                     n_bootstraps = 50L) {
  if (!requireNamespace("EpiForsk", quietly = TRUE)) return(NULL)
  if (length(unique(W_eval)) < 2 || length(unique(Y_eval)) < 2) return(NULL)

  tmp_cf <- grf::causal_forest(X_eval, Y_eval, W_eval,
                                num.trees = 500, seed = 2026)
  y_hat <- tmp_cf$Y.hat
  w_hat <- tmp_cf$W.hat

  p_0 <- pmin(pmax(y_hat - w_hat * tau_event, 1e-6), 1 - 1e-6)
  p_1 <- pmin(pmax(y_hat + (1 - w_hat) * tau_event, 1e-6), 1 - 1e-6)

  out <- EpiForsk::CForBenefit(
    forest = tmp_cf, Y = Y_eval, W = W_eval,
    X = as.data.frame(X_eval), tau_hat = tau_event,
    p_0 = p_0, p_1 = p_1,
    CI = "bootstrap", n_bootstraps = as.integer(n_bootstraps),
    verbose = FALSE, match_method = "nearest",
    match_distance = "mahalanobis"
  )
  list(
    c_for_benefit = out$c_for_benefit,
    ci_lower = if (!is.null(out$lower_CI)) out$lower_CI else NA_real_,
    ci_upper = if (!is.null(out$upper_CI)) out$upper_CI else NA_real_
  )
}

# --- Survival benefit concordance ---

calculate_survival_benefit_cindex <- function(time, event, w, tau_survival) {
  calc_c <- function(idx, risk_score) {
    if (sum(idx) < 20 || length(unique(event[idx])) < 2) return(NA_real_)
    fit <- tryCatch(
      survival::concordance(survival::Surv(time[idx], event[idx]) ~ risk_score[idx]),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    as.numeric(fit$concordance)
  }
  list(
    cindex_treated = calc_c(w == 1, tau_survival),
    cindex_control = calc_c(w == 0, -tau_survival)
  )
}

calculate_harrell_cindex_benefit <- function(time, event, w, tau_survival) {
  if (length(time) < 20 || length(unique(event)) < 2 ||
      length(unique(w)) < 2) return(NA_real_)
  score <- ifelse(w == 1, tau_survival, -tau_survival)
  fit <- tryCatch(
    survival::concordance(survival::Surv(time, event) ~ score),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NA_real_)
  as.numeric(fit$concordance)
}

# --- Full model evaluation ---

evaluate_model <- function(state, idx, tau_event, method_name,
                            dataset_name, cfb_bootstraps = 200L) {
  cfg <- state$cfg
  df  <- state$df
  if (length(idx) == 0) return(NULL)

  y_obs     <- state$y_horizon[idx]
  w_obs     <- as.numeric(df[[cfg$w_var]][idx])
  x_eval    <- state$X_all[idx, , drop = FALSE]
  time_obs  <- as.numeric(df[[cfg$y_time_var]][idx])
  event_obs <- as.numeric(df[[cfg$y_event_var]][idx])
  tau_survival <- -as.numeric(tau_event)

  # Qini
  qini_res <- evaluate_qini_dataset(
    dataset_name = paste0(method_name, "_", dataset_name),
    ite_pred = as.numeric(tau_event),
    y_obs = y_obs, w_obs = w_obs, cfg = cfg,
    X_eval = x_eval, time_obs = time_obs, event_obs = event_obs
  )

  # Survival metrics
  cidx <- calculate_survival_benefit_cindex(time_obs, event_obs, w_obs,
                                             tau_survival)
  harrell <- calculate_harrell_cindex_benefit(time_obs, event_obs, w_obs,
                                              tau_survival)

  # C-for-benefit
  cfb <- calculate_c_for_benefit(
    tau_event = as.numeric(tau_event), X_eval = x_eval,
    Y_eval = y_obs, W_eval = w_obs, n_bootstraps = cfb_bootstraps
  )

  summary_row <- qini_res$summary_row
  summary_row$method  <- method_name
  summary_row$dataset <- dataset_name
  summary_row$harrell_c <- harrell
  summary_row$benefit_c_treated <- cidx$cindex_treated
  summary_row$benefit_c_control <- cidx$cindex_control
  summary_row$c_for_benefit <- if (!is.null(cfb)) cfb$c_for_benefit else NA_real_
  summary_row$c_for_benefit_ci_lower <- if (!is.null(cfb)) cfb$ci_lower else NA_real_
  summary_row$c_for_benefit_ci_upper <- if (!is.null(cfb)) cfb$ci_upper else NA_real_

  pred_df <- data.frame(
    row_index = idx, method = method_name, dataset = dataset_name,
    tau_event = as.numeric(tau_event),
    tau_survival = tau_survival,
    stringsAsFactors = FALSE
  )

  list(
    summary_row = summary_row,
    predictions = pred_df,
    qini_curve = qini_res$qini_curve,
    qini_bootstrap_draws = qini_res$bootstrap_draws
  )
}

# --- Evaluate all 3 models on all datasets ---

evaluate_all_models <- function(state, r_results, uf_results) {
  cfg <- state$cfg
  df  <- state$df

  # Optionally apply post-hoc calibration
  tau_lasso_all   <- r_results$tau_lasso_all
  tau_boost_all   <- r_results$tau_boost_all
  tau_upliftrf_all <- uf_results$tau_upliftrf_all

  calib_params <- list()
  if (isTRUE(cfg$use_calibration)) {
    time_train  <- as.numeric(df[[cfg$y_time_var]][state$train_idx])
    event_train <- as.numeric(df[[cfg$y_event_var]][state$train_idx])
    w_train     <- as.numeric(df[[cfg$w_var]][state$train_idx])

    log_section("Calibrating R-Lasso")
    cal_lasso <- calibrate_tau_event(
      tau_lasso_all[state$train_idx], tau_lasso_all,
      time_train, event_train, w_train, cfg$horizon,
      cfg$calib_bins, cfg$calib_min_group
    )
    tau_lasso_all <- cal_lasso$tau_event_all
    calib_params$R_lasso <- cal_lasso$calib_fit
    cat(sprintf("  a=%.6f, b=%.6f\n", cal_lasso$calib_fit$a,
                cal_lasso$calib_fit$b))

    log_section("Calibrating R-Boost")
    cal_boost <- calibrate_tau_event(
      tau_boost_all[state$train_idx], tau_boost_all,
      time_train, event_train, w_train, cfg$horizon,
      cfg$calib_bins, cfg$calib_min_group
    )
    tau_boost_all <- cal_boost$tau_event_all
    calib_params$R_boost <- cal_boost$calib_fit
    cat(sprintf("  a=%.6f, b=%.6f\n", cal_boost$calib_fit$a,
                cal_boost$calib_fit$b))

    log_section("Calibrating UpliftRF")
    cal_uplift <- calibrate_tau_event(
      tau_upliftrf_all[state$train_idx], tau_upliftrf_all,
      time_train, event_train, w_train, cfg$horizon,
      cfg$calib_bins, cfg$calib_min_group
    )
    tau_upliftrf_all <- cal_uplift$tau_event_all
    calib_params$UpliftRF <- cal_uplift$calib_fit
    cat(sprintf("  a=%.6f, b=%.6f\n", cal_uplift$calib_fit$a,
                cal_uplift$calib_fit$b))
  }

  # Define evaluation sets
  eval_sets <- list(Test = state$test_idx, External = state$external_idx)

  # Define models
  models <- list(
    list(name = "UpliftRF", tau_all = tau_upliftrf_all),
    list(name = "R_lasso",  tau_all = tau_lasso_all),
    list(name = "R_boost",  tau_all = tau_boost_all)
  )

  summary_rows <- list()
  pred_rows    <- list()
  qini_curves  <- list()
  boot_rows    <- list()

  for (m in models) {
    for (dname in names(eval_sets)) {
      idx <- eval_sets[[dname]]
      if (length(idx) == 0) next

      log_section(sprintf("Evaluating %s on %s (n=%d)", m$name, dname,
                          length(idx)))
      res <- evaluate_model(
        state, idx, m$tau_all[idx], m$name, dname,
        cfb_bootstraps = cfg$cfb_bootstraps
      )
      if (is.null(res)) next

      summary_rows[[length(summary_rows) + 1L]] <- res$summary_row
      pred_rows[[length(pred_rows) + 1L]]       <- res$predictions
      if (!is.null(res$qini_curve)) {
        qini_curves[[paste0(m$name, "_", dname)]] <- res$qini_curve
      }
      if (!is.null(res$qini_bootstrap_draws) &&
          nrow(res$qini_bootstrap_draws) > 0) {
        bd <- res$qini_bootstrap_draws
        bd$method <- m$name
        boot_rows[[length(boot_rows) + 1L]] <- bd
      }
    }
  }

  summary_df <- dplyr::bind_rows(summary_rows)
  pred_df    <- dplyr::bind_rows(pred_rows)
  boot_df    <- if (length(boot_rows) > 0) dplyr::bind_rows(boot_rows) else
    data.frame()

  list(
    summary = summary_df,
    predictions = pred_df,
    qini_curves = qini_curves,
    bootstrap_draws = boot_df,
    calibration_params = calib_params,
    tau_calibrated = list(
      lasso = tau_lasso_all,
      boost = tau_boost_all,
      upliftrf = tau_upliftrf_all
    )
  )
}
