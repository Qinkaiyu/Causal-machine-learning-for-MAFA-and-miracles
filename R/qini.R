###############################################################################
# qini.R - Custom Qini curve implementation
#
# The Qini curve measures the cumulative gain from ITE-based treatment
# targeting compared to uniform (random) allocation.
#
# Algorithm:
#   1. Sort individuals by predicted ITE in ascending order
#   2. For each percentile k:
#      - Compute ATE in the bottom-k subgroup
#      - Qini value = k * ATE_k - k * overall_ATE
#   3. Integrate under the curve (trapezoidal rule) for a scalar score
#
# This implementation is custom and does not rely on external Qini packages.
#
# Also includes:
#   - Bootstrap confidence intervals (200 replications by default)
#   - AIPW-based Qini variant using doubly-robust effect scores
#   - Interaction significance tests (LRT, Wald, Score)
###############################################################################

# --- Core Qini calculation ---

calculate_qini <- function(ite_pred, y_obs, w_obs,
                            overall_ate = NULL, n_points = 100) {
  valid <- !is.na(ite_pred) & !is.na(y_obs) & !is.na(w_obs)
  ite_pred <- as.numeric(ite_pred[valid])
  y_obs    <- as.numeric(y_obs[valid])
  w_obs    <- as.numeric(w_obs[valid])

  n <- length(ite_pred)
  if (n < 2 || length(unique(w_obs)) < 2) {
    return(list(qini = NA_real_, qini_curve = NULL, ate = NA_real_))
  }

  if (is.null(overall_ate)) {
    overall_ate <- mean(y_obs[w_obs == 0]) - mean(y_obs[w_obs == 1])
  }

  # Sort by ITE ascending (lowest predicted benefit first)
  ord <- order(ite_pred, decreasing = FALSE)
  y_sorted <- y_obs[ord]
  w_sorted <- w_obs[ord]

  percentiles <- seq(0, 1, length.out = n_points + 1)
  qini_values <- numeric(n_points + 1)

  for (i in seq_along(percentiles)) {
    k <- floor(n * percentiles[i])
    if (k == 0) {
      qini_values[i] <- 0
      next
    }
    y_k <- y_sorted[1:k]
    w_k <- w_sorted[1:k]
    if (sum(w_k == 1) == 0 || sum(w_k == 0) == 0) {
      qini_values[i] <- NA_real_
      next
    }
    ate_k <- mean(y_k[w_k == 0]) - mean(y_k[w_k == 1])
    qini_values[i] <- k * ate_k - k * overall_ate
  }

  # Interpolate NA values with last non-NA
  for (i in seq_along(qini_values)) {
    if (is.na(qini_values[i])) {
      prev <- qini_values[seq_len(i - 1)]
      qini_values[i] <- if (length(prev) > 0 && any(!is.na(prev))) {
        tail(prev[!is.na(prev)], 1)
      } else {
        0
      }
    }
  }

  # Trapezoidal integration
  qini_area <- 0
  for (i in 2:length(percentiles)) {
    dx <- percentiles[i] - percentiles[i - 1]
    qini_area <- qini_area + dx * (qini_values[i] + qini_values[i - 1]) / 2
  }

  list(
    qini = qini_area / n,
    ate = overall_ate,
    qini_curve = data.frame(percentile = percentiles, qini_value = qini_values)
  )
}

# --- Binary ATE helper ---

calculate_binary_ate <- function(y_obs, w_obs) {
  ok <- !is.na(y_obs) & !is.na(w_obs)
  y_obs <- as.numeric(y_obs[ok])
  w_obs <- as.numeric(w_obs[ok])
  if (length(y_obs) < 2 || sum(w_obs == 0) == 0 || sum(w_obs == 1) == 0) {
    return(NA_real_)
  }
  mean(y_obs[w_obs == 0]) - mean(y_obs[w_obs == 1])
}

# --- AIPW effect scores ---

safe_scale_numeric <- function(x) {
  x <- as.numeric(x)
  s <- stats::sd(x, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) return(rep(0, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

fit_glm_prob <- function(y, x_df, new_df, fallback_mean = NULL) {
  if (is.null(fallback_mean)) fallback_mean <- mean(as.numeric(y), na.rm = TRUE)
  fallback <- pmin(pmax(rep(fallback_mean, nrow(new_df)), 1e-6), 1 - 1e-6)

  ok <- complete.cases(x_df) & !is.na(y)
  x_df <- x_df[ok, , drop = FALSE]
  y <- as.numeric(y[ok])
  if (nrow(x_df) < 10 || length(unique(y)) < 2) return(fallback)

  fit_df <- data.frame(y = y, x_df, check.names = FALSE)
  form <- stats::reformulate(termlabels = colnames(x_df), response = "y")
  fit <- tryCatch(
    suppressWarnings(stats::glm(form, data = fit_df, family = stats::binomial())),
    error = function(e) NULL
  )
  if (is.null(fit)) return(fallback)
  pred <- tryCatch(
    suppressWarnings(stats::predict(fit, newdata = new_df, type = "response")),
    error = function(e) rep(mean(y), nrow(new_df))
  )
  pmin(pmax(as.numeric(pred), 1e-6), 1 - 1e-6)
}

estimate_aipw_effect_scores <- function(X_eval, y_eval, w_eval, trim = 0.02) {
  ok <- complete.cases(X_eval) & !is.na(y_eval) & !is.na(w_eval)
  X_eval <- X_eval[ok, , drop = FALSE]
  y_eval <- as.numeric(y_eval[ok])
  w_eval <- as.numeric(w_eval[ok])

  if (nrow(X_eval) < 30 || length(unique(w_eval)) < 2) {
    return(list(scores = rep(NA_real_, length(y_eval)), ate = NA_real_))
  }

  x_df <- as.data.frame(X_eval)
  colnames(x_df) <- make.names(colnames(x_df), unique = TRUE)
  trim <- min(0.20, max(0.001, trim))

  e_hat <- fit_glm_prob(w_eval, x_df, x_df, fallback_mean = mean(w_eval))
  e_hat <- pmin(pmax(e_hat, trim), 1 - trim)

  mu0_hat <- fit_glm_prob(y_eval[w_eval == 0], x_df[w_eval == 0, , drop = FALSE],
                           x_df, fallback_mean = mean(y_eval[w_eval == 0]))
  mu1_hat <- fit_glm_prob(y_eval[w_eval == 1], x_df[w_eval == 1, , drop = FALSE],
                           x_df, fallback_mean = mean(y_eval[w_eval == 1]))

  scores <- mu0_hat - mu1_hat +
    ((1 - w_eval) * (y_eval - mu0_hat) / (1 - e_hat)) -
    (w_eval * (y_eval - mu1_hat) / e_hat)

  list(scores = as.numeric(scores), ate = mean(scores, na.rm = TRUE), ok = ok)
}

calculate_qini_from_effect_scores <- function(ite_pred, effect_scores,
                                               overall_ate = NULL,
                                               n_points = 100) {
  valid <- !is.na(ite_pred) & !is.na(effect_scores)
  ite_pred <- as.numeric(ite_pred[valid])
  effect_scores <- as.numeric(effect_scores[valid])
  n <- length(ite_pred)
  if (n < 2) return(list(qini = NA_real_, qini_curve = NULL, ate = NA_real_))

  if (is.null(overall_ate)) overall_ate <- mean(effect_scores)

  ord <- order(ite_pred, decreasing = FALSE)
  effect_sorted <- effect_scores[ord]
  percentiles <- seq(0, 1, length.out = n_points + 1)
  qini_values <- numeric(n_points + 1)

  for (i in seq_along(percentiles)) {
    k <- floor(n * percentiles[i])
    if (k == 0) { qini_values[i] <- 0; next }
    effect_k <- mean(effect_sorted[1:k])
    qini_values[i] <- k * effect_k - k * overall_ate
  }

  qini_area <- 0
  for (i in 2:length(percentiles)) {
    dx <- percentiles[i] - percentiles[i - 1]
    qini_area <- qini_area + dx * (qini_values[i] + qini_values[i - 1]) / 2
  }

  list(
    qini = qini_area / n, ate = overall_ate,
    qini_curve = data.frame(percentile = percentiles, qini_value = qini_values)
  )
}

# --- Bootstrap Qini confidence intervals ---

bootstrap_qini_metrics <- function(ite_pred, y_obs, w_obs,
                                    X_eval = NULL, B = 200L,
                                    seed = 2026L, trim = 0.02,
                                    overall_ate = NULL,
                                    include_aipw = FALSE) {
  valid <- !is.na(ite_pred) & !is.na(y_obs) & !is.na(w_obs)
  if (!is.null(X_eval)) valid <- valid & complete.cases(X_eval)

  ite_pred <- as.numeric(ite_pred[valid])
  y_obs    <- as.numeric(y_obs[valid])
  w_obs    <- as.numeric(w_obs[valid])
  X_eval   <- if (is.null(X_eval)) NULL else X_eval[valid, , drop = FALSE]
  n <- length(ite_pred)

  if (n < 30 || length(unique(w_obs)) < 2) return(NULL)
  if (is.null(overall_ate)) overall_ate <- calculate_binary_ate(y_obs, w_obs)

  set.seed(seed)
  draw_rows <- vector("list", B)
  for (b in seq_len(B)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    qini_std <- calculate_qini(ite_pred[idx], y_obs[idx], w_obs[idx],
                                overall_ate)$qini
    qini_aipw <- NA_real_
    if (isTRUE(include_aipw) && !is.null(X_eval)) {
      aipw_fit <- estimate_aipw_effect_scores(X_eval[idx, , drop = FALSE],
                                               y_obs[idx], w_obs[idx], trim)
      qini_aipw <- calculate_qini_from_effect_scores(
        ite_pred[idx], aipw_fit$scores, aipw_fit$ate
      )$qini
    }
    draw_rows[[b]] <- data.frame(
      bootstrap_id = b, qini = qini_std, qini_aipw = qini_aipw,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(draw_rows)
}

summarise_bootstrap_metric <- function(draws, metric_name) {
  if (is.null(draws) || !(metric_name %in% colnames(draws))) {
    return(list(n = 0L, ci_lower = NA_real_, ci_upper = NA_real_,
                p_value = NA_real_))
  }
  vals <- as.numeric(draws[[metric_name]])
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0) {
    return(list(n = 0L, ci_lower = NA_real_, ci_upper = NA_real_,
                p_value = NA_real_))
  }
  list(
    n = length(vals),
    ci_lower = as.numeric(stats::quantile(vals, 0.025, na.rm = TRUE, type = 6)),
    ci_upper = as.numeric(stats::quantile(vals, 0.975, na.rm = TRUE, type = 6)),
    p_value = min(1, 2 * min(mean(vals <= 0), mean(vals >= 0)))
  )
}

# --- Interaction significance tests ---

run_interaction_significance_tests <- function(score, y_obs, w_obs,
                                                time_obs = NULL,
                                                event_obs = NULL) {
  out_rows <- list()
  ok <- !is.na(score) & !is.na(y_obs) & !is.na(w_obs)
  d <- data.frame(
    y = as.numeric(y_obs[ok]),
    w = as.numeric(w_obs[ok]),
    score = safe_scale_numeric(score[ok])
  )

  # Logistic interaction test
  if (nrow(d) >= 30 && length(unique(d$w)) == 2 && length(unique(d$y)) == 2) {
    fit_red <- tryCatch(
      suppressWarnings(stats::glm(y ~ w + score, data = d,
                                   family = stats::binomial())),
      error = function(e) NULL)
    fit_full <- tryCatch(
      suppressWarnings(stats::glm(y ~ w * score, data = d,
                                   family = stats::binomial())),
      error = function(e) NULL)

    lrt_p <- score_p <- wald_p <- NA_real_
    if (!is.null(fit_red) && !is.null(fit_full)) {
      an <- tryCatch(anova(fit_red, fit_full, test = "Chisq"),
                     error = function(e) NULL)
      if (!is.null(an) && nrow(an) >= 2 && "Pr(>Chi)" %in% colnames(an)) {
        lrt_p <- as.numeric(an$`Pr(>Chi)`[2])
      }
      dr <- tryCatch(drop1(fit_full, scope = ~ w:score, test = "Rao"),
                     error = function(e) NULL)
      if (!is.null(dr) && "Pr(>Chi)" %in% colnames(dr) &&
          "w:score" %in% rownames(dr)) {
        score_p <- as.numeric(dr["w:score", "Pr(>Chi)"])
      }
      cf <- tryCatch(summary(fit_full)$coefficients, error = function(e) NULL)
      if (!is.null(cf) && "w:score" %in% rownames(cf) &&
          "Pr(>|z|)" %in% colnames(cf)) {
        wald_p <- as.numeric(cf["w:score", "Pr(>|z|)"])
      }
    }
    out_rows[[length(out_rows) + 1L]] <- data.frame(
      model = "logistic_event_by_horizon",
      test = c("score_test", "likelihood_ratio_test", "wald_test"),
      p_value = c(score_p, lrt_p, wald_p),
      stringsAsFactors = FALSE
    )
  }

  # Cox interaction test
  if (!is.null(time_obs) && !is.null(event_obs)) {
    ok_surv <- ok & !is.na(time_obs) & !is.na(event_obs)
    d_surv <- data.frame(
      time = as.numeric(time_obs[ok_surv]),
      event = as.numeric(event_obs[ok_surv]),
      w = as.numeric(w_obs[ok_surv]),
      score = safe_scale_numeric(score[ok_surv])
    )
    if (nrow(d_surv) >= 30 && length(unique(d_surv$w)) == 2 &&
        length(unique(d_surv$event)) == 2) {
      fit_red <- tryCatch(
        survival::coxph(survival::Surv(time, event) ~ w + score, data = d_surv),
        error = function(e) NULL)
      fit_full <- tryCatch(
        survival::coxph(survival::Surv(time, event) ~ w * score, data = d_surv),
        error = function(e) NULL)

      lrt_p <- wald_p <- NA_real_
      if (!is.null(fit_red) && !is.null(fit_full)) {
        an <- tryCatch(anova(fit_red, fit_full, test = "Chisq"),
                       error = function(e) NULL)
        if (!is.null(an) && nrow(an) >= 2) {
          p_col <- intersect(c("P(>|Chi|)", "Pr(>|Chi|)", "Pr(>Chi)"),
                             colnames(an))
          if (length(p_col) > 0) lrt_p <- as.numeric(an[[p_col[1]]][2])
        }
        cf <- tryCatch(summary(fit_full)$coefficients, error = function(e) NULL)
        if (!is.null(cf) && "w:score" %in% rownames(cf) &&
            "Pr(>|z|)" %in% colnames(cf)) {
          wald_p <- as.numeric(cf["w:score", "Pr(>|z|)"])
        }
      }
      out_rows[[length(out_rows) + 1L]] <- data.frame(
        model = "cox_time_to_event",
        test = c("likelihood_ratio_test", "wald_test"),
        p_value = c(lrt_p, wald_p),
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(out_rows) == 0) return(NULL)
  dplyr::bind_rows(out_rows)
}

# --- Full Qini evaluation for a dataset ---

evaluate_qini_dataset <- function(dataset_name, ite_pred, y_obs, w_obs,
                                   cfg, X_eval = NULL,
                                   time_obs = NULL, event_obs = NULL,
                                   overall_ate = NULL) {
  qini_std <- calculate_qini(ite_pred, y_obs, w_obs, overall_ate)

  # AIPW Qini variant
  aipw_qini <- list(qini = NA_real_, ate = NA_real_, qini_curve = NULL)
  if (isTRUE(cfg$enable_qini_aipw) && !is.null(X_eval)) {
    aipw_fit <- estimate_aipw_effect_scores(X_eval, y_obs, w_obs, trim = 0.02)
    aipw_qini <- calculate_qini_from_effect_scores(
      ite_pred[aipw_fit$ok], aipw_fit$scores, aipw_fit$ate
    )
  }

  # Bootstrap CI
  bootstrap_draws <- NULL
  if (isTRUE(cfg$enable_qini_bootstrap)) {
    bootstrap_draws <- bootstrap_qini_metrics(
      ite_pred, y_obs, w_obs,
      X_eval = if (isTRUE(cfg$enable_qini_aipw)) X_eval else NULL,
      B = cfg$qini_bootstraps, seed = cfg$qini_bootstrap_seed,
      trim = 0.02, overall_ate = qini_std$ate,
      include_aipw = isTRUE(cfg$enable_qini_aipw)
    )
    if (!is.null(bootstrap_draws) && nrow(bootstrap_draws) > 0) {
      bootstrap_draws$dataset <- dataset_name
    }
  }

  boot_std  <- summarise_bootstrap_metric(bootstrap_draws, "qini")
  boot_aipw <- summarise_bootstrap_metric(bootstrap_draws, "qini_aipw")

  summary_row <- data.frame(
    dataset = dataset_name,
    n = sum(!is.na(ite_pred) & !is.na(y_obs) & !is.na(w_obs)),
    qini = qini_std$qini,
    qini_ate = qini_std$ate,
    qini_ci_lower = boot_std$ci_lower,
    qini_ci_upper = boot_std$ci_upper,
    qini_bootstrap_p = boot_std$p_value,
    qini_bootstrap_reps = boot_std$n,
    qini_aipw = aipw_qini$qini,
    qini_aipw_ate = aipw_qini$ate,
    qini_aipw_ci_lower = boot_aipw$ci_lower,
    qini_aipw_ci_upper = boot_aipw$ci_upper,
    qini_aipw_bootstrap_p = boot_aipw$p_value,
    qini_aipw_bootstrap_reps = boot_aipw$n,
    stringsAsFactors = FALSE
  )

  list(
    dataset = dataset_name,
    qini = qini_std$qini,
    qini_curve = qini_std$qini_curve,
    ate = qini_std$ate,
    aipw_qini = aipw_qini$qini,
    aipw_qini_curve = aipw_qini$qini_curve,
    bootstrap_draws = bootstrap_draws,
    summary_row = summary_row
  )
}
