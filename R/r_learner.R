###############################################################################
# r_learner.R - R-learner framework for heterogeneous treatment effects
#
# Implements the R-learner (Nie & Wager, 2021) with two second-stage models:
#   - R-Lasso: Lasso regression via glmnet (alpha=1)
#   - R-Boost: Gradient boosting via xgboost
#
# The R-learner estimates the conditional average treatment effect (CATE)
# tau(x) by minimizing a loss involving cross-fitted nuisance parameters:
#   - e(x): propensity score (probability of treatment)
#   - m(x): marginal outcome model
#
# Pseudo-outcome: z_i = (Y_i - m(X_i)) / (W_i - e(X_i))
# Weights: w_i = (W_i - e(X_i))^2
###############################################################################

# --- Internal helpers ---

clip_probs <- function(x, lower = 1e-6, upper = 1 - 1e-6) {
  pmin(pmax(as.numeric(x), lower), upper)
}

make_stratified_folds <- function(y, w, k = 5L, seed = 2026L) {
  n <- length(y)
  strata <- paste0(as.integer(w), "_", as.integer(y))
  fold_id <- rep(NA_integer_, n)
  set.seed(seed)
  for (s in unique(strata)) {
    idx <- which(strata == s)
    idx <- sample(idx, length(idx), replace = FALSE)
    fold_id[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold_id
}

fit_binomial_glmnet <- function(x_train, y_train) {
  y_train <- as.numeric(y_train)
  if (nrow(x_train) < 30 || length(unique(y_train)) < 2) return(NULL)
  tryCatch(
    glmnet::cv.glmnet(x = x_train, y = y_train, family = "binomial",
                       alpha = 1, nfolds = 5),
    error = function(e) NULL
  )
}

predict_binomial_glmnet <- function(model, newx, fallback_prob) {
  if (is.null(model)) return(rep(fallback_prob, nrow(newx)))
  pred <- tryCatch(
    as.numeric(predict(model, newx = newx, s = "lambda.min", type = "response")),
    error = function(e) rep(fallback_prob, nrow(newx))
  )
  clip_probs(pred)
}

# --- Cross-fitted nuisance estimation ---

crossfit_nuisance <- function(X_train, y_train, w_train, k = 5L,
                               seed = 2026L, e_trim = 0.02) {
  n <- nrow(X_train)
  fold_id <- make_stratified_folds(y_train, w_train, k = k, seed = seed)
  m_hat <- rep(NA_real_, n)
  e_hat <- rep(NA_real_, n)

  for (fold in sort(unique(fold_id))) {
    idx_valid <- which(fold_id == fold)
    idx_fit   <- which(fold_id != fold)
    if (length(idx_valid) == 0 || length(idx_fit) < 30) next

    y_fit <- as.numeric(y_train[idx_fit])
    w_fit <- as.numeric(w_train[idx_fit])
    x_fit   <- X_train[idx_fit, , drop = FALSE]
    x_valid <- X_train[idx_valid, , drop = FALSE]

    y_model <- fit_binomial_glmnet(x_fit, y_fit)
    w_model <- fit_binomial_glmnet(x_fit, w_fit)

    m_hat[idx_valid] <- predict_binomial_glmnet(y_model, x_valid, mean(y_fit))
    e_hat[idx_valid] <- predict_binomial_glmnet(w_model, x_valid, mean(w_fit))
  }

  m_hat[!is.finite(m_hat)] <- mean(y_train, na.rm = TRUE)
  e_hat[!is.finite(e_hat)] <- mean(w_train, na.rm = TRUE)
  e_hat <- clip_probs(e_hat, lower = e_trim, upper = 1 - e_trim)
  m_hat <- clip_probs(m_hat)

  list(m_hat = m_hat, e_hat = e_hat, fold_id = fold_id)
}

# --- Pseudo-outcome construction ---

build_pseudo_outcome <- function(y_train, w_train, m_hat, e_hat,
                                  min_residual = 1e-3) {
  rw <- as.numeric(w_train) - e_hat
  ry <- as.numeric(y_train) - m_hat
  denom <- ifelse(abs(rw) < min_residual, NA_real_, rw)
  z <- ry / denom
  weights <- rw^2
  ok <- is.finite(z) & is.finite(weights) & weights > 0
  list(z = z, weights = weights, ok = ok)
}

# --- R-Lasso: tau estimation via Lasso ---

fit_tau_lasso <- function(X_train, z, weights) {
  cv_fit <- glmnet::cv.glmnet(
    x = X_train, y = z, weights = weights,
    family = "gaussian", alpha = 1, nfolds = 5
  )
  list(cv_fit = cv_fit, lambda = cv_fit$lambda.min)
}

predict_tau_lasso <- function(model, X_new) {
  as.numeric(predict(model$cv_fit, newx = X_new,
                     s = model$lambda, type = "response"))
}

# --- R-Boost: tau estimation via XGBoost ---

fit_tau_boost <- function(X_train, z, weights, cfg) {
  dtrain <- xgboost::xgb.DMatrix(data = X_train, label = z, weight = weights)
  params <- list(
    objective = "reg:squarederror",
    eta = cfg$rboost_eta,
    max_depth = cfg$rboost_max_depth,
    subsample = cfg$rboost_subsample,
    colsample_bytree = cfg$rboost_colsample,
    min_child_weight = cfg$rboost_min_child_weight,
    eval_metric = "rmse",
    nthread = 1L
  )
  set.seed(cfg$seed)
  xgboost::xgb.train(
    params = params, data = dtrain,
    nrounds = cfg$rboost_nrounds, verbose = 0
  )
}

predict_tau_boost <- function(model, X_new) {
  as.numeric(predict(model, newdata = xgboost::xgb.DMatrix(X_new)))
}

# --- Main function: train R-Lasso and R-Boost ---

train_r_learner <- function(state) {
  cfg <- state$cfg
  df  <- state$df
  train_idx <- state$train_idx

  X_train <- state$X_all[train_idx, , drop = FALSE]
  y_train <- state$y_horizon[train_idx]
  w_train <- as.numeric(df[[cfg$w_var]][train_idx])

  log_section("R-learner: Cross-fitting Nuisance Parameters")
  nuisance <- crossfit_nuisance(
    X_train = X_train, y_train = y_train, w_train = w_train,
    k = cfg$rlearner_cf_folds, seed = cfg$seed, e_trim = cfg$rlearner_e_trim
  )

  log_section("R-learner: Building Pseudo-outcome")
  pseudo <- build_pseudo_outcome(
    y_train = y_train, w_train = w_train,
    m_hat = nuisance$m_hat, e_hat = nuisance$e_hat,
    min_residual = cfg$rlearner_min_residual
  )

  ok <- pseudo$ok
  X_tau  <- X_train[ok, , drop = FALSE]
  z_tau  <- pseudo$z[ok]
  wt_tau <- pseudo$weights[ok]
  cat("Pseudo-outcome: kept", nrow(X_tau), "/", length(train_idx), "rows\n")
  if (nrow(X_tau) < 50) {
    stop("Too few valid pseudo-outcome rows for R-learner fitting.")
  }

  log_section("R-learner: Fitting R-Lasso")
  lasso_model <- fit_tau_lasso(X_tau, z_tau, wt_tau)

  log_section("R-learner: Fitting R-Boost")
  boost_model <- fit_tau_boost(X_tau, z_tau, wt_tau, cfg)

  # Predict on all data
  tau_lasso_all <- predict_tau_lasso(lasso_model, state$X_all)
  tau_boost_all <- predict_tau_boost(boost_model, state$X_all)

  list(
    lasso_model = lasso_model,
    boost_model = boost_model,
    tau_lasso_all = tau_lasso_all,
    tau_boost_all = tau_boost_all,
    nuisance = nuisance,
    diagnostics = data.frame(
      train_n = length(train_idx),
      pseudo_n = nrow(X_tau),
      pseudo_kept_pct = 100 * nrow(X_tau) / length(train_idx),
      cf_folds = cfg$rlearner_cf_folds,
      e_trim = cfg$rlearner_e_trim,
      min_residual = cfg$rlearner_min_residual,
      stringsAsFactors = FALSE
    )
  )
}
