###############################################################################
# uplift_rf.R - Uplift Random Forest via transformed outcome
#
# Implements the transformed outcome approach for uplift modeling:
#   Y_TO = Y * (W - p_treat) / (p_treat * (1 - p_treat))
#
# where Y is the observed outcome, W is the treatment indicator, and
# p_treat is the marginal treatment probability.
#
# A standard random forest is then trained on (X, Y_TO) to predict the
# individualized treatment effect tau(x).
#
# Reference: Athey & Imbens (2016), "Recursive partitioning for
#   heterogeneous causal effects"
###############################################################################

train_uplift_rf <- function(state) {
  cfg <- state$cfg
  df  <- state$df
  train_idx <- state$train_idx

  X_train <- state$X_all[train_idx, , drop = FALSE]
  y_train <- as.numeric(state$y_horizon[train_idx])
  w_train <- as.numeric(df[[cfg$w_var]][train_idx])

  log_section("UpliftRF: Training")

  ntree   <- cfg$upliftrf_ntree
  nodesize <- cfg$upliftrf_nodesize
  mtry <- if (is.na(cfg$upliftrf_mtry)) floor(sqrt(ncol(X_train))) else cfg$upliftrf_mtry
  mtry <- max(1L, mtry)

  # Marginal treatment probability (clipped to [0.05, 0.95])
  p_treat <- mean(w_train, na.rm = TRUE)
  p_treat <- min(0.95, max(0.05, p_treat))

  # Transformed outcome
  y_to <- y_train * ((w_train - p_treat) / (p_treat * (1 - p_treat)))

  cat(sprintf("p_treat=%.4f, ntree=%d, mtry=%d, nodesize=%d\n",
              p_treat, ntree, mtry, nodesize))

  rf <- randomForest::randomForest(
    x = X_train, y = y_to,
    ntree = ntree, mtry = mtry, nodesize = nodesize,
    importance = TRUE
  )

  # Predict on all data
  tau_upliftrf_all <- as.numeric(stats::predict(rf, newdata = state$X_all))

  list(
    model = rf,
    tau_upliftrf_all = tau_upliftrf_all,
    p_treat = p_treat,
    diagnostics = data.frame(
      method = "UpliftRF",
      p_treat = p_treat,
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      stringsAsFactors = FALSE
    )
  )
}
