###############################################################################
# plots.R - Visualization functions
#
# Generates:
#   - Qini curves for all models (overlaid)
#   - ITE distribution histograms
#   - Kaplan-Meier curves by benefit/harm groups
#   - Model comparison bar charts
###############################################################################

# --- Qini curve overlay plot ---

plot_qini_curves <- function(qini_curves, results_summary, cfg) {
  model_colors <- c(
    "UpliftRF" = "#e41a1c",
    "R_lasso"  = "#377eb8",
    "R_boost"  = "#4daf4a"
  )

  for (dname in c("Test", "External")) {
    curve_list <- list()
    for (mname in c("UpliftRF", "R_lasso", "R_boost")) {
      key <- paste0(mname, "_", dname)
      if (key %in% names(qini_curves) && !is.null(qini_curves[[key]])) {
        qc <- qini_curves[[key]]
        qc$method <- mname
        # Get Qini score for label
        sub <- results_summary[results_summary$method == mname &
                                results_summary$dataset == dname, , drop = FALSE]
        qini_val <- if (nrow(sub) > 0) round(sub$qini[1], 6) else NA
        qc$label <- paste0(mname, " (Qini=", qini_val, ")")
        curve_list[[key]] <- qc
      }
    }
    if (length(curve_list) == 0) next

    plot_df <- dplyr::bind_rows(curve_list)

    p <- ggplot(plot_df, aes(x = percentile, y = qini_value,
                              color = method, linetype = method)) +
      geom_line(linewidth = 1.1) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "#555555") +
      scale_color_manual(
        values = model_colors,
        labels = setNames(
          sapply(unique(plot_df$label), identity),
          sapply(unique(plot_df$method), identity)
        )
      ) +
      scale_linetype_manual(
        values = c("UpliftRF" = "solid", "R_lasso" = "dashed",
                   "R_boost" = "dotdash"),
        guide = "none"
      ) +
      labs(
        title = paste0("Qini Curves - ", dname, " Set"),
        x = "Population percentile (sorted by predicted ITE)",
        y = "Cumulative gain over random allocation",
        color = "Model"
      ) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "bottom")

    save_output_plot(p, paste0("figure_qini_curves_", tolower(dname), ".png"),
                     cfg, width = 9, height = 6, dpi = 300)
  }
}

# --- ITE distribution plot ---

plot_ite_distributions <- function(results, state, cfg) {
  pred_df <- results$predictions
  if (is.null(pred_df) || nrow(pred_df) == 0) return(invisible(NULL))

  plot_df <- pred_df[pred_df$dataset %in% c("Test", "External"), , drop = FALSE]
  plot_df$ite_pct <- 100 * plot_df$tau_survival

  fill_values <- c(
    "Benefit from intervention" = "#4F81BD",
    "Harm from intervention"    = "#C0504D"
  )
  plot_df$ite_group <- ifelse(plot_df$tau_survival >= 0,
                               "Benefit from intervention",
                               "Harm from intervention")
  plot_df$ite_group <- factor(plot_df$ite_group,
                               levels = c("Benefit from intervention",
                                          "Harm from intervention"))

  for (mname in c("UpliftRF", "R_lasso", "R_boost")) {
    sub <- plot_df[plot_df$method == mname, , drop = FALSE]
    if (nrow(sub) == 0) next

    bin_width <- max(0.1, diff(range(sub$ite_pct, na.rm = TRUE)) / 30)

    p <- ggplot(sub, aes(x = ite_pct, fill = ite_group)) +
      geom_histogram(binwidth = bin_width, color = NA, alpha = 0.85) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "#555555",
                 linewidth = 0.8) +
      facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
      scale_fill_manual(values = fill_values, drop = FALSE) +
      labs(
        title = paste0("ITE Distribution - ", mname),
        x = "Predicted treatment effect (%) [tau_survival * 100]",
        y = "Number of patients",
        fill = NULL
      ) +
      theme_minimal(base_size = 12)

    save_output_plot(p, paste0("figure_ite_distribution_",
                               tolower(mname), ".png"),
                     cfg, width = 10, height = 7, dpi = 300)
  }
}

# --- Kaplan-Meier curve by treatment ---

plot_km_by_treatment <- function(state, cfg) {
  df <- state$df

  for (dname in c("Test", "External")) {
    idx <- if (dname == "Test") state$test_idx else state$external_idx
    if (length(idx) == 0) next

    d <- df[idx, c(cfg$y_time_var, cfg$y_event_var, cfg$w_var), drop = FALSE]
    d <- d[complete.cases(d), , drop = FALSE]
    if (nrow(d) < 10) next

    d$grp <- factor(d[[cfg$w_var]], levels = c(0, 1),
                    labels = c("Control", "Intervention"))
    fit <- survival::survfit(
      survival::Surv(d[[cfg$y_time_var]], d[[cfg$y_event_var]]) ~ grp,
      data = d
    )
    lr <- survival::survdiff(
      survival::Surv(d[[cfg$y_time_var]], d[[cfg$y_event_var]]) ~ grp,
      data = d
    )
    p_val <- 1 - pchisq(lr$chisq, df = length(lr$n) - 1)

    out_path <- output_path(cfg, paste0("figure_km_", tolower(dname),
                                         "_by_treatment.png"))
    png(out_path, width = 1800, height = 1300, res = 220)
    plot(fit, col = c("#1f78b4", "#e31a1c"), lwd = 2, mark.time = TRUE,
         conf.int = FALSE, xlab = "Follow-up time (years)",
         ylab = "Survival probability",
         main = paste0("Kaplan-Meier by Treatment (", dname, ")"))
    legend("bottomleft",
           legend = c(paste0("Control (n=", sum(d$grp == "Control"), ")"),
                      paste0("Intervention (n=", sum(d$grp == "Intervention"),
                             ")")),
           col = c("#1f78b4", "#e31a1c"), lwd = 2, bty = "n")
    mtext(paste0("Log-rank p = ", formatC(p_val, format = "g", digits = 4)),
          side = 3, line = -1.2, adj = 1, cex = 0.9)
    dev.off()
    cat("Saved:", out_path, "\n")
  }
}

# --- Model comparison summary bar chart ---

plot_model_comparison <- function(results_summary, cfg) {
  sub <- results_summary[results_summary$dataset %in% c("Test", "External"),
                          , drop = FALSE]
  if (nrow(sub) == 0) return(invisible(NULL))

  model_colors <- c(
    "UpliftRF" = "#e41a1c",
    "R_lasso"  = "#377eb8",
    "R_boost"  = "#4daf4a"
  )

  # Qini comparison
  p_qini <- ggplot(sub, aes(x = method, y = qini, fill = method)) +
    geom_bar(stat = "identity", alpha = 0.85) +
    geom_errorbar(aes(ymin = qini_ci_lower, ymax = qini_ci_upper),
                  width = 0.2, na.rm = TRUE) +
    facet_wrap(~ dataset) +
    scale_fill_manual(values = model_colors) +
    labs(title = "Qini Score Comparison",
         x = NULL, y = "Qini Score", fill = NULL) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")

  save_output_plot(p_qini, "figure_model_qini_comparison.png",
                   cfg, width = 9, height = 5, dpi = 300)

  # C-for-benefit comparison
  if ("c_for_benefit" %in% colnames(sub) && any(!is.na(sub$c_for_benefit))) {
    p_cfb <- ggplot(sub, aes(x = method, y = c_for_benefit, fill = method)) +
      geom_bar(stat = "identity", alpha = 0.85) +
      geom_errorbar(aes(ymin = c_for_benefit_ci_lower,
                        ymax = c_for_benefit_ci_upper),
                    width = 0.2, na.rm = TRUE) +
      facet_wrap(~ dataset) +
      scale_fill_manual(values = model_colors) +
      labs(title = "C-for-Benefit Comparison",
           x = NULL, y = "C-for-Benefit", fill = NULL) +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")

    save_output_plot(p_cfb, "figure_model_cfb_comparison.png",
                     cfg, width = 9, height = 5, dpi = 300)
  }
}

# --- Generate all plots ---

generate_all_plots <- function(results, state, cfg) {
  log_section("Generating Plots")

  plot_qini_curves(results$qini_curves, results$summary, cfg)
  plot_ite_distributions(results, state, cfg)
  plot_km_by_treatment(state, cfg)
  plot_model_comparison(results$summary, cfg)
}
