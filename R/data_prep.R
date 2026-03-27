###############################################################################
# data_prep.R - Data loading, imputation, and train/test splitting
#
# Loads clinical data, performs MICE imputation for missing covariates,
# and splits into train/test/external sets using odd/even cluster sites.
###############################################################################

load_and_prepare_data <- function(cfg) {
  # --- Helper: read CSV or Excel ---
  read_tabular_data <- function(path) {
    ext <- tolower(tools::file_ext(path))
    if (ext %in% c("xlsx", "xls")) {
      return(as.data.frame(readxl::read_excel(path)))
    }
    if (ext == "csv") {
      return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE))
    }
    stop(sprintf("Unsupported data file extension: %s", ext))
  }

  bind_rows_fill <- function(a, b) {
    miss_a <- setdiff(colnames(b), colnames(a))
    miss_b <- setdiff(colnames(a), colnames(b))
    for (col in miss_a) a[[col]] <- NA
    for (col in miss_b) b[[col]] <- NA
    b <- b[, colnames(a), drop = FALSE]
    rbind(a, b)
  }

  log_section("Data Loading")

  if (!file.exists(cfg$data_file)) {
    stop("Data file not found: ", cfg$data_file)
  }
  df <- read_tabular_data(cfg$data_file)
  df$.external_override <- FALSE
  cat("Using data file:", cfg$data_file, "\n")

  # --- Load external validation data if provided ---
  external_used <- NA_character_
  if (is.character(cfg$external_data_file) && !is.na(cfg$external_data_file) &&
      nzchar(cfg$external_data_file) && file.exists(cfg$external_data_file)) {
    ext_df <- read_tabular_data(cfg$external_data_file)
    ext_df$.external_override <- TRUE
    ext_df$mafa <- 0
    if (!("cluster_site" %in% colnames(ext_df))) {
      ext_df$cluster_site <- seq_len(nrow(ext_df)) + 1000000L
    }
    if ("mafa" %in% colnames(df)) {
      df <- df[is.na(df$mafa) | df$mafa == 1, , drop = FALSE]
    }
    df <- bind_rows_fill(df, ext_df)
    external_used <- cfg$external_data_file
    cat("Using external validation file:", cfg$external_data_file, "\n")
  }

  cat("Data shape:", nrow(df), "x", ncol(df), "\n")

  # --- Validate required columns ---
  required_cols <- unique(c(
    cfg$feature_vars, cfg$y_time_var, cfg$y_event_var, cfg$w_var,
    "mafa", cfg$split_col, "age"
  ))
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # --- Age filter ---
  before_n <- nrow(df)
  if (cfg$apply_age_filter) {
    df <- df[!is.na(df$age) & df$age >= cfg$age_min, , drop = FALSE]
    cat("Filtered age >=", cfg$age_min, ":", before_n, "->", nrow(df), "\n")
  } else {
    df <- df[!is.na(df$age), , drop = FALSE]
    cat("No lower age filter applied:", before_n, "->", nrow(df), "\n")
  }

  # --- Coerce outcome and treatment to 0/1 ---
  df[[cfg$y_event_var]] <- coerce_binary01(df[[cfg$y_event_var]], cfg$y_event_var)
  df[[cfg$w_var]] <- coerce_binary01(df[[cfg$w_var]], cfg$w_var)

  if (any(df[[cfg$y_time_var]] < 0, na.rm = TRUE)) {
    stop(sprintf("%s must be non-negative.", cfg$y_time_var))
  }

  # --- MICE imputation for covariates ---
  log_section("Covariate Imputation")
  feature_vars <- intersect(cfg$feature_vars, colnames(df))
  if (length(feature_vars) == 0) {
    stop("No feature variables found in data.")
  }

  imp <- mice(df[, feature_vars, drop = FALSE],
              m = 1, method = "pmm", maxit = 5, printFlag = FALSE)
  df_imputed <- complete(imp, 1)
  for (col in feature_vars) {
    df[[col]] <- df_imputed[[col]]
  }
  X_all <- as.matrix(df[, feature_vars, drop = FALSE])

  # --- Train / Test / External split ---
  log_section("Train/Test/External Split")
  eligible_idx <- which(df$mafa == 1 & !is.na(df$mafa))
  external_idx <- if (any(df$.external_override %in% TRUE, na.rm = TRUE)) {
    which(df$.external_override %in% TRUE)
  } else {
    which(df$mafa == 0 & !is.na(df$mafa))
  }

  if (length(eligible_idx) < 10) {
    stop("Too few mafa==1 patients for train/test split.")
  }

  # Odd/even split by cluster_site
  split_col <- cfg$split_col
  cluster_int <- suppressWarnings(as.integer(df[[split_col]][eligible_idx]))
  if (any(is.na(cluster_int))) {
    stop(split_col, " cannot be safely converted to integer.")
  }
  train_idx <- eligible_idx[cluster_int %% 2 == 0]
  test_idx  <- eligible_idx[cluster_int %% 2 == 1]
  if (length(train_idx) == 0 || length(test_idx) == 0) {
    stop("Even/odd cluster split produced an empty train or test set.")
  }
  cat("Using odd/even split on", split_col, "\n")

  # Drop rows with missing outcome/treatment
  core_cols <- c(cfg$y_time_var, cfg$y_event_var, cfg$w_var)
  train_ok <- complete.cases(df[train_idx, core_cols, drop = FALSE])
  test_ok  <- complete.cases(df[test_idx, core_cols, drop = FALSE])
  if (sum(!train_ok) > 0 || sum(!test_ok) > 0) {
    cat("Dropping missing Y/D/W rows: train", sum(!train_ok), ", test", sum(!test_ok), "\n")
  }
  train_idx <- train_idx[train_ok]
  test_idx  <- test_idx[test_ok]
  external_idx <- external_idx[complete.cases(df[external_idx, core_cols, drop = FALSE])]

  # --- Compute binary event indicator at horizon ---
  y_horizon <- as.integer(
    df[[cfg$y_event_var]] == 1 & df[[cfg$y_time_var]] <= cfg$horizon
  )

  cat("Train:", length(train_idx),
      " Test:", length(test_idx),
      " External:", length(external_idx), "\n")

  list(
    cfg = cfg,
    df = df,
    X_all = X_all,
    feature_vars = feature_vars,
    y_horizon = y_horizon,
    train_idx = train_idx,
    test_idx = test_idx,
    external_idx = external_idx,
    overall_idx = c(train_idx, test_idx)
  )
}
