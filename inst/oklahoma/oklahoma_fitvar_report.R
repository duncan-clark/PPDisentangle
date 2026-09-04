# Helpers for Fit-D initialization-variability summaries in the Oklahoma report.
# Kept separate from oklahoma_report_support.R so tests do not pull geometry.

oklahoma_fitvar_d_summary <- function(results) {
  empty <- data.frame(
    rep = integer(),
    success = logical(),
    mc_total_saved_mean = numeric(),
    mc_total_saved_mean_all_or_nothing = numeric(),
    mc_total_saved_mean_observed = numeric(),
    raw_total_saved = numeric(),
    n_relabel = integer(),
    n_label_init_flips = integer(),
    stringsAsFactors = FALSE
  )
  fv <- results$fit_variability
  if (is.null(fv) || is.null(fv$replicate_summary)) return(empty)
  df <- as.data.frame(fv$replicate_summary)
  if (nrow(df) < 1L) return(empty)
  model <- toupper(as.character(df$model))
  df <- df[model %in% c("F", "D"), , drop = FALSE]
  if (nrow(df) < 1L) return(empty)
  if (!"n_label_init_flips" %in% names(df)) df$n_label_init_flips <- NA_integer_
  if (!"n_relabel" %in% names(df)) df$n_relabel <- NA_integer_
  if (!"raw_total_saved" %in% names(df)) df$raw_total_saved <- NA_real_
  if (!"mc_total_saved_mean" %in% names(df)) df$mc_total_saved_mean <- NA_real_
  if (!"mc_total_saved_mean_all_or_nothing" %in% names(df)) {
    df$mc_total_saved_mean_all_or_nothing <- df$mc_total_saved_mean
  }
  if (!"mc_total_saved_mean_observed" %in% names(df)) {
    df$mc_total_saved_mean_observed <- NA_real_
  }
  if (!"success" %in% names(df)) df$success <- TRUE
  if (!"rep" %in% names(df)) df$rep <- seq_len(nrow(df))
  df$success <- as.logical(df$success)
  df$mc_total_saved_mean <- suppressWarnings(as.numeric(df$mc_total_saved_mean))
  df$mc_total_saved_mean_all_or_nothing <- suppressWarnings(
    as.numeric(df$mc_total_saved_mean_all_or_nothing)
  )
  df$mc_total_saved_mean_observed <- suppressWarnings(
    as.numeric(df$mc_total_saved_mean_observed)
  )
  df$raw_total_saved <- suppressWarnings(as.numeric(df$raw_total_saved))
  df$n_relabel <- suppressWarnings(as.integer(df$n_relabel))
  df$n_label_init_flips <- suppressWarnings(as.integer(df$n_label_init_flips))
  df[, c("rep", "success", "mc_total_saved_mean",
         "mc_total_saved_mean_all_or_nothing", "mc_total_saved_mean_observed",
         "raw_total_saved", "n_relabel", "n_label_init_flips"), drop = FALSE]
}

oklahoma_fitvar_d_ates <- function(results, contrast = c("primary", "all_or_nothing", "observed")) {
  contrast <- match.arg(contrast)
  df <- oklahoma_fitvar_d_summary(results)
  if (nrow(df) < 1L) return(numeric(0))
  ok <- !is.na(df$success) & df$success
  col <- switch(
    contrast,
    observed = "mc_total_saved_mean_observed",
    all_or_nothing = "mc_total_saved_mean_all_or_nothing",
    primary = "mc_total_saved_mean"
  )
  if (!col %in% names(df)) return(numeric(0))
  vals <- df[[col]][ok]
  vals[is.finite(vals)]
}

oklahoma_boot_vals_for_contrast <- function(df, contrast = c("all_or_nothing", "observed")) {
  contrast <- match.arg(contrast)
  if (is.null(df) || !is.data.frame(df) || nrow(df) < 1L) return(numeric(0))
  col <- if (identical(contrast, "observed") &&
             "ate_total_mean_observed" %in% names(df)) {
    "ate_total_mean_observed"
  } else if (identical(contrast, "all_or_nothing") &&
             "ate_total_mean_all_or_nothing" %in% names(df)) {
    "ate_total_mean_all_or_nothing"
  } else if (identical(contrast, "observed")) {
    return(numeric(0))
  } else {
    "ate_total_mean"
  }
  if (!col %in% names(df)) return(numeric(0))
  vals <- suppressWarnings(as.numeric(df[[col]]))
  vals[is.finite(vals)]
}
