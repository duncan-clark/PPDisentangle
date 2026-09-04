# Helpers for Fit-D initialization-variability summaries in the Oklahoma report.
# Kept separate from oklahoma_report_support.R so tests do not pull geometry.

oklahoma_fitvar_d_summary <- function(results) {
  empty <- data.frame(
    rep = integer(),
    success = logical(),
    mc_total_saved_mean = numeric(),
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
  if (!"success" %in% names(df)) df$success <- TRUE
  if (!"rep" %in% names(df)) df$rep <- seq_len(nrow(df))
  df$success <- as.logical(df$success)
  df$mc_total_saved_mean <- suppressWarnings(as.numeric(df$mc_total_saved_mean))
  df$raw_total_saved <- suppressWarnings(as.numeric(df$raw_total_saved))
  df$n_relabel <- suppressWarnings(as.integer(df$n_relabel))
  df$n_label_init_flips <- suppressWarnings(as.integer(df$n_label_init_flips))
  df[, c("rep", "success", "mc_total_saved_mean", "raw_total_saved",
         "n_relabel", "n_label_init_flips"), drop = FALSE]
}

oklahoma_fitvar_d_ates <- function(results) {
  df <- oklahoma_fitvar_d_summary(results)
  if (nrow(df) < 1L) return(numeric(0))
  ok <- !is.na(df$success) & df$success
  vals <- df$mc_total_saved_mean[ok]
  vals[is.finite(vals)]
}
