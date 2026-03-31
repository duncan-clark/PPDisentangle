#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
rds_path <- if (length(args) >= 1L) args[[1]] else "output/oklahoma/oklahoma_results_pre_bootstrap_job5206299.rds"
out_prefix <- if (length(args) >= 2L) args[[2]] else file.path("output", "oklahoma", "sem_mismatch_investigation")

if (!file.exists(rds_path)) {
  stop("RDS file not found: ", rds_path)
}

results <- readRDS(rds_path)
fits_named <- results$fits_named
if (is.null(fits_named)) stop("No fits_named found in results object.")

model_ids <- c("B", "D", "F", "H")

extract_selected_trace <- function(fit_obj) {
  ad <- fit_obj$adaptive
  if (is.null(ad) || is.null(ad$all_metrics) || is.null(ad$class_results)) return(NULL)
  all_metrics <- ad$all_metrics
  class_results <- ad$class_results
  if (length(all_metrics) < 1L || length(class_results) < 1L) return(NULL)

  out <- vector("list", length(all_metrics))
  offset <- 0L
  for (i in seq_along(all_metrics)) {
    metric_i <- as.numeric(all_metrics[[i]])
    if (length(metric_i) < 1L) next
    best_j <- which.max(metric_i)
    class_idx <- offset + best_j
    if (class_idx > length(class_results)) break
    cm <- class_results[[class_idx]]
    out[[i]] <- data.frame(
      iter = i,
      control_as_control = as.numeric(cm$control_as_control),
      control_as_treated = as.numeric(cm$control_as_treated),
      treated_as_control = as.numeric(cm$treated_as_control),
      treated_as_treated = as.numeric(cm$treated_as_treated)
    )
    offset <- offset + length(metric_i)
  }
  out <- Filter(Negate(is.null), out)
  if (length(out) < 1L) return(NULL)
  rbindlist(out)
}

extract_final_confusion <- function(fit_obj) {
  ad <- fit_obj$adaptive
  if (is.null(ad) || is.null(ad$adaptive_labelling)) return(NULL)
  post <- ad$adaptive_labelling
  post <- post[post$t >= 0, , drop = FALSE]
  if (nrow(post) < 1L) return(NULL)
  tab <- table(
    factor(post$location_process, levels = c("control", "treated")),
    factor(post$inferred_process, levels = c("control", "treated"))
  )
  as.list(tab)
}

summary_rows <- list()
trace_rows <- list()

for (id in model_ids) {
  obj <- fits_named[[id]]
  fit <- if (!is.null(obj)) obj$fit else NULL
  if (is.null(fit) || is.null(fit$adaptive)) next

  ad <- fit$adaptive
  final_cm <- extract_final_confusion(fit)
  selected_trace <- extract_selected_trace(fit)

  accepted <- as.numeric(ad$max_metric_flips)
  avg <- as.numeric(ad$average_flips)

  final_cc <- if (!is.null(final_cm$A1)) as.numeric(final_cm$A1) else as.numeric(final_cm[[1]])
  # Table order is control/control, treated/control, control/treated, treated/treated
  # but named list from table is robust:
  final_cc <- if (!is.null(final_cm$control.control)) as.numeric(final_cm$control.control) else as.numeric(final_cm[[1]])
  final_tc <- if (!is.null(final_cm$treated.control)) as.numeric(final_cm$treated.control) else as.numeric(final_cm[[2]])
  final_ct <- if (!is.null(final_cm$control.treated)) as.numeric(final_cm$control.treated) else as.numeric(final_cm[[3]])
  final_tt <- if (!is.null(final_cm$treated.treated)) as.numeric(final_cm$treated.treated) else as.numeric(final_cm[[4]])

  trace_last <- if (!is.null(selected_trace) && nrow(selected_trace) > 0) selected_trace[nrow(selected_trace)] else NULL

  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    model = id,
    n_iter = length(accepted),
    accepted_sum = sum(accepted, na.rm = TRUE),
    accepted_mean = mean(accepted, na.rm = TRUE),
    accepted_nonzero_iters = sum(accepted > 0, na.rm = TRUE),
    avg_sum = sum(avg, na.rm = TRUE),
    avg_mean = mean(avg, na.rm = TRUE),
    final_cc = final_cc,
    final_ct = final_ct,
    final_tc = final_tc,
    final_tt = final_tt,
    final_offdiag = final_ct + final_tc,
    trace_last_offdiag = if (!is.null(trace_last)) trace_last$control_as_treated + trace_last$treated_as_control else NA_real_,
    trace_matches_final = if (!is.null(trace_last)) {
      isTRUE(all.equal(as.numeric(trace_last$control_as_control), final_cc)) &&
        isTRUE(all.equal(as.numeric(trace_last$control_as_treated), final_ct)) &&
        isTRUE(all.equal(as.numeric(trace_last$treated_as_control), final_tc)) &&
        isTRUE(all.equal(as.numeric(trace_last$treated_as_treated), final_tt))
    } else {
      NA
    },
    adaptive_history_entries = length(fit$adaptive_history),
    adaptive_history_accepted_sum = if (length(fit$adaptive_history) > 0 && !is.null(fit$adaptive_history[[1]]$max_metric_flips)) {
      sum(as.numeric(fit$adaptive_history[[1]]$max_metric_flips), na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )

  if (!is.null(selected_trace) && nrow(selected_trace) > 0) {
    selected_trace$model <- id
    selected_trace$accepted_flips <- accepted[seq_len(min(length(accepted), nrow(selected_trace)))]
    trace_rows[[length(trace_rows) + 1L]] <- selected_trace
  }
}

summary_df <- if (length(summary_rows) > 0) rbindlist(summary_rows) else data.frame()
trace_df <- if (length(trace_rows) > 0) rbindlist(trace_rows) else data.frame()

cat("\n=== SEM Confusion vs Flips Investigation ===\n")
cat("RDS: ", normalizePath(rds_path, winslash = "/", mustWork = FALSE), "\n", sep = "")
if (nrow(summary_df) < 1L) {
  cat("No B/D/F/H SEM fits found.\n")
  quit(save = "no", status = 0)
}

print(summary_df)

cat("\nInterpretation checks:\n")
for (i in seq_len(nrow(summary_df))) {
  r <- summary_df[i, ]
  cat(sprintf(
    "Model %s: final_offdiag=%d, accepted_sum=%.0f, accepted_nonzero_iters=%d, trace_matches_final=%s\n",
    r$model, as.integer(r$final_offdiag), r$accepted_sum, as.integer(r$accepted_nonzero_iters), as.character(r$trace_matches_final)
  ))
  if (isTRUE(r$final_offdiag > r$accepted_sum)) {
    cat("  -> final off-diagonal exceeds accepted-sum in this run; indicates run likely started already relabelled.\n")
  }
}

summary_out <- paste0(out_prefix, "_summary.csv")
trace_out <- paste0(out_prefix, "_selected_trace.csv")
fwrite(summary_df, summary_out)
fwrite(trace_df, trace_out)
cat("\nWrote:\n")
cat("  - ", summary_out, "\n", sep = "")
cat("  - ", trace_out, "\n", sep = "")
