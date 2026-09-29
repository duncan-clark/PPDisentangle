# Select explicitly identified estimands from frozen application results.
# Legacy unlabelled ATE objects contain the all-or-nothing contrast only.
select_paper_ate <- function(ate_obj, contrast = "observed") {
  contrast <- match.arg(contrast, c("observed", "all_or_nothing"))
  if (is.null(ate_obj)) stop("Missing saved ATE results.", call. = FALSE)
  if (is.list(ate_obj$by_contrast) && !is.null(ate_obj$by_contrast[[contrast]])) {
    return(ate_obj$by_contrast[[contrast]])
  }
  recorded <- ate_obj$contrast
  if (!is.null(recorded) && identical(as.character(recorded)[1L], contrast)) {
    return(ate_obj)
  }
  if (identical(contrast, "all_or_nothing") && is.null(recorded) &&
      !is.null(ate_obj$all_nothing_sim)) return(ate_obj)
  stop("Saved results do not contain the requested ", contrast,
       " contrast; use the matching results deposit.", call. = FALSE)
}

select_paper_bootstrap_column <- function(column_names, contrast = "observed",
                                          recorded_contrast = NULL) {
  contrast <- match.arg(contrast, c("observed", "all_or_nothing"))
  explicit <- paste0("ate_total_mean_", contrast)
  if (explicit %in% column_names) return(explicit)
  legacy_matches <- identical(as.character(recorded_contrast)[1L], contrast) ||
    (is.null(recorded_contrast) && identical(contrast, "all_or_nothing"))
  if (legacy_matches && "ate_total_mean" %in% column_names) return("ate_total_mean")
  stop("Saved bootstrap does not identify the requested ", contrast,
       " contrast; refusing to substitute another estimand.", call. = FALSE)
}
