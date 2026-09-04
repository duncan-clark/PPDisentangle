# Fit-variability D ATE extractor used by oklahoma_report.qmd.

helper <- NULL
for (cand in c(
  "inst/oklahoma/oklahoma_fitvar_report.R",
  file.path("..", "..", "inst/oklahoma/oklahoma_fitvar_report.R")
)) {
  if (file.exists(cand)) {
    helper <- normalizePath(cand, winslash = "/", mustWork = FALSE)
    break
  }
}
skip_if(is.null(helper), "oklahoma_fitvar_report.R not found")
source(helper, local = TRUE)

test_that("fitvar D ATE helper keeps successful F rows and drops C/E", {
  results <- list(
    fit_variability = list(
      replicate_summary = data.frame(
        rep = c(1L, 2L, 1L, 3L),
        model = c("E", "F", "F", "F"),
        success = c(TRUE, TRUE, FALSE, TRUE),
        mc_total_saved_mean = c(-10, 211.4, NA_real_, 198.2),
        raw_total_saved = c(-9, 200, NA_real_, 190),
        n_relabel = c(NA_integer_, 288L, NA_integer_, 301L),
        n_label_init_flips = c(NA_integer_, 285L, 12L, 301L),
        stringsAsFactors = FALSE
      )
    )
  )
  ates <- oklahoma_fitvar_d_ates(results)
  expect_equal(ates, c(211.4, 198.2))
  sm <- oklahoma_fitvar_d_summary(results)
  expect_equal(nrow(sm), 3L)
  expect_equal(sm$n_label_init_flips[sm$success], c(285L, 301L))
})

test_that("fitvar D ATE helper is empty when the stage was not run", {
  expect_equal(oklahoma_fitvar_d_ates(list()), numeric(0))
  expect_equal(nrow(oklahoma_fitvar_d_summary(list(fit_variability = NULL))), 0L)
})
