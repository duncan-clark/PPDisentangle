source(testthat::test_path("..", "..", "inst", "oklahoma", "paper", "paper_result_helpers.R"), local = TRUE)

test_that("publication selection distinguishes observed and all-or-nothing", {
  observed <- list(all_nothing_sim = data.frame(total_saved = 175.2))
  aon <- list(all_nothing_sim = data.frame(total_saved = 204.472))
  saved <- c(aon, list(contrast = "all_or_nothing",
                      by_contrast = list(observed = observed, all_or_nothing = aon)))
  expect_identical(select_paper_ate(saved, "observed"), observed)
  expect_identical(select_paper_ate(saved, "all_or_nothing"), aon)
  expect_error(select_paper_ate(aon, "observed"), "do not contain")
  expect_identical(select_paper_ate(aon, "all_or_nothing"), aon)
})

test_that("bootstrap selection cannot silently relabel a legacy estimand", {
  cols <- c("ate_total_mean", "ate_total_mean_observed", "ate_total_mean_all_or_nothing")
  expect_identical(select_paper_bootstrap_column(cols, "observed"), "ate_total_mean_observed")
  expect_identical(select_paper_bootstrap_column(cols, "all_or_nothing"), "ate_total_mean_all_or_nothing")
  expect_error(select_paper_bootstrap_column("ate_total_mean", "observed", "all_or_nothing"), "refusing")
  expect_identical(select_paper_bootstrap_column("ate_total_mean", "observed", "observed"), "ate_total_mean")
})
