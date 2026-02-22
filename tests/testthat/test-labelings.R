test_that("oracle labeling sets inferred_process = process", {
  df <- data.frame(x = 1:4, y = 1:4, process = c("control", "treated", "control", "treated"))
  result <- oracle_labeling(df)
  expect_equal(result$inferred_process, result$process)
})

test_that("naive labeling sets inferred_process = location_process", {
  df <- data.frame(
    x = 1:4, y = 1:4,
    process = c("control", "treated", "control", "treated"),
    location_process = c("control", "control", "treated", "treated")
  )
  result <- naive_labeling(df)
  expect_equal(result$inferred_process, result$location_process)
})

test_that("normalize_weights returns valid probability vector", {
  w <- c(1, 2, 3, 4, 5)
  nw <- normalize_weights(w)
  expect_equal(sum(nw), 1)
  expect_true(all(nw >= 0))

  w_same <- rep(3, 5)
  nw_same <- normalize_weights(w_same)
  expect_equal(nw_same, rep(0.2, 5))
})
