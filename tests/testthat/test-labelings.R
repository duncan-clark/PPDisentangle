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
  expect_equal(length(nw), length(w))

  w_same <- rep(3, 5)
  nw_same <- normalize_weights(w_same)
  expect_equal(nw_same, rep(0.2, 5))
})

test_that("normalize_weights uses log-sum-exp correctly", {
  w <- c(-100, -90, -80)
  nw <- normalize_weights(w)
  expect_equal(sum(nw), 1)
  expect_equal(length(nw), 3)
  expect_true(nw[3] > 0.99)
  expect_true(nw[1] < 0.001)
  expected <- exp(w - max(w)) / sum(exp(w - max(w)))
  expect_equal(nw, expected)
})

test_that("normalize_weights handles -Inf with zero weight", {
  w <- c(-500, -Inf, -490)
  nw <- normalize_weights(w)
  expect_equal(length(nw), 3)
  expect_equal(nw[2], 0)
  expect_equal(sum(nw), 1)
  expect_true(nw[3] > nw[1])
})

test_that("normalize_weights edge cases", {
  expect_equal(normalize_weights(numeric(0)), numeric(0))
  expect_equal(normalize_weights(42), 1)
  nw_inf <- normalize_weights(c(-Inf, -Inf))
  expect_equal(nw_inf, c(0.5, 0.5))
})
