test_that("KDE background lookups survive serialization when the image is captured", {
  skip_if_not_installed("spatstat.geom")
  win <- spatstat.geom::owin(c(0, 10), c(0, 10))
  covariate_im <- spatstat.geom::as.im(function(x, y) 1 + x / 10, W = win, dimyx = 32)

  make_lookup_captured <- function(ref_win, im) {
    force(im)
    cov_in_window <- im[ref_win, drop = FALSE]
    total_mass <- spatstat.geom::integral.im(cov_in_window)
    target_area <- spatstat.geom::area(ref_win)
    norm_factor <- target_area / total_mass
    im_local <- im
    force(norm_factor)
    force(im_local)
    function(x, y) {
      vals <- spatstat.geom::interp.im(im_local, x, y)
      vals[!is.finite(vals) | vals < 0] <- 0
      vals * norm_factor
    }
  }

  lookup <- make_lookup_captured(win, covariate_im)
  restored <- unserialize(serialize(lookup, NULL))
  rm(covariate_im)
  vals <- restored(c(1, 5, 9), c(1, 5, 9))
  expect_length(vals, 3L)
  expect_true(all(is.finite(vals)))
  expect_true(all(vals > 0))
})
