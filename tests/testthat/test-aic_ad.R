test_that("aic_ad matches manual formula for Gaussian AD fit", {
  set.seed(410)
  y <- simulate_ad(
    n_subjects = 70,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.35, 4))
  )

  fit <- fit_ad(y, order = 1, estimate_mu = TRUE)
  aic_pkg <- aic_ad(fit)

  N <- ncol(y)
  k <- N + N + (N - 1)
  aic_manual <- -2 * fit$log_l + 2 * k

  expect_equal(aic_pkg, aic_manual)
})
