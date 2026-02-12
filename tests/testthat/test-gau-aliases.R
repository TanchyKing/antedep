test_that("fit_gau is an alias of fit_ad", {
  set.seed(123)
  y <- simulate_ad(
    n_subjects = 40,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.3, 4))
  )

  fit_ad_obj <- fit_ad(y, order = 1)
  fit_gau_obj <- fit_gau(y, order = 1)

  expect_equal(fit_gau_obj$log_l, fit_ad_obj$log_l)
  expect_equal(fit_gau_obj$settings$order, fit_ad_obj$settings$order)
  expect_equal(fit_gau_obj$phi, fit_ad_obj$phi)
  expect_equal(fit_gau_obj$sigma, fit_ad_obj$sigma)
})

test_that("bic_gau is an alias of bic_ad", {
  set.seed(124)
  y <- simulate_ad(
    n_subjects = 50,
    n_time = 5,
    order = 1,
    phi = c(0, rep(0.4, 4))
  )
  fit <- fit_ad(y, order = 1)

  expect_equal(
    bic_gau(fit, n_subjects = nrow(y)),
    bic_ad(fit, n_subjects = nrow(y))
  )
})

test_that("bic_order_gau is an alias of bic_order_ad", {
  set.seed(125)
  y <- simulate_ad(
    n_subjects = 45,
    n_time = 6,
    order = 1,
    phi = c(0, rep(0.35, 5))
  )

  ord_ad <- bic_order_ad(y, max_order = 2)
  ord_gau <- bic_order_gau(y, max_order = 2)

  expect_equal(ord_gau$best_order, ord_ad$best_order)
  expect_equal(ord_gau$bic, ord_ad$bic)
  expect_equal(ord_gau$table$order, ord_ad$table$order)
})

test_that("aic_gau is an alias of aic_ad", {
  set.seed(126)
  y <- simulate_ad(n_subjects = 50, n_time = 5, order = 1, phi = c(0, rep(0.3, 4)))
  fit <- fit_ad(y, order = 1)
  expect_equal(aic_gau(fit), aic_ad(fit))
})

test_that("ci_gau is an alias of ci_ad", {
  set.seed(127)
  y <- simulate_ad(n_subjects = 70, n_time = 5, order = 1, phi = c(0, rep(0.35, 4)))
  fit <- fit_ad(y, order = 1)
  ci1 <- ci_ad(fit, parameters = "mu")
  ci2 <- ci_gau(fit, parameters = "mu")
  expect_equal(ci2$mu, ci1$mu)
})

test_that("em_gau is an alias of em_ad", {
  set.seed(128)
  y <- simulate_ad(n_subjects = 60, n_time = 5, order = 1, phi = c(0, rep(0.3, 4)))
  y[sample(length(y), 10)] <- NA

  fit1 <- em_ad(y, order = 1, max_iter = 15, tol = 1e-6, verbose = FALSE)
  fit2 <- em_gau(y, order = 1, max_iter = 15, tol = 1e-6, verbose = FALSE)
  expect_equal(fit2$log_l, fit1$log_l)
  expect_equal(fit2$em_iterations, fit1$em_iterations)
})

test_that("logL_gau is an alias of logL_ad", {
  set.seed(129)
  y <- simulate_ad(n_subjects = 40, n_time = 5, order = 1, phi = c(0, rep(0.4, 4)))
  fit <- fit_ad(y, order = 1)

  ll1 <- logL_ad(y, order = 1, mu = fit$mu, phi = fit$phi, sigma = fit$sigma)
  ll2 <- logL_gau(y, order = 1, mu = fit$mu, phi = fit$phi, sigma = fit$sigma)
  expect_equal(ll2, ll1)
})

test_that("simulate_gau is an alias of simulate_ad", {
  set.seed(130)
  y1 <- simulate_ad(n_subjects = 30, n_time = 5, order = 1, phi = c(0, rep(0.2, 4)))
  set.seed(130)
  y2 <- simulate_gau(n_subjects = 30, n_time = 5, order = 1, phi = c(0, rep(0.2, 4)))
  expect_equal(y2, y1)
})

test_that("Gaussian LRT aliases call AD implementations", {
  set.seed(131)
  y <- simulate_ad(n_subjects = 80, n_time = 6, order = 1, phi = c(0, rep(0.35, 5)))
  blocks <- rep(1:2, each = nrow(y) / 2)
  C <- diff(diag(ncol(y)))
  mu0 <- colMeans(y)

  o1 <- lrt_order_ad(y, p = 0, q = 1)
  o2 <- lrt_order_gau(y, p = 0, q = 1)
  expect_equal(o2$statistic, o1$statistic)

  s1 <- lrt_one_sample_ad(y, mu0 = mu0, p = 1)
  s2 <- lrt_one_sample_gau(y, mu0 = mu0, p = 1)
  expect_equal(s2$statistic, s1$statistic)

  t1 <- lrt_two_sample_ad(y, blocks = blocks, p = 1)
  t2 <- lrt_two_sample_gau(y, blocks = blocks, p = 1)
  expect_equal(t2$statistic, t1$statistic)

  c1 <- lrt_contrast_ad(y, C = C, p = 1)
  c2 <- lrt_contrast_gau(y, C = C, p = 1)
  expect_equal(c2$statistic, c1$statistic)

  h1 <- lrt_homogeneity_ad(y, blocks = blocks, p = 1)
  h2 <- lrt_homogeneity_gau(y, blocks = blocks, p = 1)
  expect_equal(h2$statistic, h1$statistic)
})
