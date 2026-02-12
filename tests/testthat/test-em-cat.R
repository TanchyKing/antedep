# test-em-cat.R - Tests for EM estimation in CAT models

test_that("em_cat returns closed-form fit when no missing data", {
  set.seed(1301)
  y <- simulate_cat(n_subjects = 120, n_time = 5, order = 1, n_categories = 2)

  fit_em <- em_cat(y, order = 1)
  fit_cf <- fit_cat(y, order = 1)

  expect_s3_class(fit_em, "cat_fit")
  expect_equal(as.numeric(fit_em$log_l), as.numeric(fit_cf$log_l), tolerance = 1e-10)
  expect_equal(fit_em$settings$method, "closed_form")
})


test_that("em_cat fits missing-data order 0 model", {
  set.seed(1302)
  y <- simulate_cat(n_subjects = 100, n_time = 6, order = 0, n_categories = 3)
  y[sample(length(y), 60)] <- NA

  fit <- em_cat(y, order = 0, max_iter = 80, tol = 1e-7)

  expect_s3_class(fit, "cat_fit")
  expect_true(is.finite(fit$log_l))
  expect_equal(fit$settings$order, 0)
  expect_true(all(vapply(fit$marginal, function(x) abs(sum(x) - 1) < 1e-8, logical(1))))
})


test_that("em_cat handles non-sequential block IDs and stores block levels", {
  set.seed(1303)
  y <- simulate_cat(n_subjects = 90, n_time = 5, order = 1, n_categories = 2)
  y[sample(length(y), 45)] <- NA
  blocks <- rep(c(2, 5, 2), length.out = nrow(y))

  fit <- em_cat(y, order = 1, blocks = blocks, homogeneous = FALSE, max_iter = 60)

  expect_equal(fit$settings$n_blocks, 2)
  expect_equal(sort(fit$settings$block_levels), c("2", "5"))
  expect_equal(length(fit$marginal), 2)
  expect_equal(length(fit$transition), 2)
})


test_that("em_cat order 1 does not report immediate convergence on first iteration", {
  set.seed(1304)
  y <- simulate_cat(n_subjects = 120, n_time = 5, order = 1, n_categories = 2)
  y[sample(length(y), 40)] <- NA

  fit <- em_cat(y, order = 1, max_iter = 50, tol = 1e-6)

  expect_true(fit$settings$n_iter > 1)
})


test_that("em_cat order 2 reports not implemented", {
  set.seed(1305)
  y <- simulate_cat(n_subjects = 40, n_time = 5, order = 2, n_categories = 2)
  y[sample(length(y), 10)] <- NA

  expect_error(
    em_cat(y, order = 2),
    "not yet implemented"
  )
})


test_that("em_cat final log-likelihood is consistent when max_iter is reached", {
  set.seed(1306)
  y <- simulate_cat(n_subjects = 140, n_time = 6, order = 1, n_categories = 2)
  y[sample(length(y), 90)] <- NA

  expect_warning(
    fit <- em_cat(y, order = 1, max_iter = 1, tol = 1e-12),
    "did not converge"
  )

  ll_check <- logL_cat(
    y = y,
    order = 1,
    marginal = fit$marginal,
    transition = fit$transition,
    n_categories = 2,
    na_action = "marginalize"
  )

  expect_equal(as.numeric(fit$log_l), as.numeric(ll_check), tolerance = 1e-8)
  expect_equal(
    as.numeric(fit$aic),
    as.numeric(-2 * fit$log_l + 2 * fit$n_params),
    tolerance = 1e-8
  )
  expect_equal(
    as.numeric(fit$bic),
    as.numeric(-2 * fit$log_l + log(nrow(y)) * fit$n_params),
    tolerance = 1e-8
  )
})
