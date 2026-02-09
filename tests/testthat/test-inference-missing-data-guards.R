test_that("CAT inference utilities reject missing-data inference clearly", {
  y_miss <- matrix(
    c(1, 2, 1,
      2, 1, NA,
      1, 2, 2,
      2, 1, 1),
    nrow = 4,
    byrow = TRUE
  )

  fit_missing_ci <- structure(
    list(
      settings = list(
        order = 1,
        na_action = "marginalize",
        na_action_effective = "marginalize"
      ),
      cell_counts = NULL
    ),
    class = "cat_fit"
  )
  expect_error(ci_cat(fit_missing_ci), "complete-data fits only")

  expect_error(
    lrt_order_cat(y_miss, order_null = 0, order_alt = 1),
    "complete data only"
  )
  expect_error(
    lrt_homogeneity_cat(y_miss, blocks = c(1, 1, 2, 2), order = 1),
    "complete data only"
  )
  expect_error(lrt_timeinvariance_cat(y_miss, order = 1), "complete data only")
  expect_error(lrt_stationarity_cat(y_miss, order = 1), "complete data only")
  expect_error(run_order_tests_cat(y_miss), "complete data only")
  expect_error(run_stationarity_tests_cat(y_miss, order = 1), "complete data only")

  fit_null <- structure(
    list(
      settings = list(order = 0, homogeneous = TRUE, n_blocks = 1,
                      na_action_effective = "marginalize"),
      cell_counts = NULL,
      log_l = -10,
      n_params = 2,
      aic = 24,
      bic = 25
    ),
    class = "cat_fit"
  )
  fit_alt <- structure(
    list(
      settings = list(order = 1, homogeneous = TRUE, n_blocks = 1,
                      na_action_effective = "marginalize"),
      cell_counts = NULL,
      log_l = -9,
      n_params = 4,
      aic = 26,
      bic = 27
    ),
    class = "cat_fit"
  )
  expect_error(
    lrt_order_cat(fit_null = fit_null, fit_alt = fit_alt),
    "complete data only"
  )
})

test_that("INAD inference utilities reject missing-data inference clearly", {
  y_miss <- matrix(
    c(1, 2, 0,
      2, NA, 1,
      0, 1, 2,
      1, 0, 1),
    nrow = 4,
    byrow = TRUE
  )
  blocks <- c(1, 1, 2, 2)

  fit_stub <- list(
    settings = list(order = 1, thinning = "binom", innovation = "pois"),
    theta = c(1, 1, 1),
    log_l = -10
  )

  expect_error(ci_inad(y_miss, fit_stub, blocks = blocks), "complete data only")
  expect_error(
    lrt_order_inad(y_miss, order_null = 0, order_alt = 1, blocks = blocks),
    "complete data only"
  )
  expect_error(
    lrt_stationarity_inad(y_miss, order = 1, blocks = blocks),
    "complete data only"
  )
  expect_error(
    lrt_homogeneity_inad(y_miss, blocks = blocks, order = 1),
    "complete data only"
  )
  expect_error(
    run_stationarity_tests_inad(y_miss, order = 1, blocks = blocks, verbose = FALSE),
    "complete data only"
  )
  expect_error(
    run_homogeneity_tests_inad(y_miss, blocks = blocks, order = 1),
    "complete data only"
  )
})

test_that("AD likelihood-ratio utilities reject missing-data inference clearly", {
  y_miss <- matrix(
    c(1.0, 2.0, 3.0,
      2.0, 1.0, 0.0,
      0.5, NA, 1.5,
      3.0, 2.0, 1.0),
    nrow = 4,
    byrow = TRUE
  )
  blocks <- c(1, 1, 2, 2)

  expect_error(lrt_order_ad(y_miss, p = 0, q = 1), "complete data only")
  expect_error(lrt_homogeneity_ad(y_miss, blocks = blocks, p = 1), "complete data only")
  expect_error(
    lrt_one_sample_ad(y_miss, mu0 = c(0, 0, 0), p = 1),
    "complete data only"
  )
  expect_error(lrt_two_sample_ad(y_miss, blocks = blocks, p = 1), "complete data only")
  expect_error(
    lrt_contrast_ad(y_miss, C = matrix(c(1, -1, 0), nrow = 1), p = 1),
    "complete data only"
  )
})
