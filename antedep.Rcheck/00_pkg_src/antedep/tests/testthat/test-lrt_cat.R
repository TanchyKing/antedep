# test-lrt_cat.R - Tests for CAT hypothesis testing functions

# ============================================================
# Test lrt_order_cat
# ============================================================
test_that("lrt_order_cat works for AD(0) vs AD(1)", {
  set.seed(123)
  # Simulate AD(1) data - should reject AD(0)
  marg <- list(t1 = c(0.6, 0.4))
  trans <- list(
    t2 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE),
    t3 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE),
    t4 = matrix(c(0.8, 0.2, 0.3, 0.7), 2, byrow = TRUE)
  )
  y <- simulate_cat(300, 4, order = 1, n_categories = 2,
                    marginal = marg, transition = trans)
  
  test <- lrt_order_cat(y, order_null = 0, order_alt = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_true(test$p_value >= 0 && test$p_value <= 1)
  expect_equal(test$order_null, 0)
  expect_equal(test$order_alt, 1)
  
  # For data with true AD(1), should reject AD(0)
  expect_true(test$p_value < 0.05)
})


test_that("lrt_order_cat works for AD(1) vs AD(2)", {
  set.seed(456)
  # Simulate AD(1) data - should NOT reject AD(1) in favor of AD(2)
  y <- simulate_cat(200, 5, order = 1, n_categories = 2)
  
  test <- lrt_order_cat(y, order_null = 1, order_alt = 2)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_equal(test$order_null, 1)
  expect_equal(test$order_alt, 2)
  
  # For data with true AD(1), should usually NOT reject AD(1)
  # (but with random data this isn't guaranteed, so just check it runs)
})


test_that("lrt_order_cat works with pre-fitted models", {
  set.seed(789)
  y <- simulate_cat(150, 4, order = 1, n_categories = 2)
  
  fit0 <- fit_cat(y, order = 0)
  fit1 <- fit_cat(y, order = 1)
  
  # Test with pre-fitted models
  test <- lrt_order_cat(fit_null = fit0, fit_alt = fit1)
  
  expect_s3_class(test, "cat_lrt")
  expect_equal(test$order_null, 0)
  expect_equal(test$order_alt, 1)
  
  # Test with y + one pre-fitted model
  test2 <- lrt_order_cat(y, order_null = 0, fit_alt = fit1)
  expect_equal(test$lrt_stat, test2$lrt_stat, tolerance = 1e-10)
})


test_that("lrt_order_cat validates inputs correctly", {
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  
  # Error: order_alt <= order_null
  expect_error(lrt_order_cat(y, order_null = 1, order_alt = 0))
  expect_error(lrt_order_cat(y, order_null = 1, order_alt = 1))
  
  # Error: negative order
  expect_error(lrt_order_cat(y, order_null = -1, order_alt = 1))
  
  # Error: order too high
  expect_error(lrt_order_cat(y, order_null = 2, order_alt = 3))
  
  # Error: no data and no pre-fitted models
  expect_error(lrt_order_cat(order_null = 0, order_alt = 1))
})


test_that("lrt_order_cat df calculation is correct", {
  set.seed(111)
  y <- simulate_cat(100, 5, order = 1, n_categories = 2)
  
  # AD(0) vs AD(1) with c=2, n=5
  # df = (c-1) * c^0 * (n-1) = 1 * 1 * 4 = 4
  test_01 <- lrt_order_cat(y, order_null = 0, order_alt = 1)
  expect_equal(test_01$df, 4)
  
  # AD(1) vs AD(2) with c=2, n=5
  # df = (c-1) * c^1 * (n-2) = 1 * 2 * 3 = 6
  test_12 <- lrt_order_cat(y, order_null = 1, order_alt = 2)
  expect_equal(test_12$df, 6)
})


test_that("run_order_tests_cat works correctly", {
  set.seed(222)
  # Simulate AD(1) data
  marg <- list(t1 = c(0.6, 0.4))
  trans <- list(
    t2 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t3 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t4 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE),
    t5 = matrix(c(0.85, 0.15, 0.25, 0.75), 2, byrow = TRUE)
  )
  y <- simulate_cat(400, 5, order = 1, n_categories = 2,
                    marginal = marg, transition = trans)
  
  result <- run_order_tests_cat(y, max_order = 2)
  
  expect_true(is.list(result))
  expect_true("tests" %in% names(result))
  expect_true("table" %in% names(result))
  expect_true("fits" %in% names(result))
  expect_true("selected_order" %in% names(result))
  
  expect_equal(nrow(result$table), 2)  # Two tests: 0v1, 1v2
  expect_true(result$selected_order %in% 0:2)
  
  # For true AD(1), should select order 1
  # (0v1 significant, 1v2 not significant)
  expect_equal(result$selected_order, 1)
})


# ============================================================
# Test lrt_homogeneity_cat
# ============================================================
test_that("lrt_homogeneity_cat detects heterogeneity", {
  set.seed(333)
  # Create two groups with different transition probabilities
  marg1 <- list(t1 = c(0.8, 0.2))
  marg2 <- list(t1 = c(0.3, 0.7))
  trans1 <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    t3 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE)
  )
  trans2 <- list(
    t2 = matrix(c(0.4, 0.6, 0.6, 0.4), 2, byrow = TRUE),
    t3 = matrix(c(0.4, 0.6, 0.6, 0.4), 2, byrow = TRUE)
  )
  
  y1 <- simulate_cat(150, 3, order = 1, n_categories = 2,
                     marginal = marg1, transition = trans1)
  y2 <- simulate_cat(150, 3, order = 1, n_categories = 2,
                     marginal = marg2, transition = trans2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 150), rep(2, 150))
  
  test <- lrt_homogeneity_cat(y, blocks, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$n_groups, 2)
  
  # Should reject homogeneity
  expect_true(test$p_value < 0.05)
})


test_that("lrt_homogeneity_cat accepts homogeneous data", {
  set.seed(444)
  # Create two groups with SAME parameters
  y1 <- simulate_cat(100, 3, order = 1, n_categories = 2)
  y2 <- simulate_cat(100, 3, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 100), rep(2, 100))
  
  test <- lrt_homogeneity_cat(y, blocks, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  
  # Should NOT reject homogeneity (p > 0.05 usually, but not guaranteed)
  # Just check that test runs and gives reasonable output
  expect_true(test$p_value >= 0 && test$p_value <= 1)
})


test_that("lrt_homogeneity_cat works with pre-fitted models", {
  set.seed(555)
  y1 <- simulate_cat(80, 3, order = 1, n_categories = 2)
  y2 <- simulate_cat(80, 3, order = 1, n_categories = 2)
  y <- rbind(y1, y2)
  blocks <- c(rep(1, 80), rep(2, 80))
  
  fit_homo <- fit_cat(y, order = 1, blocks = blocks, homogeneous = TRUE)
  fit_hetero <- fit_cat(y, order = 1, blocks = blocks, homogeneous = FALSE)
  
  test <- lrt_homogeneity_cat(fit_null = fit_homo, fit_alt = fit_hetero)
  
  expect_s3_class(test, "cat_lrt")
  expect_equal(test$n_groups, 2)
})


test_that("lrt_homogeneity_cat validates inputs", {
  y <- simulate_cat(50, 3, order = 1, n_categories = 2)
  
  # Error: no blocks provided
  expect_error(lrt_homogeneity_cat(y, blocks = NULL))
  
  # Error: no y and no pre-fitted models
  expect_error(lrt_homogeneity_cat(blocks = c(1, 1, 2, 2)))
})


# ============================================================
# Test lrt_timeinvariance_cat
# ============================================================
test_that("lrt_timeinvariance_cat works", {
  set.seed(666)
  # Simulate with time-invariant transitions (default)
  y <- simulate_cat(200, 5, order = 1, n_categories = 2)
  
  test <- lrt_timeinvariance_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  expect_true(test$lrt_stat >= 0)
  expect_true(test$df > 0)
  expect_equal(test$order, 1)
  
  # For time-invariant data, should NOT reject time-invariance
  # (but this is statistical, so just check it runs)
  expect_true(test$p_value >= 0 && test$p_value <= 1)
})


test_that("lrt_timeinvariance_cat detects time-varying transitions", {
  set.seed(777)
  # Simulate with time-VARYING transitions
  marg <- list(t1 = c(0.5, 0.5))
  trans_varying <- list(
    t2 = matrix(c(0.9, 0.1, 0.1, 0.9), 2, byrow = TRUE),
    t3 = matrix(c(0.5, 0.5, 0.5, 0.5), 2, byrow = TRUE),
    t4 = matrix(c(0.2, 0.8, 0.8, 0.2), 2, byrow = TRUE),
    t5 = matrix(c(0.7, 0.3, 0.3, 0.7), 2, byrow = TRUE)
  )
  
  y <- simulate_cat(400, 5, order = 1, n_categories = 2,
                    marginal = marg, transition = trans_varying)
  
  test <- lrt_timeinvariance_cat(y, order = 1)
  
  expect_s3_class(test, "cat_lrt")
  
  # Should reject time-invariance
  expect_true(test$p_value < 0.05)
})


test_that("lrt_timeinvariance_cat validates inputs", {
  y <- simulate_cat(50, 4, order = 1, n_categories = 2)
  
  # Error: order < 1
  expect_error(lrt_timeinvariance_cat(y, order = 0))
})


# ============================================================
# Test print methods
# ============================================================
test_that("print.cat_lrt works", {
  set.seed(888)
  y <- simulate_cat(100, 4, order = 1, n_categories = 2)
  
  test <- lrt_order_cat(y, order_null = 0, order_alt = 1)
  
  # Just verify it doesn't error
  expect_output(print(test), "Likelihood Ratio Test")
})
