testthat::test_that("fit_ad returns correctly shaped outputs", {
    testthat::skip_if_not(exists("fit_ad"))
    testthat::skip_if_not(exists("simulate_ad"))
    testthat::skip_if_not(exists("logL_ad"))

    set.seed(1)
    y <- simulate_ad(
        n_subjects = 30,
        n_time = 6,
        order = 1,
        mu = seq(0, 0.5, length.out = 6),
        phi = 0.5,
        sigma = 1
    )

    fit <- fit_ad(y, order = 1, estimate_mu = FALSE)

    testthat::expect_true(is.list(fit))
    testthat::expect_true(is.numeric(fit$mu))
    testthat::expect_equal(length(fit$mu), ncol(y))

    testthat::expect_true(is.numeric(fit$sigma))
    testthat::expect_equal(length(fit$sigma), ncol(y))
    testthat::expect_true(all(is.finite(fit$sigma)))
    testthat::expect_true(all(fit$sigma > 0))

    testthat::expect_true(is.numeric(fit$phi))
    testthat::expect_equal(length(fit$phi), ncol(y))
    testthat::expect_equal(fit$phi[1], 0)

    testthat::expect_true(is.finite(fit$log_l))
    testthat::expect_true(is.finite(fit$aic))
    testthat::expect_true(is.finite(fit$bic))
})

testthat::test_that("fit_ad log_l equals logL_ad at returned parameters", {
    testthat::skip_if_not(exists("fit_ad"))
    testthat::skip_if_not(exists("simulate_ad"))
    testthat::skip_if_not(exists("logL_ad"))

    set.seed(2)
    y <- simulate_ad(
        n_subjects = 40,
        n_time = 5,
        order = 0,
        mu = 0.2,
        sigma = 1.3
    )

    fit <- fit_ad(y, order = 0, estimate_mu = TRUE)

    ll <- logL_ad(
        y = y,
        order = 0,
        mu = fit$mu,
        sigma = fit$sigma
    )

    testthat::expect_equal(fit$log_l, ll, tolerance = 1e-6)
})

testthat::test_that("fit_ad respects structural zeros for order 2", {
    testthat::skip_if_not(exists("fit_ad"))
    testthat::skip_if_not(exists("simulate_ad"))

    set.seed(3)
    phi_true <- matrix(0, nrow = 2, ncol = 6)
    phi_true[1, 2:6] <- 0.4
    phi_true[2, 3:6] <- -0.2

    y <- simulate_ad(
        n_subjects = 50,
        n_time = 6,
        order = 2,
        mu = 0,
        phi = phi_true,
        sigma = 1
    )

    fit <- fit_ad(y, order = 2, estimate_mu = FALSE)

    testthat::expect_true(is.matrix(fit$phi))
    testthat::expect_equal(dim(fit$phi), c(2, 6))

    testthat::expect_equal(fit$phi[1, 1], 0)
    testthat::expect_equal(fit$phi[2, 1], 0)
    testthat::expect_equal(fit$phi[2, 2], 0)
})

testthat::test_that("fit_ad handles blocks and returns tau with tau[1] = 0", {
    testthat::skip_if_not(exists("fit_ad"))
    testthat::skip_if_not(exists("simulate_ad"))

    set.seed(4)
    blocks <- rep(1:2, each = 30)

    y <- simulate_ad(
        n_subjects = 60,
        n_time = 5,
        order = 1,
        mu = 0,
        phi = 0.3,
        sigma = 1,
        blocks = blocks,
        tau = 0.7
    )

    fit <- fit_ad(
        y,
        order = 1,
        blocks = blocks,
        estimate_mu = TRUE,
        init_tau = 0.2
    )

    testthat::expect_true(is.numeric(fit$tau))
    testthat::expect_equal(length(fit$tau), 2)
    testthat::expect_equal(fit$tau[1], 0)
})

testthat::test_that("fit_ad estimates are closer to truth than a very wrong phi", {
    testthat::skip_if_not(exists("fit_ad"))
    testthat::skip_if_not(exists("simulate_ad"))
    testthat::skip_if_not(exists("logL_ad"))

    set.seed(5)
    y <- simulate_ad(
        n_subjects = 120,
        n_time = 6,
        order = 1,
        mu = 0,
        phi = 0.6,
        sigma = 1
    )

    fit <- fit_ad(y, order = 1, estimate_mu = FALSE)

    ll_fit <- fit$log_l
    ll_bad <- logL_ad(
        y = y,
        order = 1,
        mu = colMeans(y),
        phi = c(0, rep(-2, 5)),
        sigma = rep(1, 6)
    )

    testthat::expect_true(ll_fit > ll_bad)
})
