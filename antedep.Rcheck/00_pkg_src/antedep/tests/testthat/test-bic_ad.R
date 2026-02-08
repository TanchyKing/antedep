# File: tests/testthat/test-bic_ad.R

test_that("bic_ad matches manual formula for order 1 without blocks", {
    skip_if_not(exists("simulate_ad"))
    skip_if_not(exists("fit_ad"))

    set.seed(1)
    y <- simulate_ad(
        n_subjects = 80,
        n_time = 5,
        order = 1,
        mu = seq(0, 0.4, length.out = 5),
        sigma = rep(1, 5),
        phi = c(0, rep(0.4, 4))
    )

    fit <- fit_ad(y, order = 1, estimate_mu = TRUE)
    b1 <- bic_ad(fit, n_subjects = nrow(y))

    N <- ncol(y)
    k <- N + N + (N - 1)

    expect_equal(b1, -2 * fit$log_l + k * log(nrow(y)))
})

test_that("bic_ad counts tau and order 2 antedependence parameters", {
    N <- 6

    fit <- list(
        log_l = -100,
        settings = list(order = 2, estimate_mu = FALSE),
        sigma = rep(1, N),
        mu = rep(0, N),
        phi = matrix(0, nrow = 2, ncol = N),
        tau = c(0, 0.1, 0.2)
    )

    b <- bic_ad(fit, n_subjects = 50)

    k <- 0
    k <- k + N
    k <- k + (2 * N - 3)
    k <- k + (3 - 1)

    expect_equal(b, -2 * fit$log_l + k * log(50))
})
