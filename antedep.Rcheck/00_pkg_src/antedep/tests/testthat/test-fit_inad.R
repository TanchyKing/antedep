test_that("fit_inad works for bolus_inad without fixed effect", {
    skip_if_not(exists("bolus_inad", where = asNamespace("antedep"), inherits = FALSE) ||
                    exists("bolus_inad", where = parent.env(environment()), inherits = TRUE))

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "bell")

    expect_type(fit, "list")
    expect_true(is.null(fit$tau))
    expect_true(is.matrix(y))
    expect_equal(length(fit$alpha), ncol(y))
    expect_equal(length(fit$theta), ncol(y))

    expect_true(is.finite(fit$log_l))

    expect_equal(fit$alpha[1], 0, tolerance = 1e-12)
    expect_true(all(is.finite(fit$alpha)))
    expect_true(all(fit$alpha >= -1e-12))
    expect_true(all(fit$alpha < 1))

    expect_true(all(is.finite(fit$theta)))
    expect_true(all(fit$theta > 0))
})

test_that("fit_inad returns log_l equal to sum(loglik_i) in no FE case", {
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 1, thinning = "binom", innovation = "bell")

    expect_true(!is.null(fit$loglik_i))
    expect_equal(fit$log_l, sum(fit$loglik_i), tolerance = 1e-8)
})

test_that("fit_inad works for order 0 without fixed effect", {
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 0, thinning = "binom", innovation = "bell")

    expect_true(is.null(fit$alpha) || length(fit$alpha) == 0)
    expect_equal(length(fit$theta), ncol(y))
    expect_true(is.finite(fit$log_l))
    expect_true(all(is.finite(fit$theta)))
    expect_true(all(fit$theta > 0))
})

test_that("fit_inad handles order 2 structural zeros without fixed effect", {
    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(y, order = 2, thinning = "binom", innovation = "bell")

    expect_true(is.matrix(fit$alpha))
    expect_equal(dim(fit$alpha), c(ncol(y), 2))

    expect_equal(fit$alpha[1, 1], 0, tolerance = 1e-12)
    expect_equal(fit$alpha[1, 2], 0, tolerance = 1e-12)
    expect_equal(fit$alpha[2, 2], 0, tolerance = 1e-12)

    expect_true(all(is.finite(fit$alpha)))
    expect_true(all(fit$alpha >= -1e-12))
    expect_true(all(fit$alpha < 1))

    expect_true(is.finite(fit$log_l))
})

test_that("fit_inad works with fixed effect and tau normalization", {
    skip_if_not_installed("nloptr")

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y
    blocks <- bolus_inad$blocks

    fit1 <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "bell",
        blocks = blocks,
        init_tau = 0.4,
        max_iter = 5
    )

    expect_true(is.finite(fit1$log_l))
    expect_true(!is.null(fit1$tau))
    expect_equal(fit1$tau[1], 0, tolerance = 1e-12)

    B <- max(as.integer(blocks))
    expect_equal(length(fit1$tau), B)

    fit2 <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "bell",
        blocks = blocks,
        init_tau = c(999, rep(0.2, B - 1)),
        max_iter = 5
    )

    expect_equal(fit2$tau[1], 0, tolerance = 1e-12)
})

test_that("fit_inad works for innovation nbinom and returns nb_inno_size", {
    skip_if_not_installed("nloptr")

    data("bolus_inad", package = "antedep", envir = environment())
    y <- bolus_inad$y

    fit <- fit_inad(
        y,
        order = 1,
        thinning = "binom",
        innovation = "nbinom",
        init_nb_inno_size = 1
    )

    expect_true(!is.null(fit$nb_inno_size))
    expect_equal(length(fit$nb_inno_size), ncol(y))
    expect_true(all(is.finite(fit$nb_inno_size)))
    expect_true(all(fit$nb_inno_size > 0))
    expect_true(is.finite(fit$log_l))
})
