test_that("simulate_inad basic structure and defaults", {
    set.seed(1)
    y <- simulate_inad(
        n_subjects = 5,
        n_time = 7
    )

    expect_true(is.matrix(y))
    expect_equal(dim(y), c(5L, 7L))
    expect_true(is.integer(y))
    expect_true(all(y >= 0))
})

test_that("order 0 ignores thinning choice", {
    set.seed(1)
    y1 <- simulate_inad(
        n_subjects = 10,
        n_time = 6,
        order = 0,
        thinning = "binom",
        innovation = "pois"
    )

    set.seed(1)
    y2 <- simulate_inad(
        n_subjects = 10,
        n_time = 6,
        order = 0,
        thinning = "nbinom",
        innovation = "pois"
    )

    expect_equal(y1, y2)
})

test_that("order 0 matches order 1 with zero alpha", {
    n_sub <- 8
    n_time <- 6

    set.seed(2)
    y0 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 0,
        thinning = "binom",
        innovation = "pois",
        theta = 5
    )

    zero_alpha <- rep(0, n_time)

    set.seed(2)
    y1 <- simulate_inad(
        n_subjects = n_sub,
        n_time = n_time,
        order = 1,
        thinning = "binom",
        innovation = "pois",
        alpha = zero_alpha,
        theta = 5
    )

    expect_equal(y0, y1)
})

test_that("order 2 runs and produces nonnegative integers", {
    set.seed(3)
    y <- simulate_inad(
        n_subjects = 4,
        n_time = 8,
        order = 2,
        thinning = "binom",
        innovation = "bell"
    )

    expect_true(is.matrix(y))
    expect_equal(dim(y), c(4L, 8L))
    expect_true(is.integer(y))
    expect_true(all(y >= 0))
})

test_that("negative binomial innovations use size and prob", {
    set.seed(4)
    y1 <- simulate_inad(
        n_subjects = 30,
        n_time = 5,
        order = 1,
        thinning = "binom",
        innovation = "nbinom",
        theta = 10,
        nb_inno_size = 0.5
    )

    set.seed(4)
    y2 <- simulate_inad(
        n_subjects = 30,
        n_time = 5,
        order = 1,
        thinning = "binom",
        innovation = "nbinom",
        theta = 5,
        nb_inno_size = 0.8
    )

    expect_true(mean(y1) > mean(y2))
})
