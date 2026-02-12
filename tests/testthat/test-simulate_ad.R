testthat::test_that("simulate_ad returns a numeric matrix with correct dimensions", {
    testthat::skip_if_not(exists("simulate_ad"))

    set.seed(1)
    y <- simulate_ad(n_subjects = 10, n_time = 5, order = 1)

    testthat::expect_true(is.matrix(y))
    testthat::expect_type(y, "double")
    testthat::expect_equal(dim(y), c(10, 5))
    testthat::expect_true(all(is.finite(y)))
})

testthat::test_that("simulate_ad is reproducible given the same seed", {
    testthat::skip_if_not(exists("simulate_ad"))

    set.seed(123)
    y1 <- simulate_ad(20, 6, order = 2)

    set.seed(123)
    y2 <- simulate_ad(20, 6, order = 2)

    testthat::expect_equal(y1, y2)

    y3 <- simulate_ad(20, 6, order = 2, seed = 999)
    y4 <- simulate_ad(20, 6, order = 2, seed = 999)
    testthat::expect_equal(y3, y4)
})

testthat::test_that("simulate_ad works for order 0, 1, 2", {
    testthat::skip_if_not(exists("simulate_ad"))

    set.seed(2)
    y0 <- simulate_ad(8, 4, order = 0)
    y1 <- simulate_ad(8, 4, order = 1)
    y2 <- simulate_ad(8, 4, order = 2)

    testthat::expect_equal(dim(y0), c(8, 4))
    testthat::expect_equal(dim(y1), c(8, 4))
    testthat::expect_equal(dim(y2), c(8, 4))
})

testthat::test_that("simulate_ad validates key inputs", {
    testthat::skip_if_not(exists("simulate_ad"))

    testthat::expect_error(simulate_ad(10, 5, order = 3))
    testthat::expect_error(simulate_ad(0, 5))
    testthat::expect_error(simulate_ad(10, 0))

    testthat::expect_error(simulate_ad(10, 5, mu = 1:3))
    testthat::expect_error(simulate_ad(10, 5, sigma = 0))
    testthat::expect_error(simulate_ad(10, 5, sigma = c(1, 1, -1, 1, 1)))

    testthat::expect_error(simulate_ad(10, 5, blocks = 1:9))
    y_nonseq_ok <- simulate_ad(10, 5, blocks = c(rep(1, 9), 0))
    testthat::expect_equal(dim(y_nonseq_ok), c(10, 5))
})

testthat::test_that("simulate_ad enforces phi shape rules for order 2", {
    testthat::skip_if_not(exists("simulate_ad"))

    testthat::expect_error(simulate_ad(10, 6, order = 2, phi = rep(0.2, 6)))
    testthat::expect_error(simulate_ad(10, 6, order = 2, phi = matrix(0, nrow = 2, ncol = 5)))

    phi_ok <- matrix(0, nrow = 2, ncol = 6)
    phi_ok[1, 2:6] <- 0.4
    phi_ok[2, 3:6] <- 0.1

    set.seed(7)
    y <- simulate_ad(10, 6, order = 2, phi = phi_ok)
    testthat::expect_equal(dim(y), c(10, 6))
})

testthat::test_that("simulate_ad handles blocks and tau as intended for order 0", {
    testthat::skip_if_not(exists("simulate_ad"))

    blocks <- rep(1:2, each = 5)

    set.seed(10)
    y_scalar_tau <- simulate_ad(
        n_subjects = 10,
        n_time = 3,
        order = 0,
        mu = 0,
        sigma = 1e-10,
        blocks = blocks,
        tau = 2
    )

    d1 <- mean(y_scalar_tau[blocks == 2, 1]) - mean(y_scalar_tau[blocks == 1, 1])
    testthat::expect_equal(d1, 2, tolerance = 1e-6)

    set.seed(11)
    y_vec_tau <- simulate_ad(
        n_subjects = 10,
        n_time = 3,
        order = 0,
        mu = 0,
        sigma = 1e-10,
        blocks = blocks,
        tau = c(999, 2)
    )

    m_block1 <- mean(y_vec_tau[blocks == 1, 1])
    m_block2 <- mean(y_vec_tau[blocks == 2, 1])

    testthat::expect_equal(m_block1, 0, tolerance = 1e-6)
    testthat::expect_equal(m_block2, 2, tolerance = 1e-6)

    testthat::expect_error(simulate_ad(
        n_subjects = 10,
        n_time = 3,
        order = 0,
        blocks = blocks,
        tau = c(0, 1, 2)
    ))

    blocks_nonseq <- c(rep(2, 5), rep(5, 5))
    set.seed(12)
    y_nonseq <- simulate_ad(
        n_subjects = 10,
        n_time = 3,
        order = 0,
        mu = 0,
        sigma = 1e-10,
        blocks = blocks_nonseq,
        tau = c(0, 2)
    )
    d_nonseq <- mean(y_nonseq[blocks_nonseq == 5, 1]) - mean(y_nonseq[blocks_nonseq == 2, 1])
    testthat::expect_equal(d_nonseq, 2, tolerance = 1e-6)
})
