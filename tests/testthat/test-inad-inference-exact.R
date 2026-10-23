test_that("chi-bar-square p-value treats chi^2_0 as a point mass at zero", {
    r <- 10
    stat <- 25.0687
    expected <- sum(choose(r, 1:r) * 2^-r * stats::pchisq(stat, 1:r, lower.tail = FALSE))
    expect_equal(antedep:::.chi_bar_pvalue(stat, r), expected, tolerance = 1e-12)
    expect_equal(antedep:::.chi_bar_pvalue(0, r), 1)
    expect_equal(antedep:::.chi_bar_pvalue(1, 0), 0)
})

test_that("parameter counts agree between fits and model tables", {
    expect_equal(antedep:::.count_params_inad(1, 12, 2, "nbinom", "nbinom"), 11 + 12 + 1 + 12)
    expect_equal(antedep:::.count_params_inad(2, 12, 2, "nbinom", "nbinom"), 21 + 12 + 1 + 12)
    cinfo <- antedep:::.parse_constraint(c("theta", "kappa"), 1, 12, "nbinom")
    expect_equal(cinfo$df, 22)
    expect_equal(antedep:::.parse_constraint("all", 1, 12, "nbinom")$df, 32)
    expect_equal(antedep:::.parse_constraint("alpha1", 2, 12, "pois")$df, 10)
    expect_error(antedep:::.parse_constraint("kappa", 1, 12, "pois"), "nbinom")
})

test_that("Bell parameter inversion is accurate to machine precision", {
    mu <- c(0.05, 0.7, 3, 25)
    th <- antedep:::.bell_theta_from_mean(mu)
    expect_lt(max(abs(th * exp(th) - mu) / mu), 1e-12)
})

test_that("Louis information equals the numerical Hessian, with blocks and constraints", {
    skip_on_cran()
    skip_if_not_installed("numDeriv")
    set.seed(42)
    n_time <- 4
    cases <- list(
        list(order = 1, thinning = "nbinom", innovation = "nbinom", constrain = "alpha"),
        list(order = 1, thinning = "pois", innovation = "bell", constrain = NULL),
        list(order = 2, thinning = "binom", innovation = "nbinom", constrain = NULL)
    )
    for (cs in cases) {
        alpha <- if (cs$order == 1) c(0, 0.3, 0.4, 0.35) else
            rbind(c(0, 0.3, 0.3, 0.35), c(0, 0, 0.2, 0.15))
        theta <- if (cs$innovation == "bell") rep(0.8, n_time) else c(3, 2, 2.5, 2)
        blocks <- rep(1:2, each = 30)
        y <- simulate_inad(60, n_time, order = cs$order, thinning = cs$thinning,
                           innovation = cs$innovation, alpha = alpha, theta = theta,
                           nb_inno_size = 2, blocks = blocks, tau = 0.6)
        fit <- fit_inad(y, cs$order, cs$thinning, cs$innovation, blocks = blocks,
                        max_iter = 50, tol = 1e-8, constrain = cs$constrain)
        li <- antedep:::.louis_full_info(y, fit, blocks)
        ni <- antedep:::.inad_numeric_info(y, fit, blocks, step = 1e-3)
        # Numerical differentiation on the log scale is unreliable for estimates
        # very close to zero, so compare the parameters away from zero.
        ok <- li$fp$tab$est > 0.01 | li$fp$tab$group == "tau"
        expect_lt(max(abs(li$info[ok, ok] - ni$info[ok, ok])) / max(abs(ni$info[ok, ok])), 1e-3)
        ci <- ci_inad(y, fit, blocks = blocks)
        expect_s3_class(ci, "inad_ci")
        expect_true(all(ci$theta$lower < ci$theta$est & ci$theta$est < ci$theta$upper))
    }
})

test_that("constrained fits re-estimate NB innovation sizes and reach the MLE", {
    skip_on_cran()
    set.seed(7)
    blocks <- rep(1:2, each = 40)
    y <- simulate_inad(80, 5, order = 1, thinning = "nbinom", innovation = "nbinom",
                       alpha = c(0, 0.4, 0.4, 0.4, 0.4), theta = c(4, 2, 3, 2, 2.5),
                       nb_inno_size = c(1, 3, 0.8, 2, 1.5), blocks = blocks, tau = 0.8)
    fit_u <- fit_inad(y, 1, "nbinom", "nbinom", blocks = blocks, max_iter = 100, tol = 1e-8)
    fit_c <- fit_inad(y, 1, "nbinom", "nbinom", blocks = blocks, max_iter = 100, tol = 1e-8,
                      constrain = "alpha")
    expect_equal(fit_c$n_params, fit_u$n_params - 3)
    expect_false(isTRUE(all.equal(fit_c$nb_inno_size, fit_u$nb_inno_size)))
    expect_lte(fit_c$log_l, fit_u$log_l + 1e-6)
    li <- antedep:::.louis_full_info(y, fit_c, blocks)
    expect_lt(max(abs(li$score)), 1e-5)
    st <- test_stationarity_inad(y, 1, "nbinom", "nbinom", blocks = blocks,
                                 constrain = "innovation", fit_unconstrained = fit_u)
    expect_equal(st$df, 8)
    expect_equal(st$table$n_params[2], fit_u$n_params - 8)
})

test_that("EM is exact for negative binomial innovations and matches direct maximization", {
    skip_on_cran()
    set.seed(3)
    blocks <- rep(1:2, each = 40)
    y <- simulate_inad(80, 4, order = 1, thinning = "pois", innovation = "nbinom",
                       alpha = c(0, 0.5, 0.4, 0.6), theta = c(3, 2, 2, 2.5),
                       nb_inno_size = 1.5, blocks = blocks, tau = 0.5)
    em <- em_inad(y, 1, "pois", "nbinom", blocks = blocks, max_iter = 5000, tol = 1e-10)
    direct <- fit_inad(y, 1, "pois", "nbinom", blocks = blocks, max_iter = 200, tol = 1e-10)
    expect_true(all(diff(em$loglik_trace) > -1e-8))
    expect_equal(em$log_l, direct$log_l, tolerance = 1e-4)
    expect_equal(em$log_l,
                 logL_inad(y, 1, "pois", "nbinom", em$alpha, em$theta, em$nb_inno_size,
                           blocks, em$tau),
                 tolerance = 1e-10)
    expect_equal(length(em$nb_inno_size), 4)
})
