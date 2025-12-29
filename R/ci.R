#' Confidence intervals for fitted INAD models
#'
#' Computes confidence intervals for selected parameters from a fitted INAD model.
#' For the fixed effect case with negative binomial thinning and Bell innovations,
#' Wald intervals for time varying alpha and theta are computed via Louis identity.
#' For block effects tau, profile likelihood intervals are computed by fixing one
#' component of tau and re maximizing the log likelihood over nuisance parameters.
#' For negative binomial innovations, Wald intervals for the innovation size
#' parameter are computed using a one dimensional observed information approximation
#' per time point, holding other parameters fixed at their fitted values.
#'
#' @param y Integer matrix with \code{n_subjects} rows and \code{n_time} columns.
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#' @param blocks Optional integer vector of length \code{n_subjects}. Required for
#'   block effect intervals. If provided, should match \code{fit$settings$blocks}.
#' @param level Confidence level between 0 and 1.
#' @param idx_time Optional integer vector of time indices for which to compute
#'   intervals. Default is all time points.
#' @param ridge Nonnegative ridge value added to the observed information matrix
#'   used for Louis based Wald intervals.
#' @param profile_maxeval Maximum number of function evaluations used in the
#'   profile likelihood refits.
#' @param profile_xtol_rel Relative tolerance used in the profile likelihood refits.
#'
#' @return An object of class \code{inad_ci}, a list with elements
#'   \code{settings}, \code{level}, \code{alpha}, \code{theta}, \code{nb_inno_size},
#'   and \code{tau}. Each non NULL interval element is a data frame with columns
#'   \code{param}, \code{est}, \code{lower}, \code{upper}, and possibly \code{se}
#'   and \code{width}.
#'
#' @examples
#' \dontrun{
#' fit <- fit_inad(y, order = 1, thinning = "nbinom", innovation = "bell", blocks = blocks)
#' ci <- ci_inad(y, fit, blocks = blocks)
#' ci$alpha
#' ci$theta
#' ci$tau
#' }
#'
#' @export
ci_inad <- function(
        y,
        fit,
        blocks = NULL,
        level = 0.95,
        idx_time = NULL,
        ridge = 0,
        profile_maxeval = 2500,
        profile_xtol_rel = 1e-6
) {
    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) stop("y must be a nonnegative integer matrix")

    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")
    if (is.null(fit$settings$thinning)) stop("fit$settings$thinning is missing")
    if (is.null(fit$settings$innovation)) stop("fit$settings$innovation is missing")
    if (is.null(fit$theta)) stop("fit$theta is missing")
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    n <- nrow(y)
    N <- ncol(y)

    if (is.null(idx_time)) idx_time <- seq_len(N)
    idx_time <- as.integer(idx_time)
    if (any(idx_time < 1) || any(idx_time > N)) stop("idx_time out of range")

    has_blocks <- !is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 0L
    if (has_blocks) {
        blocks <- as.integer(blocks)
        if (length(blocks) != n) stop("blocks must have length nrow(y)")
        if (any(blocks < 1)) stop("blocks must be positive integers starting at 1")
    }

    out <- list(settings = fit$settings, level = level)

    if (has_blocks && thinning == "nbinom" && innovation == "bell" && ord %in% c(1, 2)) {
        tmp <- .ci_alpha_theta_louis_nbt_bell_fe(
            y = y,
            fit = fit,
            blocks = blocks,
            level = level,
            idx_time = idx_time,
            ridge = ridge
        )
        out$alpha <- tmp$alpha
        out$theta <- tmp$theta
    } else {
        out$alpha <- NULL
        out$theta <- NULL
    }

    if (innovation == "nbinom") {
        out$nb_inno_size <- .ci_nb_inno_size_obsinfo_1d(
            y = y,
            fit = fit,
            blocks = if (has_blocks) blocks else NULL,
            level = level,
            idx_time = idx_time
        )
    } else {
        out$nb_inno_size <- NULL
    }

    if (has_blocks) {
        out$tau <- .ci_tau_profile_inad(
            y = y,
            fit = fit,
            blocks = blocks,
            level = level,
            maxeval = profile_maxeval,
            xtol_rel = profile_xtol_rel
        )
    } else {
        out$tau <- NULL
    }

    class(out) <- "inad_ci"
    out
}

#' Print method for INAD confidence intervals
#'
#' @param x An object of class \code{inad_ci}.
#' @param ... Unused.
#'
#' @return The input object, invisibly.
#'
#' @export
print.inad_ci <- function(x, ...) {
    cat("CI level:", x$level, "\n")
    if (!is.null(x$alpha)) cat("alpha rows:", nrow(x$alpha), "\n")
    if (!is.null(x$theta)) cat("theta rows:", nrow(x$theta), "\n")
    if (!is.null(x$nb_inno_size)) cat("nb_inno_size rows:", nrow(x$nb_inno_size), "\n")
    if (!is.null(x$tau)) cat("tau rows:", nrow(x$tau), "\n")
    invisible(x)
}

#' @keywords internal
.lse <- function(v) {
    m <- max(v)
    m + log(sum(exp(v - m)))
}

#' @keywords internal
.poster_nb_bell <- function(m, n, alpha, lambda) {
    alpha <- max(alpha, 1e-8)
    lambda <- max(lambda, 1e-8)

    k <- 0:n
    p <- 1 / (1 + alpha)
    q <- 1 - p

    lthin <- suppressWarnings(dnbinom(k, size = m, prob = p, log = TRUE))
    linno <- sapply(n - k, function(z) log(dbell(z, theta = lambda)))
    lw <- lthin + linno

    lZ <- .lse(lw)
    w <- exp(lw - lZ)

    score_nb <- -m * p + k * (p^2 / q)
    h_nb <- -m * p^2 + k * (p^3 * (2 - p) / (q^2))

    z <- n - k
    score_la <- z / lambda - exp(lambda)
    h_la <- z / (lambda^2) + exp(lambda)

    E_sNB <- sum(w * score_nb)
    E_hNB <- sum(w * h_nb)
    E_sLA <- sum(w * score_la)
    E_hLA <- sum(w * h_la)

    E_sNB2 <- sum(w * (score_nb^2))
    E_sLA2 <- sum(w * (score_la^2))
    E_sNBsLA <- sum(w * (score_nb * score_la))

    Var_sNB <- max(0, E_sNB2 - E_sNB^2)
    Var_sLA <- max(0, E_sLA2 - E_sLA^2)
    Cov_sNB_sLA <- E_sNBsLA - E_sNB * E_sLA

    c(
        E_hNB = E_hNB,
        Var_sNB = Var_sNB,
        E_hLA = E_hLA,
        Var_sLA = Var_sLA,
        Cov_sNB_sLA = Cov_sNB_sLA
    )
}

#' @keywords internal
.louis_info_i_nbt_bell_fe <- function(y, i, alpha_hat, theta_hat, tau_hat, blocks) {
    n <- nrow(y)

    if (i == 1) {
        Itt <- 0
        for (s in 1:n) {
            b <- blocks[s]
            z <- y[s, 1]
            la <- max(theta_hat + tau_hat[b], 1e-8)
            Itt <- Itt + z / (la^2) + exp(la)
        }
        mat <- matrix(Itt, 1, 1)
        dimnames(mat) <- list("theta[1]", "theta[1]")
        return(mat)
    }

    m <- y[, i - 1]
    nvec <- y[, i]

    a <- max(alpha_hat, 1e-8)
    th <- max(theta_hat, 1e-8)

    E_hNB_sum <- 0
    Var_sNB_sum <- 0
    E_hLA_theta <- 0
    Var_sLA_theta <- 0
    Cov_theta <- 0

    for (s in 1:n) {
        b <- blocks[s]
        la <- max(th + tau_hat[b], 1e-8)
        acc <- .poster_nb_bell(m = m[s], n = nvec[s], alpha = a, lambda = la)

        E_hNB_sum <- E_hNB_sum + acc["E_hNB"]
        Var_sNB_sum <- Var_sNB_sum + acc["Var_sNB"]
        E_hLA_theta <- E_hLA_theta + acc["E_hLA"]
        Var_sLA_theta <- Var_sLA_theta + acc["Var_sLA"]
        Cov_theta <- Cov_theta + acc["Cov_sNB_sLA"]
    }

    Iaa <- E_hNB_sum - Var_sNB_sum
    Itt <- E_hLA_theta - Var_sLA_theta
    Iat <- -Cov_theta

    mat <- matrix(c(Iaa, Iat, Iat, Itt), 2, 2, byrow = TRUE)
    par_names <- c(paste0("alpha[", i, "]"), paste0("theta[", i, "]"))
    dimnames(mat) <- list(par_names, par_names)
    mat
}

#' @keywords internal
.ci_wald_i_louis_nbt_bell_fe <- function(y, i, alpha_hat, theta_hat, tau_hat, blocks, level = 0.95, ridge = 0) {
    Iobs <- .louis_info_i_nbt_bell_fe(y, i, alpha_hat, theta_hat, tau_hat, blocks)
    if (ridge > 0) Iobs <- Iobs + diag(ridge, nrow(Iobs))

    V <- tryCatch(solve(Iobs), error = function(e) NULL)
    if (is.null(V)) stop("Singular information at i = ", i, ".")

    z <- qnorm(1 - (1 - level) / 2)

    if (i == 1) {
        est <- c(theta_hat)
        nms <- c("theta[1]")
    } else {
        est <- c(alpha_hat, theta_hat)
        nms <- c(paste0("alpha[", i, "]"), paste0("theta[", i, "]"))
    }

    se <- sqrt(pmax(diag(V), 0))
    L <- est - z * se
    U <- est + z * se

    if (i != 1) {
        L[1] <- max(L[1], 0)
        U[1] <- max(U[1], 0)
    }

    data.frame(
        param = nms,
        est = est,
        se = se,
        lower = L,
        upper = U,
        width = U - L,
        level = level,
        row.names = NULL
    )
}

#' @keywords internal
.ci_alpha_theta_louis_nbt_bell_fe <- function(y, fit, blocks, level, idx_time, ridge) {
    alpha_rows <- list()
    theta_rows <- list()

    for (i in idx_time) {
        if (i == 1) {
            tmp <- .ci_wald_i_louis_nbt_bell_fe(
                y = y,
                i = 1,
                alpha_hat = NA,
                theta_hat = fit$theta[1],
                tau_hat = fit$tau,
                blocks = blocks,
                level = level,
                ridge = ridge
            )
            theta_rows[[length(theta_rows) + 1L]] <- tmp
        } else {
            alpha_i <- if (is.matrix(fit$alpha)) fit$alpha[i, 1] else fit$alpha[i]
            tmp <- .ci_wald_i_louis_nbt_bell_fe(
                y = y,
                i = i,
                alpha_hat = alpha_i,
                theta_hat = fit$theta[i],
                tau_hat = fit$tau,
                blocks = blocks,
                level = level,
                ridge = ridge
            )
            alpha_rows[[length(alpha_rows) + 1L]] <- tmp[tmp$param == paste0("alpha[", i, "]"), , drop = FALSE]
            theta_rows[[length(theta_rows) + 1L]] <- tmp[tmp$param == paste0("theta[", i, "]"), , drop = FALSE]
        }
    }

    list(
        alpha = if (length(alpha_rows) > 0) do.call(rbind, alpha_rows) else NULL,
        theta = if (length(theta_rows) > 0) do.call(rbind, theta_rows) else NULL
    )
}

#' @keywords internal
.ci_nb_inno_size_obsinfo_1d <- function(y, fit, blocks, level, idx_time) {
    eps_pos <- 1e-10
    z <- qnorm(1 - (1 - level) / 2)

    if (is.null(fit$nb_inno_size)) return(NULL)

    nb <- as.numeric(fit$nb_inno_size)
    if (length(nb) == 1L) nb <- rep(nb, ncol(y))
    if (length(nb) != ncol(y)) stop("fit$nb_inno_size must be length 1 or ncol(y)")

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    if (innovation != "nbinom") return(NULL)

    res <- vector("list", length(idx_time))

    for (jj in seq_along(idx_time)) {
        i <- idx_time[jj]
        s0 <- max(nb[i], eps_pos)

        ll_at <- function(sz) {
            sz <- max(sz, eps_pos)
            nb_try <- nb
            nb_try[i] <- sz

            ll <- tryCatch(
                logL_inad(
                    y = y,
                    order = ord,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = fit$alpha,
                    theta = fit$theta,
                    nb_inno_size = nb_try,
                    blocks = blocks,
                    tau = if (!is.null(blocks)) fit$tau else NULL
                ),
                error = function(e) NA_real_
            )
            ll
        }

        h <- max(1e-4, 1e-3 * s0)
        f0 <- ll_at(s0)
        fp <- ll_at(s0 + h)
        fm <- ll_at(max(eps_pos, s0 - h))

        if (!is.finite(f0) || !is.finite(fp) || !is.finite(fm)) {
            res[[jj]] <- data.frame(
                param = paste0("nb_inno_size[", i, "]"),
                est = s0,
                se = NA_real_,
                lower = NA_real_,
                upper = NA_real_,
                width = NA_real_,
                level = level,
                row.names = NULL
            )
            next
        }

        d2 <- (fp - 2 * f0 + fm) / (h^2)
        I <- max(0, -d2)

        if (I <= 0 || !is.finite(I)) {
            se <- NA_real_
            L <- NA_real_
            U <- NA_real_
        } else {
            se <- sqrt(1 / I)
            L <- max(0, s0 - z * se)
            U <- max(0, s0 + z * se)
        }

        res[[jj]] <- data.frame(
            param = paste0("nb_inno_size[", i, "]"),
            est = s0,
            se = se,
            lower = L,
            upper = U,
            width = if (is.finite(L) && is.finite(U)) (U - L) else NA_real_,
            level = level,
            row.names = NULL
        )
    }

    do.call(rbind, res)
}

#' @keywords internal
.ci_tau_profile_inad <- function(
        y,
        fit,
        blocks,
        level,
        maxeval,
        xtol_rel,
        expand = 2,
        max_bracket_iter = 20
) {
    if (!requireNamespace("nloptr", quietly = TRUE)) stop("nloptr is required for tau profiling")

    if (is.null(fit$tau) || length(fit$tau) < 2) stop("fit$tau must exist with at least two blocks")
    if (is.null(fit$alpha)) stop("fit$alpha is missing for profiling refits")
    if (is.null(fit$theta)) stop("fit$theta is missing for profiling refits")

    n <- nrow(y)
    N <- ncol(y)
    B <- max(blocks)

    if (length(fit$tau) != B) stop("length(fit$tau) must equal max(blocks)")
    if (length(fit$theta) != N) stop("length(fit$theta) must equal ncol(y)")

    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation

    eps_pos <- 1e-10
    eps_alpha <- 1e-10
    alpha_ub <- 1 - eps_alpha

    lambda_ok <- function(theta_vec, tau_vec) {
        lam <- outer(rep(1, n), theta_vec) + matrix(tau_vec[blocks], n, N)
        all(is.finite(lam)) && all(lam > eps_pos)
    }

    pack_par <- function(alpha, theta, tau, nb_inno_size, ord, B, N, innovation, tau_fix_idx) {
        out <- numeric(0)
        tau_free_idx <- 2:B
        if (length(tau_fix_idx) > 0) tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
        if (length(tau_free_idx) > 0) out <- c(out, tau[tau_free_idx])

        if (ord == 1) {
            if (N >= 2) out <- c(out, alpha[2:N])
        }
        if (ord == 2) {
            if (N >= 2) out <- c(out, alpha[2:N, 1])
            if (N >= 3) out <- c(out, alpha[3:N, 2])
        }

        out <- c(out, theta)

        if (innovation == "nbinom") out <- c(out, nb_inno_size)

        out
    }

    unpack_par <- function(par, alpha0, theta0, tau0, nb0, ord, B, N, innovation, tau_fix_idx, tau_fix_val) {
        k <- 0

        tau <- tau0
        tau[1] <- 0
        tau_free_idx <- 2:B
        if (length(tau_fix_idx) > 0) tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
        if (length(tau_free_idx) > 0) {
            tau[tau_free_idx] <- par[(k + 1):(k + length(tau_free_idx))]
            k <- k + length(tau_free_idx)
        }
        if (length(tau_fix_idx) > 0) {
            tau[tau_fix_idx] <- tau_fix_val
            tau[1] <- 0
        }

        alpha <- alpha0
        if (ord == 1) {
            if (N >= 2) {
                alpha[2:N] <- par[(k + 1):(k + (N - 1))]
                k <- k + (N - 1)
            }
            alpha[1] <- 0
            alpha <- pmin(pmax(alpha, 0), alpha_ub)
            alpha[1] <- 0
        }
        if (ord == 2) {
            if (N >= 2) {
                alpha[2:N, 1] <- par[(k + 1):(k + (N - 1))]
                k <- k + (N - 1)
            }
            if (N >= 3) {
                alpha[3:N, 2] <- par[(k + 1):(k + (N - 2))]
                k <- k + (N - 2)
            }
            alpha[1, 1] <- 0
            alpha[1, 2] <- 0
            if (N >= 2) alpha[2, 2] <- 0
            alpha <- pmin(pmax(alpha, 0), alpha_ub)
            alpha[1, 1] <- 0
            alpha[1, 2] <- 0
            if (N >= 2) alpha[2, 2] <- 0
        }

        theta <- par[(k + 1):(k + N)]
        k <- k + N
        theta <- pmax(theta, eps_pos)

        nb_inno_size <- nb0
        if (innovation == "nbinom") {
            nb_inno_size <- par[(k + 1):(k + N)]
            k <- k + N
            nb_inno_size <- pmax(nb_inno_size, eps_pos)
        }

        list(alpha = alpha, theta = theta, tau = tau, nb_inno_size = nb_inno_size)
    }

    refit_tau_fixed <- function(tau_fix_idx, tau_fix_val) {
        alpha0 <- fit$alpha
        theta0 <- pmax(as.numeric(fit$theta), eps_pos)

        tau0 <- as.numeric(fit$tau)
        tau0[1] <- 0

        nb0 <- fit$nb_inno_size
        if (innovation == "nbinom") {
            if (is.null(nb0)) nb0 <- rep(1, N)
            nb0 <- as.numeric(nb0)
            if (length(nb0) == 1L) nb0 <- rep(nb0, N)
            if (length(nb0) != N) stop("fit$nb_inno_size must be length 1 or ncol(y)")
            nb0 <- pmax(nb0, eps_pos)
        } else {
            nb0 <- NULL
        }

        par0 <- pack_par(alpha0, theta0, tau0, nb0, ord, B, N, innovation, tau_fix_idx)

        obj <- function(par) {
            cur <- unpack_par(
                par = par,
                alpha0 = alpha0,
                theta0 = theta0,
                tau0 = tau0,
                nb0 = nb0,
                ord = ord,
                B = B,
                N = N,
                innovation = innovation,
                tau_fix_idx = tau_fix_idx,
                tau_fix_val = tau_fix_val
            )
            if (!lambda_ok(cur$theta, cur$tau)) return(1e8)

            ll <- tryCatch(
                logL_inad(
                    y = y,
                    order = ord,
                    thinning = thinning,
                    innovation = innovation,
                    alpha = cur$alpha,
                    theta = cur$theta,
                    nb_inno_size = cur$nb_inno_size,
                    blocks = blocks,
                    tau = cur$tau
                ),
                error = function(e) NA_real_
            )
            if (!is.finite(ll)) return(1e8)
            -ll
        }

        opt <- nloptr::nloptr(
            x0 = par0,
            eval_f = obj,
            opts = list(
                algorithm = "NLOPT_LN_BOBYQA",
                xtol_rel = xtol_rel,
                maxeval = maxeval
            )
        )

        -opt$objective
    }

    prof_target <- function(b, tval) {
        ll_prof <- refit_tau_fixed(tau_fix_idx = b, tau_fix_val = tval)
        2 * (fit$log_l - ll_prof) - qchisq(level, df = 1)
    }

    bracket_one_side <- function(b, direction, step0) {
        tau_mle <- fit$tau[b]
        lower_bound <- -min(fit$theta)
        upper_bound <- max(abs(tau_mle) + 1, 1)

        x1 <- tau_mle
        step <- step0 * direction

        for (k in 1:max_bracket_iter) {
            x2 <- x1 + step
            x2 <- max(lower_bound, min(upper_bound, x2))
            f2 <- prof_target(b, x2)
            if (is.finite(f2) && f2 >= 0) return(sort(c(x1, x2)))
            if (x2 == lower_bound || x2 == upper_bound) break
            x1 <- x2
            step <- step * expand
        }

        NULL
    }

    tau_ci <- vector("list", length = B - 1)
    for (b in 2:B) {
        tau_mle <- fit$tau[b]
        lower_bound <- -min(fit$theta)
        upper_bound <- max(abs(tau_mle) + 1, 1)
        step0 <- (upper_bound - lower_bound) / 40
        if (!is.finite(step0) || step0 <= 0) step0 <- 0.1

        f_mle <- prof_target(b, tau_mle)
        if (!is.finite(f_mle) || f_mle > 1e-6) stop("Profile target at MLE not near zero for tau[", b, "]")

        left_int <- bracket_one_side(b, direction = -1, step0 = step0)
        right_int <- bracket_one_side(b, direction = 1, step0 = step0)

        left <- NA_real_
        right <- NA_real_

        if (!is.null(left_int)) left <- uniroot(function(x) prof_target(b, x), interval = left_int)$root
        if (!is.null(right_int)) right <- uniroot(function(x) prof_target(b, x), interval = right_int)$root

        tau_ci[[b - 1]] <- data.frame(
            param = paste0("tau[", b, "]"),
            est = tau_mle,
            lower = left,
            upper = right,
            level = level,
            row.names = NULL
        )
    }

    do.call(rbind, tau_ci)
}
