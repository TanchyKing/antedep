#' INAD log likelihood (full data)
#'
#' If blocks is NULL, this computes the log likelihood as the sum of per time
#' contributions from logL_inad_i for computational convenience.
#'
#' @param y Integer matrix n_sub by n_time.
#' @param order Integer in {0, 1, 2}.
#' @param thinning One of "binom", "pois", "nbinom".
#' @param innovation One of "pois", "bell", "nbinom".
#' @param alpha Thinning parameters. For order 1, numeric length 1 or n_time.
#'   For order 2, either a matrix n_time by 2 or a list(alpha1, alpha2).
#' @param theta Innovation mean parameters. Numeric length 1 or n_time.
#' @param nb_inno_size Size parameter for innovation "nbinom". Numeric length 1 or n_time.
#' @param blocks Optional integer vector of length n_sub. If NULL, no fixed effect.
#' @param tau Optional numeric vector. Only used if blocks is not NULL.
#'
#' @return A scalar log likelihood.
#' @export
logL_inad <- function(
        y,
        order = 1,
        thinning = c("binom", "pois", "nbinom"),
        innovation = c("pois", "bell", "nbinom"),
        alpha,
        theta,
        nb_inno_size = NULL,
        blocks = NULL,
        tau = 0
) {
    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)

    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) return(-Inf)

    n <- nrow(y)
    N <- ncol(y)
    if (!(order %in% c(0, 1, 2))) stop("order must be 0, 1, or 2")

    if (is.null(blocks)) {
        return(sum(vapply(
            seq_len(N),
            function(i) logL_inad_i(
                y = y,
                i = i,
                order = order,
                thinning = thinning,
                innovation = innovation,
                alpha = alpha,
                theta = theta,
                nb_inno_size = nb_inno_size
            ),
            numeric(1)
        )))
    }

    if (length(blocks) != n) stop("length(blocks) must equal nrow(y)")
    blocks <- as.integer(blocks)

    B <- max(blocks)
    tau <- as.numeric(tau)

    if (B <= 1L) {
        tau <- 0
    } else if (length(tau) == 1L) {
        tau <- c(0, rep(tau, B - 1L))
    } else {
        if (length(tau) < B) stop("tau must be length 1 or at least max(blocks)")
        tau <- tau[seq_len(B)]
        tau[1] <- 0
    }

    theta <- as.numeric(theta)
    if (length(theta) == 1L) theta <- rep(theta, N)
    if (length(theta) != N) stop("theta must be length 1 or ncol(y)")

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) stop("nb_inno_size is required for innovation='nbinom'")
        nb_inno_size <- as.numeric(nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("nb_inno_size must be length 1 or ncol(y)")
        if (any(nb_inno_size <= 0)) return(-Inf)
    }

    unpack_alpha <- function(order, alpha, N) {
        if (order == 0) return(list(a1 = NULL, a2 = NULL))

        if (order == 1) {
            a <- as.numeric(alpha)
            if (length(a) == 1L) a <- rep(a, N)
            if (length(a) != N) stop("alpha must be length 1 or ncol(y) for order=1")
            return(list(a1 = a, a2 = NULL))
        }

        if (is.matrix(alpha) && ncol(alpha) >= 2) {
            if (nrow(alpha) == 1L) alpha <- alpha[rep.int(1L, N), , drop = FALSE]
            if (nrow(alpha) != N) stop("for order=2, alpha matrix must have nrow 1 or ncol(y)")
            return(list(a1 = as.numeric(alpha[, 1]), a2 = as.numeric(alpha[, 2])))
        }

        if (is.list(alpha) && length(alpha) >= 2) {
            a1 <- as.numeric(alpha[[1]])
            a2 <- as.numeric(alpha[[2]])
            if (length(a1) == 1L) a1 <- rep(a1, N)
            if (length(a2) == 1L) a2 <- rep(a2, N)
            if (length(a1) != N || length(a2) != N) stop("alpha[[1]] and alpha[[2]] must be length 1 or ncol(y)")
            return(list(a1 = a1, a2 = a2))
        }

        stop("for order=2, alpha must be a matrix with 2 columns or a list(alpha1, alpha2)")
    }

    ap <- unpack_alpha(order, alpha, N)
    alpha1 <- ap$a1
    alpha2 <- ap$a2

    log_sum_exp <- function(x) {
        m <- max(x)
        if (!is.finite(m)) return(-Inf)
        m + log(sum(exp(x - m)))
    }

    thin_vec <- function(k_vals, yprev, a, thinning) {
        if (!is.finite(a) || a < 0) return(rep(0, length(k_vals)))

        if (thinning == "binom") {
            if (a > 1) return(rep(0, length(k_vals)))
            return(dbinom(k_vals, size = yprev, prob = a))
        }

        if (thinning == "pois") {
            return(dpois(k_vals, lambda = a * yprev))
        }

        p_nb <- 1 / (1 + a)
        dnbinom(k_vals, size = yprev, prob = p_nb)
    }

    innov_vec <- function(u_vals, lam, innovation, sz) {
        if (!is.finite(lam) || lam <= 0) return(rep(0, length(u_vals)))

        if (innovation == "pois") return(dpois(u_vals, lambda = lam))
        if (innovation == "bell") return(dbell(u_vals, theta = lam))
        dnbinom(u_vals, size = sz, mu = lam)
    }

    loglik <- 0

    for (s in 1:n) {
        b <- blocks[s]
        for (i in 1:N) {
            lam <- theta[i] + tau[b]
            if (!is.finite(lam) || lam <= 0) return(-Inf)

            y_si <- y[s, i]

            if (i == 1 || order == 0) {
                p0 <- innov_vec(y_si, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)
                if (!is.finite(p0) || p0 <= 0) return(-Inf)
                loglik <- loglik + log(p0)
                next
            }

            if (order == 1 || i == 2) {
                y_prev <- y[s, i - 1]
                k_vals <- 0:y_si

                thin_vals <- thin_vec(k_vals, y_prev, alpha1[i], thinning)
                innov_vals <- innov_vec(y_si - k_vals, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

                conv <- sum(thin_vals * innov_vals)
                if (!is.finite(conv) || conv <= 0) return(-Inf)

                loglik <- loglik + log(conv)
                next
            }

            y1 <- y[s, i - 1]
            y2 <- y[s, i - 2]

            k1_vals <- 0:y_si
            log_terms <- rep(-Inf, length(k1_vals))

            for (idx in seq_along(k1_vals)) {
                k1 <- k1_vals[idx]
                thin1 <- thin_vec(k1, y1, alpha1[i], thinning)
                if (!is.finite(thin1) || thin1 <= 0) next

                remain <- y_si - k1
                k2_vals <- 0:remain

                thin2_vals <- thin_vec(k2_vals, y2, alpha2[i], thinning)
                innov_vals <- innov_vec(remain - k2_vals, lam, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

                inner_sum <- sum(thin2_vals * innov_vals)
                if (!is.finite(inner_sum) || inner_sum <= 0) next

                log_terms[idx] <- log(thin1) + log(inner_sum)
            }

            log_conv <- log_sum_exp(log_terms)
            if (!is.finite(log_conv)) return(-Inf)
            loglik <- loglik + log_conv
        }
    }

    loglik
}

#' INAD log likelihood contribution at time i (no fixed effect)
#'
#' Returns the time i contribution, summed over subjects, for the no fixed effect model.
#'
#' @param y Integer matrix n_sub by n_time.
#' @param i Time index in 1..ncol(y).
#' @param order Integer in {0, 1, 2}.
#' @param thinning One of "binom", "pois", "nbinom".
#' @param innovation One of "pois", "bell", "nbinom".
#' @param alpha Thinning parameters. For order 1, numeric length 1 or n_time.
#'   For order 2, either a matrix n_time by 2 or a list(alpha1, alpha2).
#' @param theta Innovation mean parameter at time i, or a vector length 1 or n_time.
#' @param nb_inno_size Size parameter for innovation "nbinom". Numeric length 1 or n_time.
#'
#' @return A scalar log likelihood contribution for time i.
#' @export
logL_inad_i <- function(
        y,
        i,
        order = 1,
        thinning = c("binom", "pois", "nbinom"),
        innovation = c("pois", "bell", "nbinom"),
        alpha,
        theta,
        nb_inno_size = NULL
) {
    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)

    if (!is.matrix(y)) y <- as.matrix(y)
    if (any(!is.finite(y)) || any(y < 0) || any(y != floor(y))) return(-Inf)

    n <- nrow(y)
    N <- ncol(y)

    if (!(order %in% c(0, 1, 2))) stop("order must be 0, 1, or 2")
    if (length(i) != 1L || i < 1L || i > N) stop("i must be in 1..ncol(y)")

    theta <- as.numeric(theta)
    if (length(theta) == 1L) theta <- rep(theta, N)
    if (length(theta) != N) stop("theta must be length 1 or ncol(y)")

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) stop("nb_inno_size is required for innovation='nbinom'")
        nb_inno_size <- as.numeric(nb_inno_size)
        if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, N)
        if (length(nb_inno_size) != N) stop("nb_inno_size must be length 1 or ncol(y)")
        if (any(nb_inno_size <= 0)) return(-Inf)
    }

    unpack_alpha <- function(order, alpha, N) {
        if (order == 0) return(list(a1 = NULL, a2 = NULL))

        if (order == 1) {
            a <- as.numeric(alpha)
            if (length(a) == 1L) a <- rep(a, N)
            if (length(a) != N) stop("alpha must be length 1 or ncol(y) for order=1")
            return(list(a1 = a, a2 = NULL))
        }

        if (is.matrix(alpha) && ncol(alpha) >= 2) {
            if (nrow(alpha) == 1L) alpha <- alpha[rep.int(1L, N), , drop = FALSE]
            if (nrow(alpha) != N) stop("for order=2, alpha matrix must have nrow 1 or ncol(y)")
            return(list(a1 = as.numeric(alpha[, 1]), a2 = as.numeric(alpha[, 2])))
        }

        if (is.list(alpha) && length(alpha) >= 2) {
            a1 <- as.numeric(alpha[[1]])
            a2 <- as.numeric(alpha[[2]])
            if (length(a1) == 1L) a1 <- rep(a1, N)
            if (length(a2) == 1L) a2 <- rep(a2, N)
            if (length(a1) != N || length(a2) != N) stop("alpha[[1]] and alpha[[2]] must be length 1 or ncol(y)")
            return(list(a1 = a1, a2 = a2))
        }

        stop("for order=2, alpha must be a matrix with 2 columns or a list(alpha1, alpha2)")
    }

    ap <- unpack_alpha(order, alpha, N)
    alpha1 <- ap$a1
    alpha2 <- ap$a2

    log_sum_exp <- function(x) {
        m <- max(x)
        if (!is.finite(m)) return(-Inf)
        m + log(sum(exp(x - m)))
    }

    thin_vec <- function(k_vals, yprev, a, thinning) {
        if (!is.finite(a) || a < 0) return(rep(0, length(k_vals)))

        if (thinning == "binom") {
            if (a > 1) return(rep(0, length(k_vals)))
            return(dbinom(k_vals, size = yprev, prob = a))
        }

        if (thinning == "pois") {
            return(dpois(k_vals, lambda = a * yprev))
        }

        p_nb <- 1 / (1 + a)
        dnbinom(k_vals, size = yprev, prob = p_nb)
    }

    innov_vec <- function(u_vals, lam, innovation, sz) {
        if (!is.finite(lam) || lam <= 0) return(rep(0, length(u_vals)))

        if (innovation == "pois") return(dpois(u_vals, lambda = lam))
        if (innovation == "bell") return(dbell(u_vals, theta = lam))
        dnbinom(u_vals, size = sz, mu = lam)
    }

    loglik_i <- 0
    th_i <- theta[i]
    if (!is.finite(th_i) || th_i <= 0) return(-Inf)

    for (s in 1:n) {
        y_si <- y[s, i]

        if (i == 1 || order == 0) {
            p0 <- innov_vec(y_si, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)
            if (!is.finite(p0) || p0 <= 0) return(-Inf)
            loglik_i <- loglik_i + log(p0)
            next
        }

        if (order == 1 || i == 2) {
            y_prev <- y[s, i - 1]
            k_vals <- 0:y_si

            thin_vals <- thin_vec(k_vals, y_prev, alpha1[i], thinning)
            innov_vals <- innov_vec(y_si - k_vals, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

            conv <- sum(thin_vals * innov_vals)
            if (!is.finite(conv) || conv <= 0) return(-Inf)

            loglik_i <- loglik_i + log(conv)
            next
        }

        y1 <- y[s, i - 1]
        y2 <- y[s, i - 2]

        k1_vals <- 0:y_si
        log_terms <- rep(-Inf, length(k1_vals))

        for (idx in seq_along(k1_vals)) {
            k1 <- k1_vals[idx]
            thin1 <- thin_vec(k1, y1, alpha1[i], thinning)
            if (!is.finite(thin1) || thin1 <= 0) next

            remain <- y_si - k1
            k2_vals <- 0:remain

            thin2_vals <- thin_vec(k2_vals, y2, alpha2[i], thinning)
            innov_vals <- innov_vec(remain - k2_vals, th_i, innovation, if (innovation == "nbinom") nb_inno_size[i] else NA_real_)

            inner_sum <- sum(thin2_vals * innov_vals)
            if (!is.finite(inner_sum) || inner_sum <= 0) next

            log_terms[idx] <- log(thin1) + log(inner_sum)
        }

        log_conv <- log_sum_exp(log_terms)
        if (!is.finite(log_conv)) return(-Inf)
        loglik_i <- loglik_i + log_conv
    }

    loglik_i
}
