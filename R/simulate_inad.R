#' Simulate integer valued antedependence series
#'
#' Generate longitudinal count data from an INAD model using a thinning
#' operator and an innovation distribution.
#'
#' Time 1 observations are generated from the innovation distribution alone.
#' For times 2 to `n_time`, counts are generated as thinning of previous
#' counts plus independent innovations. When `order = 0`, all time points
#' are generated from the innovation distribution and the thinning operator
#' and `alpha` are ignored.
#'
#' @param n_subjects number of subjects
#' @param n_time number of time points
#' @param order antedependence order, 0, 1 or 2
#' @param thinning thinning operator, one of `"binom"`, `"pois"`, `"nbinom"`
#' @param innovation innovation distribution, one of `"pois"`, `"bell"`,
#'   `"nbinom"`
#' @param alpha thinning parameter or vector or matrix; if `NULL`, defaults
#'   are used depending on the order
#' @param theta innovation parameter or vector; if `NULL`, defaults are used
#'   depending on the innovation type
#' @param nb_inno_size probability parameter for negative binomial innovations
#'   when `innovation = "nbinom"`; if `NULL`, defaults to 0.5
#'
#' @return integer matrix of counts with dimension `n_subjects` by `n_time`
#' @export
simulate_inad <- function(n_subjects,
                          n_time,
                          order = 1L,
                          thinning = c("binom", "pois", "nbinom"),
                          innovation = c("pois", "bell", "nbinom"),
                          alpha = NULL,
                          theta = NULL,
                          nb_inno_size = NULL) {

    thinning <- match.arg(thinning)
    innovation <- match.arg(innovation)

    if (!order %in% c(0L, 1L, 2L)) {
        stop("order must be 0, 1 or 2")
    }

    if (length(n_subjects) != 1L || n_subjects <= 0) {
        stop("n_subjects must be a positive integer")
    }
    if (length(n_time) != 1L || n_time <= 0) {
        stop("n_time must be a positive integer")
    }

    n_subjects <- as.integer(n_subjects)
    n_time <- as.integer(n_time)
    order <- as.integer(order)

    if (is.null(theta)) {
        if (innovation == "pois") {
            theta <- rep(10, n_time)
        } else if (innovation == "bell") {
            theta <- rep(2, n_time)
        } else if (innovation == "nbinom") {
            theta <- rep(10, n_time)
        }
    } else if (length(theta) == 1L) {
        theta <- rep(theta, n_time)
    }

    if (length(theta) != n_time) {
        stop("theta must be NULL, scalar, or length n_time")
    }

    if (innovation == "nbinom") {
        if (is.null(nb_inno_size)) {
            nb_inno_size <- rep(0.5, n_time)
        } else if (length(nb_inno_size) == 1L) {
            nb_inno_size <- rep(nb_inno_size, n_time)
        }
        if (length(nb_inno_size) != n_time) {
            stop("nb_inno_size must be NULL, scalar, or length n_time")
        }
        if (any(nb_inno_size <= 0 | nb_inno_size >= 1)) {
            stop("nb_inno_size must lie in (0, 1)")
        }
    }

    if (is.null(alpha)) {
        if (order == 0L) {
            alpha <- 0
        } else if (order == 1L) {
            alpha <- rep(0.5, n_time)
            alpha[1L] <- 0
        } else {
            alpha <- matrix(0, nrow = 2L, ncol = n_time)
            alpha[1L, 2:n_time] <- 0.5
            if (n_time >= 3L) {
                alpha[2L, 3:n_time] <- 0.25
            }
        }
    } else {
        if (order == 0L) {
            if (length(alpha) != 1L || alpha != 0) {
                warning("order = 0, alpha is ignored and treated as 0")
            }
        } else if (order == 1L) {
            if (length(alpha) == 1L) {
                alpha <- rep(alpha, n_time)
            }
            if (length(alpha) != n_time) {
                stop("for order 1, alpha must be scalar or length n_time")
            }
        } else {
            if (!is.matrix(alpha) || nrow(alpha) != 2L || ncol(alpha) != n_time) {
                stop("for order 2, alpha must be a 2 by n_time matrix")
            }
        }
    }

    if (order > 0L) {
        if (thinning %in% c("binom", "nbinom")) {
            if (order == 1L) {
                if (any(alpha < 0 | alpha >= 1)) {
                    stop("alpha must lie in [0, 1) for binom and nbinom thinning")
                }
            } else {
                if (any(alpha < 0 | alpha >= 1)) {
                    stop("alpha must lie in [0, 1) for binom and nbinom thinning")
                }
            }
        } else if (thinning == "pois") {
            if (order == 1L) {
                if (any(alpha < 0)) {
                    stop("alpha must be nonnegative for pois thinning")
                }
            } else {
                if (any(alpha < 0)) {
                    stop("alpha must be nonnegative for pois thinning")
                }
            }
        }
    }

    y <- matrix(0L, nrow = n_subjects, ncol = n_time)

    draw_innovation <- function(t_index) {
        if (innovation == "pois") {
            stats::rpois(n_subjects, lambda = theta[t_index])
        } else if (innovation == "bell") {
            rbell(n_subjects, theta = theta[t_index])
        } else {
            stats::rnbinom(
                n_subjects,
                size = theta[t_index],
                prob = nb_inno_size[t_index]
            )
        }
    }

    if (order == 0L) {
        for (t in seq_len(n_time)) {
            y[, t] <- as.integer(draw_innovation(t))
        }
        return(y)
    }

    y[, 1L] <- as.integer(draw_innovation(1L))

    if (n_time >= 2L) {
        for (t in 2L:n_time) {
            thinned_total <- integer(n_subjects)

            if (order >= 1L) {
                alpha_t1 <- if (order == 1L) alpha[t] else alpha[1L, t]
                if (alpha_t1 > 0) {
                    y_prev1 <- y[, t - 1L]
                    if (thinning == "binom") {
                        thinned1 <- stats::rbinom(
                            n = n_subjects,
                            size = y_prev1,
                            prob = alpha_t1
                        )
                    } else if (thinning == "pois") {
                        thinned1 <- stats::rpois(
                            n = n_subjects,
                            lambda = alpha_t1 * y_prev1
                        )
                    } else {
                        prob_nb1 <- 1 / (1 + alpha_t1)
                        thinned1 <- stats::rnbinom(
                            n = n_subjects,
                            size = y_prev1,
                            prob = prob_nb1
                        )
                    }
                    thinned_total <- thinned_total + thinned1
                }
            }

            if (order >= 2L && t > 2L) {
                alpha_t2 <- alpha[2L, t]
                if (alpha_t2 > 0) {
                    y_prev2 <- y[, t - 2L]
                    if (thinning == "binom") {
                        thinned2 <- stats::rbinom(
                            n = n_subjects,
                            size = y_prev2,
                            prob = alpha_t2
                        )
                    } else if (thinning == "pois") {
                        thinned2 <- stats::rpois(
                            n = n_subjects,
                            lambda = alpha_t2 * y_prev2
                        )
                    } else {
                        prob_nb2 <- 1 / (1 + alpha_t2)
                        thinned2 <- stats::rnbinom(
                            n = n_subjects,
                            size = y_prev2,
                            prob = prob_nb2
                        )
                    }
                    thinned_total <- thinned_total + thinned2
                }
            }

            eps_t <- draw_innovation(t)
            y[, t] <- as.integer(thinned_total + eps_t)
        }
    }

    y
}
