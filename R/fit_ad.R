#' Fit Gaussian antedependence model
#'
#' Fits an AD(0), AD(1), or AD(2) model for Gaussian longitudinal data.
#' Missing values can be handled by complete-case deletion or by EM.
#'
#' @param y Numeric matrix (n_subjects x n_time). May contain NA.
#' @param order Integer 0, 1, or 2.
#' @param blocks Optional vector of block membership (length n_subjects).
#' @param na_action One of "fail", "complete", or "em".
#' @param estimate_mu Logical, whether to estimate mu (default TRUE).
#' @param em_max_iter Maximum EM iterations (only used when na_action = "em").
#' @param em_tol EM convergence tolerance (only used when na_action = "em").
#' @param em_verbose Logical, print EM progress (only used when na_action = "em").
#' @param ... Passed through to the EM fitter.
#'
#' @return A list with components including mu, phi, sigma, tau, log_l, n_obs, n_missing.
#' @export
fit_ad <- function(y, order = 1, blocks = NULL,
                   na_action = c("fail", "complete", "em"),
                   estimate_mu = TRUE,
                   em_max_iter = 100, em_tol = 1e-6, em_verbose = FALSE, ...) {
    na_action <- match.arg(na_action)

    if (!is.matrix(y)) y <- as.matrix(y)
    storage.mode(y) <- "double"

    if (!order %in% c(0, 1, 2)) {
        stop("order must be 0, 1, or 2.", call. = FALSE)
    }

    if (!is.null(blocks)) {
        if (length(blocks) != nrow(y)) {
            stop("blocks must have length n_subjects.", call. = FALSE)
        }
        blocks <- as.integer(factor(blocks))
    }

    has_missing <- any(is.na(y))

    if (has_missing) {
        if (na_action == "fail") {
            stop("y contains NA values. Use na_action = 'complete' or 'em'.", call. = FALSE)
        }

        if (na_action == "complete") {
            cc <- .extract_complete_cases(y, blocks)
            y <- cc$y
            blocks <- cc$blocks
            has_missing <- FALSE
        }

        if (na_action == "em") {
            return(.fit_ad_em(
                y = y,
                order = order,
                blocks = blocks,
                estimate_mu = estimate_mu,
                max_iter = em_max_iter,
                tol = em_tol,
                verbose = em_verbose,
                ...
            ))
        }
    }

    if (!is.matrix(y)) y <- as.matrix(y)
    n_subjects <- nrow(y)
    n_time <- ncol(y)

    if (n_subjects < 2) {
        stop("Need at least 2 subjects to fit the model.", call. = FALSE)
    }

    if (!is.null(blocks)) {
        k <- length(unique(blocks))
        tau <- numeric(k)
        tau[1] <- 0

        idx_ref <- which(blocks == 1)
        if (length(idx_ref) == 0) {
            stop("Reference block (level 1) has no subjects.", call. = FALSE)
        }

        mu <- if (isTRUE(estimate_mu)) {
            colMeans(y[idx_ref, , drop = FALSE])
        } else {
            rep(0, n_time)
        }

        if (k > 1) {
            for (b in 2:k) {
                idx_b <- which(blocks == b)
                if (length(idx_b) == 0) next
                diffs <- y[idx_b, , drop = FALSE] - matrix(mu, nrow = length(idx_b), ncol = n_time, byrow = TRUE)
                tau[b] <- mean(diffs)
            }
        }

        mu_mat <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
        tau_mat <- matrix(tau[blocks], nrow = n_subjects, ncol = n_time)
        y_center <- y - mu_mat - tau_mat
    } else {
        tau <- 0
        mu <- if (isTRUE(estimate_mu)) colMeans(y) else rep(0, n_time)
        y_center <- y - matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
    }

    if (order == 0) {
        phi <- numeric(0)
        sigma <- apply(y_center, 2, function(z) sqrt(max(mean(z^2), 1e-8)))
    }

    if (order == 1) {
        phi <- numeric(n_time)
        sigma <- numeric(n_time)
        phi[1] <- 0

        sigma[1] <- sqrt(max(mean(y_center[, 1]^2), 1e-8))

        for (t in 2:n_time) {
            x <- y_center[, t - 1]
            yy <- y_center[, t]
            denom <- sum(x^2)

            if (!is.finite(denom) || denom <= 1e-12) {
                phi[t] <- 0
                resid <- yy
            } else {
                phi[t] <- sum(x * yy) / denom
                resid <- yy - phi[t] * x
            }

            sigma[t] <- sqrt(max(mean(resid^2), 1e-8))
        }
    }

    if (order == 2) {
        phi <- matrix(0, nrow = 2, ncol = n_time)
        sigma <- numeric(n_time)

        sigma[1] <- sqrt(max(mean(y_center[, 1]^2), 1e-8))

        if (n_time >= 2) {
            x <- y_center[, 1]
            yy <- y_center[, 2]
            denom <- sum(x^2)

            if (!is.finite(denom) || denom <= 1e-12) {
                phi[1, 2] <- 0
                resid <- yy
            } else {
                phi[1, 2] <- sum(x * yy) / denom
                resid <- yy - phi[1, 2] * x
            }
            sigma[2] <- sqrt(max(mean(resid^2), 1e-8))
        }

        if (n_time >= 3) {
            for (t in 3:n_time) {
                X <- cbind(y_center[, t - 1], y_center[, t - 2])
                yy <- y_center[, t]

                fit <- tryCatch(
                    stats::lm.fit(x = X, y = yy),
                    error = function(e) NULL
                )

                if (is.null(fit) || any(!is.finite(fit$coefficients))) {
                    phi[, t] <- c(0, 0)
                    resid <- yy
                } else {
                    phi[, t] <- fit$coefficients
                    resid <- fit$residuals
                }

                sigma[t] <- sqrt(max(mean(resid^2), 1e-8))
            }
        }
    }

    log_l <- logL_ad(
        y = y,
        order = order,
        mu = mu,
        phi = phi,
        sigma = sigma,
        blocks = blocks,
        tau = tau,
        na_action = "fail"
    )

    n_blocks <- if (is.null(blocks)) 1 else length(unique(blocks))
    n_params <- .count_params_ad(order, n_time, n_blocks)

    result <- list(
        mu = mu,
        phi = phi,
        sigma = sigma,
        tau = tau,
        log_l = log_l,
        aic = -2 * log_l + 2 * n_params,
        bic = -2 * log_l + n_params * log(n_subjects),
        n_obs = sum(!is.na(y)),
        n_missing = sum(is.na(y)),
        pct_missing = mean(is.na(y)) * 100,
        missing_pattern = "complete",
        settings = list(
            order = order,
            n_time = n_time,
            n_subjects = n_subjects,
            blocks = blocks,
            estimate_mu = estimate_mu,
            na_action = na_action
        )
    )

    class(result) <- c("ad_fit", "list")
    result
}
