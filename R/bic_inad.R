#' Bayesian information criterion for fitted INAD models
#'
#' Computes BIC using the fitted log likelihood and a parameter count that
#' respects structural zeros and identifiability constraints.
#'
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#' @param n_subjects Number of subjects, typically \code{nrow(y)}.
#'
#' @return A numeric scalar BIC value.
#'
#' @details
#' The BIC is computed as:
#' \deqn{BIC = -2 \times \ell + k \times \log(N)}
#' where \eqn{\ell} is the log-likelihood, \eqn{k} is the number of free parameters,
#' and \eqn{N} is the number of subjects.
#' @export
bic_inad <- function(fit, n_subjects) {
    if (missing(n_subjects) || length(n_subjects) != 1 || !is.finite(n_subjects) || n_subjects <= 0) {
        stop("n_subjects must be a positive finite scalar")
    }

    k <- .count_params_inad_fit(fit)
    -2 * fit$log_l + k * log(n_subjects)
}

#' Akaike information criterion for fitted INAD models
#'
#' Computes AIC using the fitted log likelihood and a parameter count that
#' respects structural zeros and identifiability constraints.
#'
#' @param fit A fitted model object returned by \code{\link{fit_inad}}.
#'
#' @return A numeric scalar AIC value.
#'
#' @details
#' The AIC is computed as:
#' \deqn{AIC = -2 \times \ell + 2k}
#' where \eqn{\ell} is the log-likelihood and \eqn{k} is the number of free
#' parameters.
#' @export
aic_inad <- function(fit) {
    k <- .count_params_inad_fit(fit)
    -2 * fit$log_l + 2 * k
}

#' @keywords internal
.count_params_inad_fit <- function(fit) {
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")
    if (is.null(fit$settings$innovation)) stop("fit$settings$innovation is missing")

    ord <- fit$settings$order
    innovation <- fit$settings$innovation

    N <- length(fit$theta)
    if (N < 1) stop("fit$theta must be nonempty")

    k <- 0
    k <- k + N

    if (ord == 1) k <- k + (N - 1)
    if (ord == 2) k <- k + (2 * N - 3)

    if (!is.null(fit$tau)) {
        B <- length(fit$tau)
        if (B >= 1) k <- k + (B - 1)
    }

    if (!is.null(fit$nb_inno_size) && innovation == "nbinom") {
        if (length(fit$nb_inno_size) == 1L) k <- k + 1 else k <- k + N
    }

    k
}
