# File: R/bic_ad.R

#' Bayesian information criterion for fitted antedependence models
#'
#' Computes BIC using the fitted log likelihood and a parameter count that
#' respects identifiability constraints for the antedependence parameters.
#'
#' @param fit A fitted model object returned by \code{\link{fit_ad}}.
#' @param n_subjects Number of subjects, typically \code{nrow(y)}.
#'
#' @return A numeric scalar BIC value.
#' @export
bic_ad <- function(fit, n_subjects) {
    if (is.null(fit$log_l) || !is.finite(fit$log_l)) stop("fit$log_l must be finite")
    if (is.null(fit$settings$order)) stop("fit$settings$order is missing")

    if (missing(n_subjects) || length(n_subjects) != 1 || !is.finite(n_subjects) || n_subjects <= 0) {
        stop("n_subjects must be a positive finite scalar")
    }

    ord <- fit$settings$order
    if (!(ord %in% c(0, 1, 2))) stop("fit$settings$order must be 0, 1, or 2")

    N <- NULL
    if (!is.null(fit$sigma)) N <- length(as.numeric(fit$sigma))
    if (is.null(N) && !is.null(fit$mu)) N <- length(as.numeric(fit$mu))
    if (is.null(N) || N < 1) stop("Cannot infer n_time from fit; need fit$sigma or fit$mu")

    estimate_mu <- fit$settings$estimate_mu
    if (is.null(estimate_mu)) estimate_mu <- !is.null(fit$mu)
    estimate_mu <- isTRUE(estimate_mu)

    k <- 0

    if (estimate_mu) k <- k + N

    k <- k + N

    if (ord == 1) k <- k + (N - 1)
    if (ord == 2) k <- k + (2 * N - 3)

    if (!is.null(fit$tau)) {
        B <- length(as.numeric(fit$tau))
        if (B >= 1) k <- k + (B - 1)
    }

    -2 * fit$log_l + k * log(n_subjects)
}
