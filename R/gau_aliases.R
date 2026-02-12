#' Fit Gaussian AD model (alias of fit_ad)
#'
#' Convenience alias for \code{\link{fit_ad}}.
#'
#' @param ... Passed through to \code{\link{fit_ad}}.
#' @return See \code{\link{fit_ad}}.
#' @seealso \code{\link{fit_ad}}, \code{\link{fit_cat}}, \code{\link{fit_inad}}
#' @export
fit_gau <- function(...) {
  fit_ad(...)
}

#' Compute AIC for Gaussian AD model (alias of aic_ad)
#'
#' Convenience alias for \code{\link{aic_ad}}.
#'
#' @param ... Passed through to \code{\link{aic_ad}}.
#' @return See \code{\link{aic_ad}}.
#' @seealso \code{\link{aic_ad}}, \code{\link{aic_cat}}, \code{\link{aic_inad}}
#' @export
aic_gau <- function(...) {
  aic_ad(...)
}

#' Compute BIC for Gaussian AD model (alias of bic_ad)
#'
#' Convenience alias for \code{\link{bic_ad}}.
#'
#' @param ... Passed through to \code{\link{bic_ad}}.
#' @return See \code{\link{bic_ad}}.
#' @seealso \code{\link{bic_ad}}, \code{\link{bic_cat}}, \code{\link{bic_inad}}
#' @export
bic_gau <- function(...) {
  bic_ad(...)
}

#' Compare Gaussian AD orders by BIC (alias of bic_order_ad)
#'
#' Convenience alias for \code{\link{bic_order_ad}}.
#'
#' @param ... Passed through to \code{\link{bic_order_ad}}.
#' @return See \code{\link{bic_order_ad}}.
#' @seealso \code{\link{bic_order_ad}}, \code{\link{bic_order_cat}}, \code{\link{bic_order_inad}}
#' @export
bic_order_gau <- function(...) {
  bic_order_ad(...)
}

#' Confidence intervals for Gaussian AD model (alias of ci_ad)
#'
#' Convenience alias for \code{\link{ci_ad}}.
#'
#' @param ... Passed through to \code{\link{ci_ad}}.
#' @return See \code{\link{ci_ad}}.
#' @seealso \code{\link{ci_ad}}, \code{\link{ci_cat}}, \code{\link{ci_inad}}
#' @export
ci_gau <- function(...) {
  ci_ad(...)
}

#' EM algorithm for Gaussian AD model (alias of em_ad)
#'
#' Convenience alias for \code{\link{em_ad}}.
#'
#' @param ... Passed through to \code{\link{em_ad}}.
#' @return See \code{\link{em_ad}}.
#' @seealso \code{\link{em_ad}}, \code{\link{em_inad}}
#' @export
em_gau <- function(...) {
  em_ad(...)
}

#' Log-likelihood for Gaussian AD model (alias of logL_ad)
#'
#' Convenience alias for \code{\link{logL_ad}}.
#'
#' @param ... Passed through to \code{\link{logL_ad}}.
#' @return See \code{\link{logL_ad}}.
#' @seealso \code{\link{logL_ad}}, \code{\link{logL_cat}}, \code{\link{logL_inad}}
#' @export
logL_gau <- function(...) {
  logL_ad(...)
}

#' Simulate Gaussian AD data (alias of simulate_ad)
#'
#' Convenience alias for \code{\link{simulate_ad}}.
#'
#' @param ... Passed through to \code{\link{simulate_ad}}.
#' @return See \code{\link{simulate_ad}}.
#' @seealso \code{\link{simulate_ad}}, \code{\link{simulate_cat}}, \code{\link{simulate_inad}}
#' @export
simulate_gau <- function(...) {
  simulate_ad(...)
}

#' LRT for Gaussian AD order (alias of lrt_order_ad)
#'
#' Convenience alias for \code{\link{lrt_order_ad}}.
#'
#' @param ... Passed through to \code{\link{lrt_order_ad}}.
#' @return See \code{\link{lrt_order_ad}}.
#' @seealso \code{\link{lrt_order_ad}}, \code{\link{lrt_order_cat}}, \code{\link{lrt_order_inad}}
#' @export
lrt_order_gau <- function(...) {
  lrt_order_ad(...)
}

#' One-sample mean test for Gaussian AD (alias of lrt_one_sample_ad)
#'
#' Convenience alias for \code{\link{lrt_one_sample_ad}}.
#'
#' @param ... Passed through to \code{\link{lrt_one_sample_ad}}.
#' @return See \code{\link{lrt_one_sample_ad}}.
#' @seealso \code{\link{lrt_one_sample_ad}}, \code{\link{lrt_two_sample_ad}}
#' @export
lrt_one_sample_gau <- function(...) {
  lrt_one_sample_ad(...)
}

#' Two-sample mean profile test for Gaussian AD (alias of lrt_two_sample_ad)
#'
#' Convenience alias for \code{\link{lrt_two_sample_ad}}.
#'
#' @param ... Passed through to \code{\link{lrt_two_sample_ad}}.
#' @return See \code{\link{lrt_two_sample_ad}}.
#' @seealso \code{\link{lrt_two_sample_ad}}, \code{\link{lrt_one_sample_ad}}
#' @export
lrt_two_sample_gau <- function(...) {
  lrt_two_sample_ad(...)
}

#' Linear contrast test for Gaussian AD means (alias of lrt_contrast_ad)
#'
#' Convenience alias for \code{\link{lrt_contrast_ad}}.
#'
#' @param ... Passed through to \code{\link{lrt_contrast_ad}}.
#' @return See \code{\link{lrt_contrast_ad}}.
#' @seealso \code{\link{lrt_contrast_ad}}
#' @export
lrt_contrast_gau <- function(...) {
  lrt_contrast_ad(...)
}

#' Covariance homogeneity test for Gaussian AD (alias of lrt_homogeneity_ad)
#'
#' Convenience alias for \code{\link{lrt_homogeneity_ad}}.
#'
#' @param ... Passed through to \code{\link{lrt_homogeneity_ad}}.
#' @return See \code{\link{lrt_homogeneity_ad}}.
#' @seealso \code{\link{lrt_homogeneity_ad}}, \code{\link{lrt_homogeneity_cat}}, \code{\link{lrt_homogeneity_inad}}
#' @export
lrt_homogeneity_gau <- function(...) {
  lrt_homogeneity_ad(...)
}
