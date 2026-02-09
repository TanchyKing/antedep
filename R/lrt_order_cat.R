# lrt_order_cat.R - Likelihood ratio tests for order in categorical AD models

.cat_fit_uses_missing_likelihood <- function(fit) {
  na_effective <- fit$settings$na_action_effective
  na_raw <- fit$settings$na_action
  identical(na_effective, "marginalize") ||
    identical(na_raw, "marginalize") ||
    is.null(fit$cell_counts)
}

.stop_cat_missing_inference <- function(fn_name) {
  stop(
    paste0(
      fn_name,
      " currently supports complete data only. Missing-data CAT likelihood-ratio inference is not implemented yet."
    ),
    call. = FALSE
  )
}

#' Likelihood ratio test for antedependence order (categorical data)
#'
#' Tests whether a higher-order AD model provides significantly better fit
#' than a lower-order model for categorical longitudinal data.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c. Can be NULL if both fit_null and
#'   fit_alt are provided.
#' @param order_null Order under the null hypothesis (default 0).
#' @param order_alt Order under the alternative hypothesis (default 1).
#'   Must be greater than order_null.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership.
#' @param homogeneous Logical. If TRUE (default), parameters are shared across
#'   all groups.
#' @param n_categories Number of categories. If NULL, inferred from data.
#' @param fit_null Optional pre-fitted model under null hypothesis (class "cat_fit").
#'   If provided, y is not required for fitting under H0.
#' @param fit_alt Optional pre-fitted model under alternative hypothesis.
#'   If provided, y is not required for fitting under H1.
#'
#' @return A list of class \code{"cat_lrt"} containing:
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{fit_null}{Fitted model under H0}
#'   \item{fit_alt}{Fitted model under H1}
#'   \item{order_null}{Order under null}
#'   \item{order_alt}{Order under alternative}
#'   \item{table}{Summary data frame}
#'
#' @details
#' The likelihood ratio test statistic is:
#' \deqn{\lambda = -2[\ell_0 - \ell_1]}
#' where \eqn{\ell_0} and \eqn{\ell_1} are the maximized log-likelihoods under
#' the null and alternative hypotheses.
#'
#' Under H0, \eqn{\lambda} follows a chi-square distribution with degrees of
#' freedom equal to the difference in the number of free parameters.
#'
#' For testing AD(p) vs AD(p+1), the degrees of freedom are:
#' \deqn{df = (c-1) \times c^p \times (n - p - 1)}
#' where c is the number of categories and n is the number of time points.
#'
#' @examples
#' \dontrun{
#' # Simulate AD(1) data
#' set.seed(123)
#' y <- simulate_cat(200, 6, order = 1, n_categories = 2)
#'
#' # Test AD(0) vs AD(1)
#' test_01 <- lrt_order_cat(y, order_null = 0, order_alt = 1)
#' print(test_01$table)
#'
#' # Test AD(1) vs AD(2)
#' test_12 <- lrt_order_cat(y, order_null = 1, order_alt = 2)
#' print(test_12$table)
#'
#' # Using pre-fitted models
#' fit0 <- fit_cat(y, order = 0)
#' fit1 <- fit_cat(y, order = 1)
#' test_prefitted <- lrt_order_cat(fit_null = fit0, fit_alt = fit1)
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @seealso \code{\link{fit_cat}}, \code{\link{bic_order_cat}}
#'
#' @export
lrt_order_cat <- function(y = NULL, order_null = 0, order_alt = 1,
                          blocks = NULL, homogeneous = TRUE, n_categories = NULL,
                          fit_null = NULL, fit_alt = NULL) {
  
 # Validate that we have either y or pre-fitted models
  if (is.null(y) && (is.null(fit_null) || is.null(fit_alt))) {
    stop("Either y must be provided, or both fit_null and fit_alt must be provided")
  }
  use_missing <- !is.null(y) && anyNA(y)
  na_action_fit <- if (use_missing) "marginalize" else "fail"
  
  # Validate orders
  if (order_alt <= order_null) {
    stop("order_alt must be greater than order_null")
  }
  if (order_null < 0) {
    stop("order_null must be non-negative")
  }
  if (order_alt > 2) {
    stop("order_alt > 2 not currently supported")
  }
  
  # Fit models if not provided
  if (is.null(fit_null)) {
    fit_null <- fit_cat(y, order = order_null, blocks = blocks,
                        homogeneous = homogeneous, n_categories = n_categories,
                        na_action = na_action_fit)
  } else {
    if (!inherits(fit_null, "cat_fit")) {
      stop("fit_null must be a cat_fit object")
    }
    order_null <- fit_null$settings$order
  }
  
  if (is.null(fit_alt)) {
    fit_alt <- fit_cat(y, order = order_alt, blocks = blocks,
                       homogeneous = homogeneous, n_categories = n_categories,
                       na_action = na_action_fit)
  } else {
    if (!inherits(fit_alt, "cat_fit")) {
      stop("fit_alt must be a cat_fit object")
    }
    order_alt <- fit_alt$settings$order
  }
  
  # Verify order relationship from fitted models
  if (fit_alt$settings$order <= fit_null$settings$order) {
    stop("Alternative model order must be greater than null model order")
  }
  
  # Extract info
  log_l_null <- fit_null$log_l
  log_l_alt <- fit_alt$log_l
  n_params_null <- fit_null$n_params
  n_params_alt <- fit_alt$n_params
  
  # Compute LRT statistic
  lrt_stat <- -2 * (log_l_null - log_l_alt)
  
  # Degrees of freedom
  df <- n_params_alt - n_params_null
  
  # P-value from chi-square
  p_value <- stats::pchisq(lrt_stat, df = df, lower.tail = FALSE)
  
  # Build summary table
  table_df <- data.frame(
    model = c("H0 (null)", "H1 (alt)"),
    order = c(fit_null$settings$order, fit_alt$settings$order),
    log_l = c(log_l_null, log_l_alt),
    n_params = c(n_params_null, n_params_alt),
    aic = c(fit_null$aic, fit_alt$aic),
    bic = c(fit_null$bic, fit_alt$bic)
  )
  
  # Assemble output
  out <- list(
    lrt_stat = lrt_stat,
    df = df,
    p_value = p_value,
    fit_null = fit_null,
    fit_alt = fit_alt,
    order_null = fit_null$settings$order,
    order_alt = fit_alt$settings$order,
    table = table_df
  )
  
  class(out) <- "cat_lrt"
  out
}


#' Print method for cat_lrt objects
#'
#' @param x A cat_lrt object
#' @param ... Additional arguments (ignored)
#'
#' @export
print.cat_lrt <- function(x, ...) {
  cat("Likelihood Ratio Test for Categorical AD Order\n")
  cat("===============================================\n\n")
  
  cat("H0: AD(", x$order_null, ")\n", sep = "")
  cat("H1: AD(", x$order_alt, ")\n\n", sep = "")
  
  cat("Test statistic:", round(x$lrt_stat, 4), "\n")
  cat("Degrees of freedom:", x$df, "\n")
  cat("P-value:", format.pval(x$p_value, digits = 4), "\n\n")
  
  cat("Model comparison:\n")
  print(x$table, row.names = FALSE)
  
  invisible(x)
}


#' Run all pairwise order tests
#'
#' Performs sequential likelihood ratio tests for AD orders 0 vs 1, 1 vs 2, etc.
#'
#' @param y Integer matrix of categorical data (n_subjects x n_time).
#' @param max_order Maximum order to test. Default is 2.
#' @param blocks Optional block membership vector.
#' @param homogeneous Whether to use homogeneous parameters across blocks.
#' @param n_categories Number of categories (inferred if NULL).
#'
#' @return A list containing:
#'   \item{tests}{List of lrt_order_cat results for each comparison}
#'   \item{table}{Summary data frame with all comparisons}
#'   \item{fits}{List of all fitted models}
#'   \item{selected_order}{Recommended order based on sequential testing at alpha=0.05}
#'
#' @details
#' This function performs forward selection: starting from order 0, it tests
#' whether increasing the order provides significant improvement. The selected
#' order is the highest order where the test was significant (at alpha = 0.05).
#'
#' @examples
#' \dontrun{
#' y <- simulate_cat(200, 6, order = 1, n_categories = 2)
#' result <- run_order_tests_cat(y, max_order = 2)
#' print(result$table)
#' cat("Selected order:", result$selected_order, "\n")
#' }
#'
#' @export
run_order_tests_cat <- function(y, max_order = 2, blocks = NULL, 
                                 homogeneous = TRUE, n_categories = NULL) {
  use_missing <- anyNA(y)
  
  # Validate data
  validated <- .validate_y_cat(y, n_categories, allow_na = use_missing)
  y <- validated$y
  n_categories <- validated$n_categories
  n_time <- ncol(y)
  
  # Adjust max_order if needed
  max_possible <- n_time - 1
  if (max_order > max_possible) {
    warning("max_order reduced to ", max_possible, " (n_time - 1)")
    max_order <- max_possible
  }
  if (max_order > 2) {
    warning("max_order > 2 not supported; setting to 2")
    max_order <- 2
  }
  
  # Fit all models first
  orders <- 0:max_order
  fits <- vector("list", length(orders))
  names(fits) <- paste0("order_", orders)
  
  for (i in seq_along(orders)) {
    fits[[i]] <- fit_cat(y, order = orders[i], blocks = blocks,
                         homogeneous = homogeneous, n_categories = n_categories,
                         na_action = if (use_missing) "marginalize" else "fail")
  }
  
  # Run pairwise tests
  n_tests <- max_order
  tests <- vector("list", n_tests)
  names(tests) <- paste0("order_", 0:(max_order-1), "_vs_", 1:max_order)
  
  test_results <- data.frame(
    comparison = character(n_tests),
    order_null = integer(n_tests),
    order_alt = integer(n_tests),
    lrt_stat = numeric(n_tests),
    df = integer(n_tests),
    p_value = numeric(n_tests),
    significant = logical(n_tests),
    stringsAsFactors = FALSE
  )
  
  for (i in seq_len(n_tests)) {
    p_null <- orders[i]
    p_alt <- orders[i + 1]
    
    tests[[i]] <- lrt_order_cat(fit_null = fits[[i]], fit_alt = fits[[i + 1]])
    
    test_results$comparison[i] <- paste0("AD(", p_null, ") vs AD(", p_alt, ")")
    test_results$order_null[i] <- p_null
    test_results$order_alt[i] <- p_alt
    test_results$lrt_stat[i] <- tests[[i]]$lrt_stat
    test_results$df[i] <- tests[[i]]$df
    test_results$p_value[i] <- tests[[i]]$p_value
    test_results$significant[i] <- tests[[i]]$p_value < 0.05
  }
  
  # Determine selected order via forward selection
  # Start at order 0, increase while tests are significant
  selected_order <- 0
  for (i in seq_len(n_tests)) {
    if (test_results$significant[i]) {
      selected_order <- test_results$order_alt[i]
    } else {
      break
    }
  }
  
  list(
    tests = tests,
    table = test_results,
    fits = fits,
    selected_order = selected_order
  )
}
