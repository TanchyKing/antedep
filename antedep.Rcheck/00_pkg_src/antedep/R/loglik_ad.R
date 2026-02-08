#' Log-likelihood for Gaussian AD models (with missing data support)
#'
#' Computes the log-likelihood for Gaussian antedependence models of order 0, 1, or 2.
#' Supports missing data under MAR assumption via na_action parameter.
#'
#' @param y Numeric matrix with n_subjects rows and n_time columns. May contain NA.
#' @param order Antedependence order, one of 0, 1, or 2.
#' @param mu Mean vector (length n_time).
#' @param phi Dependence coefficient(s). For order 1: vector of length n_time-1.
#'   For order 2: matrix with 2 columns or vector of length 2*(n_time-2).
#' @param sigma Innovation standard deviations (length n_time).
#' @param blocks Integer vector of block membership (length n_subjects), or NULL.
#' @param tau Block effects, first element constrained to zero
#' @param na_action How to handle missing values:
#'   \itemize{
#'     \item marginalize: Compute observed-data likelihood (default)
#'     \item complete: Use only complete cases
#'     \item fail: Error if any NA present
#'   }
#'
#' @return Scalar log-likelihood value.
#'
#' @details
#' For complete data (no NA), all three na_action options give the same result.
#' 
#' For missing data:
#' \itemize{
#'   \item marginalize: Uses MVN marginalization to compute P(Y_obs). This is
#'     the correct observed-data likelihood for MAR missing data.
#'   \item complete: Removes subjects with any missing values. May lose information.
#'   \item fail: Stops with error. Useful to ensure no missing data present.
#' }
#'
#' @export
#' @importFrom stats dnorm
logL_ad <- function(y, order, mu, phi, sigma, blocks = NULL, tau = 0,
                    na_action = c("marginalize", "complete", "fail")) {
  
  na_action <- match.arg(na_action)
  
  # Input validation
  if (!is.matrix(y)) y <- as.matrix(y)
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Check for missing data
  has_missing <- any(is.na(y))
  
  if (has_missing) {
    if (na_action == "fail") {
      stop("y contains NA values. Use na_action = 'marginalize' or 'complete'",
           call. = FALSE)
    }
    
    if (na_action == "complete") {
      # Extract complete cases
      cc <- .extract_complete_cases(y, blocks, warn = FALSE)
      y <- cc$y
      blocks <- cc$blocks
      n_subjects <- nrow(y)
      has_missing <- FALSE  # After filtering
    }
    
    if (na_action == "marginalize" && has_missing) {
      # Use MVN marginalization
      return(.logL_ad_missing(y, order, mu, phi, sigma, blocks, tau))
    }
  }
  
  # ==== Complete data likelihood (original code) ====
  
  # Validate parameters
  if (length(mu) != n_time) stop("mu must have length n_time")
  if (length(sigma) != n_time) stop("sigma must have length n_time")
  if (any(sigma <= 0)) stop("sigma must be positive")
  
  if (order == 1) {
    if (length(phi) != n_time - 1) stop("phi must have length n_time - 1 for order 1")
  } else if (order == 2) {
    if (is.matrix(phi)) {
      if (nrow(phi) != n_time - 2 || ncol(phi) != 2) {
        stop("phi must be (n_time - 2) x 2 matrix for order 2")
      }
    } else {
      if (length(phi) != 2 * (n_time - 2)) {
        stop("phi must have length 2*(n_time - 2) for order 2")
      }
      # Convert to matrix
      phi <- matrix(phi, ncol = 2, byrow = TRUE)
    }
  }
  
  # Handle blocks
  if (!is.null(blocks)) {
    if (length(blocks) != n_subjects) stop("blocks must have length n_subjects")
    n_blocks <- length(unique(blocks))
    if (length(tau) != n_blocks) stop("tau must have length equal to number of blocks")
    if (abs(tau[1]) > 1e-10) {
      warning("tau[1] should be 0 for identifiability. Setting tau[1] = 0.")
      tau[1] <- 0
    }
  } else {
    tau <- 0
  }
  
  # Compute log-likelihood
  log_lik <- 0
  
  for (s in 1:n_subjects) {
    y_s <- y[s, ]
    
    # Adjust mean for block effect
    if (!is.null(blocks)) {
      mu_s <- mu + tau[blocks[s]]
    } else {
      mu_s <- mu
    }
    
    # Time 1 contribution
    log_lik <- log_lik + dnorm(y_s[1], mean = mu_s[1], sd = sigma[1], log = TRUE)
    
    if (order == 0) {
      # Independence: all times independent
      for (t in 2:n_time) {
        log_lik <- log_lik + dnorm(y_s[t], mean = mu_s[t], sd = sigma[t], log = TRUE)
      }
      
    } else if (order == 1) {
      # Order 1: Y_t | Y_{t-1} ~ N(mu_t + phi_t * (Y_{t-1} - mu_{t-1}), sigma_t^2)
      for (t in 2:n_time) {
        cond_mean <- mu_s[t] + phi[t-1] * (y_s[t-1] - mu_s[t-1])
        log_lik <- log_lik + dnorm(y_s[t], mean = cond_mean, sd = sigma[t], log = TRUE)
      }
      
    } else if (order == 2) {
      # Order 2
      # Time 2
      if (n_time >= 2) {
        cond_mean_2 <- mu_s[2] + phi[1, 1] * (y_s[1] - mu_s[1])
        log_lik <- log_lik + dnorm(y_s[2], mean = cond_mean_2, sd = sigma[2], log = TRUE)
      }
      
      # Times 3 to n
      for (t in 3:n_time) {
        cond_mean <- mu_s[t] + 
          phi[t-2, 1] * (y_s[t-1] - mu_s[t-1]) +
          phi[t-2, 2] * (y_s[t-2] - mu_s[t-2])
        log_lik <- log_lik + dnorm(y_s[t], mean = cond_mean, sd = sigma[t], log = TRUE)
      }
    }
  }
  
  log_lik
}

#' Count parameters for AD model
#'
#' @param order Antedependence order
#' @param n_time Number of time points
#' @param n_blocks Number of blocks (default 1)
#'
#' @return Number of free parameters
#'
#' @keywords internal
.count_params_ad <- function(order, n_time, n_blocks = 1) {
  # mu: n_time parameters
  n_mu <- n_time
  
  # phi: depends on order
  if (order == 0) {
    n_phi <- 0
  } else if (order == 1) {
    n_phi <- n_time - 1
  } else if (order == 2) {
    n_phi <- 2 * (n_time - 2)
  } else {
    stop("order must be 0, 1, or 2")
  }
  
  # sigma: n_time parameters
  n_sigma <- n_time
  
  # tau: (n_blocks - 1) parameters (tau[1] = 0)
  n_tau <- max(0, n_blocks - 1)
  
  n_mu + n_phi + n_sigma + n_tau
}
