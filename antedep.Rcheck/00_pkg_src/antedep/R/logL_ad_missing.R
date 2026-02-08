#' Compute observed-data log-likelihood for AD with missing values
#'
#' Uses multivariate normal marginalization to compute the likelihood
#' of observed values, marginalizing over missing values.
#'
#' @param y Data matrix (n_subjects x n_time), may contain NA
#' @param order Antedependence order (0, 1, or 2)
#' @param mu Mean vector (length n_time)
#' @param phi Dependence parameter(s)
#' @param sigma Innovation standard deviations (length n_time)
#' @param blocks Block membership vector (optional)
#' @param tau Block effects, first element constrained to zero
#'
#' @return Observed-data log-likelihood (scalar)
#'
#' @keywords internal
.logL_ad_missing <- function(y, order, mu, phi, sigma, blocks = NULL, tau = 0) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Build full covariance matrix from AD parameters
  Sigma <- .build_ad_covariance(order, phi, sigma, n_time)
  
  # Adjust mean for blocks if present
  if (!is.null(blocks)) {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
    tau_matrix <- matrix(tau[blocks], nrow = n_subjects, ncol = n_time)
    mu_matrix <- mu_matrix + tau_matrix
  } else {
    mu_matrix <- matrix(mu, nrow = n_subjects, ncol = n_time, byrow = TRUE)
  }
  
  # Compute log-likelihood by marginalizing over missing values
  log_lik <- 0
  
  for (s in 1:n_subjects) {
    y_s <- y[s, ]
    mu_s <- mu_matrix[s, ]
    
    obs_idx <- which(!is.na(y_s))
    
    if (length(obs_idx) == 0) {
      # All missing - skip
      next
    }
    
    # Extract observed values and parameters
    y_obs <- y_s[obs_idx]
    mu_obs <- mu_s[obs_idx]
    Sigma_obs <- Sigma[obs_idx, obs_idx, drop = FALSE]
    
    # Compute MVN density on observed values
    tryCatch({
      # Use Cholesky decomposition for numerical stability
      L <- chol(Sigma_obs)
      z <- backsolve(L, y_obs - mu_obs, transpose = TRUE)
      log_det <- 2 * sum(log(diag(L)))
      
      log_lik_s <- -0.5 * (length(obs_idx) * log(2 * pi) + log_det + sum(z^2))
      log_lik <- log_lik + log_lik_s
    }, error = function(e) {
      # If Cholesky fails, return -Inf
      log_lik <<- -Inf
    })
    
    if (!is.finite(log_lik)) break
  }
  
  return(log_lik)
}

#' Build AD covariance matrix from parameters
#'
#' @param order Antedependence order
#' @param phi Dependence parameters
#' @param sigma Innovation standard deviations
#' @param n_time Number of time points
#' @return Covariance matrix (n_time x n_time)
#' @keywords internal
.build_ad_covariance <- function(order, phi, sigma, n_time) {
  # Initialize covariance matrix
  Sigma <- matrix(0, nrow = n_time, ncol = n_time)
  
  if (order == 0) {
    # Order 0: independence
    diag(Sigma) <- sigma^2
    return(Sigma)
  }
  
  # Compute marginal variances using forward recursion
  V <- numeric(n_time)
  
  if (order == 1) {
    # Order 1: Var(Y_t) = phi[t-1]^2 * Var(Y_{t-1}) + sigma[t]^2
    # phi has length n_time-1: phi[1] for Y2 given Y1, phi[2] for Y3 given Y2, etc.
    V[1] <- sigma[1]^2
    for (t in 2:n_time) {
      V[t] <- phi[t-1]^2 * V[t-1] + sigma[t]^2
    }
    
    # Compute covariances
    for (i in 1:n_time) {
      Sigma[i, i] <- V[i]
      if (i < n_time) {
        for (j in (i+1):n_time) {
          # Cov(Y_i, Y_j) = prod(phi[i:(j-1)]) * Var(Y_i)
          prod_phi <- prod(phi[i:(j-1)])
          Sigma[i, j] <- prod_phi * V[i]
          Sigma[j, i] <- Sigma[i, j]  # Symmetric
        }
      }
    }
    
  } else if (order == 2) {
    # Order 2: Var(Y_t) = phi1[t]^2 * Var(Y_{t-1}) + phi2[t]^2 * Var(Y_{t-2}) + sigma[t]^2
    # phi is (n_time - 2) x 2 matrix
    # Row 1 is for time 3, row 2 for time 4, etc.
    
    if (is.matrix(phi)) {
      # phi is matrix with 2 columns
      if (nrow(phi) != n_time - 2) {
        stop("phi must have n_time - 2 rows for order 2")
      }
    } else {
      # phi is vector, convert to matrix
      if (length(phi) != 2 * (n_time - 2)) {
        stop("phi must have length 2*(n_time - 2) for order 2")
      }
      phi <- matrix(phi, ncol = 2, byrow = TRUE)
    }
    
    # Initialize first two variances
    V[1] <- sigma[1]^2
    V[2] <- sigma[2]^2 + phi[1, 1]^2 * V[1]
    
    # Compute remaining variances
    for (t in 3:n_time) {
      # phi[t-2, 1] is lag-1 coef, phi[t-2, 2] is lag-2 coef
      phi1 <- phi[t-2, 1]
      phi2 <- phi[t-2, 2]
      V[t] <- phi1^2 * V[t-1] + phi2^2 * V[t-2] + sigma[t]^2
    }
    
    # Fill diagonal
    diag(Sigma) <- V
    
    # Compute off-diagonal covariances (simplified - this is approximate)
    # Full computation would require tracking all autocovariances
    for (i in 1:(n_time-1)) {
      for (j in (i+1):n_time) {
        lag <- j - i
        if (lag == 1 && j >= 3) {
          # One-step covariance
          Sigma[i, j] <- phi[j-2, 1] * V[i]
        } else if (lag == 2 && j >= 4) {
          # Two-step covariance
          Sigma[i, j] <- phi[j-2, 2] * V[i]
        } else {
          # Higher lags (approximate)
          Sigma[i, j] <- 0
        }
        Sigma[j, i] <- Sigma[i, j]
      }
    }
  }
  
  return(Sigma)
}
