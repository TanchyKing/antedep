# lrt_stationarity_cat.R - Likelihood ratio test for stationarity in categorical AD

#' Likelihood ratio test for stationarity (categorical data)
#'
#' Tests whether a categorical antedependence process is strictly stationary,
#' meaning both the marginal distribution and transition probabilities are
#' constant over time.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c.
#' @param order Antedependence order p. Default is 1.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership.
#' @param homogeneous Logical. If TRUE (default), parameters are shared across
#'   all groups.
#' @param n_categories Number of categories. If NULL, inferred from data.
#'
#' @return A list of class \code{"cat_lrt"} containing:
#'   \item{lrt_stat}{Likelihood ratio test statistic}
#'   \item{df}{Degrees of freedom}
#'   \item{p_value}{P-value from chi-square distribution}
#'   \item{fit_null}{Fitted stationary model (H0)}
#'   \item{fit_alt}{Fitted non-stationary model (H1)}
#'   \item{table}{Summary data frame}
#'
#' @details
#' Strict stationarity requires:
#' \enumerate{
#'   \item The marginal distribution P(Yk) is constant for all k
#'   \item The transition probabilities P(Yk | Y(k-p), ..., Y(k-1)) are
#'     constant for all k > p
#' }
#'
#' This is stronger than time-invariance, which only requires condition 2.
#'
#' The null hypothesis is tested against the general (non-stationary) AD(p)
#' model. The degrees of freedom are:
#' \deqn{df = (c-1)(n-1) + (c-1)c^p(n-p-1)}
#' for order p >= 1, which accounts for both marginal and transition constraints.
#'
#' @examples
#' \dontrun{
#' # Simulate stationary AD(1) data
#' set.seed(123)
#' y <- simulate_cat(200, 6, order = 1, n_categories = 2)
#'
#' # Test stationarity
#' test <- lrt_stationarity_cat(y, order = 1)
#' print(test)
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @seealso \code{\link{lrt_timeinvariance_cat}}, \code{\link{lrt_order_cat}}
#'
#' @export
lrt_stationarity_cat <- function(y, order = 1, blocks = NULL,
                                  homogeneous = TRUE, n_categories = NULL) {
  if (anyNA(y)) {
    .stop_cat_missing_inference("lrt_stationarity_cat")
  }
  
  # Validate data
 validated <- .validate_y_cat(y, n_categories)
  y <- validated$y
  c <- validated$n_categories
  n_time <- ncol(y)
  n_subjects <- nrow(y)
  p <- as.integer(order)
  
  if (p < 0) {
    stop("order must be non-negative")
  }
  if (p >= n_time) {
    stop("order must be less than n_time")
  }
  
  # Fit unconstrained (non-stationary) model - this is the alternative
  fit_alt <- fit_cat(y, order = p, blocks = blocks,
                     homogeneous = homogeneous, n_categories = c)
  
  # Fit constrained (stationary) model
  fit_null <- .fit_cat_stationary(y, p, c, blocks, homogeneous)
  
  # Compute LRT statistic
  lrt_stat <- -2 * (fit_null$log_l - fit_alt$log_l)
  
  # Degrees of freedom
  df <- fit_alt$n_params - fit_null$n_params
  
  # P-value from chi-square
  p_value <- stats::pchisq(lrt_stat, df = df, lower.tail = FALSE)
  
  # Build summary table
  table_df <- data.frame(
    model = c("Stationary (H0)", "Non-stationary (H1)"),
    log_l = c(fit_null$log_l, fit_alt$log_l),
    n_params = c(fit_null$n_params, fit_alt$n_params),
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
    order = p,
    table = table_df
  )
  
  class(out) <- "cat_lrt"
  out
}


#' Fit stationary categorical AD model
#'
#' Internal function to fit a model where both marginal and transition
#' probabilities are constant over time.
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param blocks Block vector (or NULL)
#' @param homogeneous Whether to pool across blocks
#'
#' @return A cat_fit-like object with stationary constraints
#'
#' @keywords internal
.fit_cat_stationary <- function(y, p, c, blocks = NULL, homogeneous = TRUE) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  
  # Process blocks
  if (is.null(blocks)) {
    blocks <- rep(1L, n_subjects)
  }
  n_blocks <- max(blocks)
  
  if (homogeneous) {
    n_pops <- 1
  } else {
    n_pops <- n_blocks
  }
  
  # For each population, fit stationary model
  if (n_pops == 1) {
    result <- .fit_cat_stationary_single(y, p, c, n_time, subject_mask = NULL)
    marginal <- result$marginal
    transition <- result$transition
    log_l <- result$log_l
    n_params <- result$n_params
  } else {
    marginal <- vector("list", n_pops)
    transition <- vector("list", n_pops)
    log_l <- 0
    n_params <- 0
    
    for (g in seq_len(n_pops)) {
      mask <- (blocks == g)
      result_g <- .fit_cat_stationary_single(y, p, c, n_time, subject_mask = mask)
      marginal[[g]] <- result_g$marginal
      transition[[g]] <- result_g$transition
      log_l <- log_l + result_g$log_l
      n_params <- n_params + result_g$n_params
    }
    names(marginal) <- paste0("block_", seq_len(n_pops))
    names(transition) <- paste0("block_", seq_len(n_pops))
  }
  
  # Compute AIC and BIC
  aic <- -2 * log_l + 2 * n_params
  bic <- -2 * log_l + n_params * log(n_subjects)
  
  out <- list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    aic = aic,
    bic = bic,
    n_params = n_params,
    convergence = 0L,
    settings = list(
      order = p,
      n_categories = c,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = if (n_blocks > 1) blocks else NULL,
      homogeneous = homogeneous,
      n_blocks = n_blocks,
      stationary = TRUE
    )
  )
  
  class(out) <- "cat_fit"
  out
}


#' Fit stationary model for single population
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param subject_mask Logical mask for subjects (NULL = all)
#'
#' @return List with marginal, transition, log_l, n_params
#'
#' @keywords internal
.fit_cat_stationary_single <- function(y, p, c, n_time, subject_mask = NULL) {
  
  # Apply mask
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  N <- nrow(y_sub)
  
  log_l <- 0
  
  if (p == 0) {
    # Independence + stationarity = same marginal at all time points
    # Pool all observations across all time points
    pooled_counts <- rep(0L, c)
    for (k in seq_len(n_time)) {
      for (cat in seq_len(c)) {
        pooled_counts[cat] <- pooled_counts[cat] + sum(y_sub[, k] == cat)
      }
    }
    
    # Estimate common marginal
    total <- sum(pooled_counts)
    pooled_probs <- pooled_counts / total
    
    marginal <- list(stationary = pooled_probs)
    names(marginal$stationary) <- paste0("cat_", 1:c)
    transition <- NULL
    
    # Log-likelihood: evaluate at each time point with pooled probs
    for (k in seq_len(n_time)) {
      counts_k <- .count_cells_table_cat(y_sub, k, c, subject_mask = NULL)
      log_l <- log_l + .loglik_contribution(counts_k, pooled_probs)
    }
    
    # Number of parameters: just (c-1) for the common marginal
    n_params <- c - 1
    
  } else {
    # AD(p) with p >= 1
    # Stationarity means:
    # 1. Marginal P(Y_k) is constant for all k
    # 2. Transition P(Y_k | Y_{k-p},...,Y_{k-1}) is constant for k > p
    
    # Step 1: Estimate stationary marginal by pooling across all time points
    pooled_marginal_counts <- rep(0L, c)
    for (k in seq_len(n_time)) {
      for (cat in seq_len(c)) {
        pooled_marginal_counts[cat] <- pooled_marginal_counts[cat] + sum(y_sub[, k] == cat)
      }
    }
    stationary_marginal <- pooled_marginal_counts / sum(pooled_marginal_counts)
    
    # Step 2: Estimate stationary transition by pooling across k = p+1 to n
    pooled_trans_counts <- array(0L, dim = rep(c, p + 1))
    for (k in (p + 1):n_time) {
      time_indices <- (k - p):k
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      pooled_trans_counts <- pooled_trans_counts + counts_k
    }
    stationary_trans <- .counts_to_transition_probs(pooled_trans_counts)
    
    # Store parameters
    marginal <- list(stationary = stationary_marginal)
    names(marginal$stationary) <- paste0("cat_", 1:c)
    transition <- list(stationary = stationary_trans)
    
    # Step 3: Compute log-likelihood
    # For k = 1: use stationary marginal
    counts_1 <- .count_cells_table_cat(y_sub, 1, c, subject_mask = NULL)
    log_l <- log_l + .loglik_contribution(counts_1, stationary_marginal)
    
    # For k = 2 to p: need joint probabilities under stationarity
    # Under stationarity, P(Y_1,...,Y_k) = P(Y_1) * prod_{j=2}^{k} P(Y_j|Y_{j-1},...,Y_{j-min(j-1,p)})
    # For simplicity, we use the stationary marginal for the marginal contribution
    # and compute the conditional part
    if (p >= 2 && n_time >= 2) {
      for (k in 2:min(p, n_time)) {
        # For these early time points, we need to handle carefully
        # Use stationary marginal for the marginal part
        counts_k <- .count_cells_table_cat(y_sub, k, c, subject_mask = NULL)
        log_l <- log_l + .loglik_contribution(counts_k, stationary_marginal)
      }
    }
    
    # For k = p+1 to n: use stationary transitions
    for (k in (p + 1):n_time) {
      time_indices <- (k - p):k
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      log_l <- log_l + .loglik_contribution(counts_k, stationary_trans)
    }
    
    # Number of parameters:
    # - Stationary marginal: (c-1)
    # - Stationary transition: (c-1) * c^p
    n_params <- (c - 1) + (c - 1) * c^p
  }
  
  list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    n_params = n_params
  )
}


#' Run all stationarity-related tests for categorical AD
#'
#' Performs tests for time-invariance and strict stationarity.
#'
#' @param y Integer matrix of categorical data (n_subjects x n_time).
#' @param order Antedependence order p. Default is 1.
#' @param blocks Optional block membership vector.
#' @param homogeneous Whether to use homogeneous parameters across blocks.
#' @param n_categories Number of categories (inferred if NULL).
#'
#' @return A list containing:
#'   \item{time_invariance}{Result of lrt_timeinvariance_cat}
#'   \item{stationarity}{Result of lrt_stationarity_cat}
#'   \item{table}{Summary data frame}
#'
#' @examples
#' \dontrun{
#' y <- simulate_cat(200, 6, order = 1, n_categories = 2)
#' result <- run_stationarity_tests_cat(y, order = 1)
#' print(result$table)
#' }
#'
#' @export
run_stationarity_tests_cat <- function(y, order = 1, blocks = NULL,
                                        homogeneous = TRUE, n_categories = NULL) {
  if (anyNA(y)) {
    .stop_cat_missing_inference("run_stationarity_tests_cat")
  }
  
  # Run time-invariance test
  test_ti <- lrt_timeinvariance_cat(y, order = order, blocks = blocks,
                                     homogeneous = homogeneous, 
                                     n_categories = n_categories)
  
  # Run stationarity test
  test_stat <- lrt_stationarity_cat(y, order = order, blocks = blocks,
                                     homogeneous = homogeneous,
                                     n_categories = n_categories)
  
  # Build summary table
  table_df <- data.frame(
    test = c("Time-invariance", "Stationarity"),
    lrt_stat = c(test_ti$lrt_stat, test_stat$lrt_stat),
    df = c(test_ti$df, test_stat$df),
    p_value = c(test_ti$p_value, test_stat$p_value),
    significant = c(test_ti$p_value < 0.05, test_stat$p_value < 0.05),
    stringsAsFactors = FALSE
  )
  
  list(
    time_invariance = test_ti,
    stationarity = test_stat,
    table = table_df
  )
}
