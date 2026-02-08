# fit_cat.R - Maximum likelihood estimation for categorical antedependence models

#' Fit a categorical antedependence model
#'
#' Computes maximum likelihood estimates for the parameters of an AD(p) model
#' for categorical longitudinal data. The model is parameterized by transition
#' probabilities, and MLEs are obtained in closed form.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Each entry
#'   should be a category code from 1 to c, where c is the number of categories.
#' @param order Antedependence order p. Must be 0, 1, or 2. Default is 1.
#' @param blocks Optional integer vector of length n_subjects specifying group
#'   membership. If NULL, all subjects are treated as one group.
#' @param homogeneous Logical. If TRUE (default), parameters are shared across
#'   all groups (blocks are ignored for estimation). If FALSE, separate
#'   transition probabilities are estimated for each group.
#' @param n_categories Number of categories. If NULL (default), inferred from
#'   the maximum value in y.
#'
#' @return A list of class \code{"cat_fit"} containing:
#'   \item{marginal}{List of marginal/joint probabilities for initial time points}
#'   \item{transition}{List of transition probability arrays for k = p+1 to n}
#'   \item{log_l}{Log-likelihood at MLE}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{n_params}{Number of free parameters}
#'   \item{cell_counts}{List of observed cell counts}
#'   \item{convergence}{Always 0 (closed-form solution)}
#'   \item{settings}{List of model settings}
#'
#' @details
#' For AD(p), the model decomposes as:
#' \deqn{P(Y_1, \ldots, Y_n) = P(Y_1, \ldots, Y_p) \times \prod_{k=p+1}^{n} P(Y_k | Y_{k-p}, \ldots, Y_{k-1})}
#'
#' MLEs are computed as empirical proportions:
#' \itemize{
#'   \item Marginal/joint probabilities: count / N
#'   \item Transition probabilities: conditional count / marginal count
#' }
#'
#' Empty cells receive probability 0 (if denominator is also 0).
#'
#' @examples
#' \dontrun{
#' # Simulate binary AD(1) data
#' set.seed(123)
#' y <- simulate_cat(n_subjects = 100, n_time = 5, order = 1, n_categories = 2)
#'
#' # Fit model
#' fit <- fit_cat(y, order = 1)
#' print(fit)
#'
#' # Compare orders
#' fit0 <- fit_cat(y, order = 0)
#' fit1 <- fit_cat(y, order = 1)
#' fit2 <- fit_cat(y, order = 2)
#' c(AIC_0 = fit0$aic, AIC_1 = fit1$aic, AIC_2 = fit2$aic)
#' }
#'
#' @references
#' Xie, Y. and Zimmerman, D. L. (2013). Antedependence models for nonstationary
#' categorical longitudinal data with ignorable missingness: likelihood-based
#' inference. \emph{Statistics in Medicine}, 32, 3274-3289.
#'
#' @export
fit_cat <- function(y, order = 1, blocks = NULL, homogeneous = TRUE, 
                    n_categories = NULL) {
  

  # Validate inputs
  validated <- .validate_y_cat(y, n_categories)
  y <- validated$y
  c <- validated$n_categories
  
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  p <- as.integer(order)
  
  # Validate order
  if (p < 0) stop("order must be non-negative")
  if (p > 2) stop("order > 2 not currently supported")
  if (p >= n_time) stop("order must be less than number of time points")
  
  # Validate and process blocks
  blocks <- .validate_blocks_cat(blocks, n_subjects)
  n_blocks <- max(blocks)
  
  # Determine number of populations to estimate
  if (homogeneous) {
    n_pops <- 1
  } else {
    n_pops <- n_blocks
  }
  
  # Count parameters
  n_params <- .count_params_cat(p, c, n_time, n_blocks, homogeneous)
  
  # Initialize storage for parameters
  if (n_pops == 1) {
    result <- .fit_cat_single_pop(y, p, c, n_time, subject_mask = NULL)
    marginal <- result$marginal
    transition <- result$transition
    cell_counts <- result$cell_counts
    log_l <- result$log_l
  } else {
    # Fit separately for each block
    marginal <- vector("list", n_pops)
    transition <- vector("list", n_pops)
    cell_counts <- vector("list", n_pops)
    log_l <- 0
    
    for (g in seq_len(n_pops)) {
      mask <- (blocks == g)
      result_g <- .fit_cat_single_pop(y, p, c, n_time, subject_mask = mask)
      marginal[[g]] <- result_g$marginal
      transition[[g]] <- result_g$transition
      cell_counts[[g]] <- result_g$cell_counts
      log_l <- log_l + result_g$log_l
    }
    names(marginal) <- paste0("block_", seq_len(n_pops))
    names(transition) <- paste0("block_", seq_len(n_pops))
    names(cell_counts) <- paste0("block_", seq_len(n_pops))
  }
  
  # Compute AIC and BIC
  aic <- -2 * log_l + 2 * n_params
  bic <- -2 * log_l + n_params * log(n_subjects)
  
  # Assemble output
  out <- list(
    marginal = marginal,
    transition = transition,
    log_l = log_l,
    aic = aic,
    bic = bic,
    n_params = n_params,
    cell_counts = cell_counts,
    convergence = 0L,
    settings = list(
      order = p,
      n_categories = c,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = if (n_blocks > 1) blocks else NULL,
      homogeneous = homogeneous,
      n_blocks = n_blocks
    )
  )
  
  class(out) <- "cat_fit"
  out
}


#' Fit CAT model for a single population
#'
#' @param y Data matrix
#' @param p Order
#' @param c Number of categories
#' @param n_time Number of time points
#' @param subject_mask Logical vector for subject selection (NULL = all)
#'
#' @return List with marginal, transition, cell_counts, log_l
#'
#' @keywords internal
.fit_cat_single_pop <- function(y, p, c, n_time, subject_mask = NULL) {
  
  # Apply mask
  if (!is.null(subject_mask)) {
    y_sub <- y[subject_mask, , drop = FALSE]
  } else {
    y_sub <- y
  }
  N <- nrow(y_sub)
  
  # Storage
  marginal <- list()
  transition <- list()
  cell_counts <- list()
  log_l <- 0
  
  if (p == 0) {
    # Independence model: just estimate marginal at each time point
    for (k in seq_len(n_time)) {
      counts_k <- .count_cells_table_cat(y_sub, k, c, subject_mask = NULL)
      probs_k <- counts_k / N
      
      marginal[[k]] <- as.numeric(probs_k)
      names(marginal[[k]]) <- paste0("cat_", 1:c)
      cell_counts[[paste0("t", k)]] <- counts_k
      
      # Log-likelihood contribution
      log_l <- log_l + .loglik_contribution(counts_k, probs_k)
    }
    names(marginal) <- paste0("t", seq_len(n_time))
    
  } else {
    # AD(p) model with p >= 1
    
    # Handle initial time points (k = 1 to p)
    # For these, we need the full joint distribution P(Y_1, ..., Y_k)
    
    for (k in seq_len(min(p, n_time))) {
      time_indices <- seq_len(k)
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      
      if (k == 1) {
        # Marginal of Y_1
        probs_k <- counts_k / N
        marginal[["t1"]] <- as.numeric(probs_k)
        names(marginal[["t1"]]) <- paste0("cat_", 1:c)
      } else {
        # P(Y_k | Y_1, ..., Y_{k-1})
        # counts_k has dimensions c^k
        # Need to compute conditional probabilities
        probs_k <- .counts_to_transition_probs(counts_k)
        marginal[[paste0("t", k, "_given_1to", k-1)]] <- probs_k
      }
      
      cell_counts[[paste0("t1_to_t", k)]] <- counts_k
      
      # Log-likelihood contribution for this time point
      # For k=1: sum N_y1 * log(pi_y1)
      # For k>1: sum N_{y1...yk} * log(pi_{yk|y1...y_{k-1}})
      if (k == 1) {
        log_l <- log_l + .loglik_contribution(counts_k, probs_k)
      } else {
        # Need joint counts and transition probs
        # LL contribution = sum over all (y1,...,yk) of N * log(pi_{yk|y1...y_{k-1}})
        log_l <- log_l + .loglik_contribution(counts_k, probs_k)
      }
    }
    
    # Handle time points k = p+1 to n (only condition on last p values)
    for (k in (p + 1):n_time) {
      # Count (Y_{k-p}, ..., Y_{k-1}, Y_k) combinations
      time_indices <- (k - p):k
      counts_k <- .count_cells_table_cat(y_sub, time_indices, c, subject_mask = NULL)
      
      # Compute transition probabilities P(Y_k | Y_{k-p}, ..., Y_{k-1})
      probs_k <- .counts_to_transition_probs(counts_k)
      
      transition[[paste0("t", k)]] <- probs_k
      cell_counts[[paste0("t", k-p, "_to_t", k)]] <- counts_k
      
      # Log-likelihood contribution
      log_l <- log_l + .loglik_contribution(counts_k, probs_k)
    }
  }
  
  list(
    marginal = marginal,
    transition = transition,
    cell_counts = cell_counts,
    log_l = log_l
  )
}
