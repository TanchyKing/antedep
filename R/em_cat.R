# em_cat.R - EM estimation for categorical antedependence models

#' EM algorithm for categorical AD model estimation
#'
#' Fits categorical antedependence models with missing outcomes using the
#' Expectation-Maximization (EM) algorithm for orders 0 and 1.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns. Values are
#'   category codes in \code{1, ..., n_categories}; \code{NA} is allowed.
#' @param order Antedependence order. Supported values are \code{0} and
#'   \code{1}. Order \code{2} is not yet implemented in \code{em_cat()}.
#' @param blocks Optional block/group vector of length \code{n_subjects}. Any
#'   coding is accepted (e.g., non-sequential integers or factor levels).
#' @param homogeneous Logical. If \code{TRUE}, a single parameter set is fitted
#'   across blocks. If \code{FALSE}, separate parameters are fitted by block.
#' @param n_categories Number of categories. If \code{NULL}, inferred from
#'   observed data.
#' @param max_iter Maximum number of EM iterations.
#' @param tol Convergence tolerance on absolute log-likelihood change.
#' @param epsilon Small positive constant used for smoothing and numerical
#'   stability.
#' @param verbose Logical; if \code{TRUE}, print EM progress.
#'
#' @return A \code{cat_fit} object with fields matching \code{\link{fit_cat}}.
#'
#' @details
#' For complete data (no missing values), this function defers to
#' \code{\link{fit_cat}} with closed-form MLEs.
#'
#' For missing data and orders 0/1, each EM iteration computes expected
#' sufficient statistics with a forward-backward E-step, then updates
#' probabilities by normalized expected counts in the M-step.
#'
#' A final E-step is run before returning so that \code{log_l}/AIC/BIC and
#' expected cell counts correspond exactly to the returned parameter values.
#'
#' @seealso \code{\link{fit_cat}}, \code{\link{logL_cat}}
#' @export
em_cat <- function(y, order = 1, blocks = NULL, homogeneous = TRUE,
                   n_categories = NULL, max_iter = 100, tol = 1e-6,
                   epsilon = 1e-8, verbose = FALSE) {
  validated <- .validate_y_cat(y, n_categories = n_categories, allow_na = TRUE)
  y <- validated$y
  c <- validated$n_categories

  n_subjects <- nrow(y)
  n_time <- ncol(y)
  p <- as.integer(order)

  if (!p %in% c(0L, 1L, 2L)) stop("order must be 0, 1, or 2")
  if (p >= n_time) stop("order must be less than number of time points")
  if (max_iter < 1) stop("max_iter must be >= 1")
  if (!is.finite(tol) || tol <= 0) stop("tol must be a positive finite number")
  if (!is.finite(epsilon) || epsilon <= 0) stop("epsilon must be a positive finite number")

  if (p == 2L) {
    stop(
      paste0(
        "em_cat with order=2 is not yet implemented.\n",
        "Use fit_cat(..., na_action = 'marginalize') for order 2 with missing data."
      )
    )
  }

  block_info <- .normalize_blocks_em_cat(blocks, n_subjects)
  blocks_id <- block_info$blocks_id
  block_levels <- block_info$block_levels
  n_blocks <- block_info$n_blocks

  # Complete data: defer to closed-form MLEs.
  if (!any(is.na(y))) {
    if (verbose) message("No missing data; using closed-form MLEs via fit_cat().")
    fit <- fit_cat(
      y = y,
      order = p,
      blocks = blocks_id,
      homogeneous = homogeneous,
      n_categories = c,
      na_action = "complete"
    )
    fit$settings$method <- "closed_form"
    fit$settings$block_levels <- block_levels
    return(fit)
  }

  init <- .initialize_parameters_cat_em(
    y = y, order = p, blocks_id = blocks_id, homogeneous = homogeneous,
    c = c, n_time = n_time, epsilon = epsilon, verbose = verbose
  )
  marginal <- init$marginal
  transition <- init$transition

  max_iter <- as.integer(max_iter)
  converged <- FALSE
  iter <- 0L
  logL_prev <- -Inf
  logL_current <- -Inf
  delta_logL <- Inf

  for (iter in seq_len(max_iter)) {
    e_result <- .e_step_cat_em(
      y = y, order = p, marginal = marginal, transition = transition,
      blocks_id = blocks_id, homogeneous = homogeneous, c = c
    )

    if (isTRUE(e_result$failed)) {
      stop(sprintf("EM iteration %d: %s", iter, e_result$message))
    }

    logL_current <- e_result$log_likelihood
    delta_logL <- logL_current - logL_prev

    if (verbose) {
      message(sprintf("Iteration %d: logL = %.6f, change = %.6f", iter, logL_current, delta_logL))
    }

    # Convergence check is done before M-step.
    if (iter > 1L && abs(delta_logL) < tol) {
      converged <- TRUE
      break
    }

    if (delta_logL < -1e-6) {
      warning(
        sprintf(
          "Likelihood decreased at iteration %d (%.6f -> %.6f); possible numerical issue.",
          iter, logL_prev, logL_current
        )
      )
    }

    logL_prev <- logL_current

    params_new <- .m_step_cat_em(
      counts = e_result$counts, order = p, c = c, n_time = n_time,
      n_blocks = n_blocks, homogeneous = homogeneous, epsilon = epsilon
    )
    marginal <- params_new$marginal
    transition <- params_new$transition
  }

  # Final E-step ensures returned fit metrics/cell counts correspond exactly
  # to returned parameter values, including max_iter exits.
  e_final <- .e_step_cat_em(
    y = y, order = p, marginal = marginal, transition = transition,
    blocks_id = blocks_id, homogeneous = homogeneous, c = c
  )
  if (isTRUE(e_final$failed)) {
    stop(sprintf("Final likelihood evaluation failed: %s", e_final$message))
  }
  logL_current <- e_final$log_likelihood
  delta_logL <- logL_current - logL_prev
  cell_counts <- .expected_counts_to_cell_counts_cat_em(
    counts = e_final$counts,
    order = p,
    n_time = n_time,
    n_blocks = n_blocks,
    homogeneous = homogeneous
  )

  if (!converged && iter == max_iter && abs(delta_logL) >= tol) {
    warning(sprintf("EM did not converge after %d iterations (change = %.6f).", max_iter, delta_logL))
  }

  n_params <- .count_params_cat(
    order = p,
    n_categories = c,
    n_time = n_time,
    n_blocks = n_blocks,
    homogeneous = homogeneous
  )

  aic <- -2 * logL_current + 2 * n_params
  bic <- -2 * logL_current + log(n_subjects) * n_params

  out <- list(
    marginal = marginal,
    transition = transition,
    log_l = as.numeric(logL_current),
    aic = as.numeric(aic),
    bic = as.numeric(bic),
    n_params = as.integer(n_params),
    cell_counts = cell_counts,
    convergence = if (converged) 0L else 1L,
    settings = list(
      order = p,
      n_categories = c,
      n_time = n_time,
      n_subjects = n_subjects,
      blocks = if (n_blocks > 1) blocks_id else NULL,
      homogeneous = homogeneous,
      n_blocks = n_blocks,
      na_action = "em",
      na_action_effective = "em",
      method = "em",
      max_iter = max_iter,
      tol = tol,
      epsilon = epsilon,
      n_iter = as.integer(iter),
      converged = converged,
      delta_log_l = as.numeric(delta_logL),
      block_levels = block_levels
    )
  )

  class(out) <- "cat_fit"
  out
}


#' @keywords internal
.normalize_blocks_em_cat <- function(blocks, n_subjects) {
  if (is.null(blocks)) {
    return(list(blocks_id = rep(1L, n_subjects), block_levels = "1", n_blocks = 1L))
  }
  if (length(blocks) != n_subjects) {
    stop("blocks must have length equal to number of subjects (", n_subjects, ")")
  }
  if (any(is.na(blocks))) {
    stop("blocks must not contain missing values")
  }

  block_fac <- factor(blocks)
  list(
    blocks_id = as.integer(block_fac),
    block_levels = levels(block_fac),
    n_blocks = nlevels(block_fac)
  )
}


#' @keywords internal
.initialize_parameters_cat_em <- function(y, order, blocks_id, homogeneous,
                                          c, n_time, epsilon, verbose) {
  complete_subjects <- stats::complete.cases(y)
  n_complete <- sum(complete_subjects)
  n_blocks <- max(blocks_id)

  if (homogeneous || n_blocks == 1L) {
    if (n_complete >= 10L) {
      if (verbose) message("Initializing from complete cases.")
      return(.initialize_from_complete_cat_em(y[complete_subjects, , drop = FALSE], order, c, n_time, epsilon))
    }
    if (verbose) message("Initializing from random positive probabilities.")
    return(.initialize_uniform_cat_em(order, c, n_time))
  }

  marginal <- vector("list", n_blocks)
  transition <- vector("list", n_blocks)

  for (b in seq_len(n_blocks)) {
    idx_block <- which(blocks_id == b)
    idx_complete <- idx_block[complete_subjects[idx_block]]

    if (length(idx_complete) >= 5L) {
      params_b <- .initialize_from_complete_cat_em(y[idx_complete, , drop = FALSE], order, c, n_time, epsilon)
    } else {
      params_b <- .initialize_uniform_cat_em(order, c, n_time)
    }
    marginal[[b]] <- params_b$marginal
    transition[[b]] <- params_b$transition
  }

  names(marginal) <- paste0("block_", seq_len(n_blocks))
  names(transition) <- paste0("block_", seq_len(n_blocks))
  list(marginal = marginal, transition = transition)
}


#' @keywords internal
.initialize_from_complete_cat_em <- function(y_complete, order, c, n_time, epsilon) {
  if (order == 0L) {
    marginal <- list()
    for (t in seq_len(n_time)) {
      counts_t <- table(factor(y_complete[, t], levels = seq_len(c)))
      probs_t <- (as.numeric(counts_t) + epsilon) / (sum(counts_t) + c * epsilon)
      names(probs_t) <- paste0("cat_", seq_len(c))
      marginal[[paste0("t", t)]] <- probs_t
    }
    return(list(marginal = marginal, transition = list()))
  }

  counts_t1 <- table(factor(y_complete[, 1], levels = seq_len(c)))
  pi1 <- (as.numeric(counts_t1) + epsilon) / (sum(counts_t1) + c * epsilon)
  names(pi1) <- paste0("cat_", seq_len(c))
  marginal <- list(t1 = pi1)

  transition <- list()
  for (t in 2:n_time) {
    counts_t <- table(
      factor(y_complete[, t - 1], levels = seq_len(c)),
      factor(y_complete[, t], levels = seq_len(c))
    )
    counts_t <- counts_t + epsilon
    trans_t <- sweep(counts_t, 1, rowSums(counts_t), "/")
    transition[[paste0("t", t)]] <- unname(trans_t)
  }

  list(marginal = marginal, transition = transition)
}


#' @keywords internal
.initialize_uniform_cat_em <- function(order, c, n_time) {
  if (order == 0L) {
    marginal <- list()
    for (t in seq_len(n_time)) {
      probs <- rep(1 / c, c)
      names(probs) <- paste0("cat_", seq_len(c))
      marginal[[paste0("t", t)]] <- probs
    }
    return(list(marginal = marginal, transition = list()))
  }

  pi1 <- rep(1 / c, c)
  names(pi1) <- paste0("cat_", seq_len(c))
  marginal <- list(t1 = pi1)

  transition <- list()
  for (t in 2:n_time) {
    trans_t <- matrix(stats::rexp(c * c, rate = 1), nrow = c, ncol = c)
    trans_t <- sweep(trans_t, 1, rowSums(trans_t), "/")
    transition[[paste0("t", t)]] <- trans_t
  }

  list(marginal = marginal, transition = transition)
}


#' @keywords internal
.forward_backward_cat_order0 <- function(y_i, marginal, c) {
  n_time <- length(y_i)
  gamma <- matrix(0, nrow = n_time, ncol = c)
  scale <- numeric(n_time)
  log_likelihood <- 0

  for (t in seq_len(n_time)) {
    marg_t <- marginal[[paste0("t", t)]]
    if (is.null(marg_t)) marg_t <- marginal[[t]]
    if (is.null(marg_t) || length(marg_t) != c) {
      return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Invalid marginal at t=%d", t)))
    }
    marg_t <- as.numeric(marg_t)
    marg_t <- marg_t / sum(marg_t)

    if (is.na(y_i[t])) {
      gamma[t, ] <- marg_t
      scale[t] <- 1
    } else {
      y_obs <- as.integer(y_i[t])
      if (y_obs < 1L || y_obs > c) {
        return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Observed value out of range at t=%d", t)))
      }
      gamma[t, y_obs] <- 1
      scale[t] <- marg_t[y_obs]
    }

    if (!is.finite(scale[t]) || scale[t] <= 0) {
      return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Zero probability at t=%d", t)))
    }

    log_likelihood <- log_likelihood + log(scale[t])
  }

  list(
    gamma = gamma,
    xi = array(0, dim = c(0, c, c)),
    scale = scale,
    log_likelihood = log_likelihood,
    failed = FALSE
  )
}


#' @keywords internal
.forward_backward_cat_order1 <- function(y_i, marginal, transition, c) {
  n_time <- length(y_i)
  alpha <- matrix(0, nrow = n_time, ncol = c)
  beta <- matrix(0, nrow = n_time, ncol = c)
  scale <- numeric(n_time)
  emission <- matrix(1, nrow = n_time, ncol = c)

  for (t in seq_len(n_time)) {
    if (!is.na(y_i[t])) {
      y_obs <- as.integer(y_i[t])
      if (y_obs < 1L || y_obs > c) {
        return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Observed value out of range at t=%d", t)))
      }
      emission[t, ] <- 0
      emission[t, y_obs] <- 1
    }
  }

  pi1 <- as.numeric(marginal$t1)
  if (length(pi1) != c) {
    return(list(failed = TRUE, log_likelihood = -Inf, message = "Invalid t1 marginal"))
  }
  pi1 <- pi1 / sum(pi1)

  alpha[1, ] <- pi1 * emission[1, ]
  scale[1] <- sum(alpha[1, ])
  if (!is.finite(scale[1]) || scale[1] <= 0) {
    return(list(failed = TRUE, log_likelihood = -Inf, message = "Zero probability at t=1"))
  }
  alpha[1, ] <- alpha[1, ] / scale[1]

  if (n_time >= 2L) {
    for (t in 2:n_time) {
      trans_t <- transition[[paste0("t", t)]]
      if (is.null(trans_t) || any(dim(trans_t) != c(c, c))) {
        return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Invalid transition at t=%d", t)))
      }

      alpha[t, ] <- as.numeric(alpha[t - 1, ] %*% trans_t) * emission[t, ]
      scale[t] <- sum(alpha[t, ])
      if (!is.finite(scale[t]) || scale[t] <= 0) {
        return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Zero probability at t=%d", t)))
      }
      alpha[t, ] <- alpha[t, ] / scale[t]
    }
  }

  beta[n_time, ] <- 1
  if (n_time >= 2L) {
    for (t in (n_time - 1L):1L) {
      trans_next <- transition[[paste0("t", t + 1L)]]
      beta[t, ] <- as.numeric(trans_next %*% (emission[t + 1L, ] * beta[t + 1L, ])) / scale[t + 1L]
    }
  }

  gamma <- alpha * beta
  gamma_den <- rowSums(gamma)
  if (any(!is.finite(gamma_den) | gamma_den <= 0)) {
    return(list(failed = TRUE, log_likelihood = -Inf, message = "Failed to normalize posterior gamma"))
  }
  gamma <- gamma / gamma_den

  xi <- array(0, dim = c(n_time - 1L, c, c))
  if (n_time >= 2L) {
    for (t in 1:(n_time - 1L)) {
      trans_next <- transition[[paste0("t", t + 1L)]]
      w <- emission[t + 1L, ] * beta[t + 1L, ]

      for (a in seq_len(c)) {
        xi[t, a, ] <- alpha[t, a] * trans_next[a, ] * w
      }

      xi_sum <- sum(xi[t, , ])
      if (!is.finite(xi_sum) || xi_sum <= 0) {
        return(list(failed = TRUE, log_likelihood = -Inf, message = sprintf("Failed to normalize posterior xi at t=%d", t)))
      }
      xi[t, , ] <- xi[t, , ] / xi_sum
    }
  }

  list(
    gamma = gamma,
    xi = xi,
    scale = scale,
    log_likelihood = sum(log(scale)),
    failed = FALSE
  )
}


#' @keywords internal
.initialize_expected_counts_cat_em <- function(order, c, n_time, n_blocks, homogeneous) {
  if (order == 0L) {
    if (homogeneous || n_blocks == 1L) {
      return(list(marginal = array(0, dim = c(n_time, c)), transition = list()))
    }
    return(list(marginal = array(0, dim = c(n_blocks, n_time, c)), transition = list()))
  }

  if (homogeneous || n_blocks == 1L) {
    return(list(
      marginal = numeric(c),
      transition = array(0, dim = c(n_time - 1L, c, c))
    ))
  }

  list(
    marginal = array(0, dim = c(n_blocks, c)),
    transition = array(0, dim = c(n_blocks, n_time - 1L, c, c))
  )
}


#' @keywords internal
.accumulate_counts_cat_em <- function(counts, fb_result, order, block_i, homogeneous, n_time) {
  if (order == 0L) {
    if (homogeneous) {
      counts$marginal <- counts$marginal + fb_result$gamma
    } else {
      counts$marginal[block_i, , ] <- counts$marginal[block_i, , ] + fb_result$gamma
    }
    return(counts)
  }

  if (homogeneous) {
    counts$marginal <- counts$marginal + fb_result$gamma[1, ]
    counts$transition <- counts$transition + fb_result$xi
  } else {
    counts$marginal[block_i, ] <- counts$marginal[block_i, ] + fb_result$gamma[1, ]
    counts$transition[block_i, , , ] <- counts$transition[block_i, , , ] + fb_result$xi
  }

  counts
}


#' @keywords internal
.e_step_cat_em <- function(y, order, marginal, transition, blocks_id, homogeneous, c) {
  n_subjects <- nrow(y)
  n_time <- ncol(y)
  n_blocks <- max(blocks_id)

  expected_counts <- .initialize_expected_counts_cat_em(order, c, n_time, n_blocks, homogeneous)
  total_log_likelihood <- 0

  for (i in seq_len(n_subjects)) {
    y_i <- y[i, ]
    block_i <- blocks_id[i]

    if (homogeneous) {
      marg_i <- marginal
      trans_i <- transition
    } else {
      marg_i <- marginal[[block_i]]
      trans_i <- transition[[block_i]]
    }

    if (order == 0L) {
      fb_result <- .forward_backward_cat_order0(y_i, marg_i, c)
    } else {
      fb_result <- .forward_backward_cat_order1(y_i, marg_i, trans_i, c)
    }

    if (isTRUE(fb_result$failed)) {
      return(list(
        failed = TRUE,
        log_likelihood = -Inf,
        failed_subject = i,
        message = sprintf("Subject %d: %s", i, fb_result$message)
      ))
    }

    total_log_likelihood <- total_log_likelihood + fb_result$log_likelihood
    expected_counts <- .accumulate_counts_cat_em(
      counts = expected_counts, fb_result = fb_result, order = order,
      block_i = block_i, homogeneous = homogeneous, n_time = n_time
    )
  }

  list(
    counts = expected_counts,
    log_likelihood = total_log_likelihood,
    failed = FALSE
  )
}


#' Convert EM expected counts to fit_cat-style cell_counts
#'
#' @param counts Expected counts from the E-step.
#' @param order Model order.
#' @param n_time Number of time points.
#' @param n_blocks Number of blocks.
#' @param homogeneous Logical; whether parameters are shared across blocks.
#'
#' @return A list shaped like \code{fit_cat(...)} \code{cell_counts}.
#'
#' @keywords internal
.expected_counts_to_cell_counts_cat_em <- function(counts, order, n_time, n_blocks, homogeneous) {
  if (order == 0L) {
    if (homogeneous || n_blocks == 1L) {
      out <- list()
      for (t in seq_len(n_time)) {
        out[[paste0("t", t)]] <- as.numeric(counts$marginal[t, ])
      }
      return(out)
    }

    out <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
      out_b <- list()
      for (t in seq_len(n_time)) {
        out_b[[paste0("t", t)]] <- as.numeric(counts$marginal[b, t, ])
      }
      out[[b]] <- out_b
    }
    names(out) <- paste0("block_", seq_len(n_blocks))
    return(out)
  }

  if (order == 1L) {
    if (homogeneous || n_blocks == 1L) {
      out <- list()
      out[["t1_to_t1"]] <- as.numeric(counts$marginal)
      for (t in 2:n_time) {
        out[[paste0("t", t - 1L, "_to_t", t)]] <- unname(counts$transition[t - 1L, , ])
      }
      return(out)
    }

    out <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
      out_b <- list()
      out_b[["t1_to_t1"]] <- as.numeric(counts$marginal[b, ])
      for (t in 2:n_time) {
        out_b[[paste0("t", t - 1L, "_to_t", t)]] <- unname(counts$transition[b, t - 1L, , ])
      }
      out[[b]] <- out_b
    }
    names(out) <- paste0("block_", seq_len(n_blocks))
    return(out)
  }

  stop("Only orders 0 and 1 are supported in em_cat")
}


#' @keywords internal
.m_step_cat_em <- function(counts, order, c, n_time, n_blocks, homogeneous, epsilon) {
  if (order == 0L) {
    if (homogeneous || n_blocks == 1L) {
      marginal <- list()
      for (t in seq_len(n_time)) {
        probs_t <- counts$marginal[t, ] + epsilon
        probs_t <- probs_t / sum(probs_t)
        names(probs_t) <- paste0("cat_", seq_len(c))
        marginal[[paste0("t", t)]] <- probs_t
      }
      return(list(marginal = marginal, transition = list()))
    }

    marginal <- vector("list", n_blocks)
    transition <- vector("list", n_blocks)
    for (b in seq_len(n_blocks)) {
      marg_b <- list()
      for (t in seq_len(n_time)) {
        probs_t <- counts$marginal[b, t, ] + epsilon
        probs_t <- probs_t / sum(probs_t)
        names(probs_t) <- paste0("cat_", seq_len(c))
        marg_b[[paste0("t", t)]] <- probs_t
      }
      marginal[[b]] <- marg_b
      transition[[b]] <- list()
    }
    names(marginal) <- paste0("block_", seq_len(n_blocks))
    names(transition) <- paste0("block_", seq_len(n_blocks))
    return(list(marginal = marginal, transition = transition))
  }

  if (homogeneous || n_blocks == 1L) {
    pi1 <- counts$marginal + epsilon
    pi1 <- pi1 / sum(pi1)
    names(pi1) <- paste0("cat_", seq_len(c))
    marginal <- list(t1 = pi1)

    transition <- list()
    for (t in 2:n_time) {
      trans_t <- counts$transition[t - 1L, , ] + epsilon
      trans_t <- sweep(trans_t, 1, rowSums(trans_t), "/")
      transition[[paste0("t", t)]] <- trans_t
    }
    return(list(marginal = marginal, transition = transition))
  }

  marginal <- vector("list", n_blocks)
  transition <- vector("list", n_blocks)
  for (b in seq_len(n_blocks)) {
    pi1 <- counts$marginal[b, ] + epsilon
    pi1 <- pi1 / sum(pi1)
    names(pi1) <- paste0("cat_", seq_len(c))
    marginal[[b]] <- list(t1 = pi1)

    trans_b <- list()
    for (t in 2:n_time) {
      trans_t <- counts$transition[b, t - 1L, , ] + epsilon
      trans_t <- sweep(trans_t, 1, rowSums(trans_t), "/")
      trans_b[[paste0("t", t)]] <- trans_t
    }
    transition[[b]] <- trans_b
  }
  names(marginal) <- paste0("block_", seq_len(n_blocks))
  names(transition) <- paste0("block_", seq_len(n_blocks))

  list(marginal = marginal, transition = transition)
}
