simulate_inad_forward <- function(fit,
                                  history,
                                  blocks = NULL,
                                  start_time = ncol(history) + 1L,
                                  h = 1L,
                                  n_sims = 1000L,
                                  seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  history <- as.matrix(history)
  storage.mode(history) <- "integer"

  n_test <- nrow(history)
  n_obs <- ncol(history)
  start_time <- as.integer(start_time)
  h <- as.integer(h)
  n_sims <- as.integer(n_sims)

  if (start_time != n_obs + 1L) {
    stop("start_time must equal ncol(history) + 1 for this POC helper")
  }
  if (h < 1L || n_sims < 1L) stop("h and n_sims must be positive")

  settings <- fit$settings
  order <- as.integer(settings$order %||% 1L)
  thinning <- settings$thinning
  innovation <- settings$innovation
  n_time <- as.integer(settings$n_time %||% length(fit$theta))

  if (order != 1L) {
    stop("simulate_inad_forward POC helper currently supports order = 1 only")
  }
  if (start_time + h - 1L > n_time) {
    stop("Requested horizon exceeds fit time grid")
  }

  tau <- fit$tau
  has_tau <- !is.null(tau) && length(tau) > 1L && any(abs(as.numeric(tau[-1L])) > 1e-12)
  if (has_tau && is.null(blocks)) {
    stop("blocks must be supplied when fit has non-trivial tau")
  }
  if (is.null(blocks)) blocks <- rep(1L, n_test)
  blocks <- as.integer(blocks)
  if (length(blocks) != n_test) stop("blocks must have length nrow(history)")
  if (is.null(tau)) tau <- 0
  tau <- as.numeric(tau)
  if (any(blocks < 1L | blocks > length(tau))) {
    stop("blocks contain levels outside fit$tau")
  }

  alpha <- as.numeric(fit$alpha)
  theta <- as.numeric(fit$theta)
  if (length(alpha) == 1L) alpha <- rep(alpha, n_time)
  if (length(theta) == 1L) theta <- rep(theta, n_time)
  if (length(alpha) != n_time || length(theta) != n_time) {
    stop("fit alpha/theta lengths do not match fit time grid")
  }

  nb_inno_size <- fit$nb_inno_size
  if (innovation == "nbinom") {
    if (is.null(nb_inno_size)) stop("nb_inno_size is required for NB innovation")
    nb_inno_size <- as.numeric(nb_inno_size)
    if (length(nb_inno_size) == 1L) nb_inno_size <- rep(nb_inno_size, n_time)
    if (length(nb_inno_size) != n_time) stop("nb_inno_size length mismatch")
  }

  bell_theta_from_mean <- getFromNamespace(".bell_theta_from_mean", "antedep")

  draw_innovation <- function(t_index) {
    eff_mean <- if (innovation == "bell") {
      theta[t_index] * exp(theta[t_index]) + tau[blocks]
    } else {
      theta[t_index] + tau[blocks]
    }
    if (any(!is.finite(eff_mean)) || any(eff_mean <= 0)) {
      stop("Non-positive effective innovation mean")
    }
    if (innovation == "pois") {
      stats::rpois(n_test, lambda = eff_mean)
    } else if (innovation == "bell") {
      eff_theta <- bell_theta_from_mean(eff_mean)
      as.integer(vapply(eff_theta, function(th) antedep::rbell(1L, theta = th), integer(1L)))
    } else {
      stats::rnbinom(n_test, size = nb_inno_size[t_index], mu = eff_mean)
    }
  }

  draw_thinning <- function(y_prev, alpha_t) {
    if (alpha_t <= 0) return(integer(n_test))
    if (thinning == "binom") {
      stats::rbinom(n_test, size = y_prev, prob = alpha_t)
    } else if (thinning == "pois") {
      stats::rpois(n_test, lambda = alpha_t * y_prev)
    } else {
      out <- integer(n_test)
      idx <- y_prev > 0L
      if (any(idx)) {
        out[idx] <- stats::rnbinom(
          n = sum(idx),
          size = y_prev[idx],
          prob = 1 / (1 + alpha_t)
        )
      }
      out
    }
  }

  out <- array(NA_integer_, dim = c(n_test, h, n_sims))
  for (ss in seq_len(n_sims)) {
    y_prev <- history[, n_obs]
    for (step in seq_len(h)) {
      t_index <- start_time + step - 1L
      y_next <- draw_thinning(y_prev, alpha[t_index]) + draw_innovation(t_index)
      out[, step, ss] <- as.integer(y_next)
      y_prev <- y_next
    }
  }
  out
}

`%||%` <- function(x, y) if (is.null(x)) y else x
