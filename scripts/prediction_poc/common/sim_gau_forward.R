ppc_gau_blocks <- function(fit, blocks, n) {
  tau <- if (is.null(fit$tau)) 0 else as.numeric(fit$tau)
  has_tau <- length(tau) > 1L && any(abs(tau[-1L]) > 1e-12)
  if (is.null(blocks)) {
    if (has_tau) stop("blocks are required when fit has nonzero tau")
    return(rep(1L, n))
  }
  if (length(blocks) != n) stop("blocks must have length nrow(history)")
  as.integer(blocks)
}

ppc_gau_mean_at <- function(fit, t_index, blocks) {
  tau <- if (is.null(fit$tau)) 0 else as.numeric(fit$tau)
  if (length(tau) == 1L) tau <- rep(tau, max(blocks))
  as.numeric(fit$mu[t_index]) + tau[blocks]
}

ppc_gau_phi_at <- function(fit, t_index) {
  phi <- fit$phi
  if (is.matrix(phi)) return(as.numeric(phi[1L, t_index]))
  as.numeric(phi[t_index])
}

ppc_gau_forward_moments <- function(fit, history, blocks = NULL,
                                    start_time = ncol(history) + 1L,
                                    h = 1L) {
  history <- as.matrix(history)
  storage.mode(history) <- "double"
  n_test <- nrow(history)
  n_obs <- ncol(history)
  h <- as.integer(h)
  start_time <- as.integer(start_time)
  if (h < 1L) stop("h must be positive")
  if (start_time != n_obs + 1L) stop("start_time must equal ncol(history) + 1")

  order <- as.integer(fit$settings$order)
  if (order != 1L) stop("Gaussian POC currently supports order = 1 only")
  n_time <- as.integer(fit$settings$n_time)
  if (start_time + h - 1L > n_time) stop("Requested horizon exceeds fit time grid")

  blocks <- ppc_gau_blocks(fit, blocks, n_test)
  mean_out <- matrix(NA_real_, nrow = n_test, ncol = h)
  sd_out <- matrix(NA_real_, nrow = n_test, ncol = h)

  prev_mean <- history[, n_obs]
  prev_var <- rep(0, n_test)
  for (step in seq_len(h)) {
    tt <- start_time + step - 1L
    phi_t <- ppc_gau_phi_at(fit, tt)
    m_t <- ppc_gau_mean_at(fit, tt, blocks)
    m_prev <- ppc_gau_mean_at(fit, tt - 1L, blocks)
    mean_t <- m_t + phi_t * (prev_mean - m_prev)
    var_t <- phi_t^2 * prev_var + as.numeric(fit$sigma[tt])^2
    mean_out[, step] <- mean_t
    sd_out[, step] <- sqrt(pmax(var_t, 0))
    prev_mean <- mean_t
    prev_var <- var_t
  }

  colnames(mean_out) <- colnames(sd_out) <- paste0("h", seq_len(h))
  list(mean = mean_out, sd = sd_out)
}

ppc_simulate_gau_forward <- function(fit, history, blocks = NULL,
                                     start_time = ncol(history) + 1L,
                                     h = 1L, n_sims = 1000L, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  history <- as.matrix(history)
  storage.mode(history) <- "double"
  n_test <- nrow(history)
  n_obs <- ncol(history)
  h <- as.integer(h)
  n_sims <- as.integer(n_sims)
  if (h < 1L || n_sims < 1L) stop("h and n_sims must be positive")

  order <- as.integer(fit$settings$order)
  if (order != 1L) stop("Gaussian POC currently supports order = 1 only")
  blocks <- ppc_gau_blocks(fit, blocks, n_test)

  out <- array(NA_real_, dim = c(n_test, h, n_sims))
  for (ss in seq_len(n_sims)) {
    y_prev <- history[, n_obs]
    for (step in seq_len(h)) {
      tt <- start_time + step - 1L
      phi_t <- ppc_gau_phi_at(fit, tt)
      m_t <- ppc_gau_mean_at(fit, tt, blocks)
      m_prev <- ppc_gau_mean_at(fit, tt - 1L, blocks)
      y_next <- m_t + phi_t * (y_prev - m_prev) +
        stats::rnorm(n_test, mean = 0, sd = as.numeric(fit$sigma[tt]))
      out[, step, ss] <- y_next
      y_prev <- y_next
    }
  }
  out
}

ppc_score_gau_normal <- function(observed, mean, sd) {
  sd <- pmax(as.numeric(sd), 1e-8)
  z <- (observed - mean) / sd
  crps <- sd * (z * (2 * stats::pnorm(z) - 1) +
                  2 * stats::dnorm(z) - 1 / sqrt(pi))
  data.frame(
    mean = as.numeric(mean),
    sd = sd,
    log_score = -stats::dnorm(observed, mean = mean, sd = sd, log = TRUE),
    crps = crps,
    squared_error = (mean - observed)^2
  )
}
