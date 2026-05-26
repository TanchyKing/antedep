ppc_score_samples <- function(samples, observed, epsilon = 1 / (length(samples) + 1)) {
  samples <- as.integer(samples)
  observed <- as.integer(observed)
  max_y <- max(c(samples, observed), na.rm = TRUE)
  support <- 0:max_y
  counts <- tabulate(samples + 1L, nbins = max_y + 1L)
  prob <- counts / length(samples)
  prob <- prob + epsilon
  prob <- prob / sum(prob)

  cdf <- cumsum(prob)
  obs_idx <- observed + 1L
  u <- stats::runif(1L)
  pit <- (if (obs_idx > 1L) cdf[obs_idx - 1L] else 0) + u * prob[obs_idx]
  indicator <- as.numeric(support >= observed)

  c(
    mean = mean(samples),
    mc_se = stats::sd(samples) / sqrt(length(samples)),
    log_score = -log(prob[obs_idx]),
    rps = sum((cdf - indicator)^2),
    pit = pit,
    cover80 = as.numeric(observed >= stats::quantile(samples, 0.10, type = 1) &&
                           observed <= stats::quantile(samples, 0.90, type = 1)),
    cover95 = as.numeric(observed >= stats::quantile(samples, 0.025, type = 1) &&
                           observed <= stats::quantile(samples, 0.975, type = 1))
  )
}

ppc_score_pmf <- function(observed, support, prob, epsilon) {
  observed <- as.integer(observed)
  prob <- as.numeric(prob)
  prob <- prob + epsilon
  prob <- prob / sum(prob)
  cdf <- cumsum(prob)
  obs_idx <- match(observed, support)
  if (is.na(obs_idx)) stop("observed value is not in support")
  u <- stats::runif(1L)
  pit <- (if (obs_idx > 1L) cdf[obs_idx - 1L] else 0) + u * prob[obs_idx]
  indicator <- as.numeric(support >= observed)
  c(log_score = -log(prob[obs_idx]), rps = sum((cdf - indicator)^2), pit = pit)
}

ppc_score_nb <- function(observed, samples, mu, size, epsilon = 1 / (length(samples) + 1)) {
  max_y <- max(c(samples, observed), na.rm = TRUE)
  support <- 0:max_y
  out <- ppc_score_nb_support(observed, support, mu, size, epsilon)
  c(mean = mu, mc_se = NA_real_, out, cover80 = NA_real_, cover95 = NA_real_)
}

ppc_score_pois <- function(observed, samples, lambda, epsilon = 1 / (length(samples) + 1)) {
  max_y <- max(c(samples, observed), na.rm = TRUE)
  support <- 0:max_y
  out <- ppc_score_pois_support(observed, support, lambda, epsilon)
  c(mean = lambda, mc_se = NA_real_, out, cover80 = NA_real_, cover95 = NA_real_)
}

ppc_score_nb_support <- function(observed, support, mu, size, epsilon) {
  prob <- stats::dnbinom(support, size = size, mu = mu)
  ppc_score_pmf(observed, support, prob, epsilon)
}

ppc_score_pois_support <- function(observed, support, lambda, epsilon) {
  prob <- stats::dpois(support, lambda = lambda)
  ppc_score_pmf(observed, support, prob, epsilon)
}
