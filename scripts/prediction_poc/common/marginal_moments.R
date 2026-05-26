ppc_moment_frame <- function(model, group, time, mean, variance, source) {
  data.frame(
    model = as.character(model),
    group = as.integer(group),
    time = as.integer(time),
    mean = as.numeric(mean),
    variance = as.numeric(variance),
    source = as.character(source),
    stringsAsFactors = FALSE
  )
}

ppc_empirical_bolus_moments <- function(y, blocks) {
  out <- list()
  for (group in sort(unique(as.integer(blocks)))) {
    yy <- y[blocks == group, , drop = FALSE]
    out[[length(out) + 1L]] <- ppc_moment_frame(
      model = "Empirical bolus",
      group = group,
      time = seq_len(ncol(yy)),
      mean = colMeans(yy),
      variance = apply(yy, 2L, stats::var),
      source = "sample moments"
    )
  }
  do.call(rbind, out)
}

ppc_marginal_moments_inad_order1 <- function(fit, group_levels, label) {
  if (!identical(as.integer(fit$settings$order), 1L)) {
    stop("Only order-1 INAD marginal moments are implemented")
  }
  if (!identical(fit$settings$thinning, "nbinom") ||
      !identical(fit$settings$innovation, "nbinom")) {
    stop("Only NBT-NBI INAD marginal moments are implemented")
  }

  theta <- as.numeric(fit$theta)
  alpha <- as.numeric(fit$alpha)
  if (length(alpha) == 1L) alpha <- c(0, rep(alpha, length(theta) - 1L))
  nb_size <- as.numeric(fit$nb_inno_size)
  if (length(nb_size) == 1L) nb_size <- rep(nb_size, length(theta))
  tau <- as.numeric(fit$tau)
  if (!length(tau)) tau <- 0

  out <- list()
  for (group in group_levels) {
    tau_g <- if (length(tau) >= group) tau[[group]] else 0
    innovation_mean <- pmax(theta + tau_g, 1e-10)
    mu <- variance <- numeric(length(theta))
    mu[[1L]] <- innovation_mean[[1L]]
    variance[[1L]] <- innovation_mean[[1L]] + innovation_mean[[1L]]^2 / nb_size[[1L]]

    if (length(theta) >= 2L) {
      for (tt in 2:length(theta)) {
        a <- alpha[[tt]]
        mu[[tt]] <- a * mu[[tt - 1L]] + innovation_mean[[tt]]
        variance[[tt]] <- a * (1 + a) * mu[[tt - 1L]] +
          a^2 * variance[[tt - 1L]] +
          innovation_mean[[tt]] + innovation_mean[[tt]]^2 / nb_size[[tt]]
      }
    }
    out[[length(out) + 1L]] <- ppc_moment_frame(
      model = label,
      group = group,
      time = seq_along(theta),
      mean = mu,
      variance = variance,
      source = "analytic NBT-NBI recursion"
    )
  }
  do.call(rbind, out)
}

ppc_fit_marginal_nb_glm_full <- function(y, blocks) {
  dat <- ppc_long_marginal(y, blocks)
  MASS::glm.nb(y ~ group + factor(time), data = dat)
}

ppc_marginal_moments_nb_glm <- function(fit, group_levels, n_time) {
  out <- list()
  for (group in group_levels) {
    nd <- data.frame(group = factor(group, levels = group_levels), time = seq_len(n_time))
    mu <- as.numeric(stats::predict(fit, newdata = nd, type = "response"))
    out[[length(out) + 1L]] <- ppc_moment_frame(
      model = "NB GLM marginal",
      group = group,
      time = seq_len(n_time),
      mean = mu,
      variance = mu + mu^2 / fit$theta,
      source = "MASS::glm.nb closed form"
    )
  }
  do.call(rbind, out)
}

ppc_long_full <- function(y, blocks) {
  dat <- ppc_long_marginal(y, blocks)
  dat$subject <- factor(dat$subject)
  dat$time_fac <- factor(dat$time, levels = seq_len(ncol(y)))
  dat
}

ppc_fit_marginal_nb_glmm_full <- function(y, blocks) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    return(structure(
      list(message = "Package glmmTMB is not installed"),
      class = c("ppc_baseline_error", "error", "condition")
    ))
  }
  dat <- ppc_long_full(y, blocks)
  tryCatch(
    glmmTMB::glmmTMB(
      y ~ group + time_fac + (1 | subject),
      family = glmmTMB::nbinom2(link = "log"),
      data = dat
    ),
    error = function(e) e
  )
}

ppc_marginal_moments_nb_glmm <- function(fit, group_levels, n_time) {
  if (inherits(fit, "error")) return(NULL)
  re_var <- as.numeric(glmmTMB::VarCorr(fit)$cond$subject[1L, 1L])
  size <- as.numeric(stats::sigma(fit))
  out <- list()
  for (group in group_levels) {
    nd <- data.frame(
      subject = factor(rep("new_subject", n_time)),
      group = factor(group, levels = group_levels),
      time_fac = factor(seq_len(n_time), levels = seq_len(n_time))
    )
    mu_cond_zero <- as.numeric(stats::predict(
      fit,
      newdata = nd,
      type = "response",
      re.form = NA,
      allow.new.levels = TRUE
    ))
    lambda_mean <- mu_cond_zero * exp(re_var / 2)
    lambda_second <- mu_cond_zero^2 * exp(2 * re_var)
    variance <- lambda_mean + lambda_second / size +
      (lambda_second - lambda_mean^2)
    out[[length(out) + 1L]] <- ppc_moment_frame(
      model = "NB GLMM marginal",
      group = group,
      time = seq_len(n_time),
      mean = lambda_mean,
      variance = variance,
      source = "glmmTMB random-intercept marginalization"
    )
  }
  do.call(rbind, out)
}

ppc_marginal_moments_tscount_mc <- function(fit, y, blocks, cfg,
                                            n_sims = 2000L,
                                            seed = 20260521L) {
  if (inherits(fit, "error")) return(NULL)
  set.seed(seed)
  group_levels <- sort(unique(as.integer(blocks)))
  out <- list()

  for (group in group_levels) {
    sims <- matrix(NA_real_, nrow = n_sims, ncol = ncol(y))
    starts <- sample(y[blocks == group, 1L], n_sims, replace = TRUE)
    sims[, 1L] <- starts
    if (ncol(y) >= 2L) {
      for (tt in 2:ncol(y)) {
        mu <- ppc_predict_tscount_mu(
          fit = fit,
          blocks = rep(group, n_sims),
          time = tt,
          y_lag1 = sims[, tt - 1L],
          cfg = cfg
        )
        if (identical(fit$family, "poisson") || !is.finite(fit$size)) {
          sims[, tt] <- stats::rpois(n_sims, lambda = mu)
        } else {
          sims[, tt] <- stats::rnbinom(n_sims, size = fit$size, mu = mu)
        }
      }
    }
    out[[length(out) + 1L]] <- ppc_moment_frame(
      model = "tscount tsglm MC",
      group = group,
      time = seq_len(ncol(y)),
      mean = colMeans(sims),
      variance = apply(sims, 2L, stats::var),
      source = paste("MC panel recursion with", n_sims, "paths")
    )
  }
  do.call(rbind, out)
}

ppc_moment_l1_distance <- function(moments) {
  empirical <- moments[moments$model == "Empirical bolus",
                       c("group", "time", "mean", "variance")]
  names(empirical)[3:4] <- c("emp_mean", "emp_variance")
  fitted <- moments[moments$model != "Empirical bolus", ]
  merged <- merge(fitted, empirical, by = c("group", "time"), all.x = TRUE)
  merged$abs_mean_error <- abs(merged$mean - merged$emp_mean)
  merged$abs_variance_error <- abs(merged$variance - merged$emp_variance)
  stats::aggregate(
    cbind(abs_mean_error, abs_variance_error) ~ model + group,
    merged,
    sum,
    na.rm = TRUE
  )
}

ppc_moment_by_time_table <- function(moments, group, metric = c("mean", "variance")) {
  metric <- match.arg(metric)
  dat <- moments[moments$group == group, c("time", "model", metric)]
  models <- unique(dat$model)
  out <- data.frame(time = sort(unique(dat$time)))
  for (model in models) {
    cur <- dat[dat$model == model, c("time", metric)]
    names(cur)[2L] <- model
    out <- merge(out, cur, by = "time", all.x = TRUE, sort = TRUE)
  }
  out
}

ppc_empirical_table3_style <- function(moments) {
  empirical <- moments[moments$model == "Empirical bolus",
                       c("group", "time", "mean", "variance")]
  out <- data.frame(time = sort(unique(empirical$time)))
  for (group in sort(unique(empirical$group))) {
    cur <- empirical[empirical$group == group, c("time", "mean", "variance")]
    names(cur)[2:3] <- paste0("group_", group, c("_mean", "_variance"))
    out <- merge(out, cur, by = "time", all.x = TRUE, sort = TRUE)
  }
  out
}

ppc_moment_relative_error_table <- function(moments, group,
                                            metric = c("mean", "variance"),
                                            scale = 100) {
  metric <- match.arg(metric)
  dat <- moments[moments$group == group, c("time", "model", metric)]
  empirical <- dat[dat$model == "Empirical bolus", c("time", metric)]
  names(empirical)[2L] <- "empirical"
  out <- empirical[order(empirical$time), ]

  fitted_models <- setdiff(unique(dat$model), "Empirical bolus")
  for (model in fitted_models) {
    cur <- dat[dat$model == model, c("time", metric)]
    names(cur)[2L] <- "fitted"
    cur <- merge(empirical, cur, by = "time", all.x = TRUE, sort = TRUE)
    rel <- scale * (cur$fitted - cur$empirical) / cur$empirical
    out[[model]] <- rel[match(out$time, cur$time)]
  }
  out
}
