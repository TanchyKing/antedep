ppc_simulate_inad_dgp <- function(cfg, rep_id) {
  alpha <- c(0, rep(cfg$alpha, cfg$n_time - 1L))
  y <- antedep::simulate_inad(
    n_subjects = cfg$n_subjects,
    n_time = cfg$n_time,
    order = cfg$order,
    thinning = cfg$thinning,
    innovation = cfg$innovation,
    alpha = alpha,
    theta = cfg$theta,
    nb_inno_size = cfg$nb_inno_size,
    blocks = cfg$blocks,
    tau = cfg$tau,
    seed = ppc_rep_seed(cfg, rep_id, 1L)
  )
  list(
    y = unclass(y),
    blocks = cfg$blocks,
    truth = list(alpha = alpha, theta = cfg$theta, nb_inno_size = cfg$nb_inno_size, tau = cfg$tau)
  )
}

ppc_fit_inad_rep <- function(y_train, blocks_train, cfg) {
  init_alpha <- c(0, rep(cfg$alpha, cfg$n_time - 1L))
  unconstrained <- antedep::fit_inad(
    y = y_train,
    order = cfg$order,
    thinning = cfg$thinning,
    innovation = cfg$innovation,
    blocks = blocks_train,
    max_iter = cfg$fit_max_iter,
    init_alpha = init_alpha,
    init_theta = cfg$theta,
    init_tau = cfg$tau,
    init_nb_inno_size = cfg$nb_inno_size
  )
  constrained <- fit_inad_alpha_constant_poc(
    y = y_train,
    blocks = blocks_train,
    fit_unconstrained = unconstrained,
    thinning = cfg$thinning,
    innovation = cfg$innovation,
    max_iter = cfg$fit_max_iter
  )
  list(unconstrained = unconstrained, constrained_alpha = constrained)
}

ppc_fit_inad_rep_default_init <- function(y_train, blocks_train, cfg) {
  antedep::fit_inad(
    y = y_train,
    order = cfg$order,
    thinning = cfg$thinning,
    innovation = cfg$innovation,
    blocks = blocks_train,
    max_iter = cfg$fit_max_iter
  )
}

ppc_fit_glm_rep <- function(y_train, blocks_train) {
  train_long <- ppc_long_lag1(y_train, blocks_train)
  nb_long <- ppc_long_marginal(y_train, blocks_train)
  nb <- tryCatch(
    MASS::glm.nb(y ~ group + factor(time), data = nb_long),
    error = function(e) e
  )
  pois <- stats::glm(
    y ~ group + factor(time) + log(y_lag1 + 1),
    family = stats::poisson(link = "log"),
    data = train_long
  )
  list(
    nb = nb,
    poisson = pois,
    group_levels = levels(train_long$group),
    nb_glm_form = "marginal"
  )
}

ppc_long_marginal <- function(y, blocks, subjects = seq_len(nrow(y))) {
  y <- as.matrix(y)
  grid <- expand.grid(
    subject_index = seq_along(subjects),
    time = seq_len(ncol(y)),
    KEEP.OUT.ATTRS = FALSE
  )
  subject_rows <- as.integer(subjects[grid$subject_index])
  data.frame(
    subject = subject_rows,
    group = factor(blocks[subject_rows], levels = sort(unique(as.integer(blocks)))),
    time = as.integer(grid$time),
    y = as.numeric(y[cbind(subject_rows, grid$time)])
  )
}

ppc_newdata <- function(blocks, time, y_lag1, group_levels) {
  data.frame(
    group = factor(blocks, levels = group_levels),
    time = time,
    y_lag1 = as.numeric(y_lag1)
  )
}

ppc_predict_nb_marginal_mean <- function(glm_fits, blocks, time) {
  if (inherits(glm_fits$nb, "error")) return(rep(NA_real_, length(blocks)))
  nd <- ppc_newdata(blocks, time, rep(0, length(blocks)), glm_fits$group_levels)
  as.numeric(stats::predict(glm_fits$nb, newdata = nd, type = "response"))
}

ppc_predict_glm_means_rolling <- function(glm_fits, y_test, blocks_test, time) {
  nd <- ppc_newdata(blocks_test, time, y_test[, time - 1L], glm_fits$group_levels)
  list(
    nb = ppc_predict_nb_marginal_mean(glm_fits, blocks_test, time),
    poisson = as.numeric(stats::predict(glm_fits$poisson, newdata = nd, type = "response"))
  )
}

ppc_predict_glm_means_recursive <- function(glm_fits, y_test, blocks_test, cfg) {
  out_nb <- matrix(NA_real_, nrow = nrow(y_test), ncol = cfg$recursive_h)
  out_pois <- matrix(NA_real_, nrow = nrow(y_test), ncol = cfg$recursive_h)
  prev_pois <- y_test[, cfg$t_split]
  for (hh in seq_len(cfg$recursive_h)) {
    tt <- cfg$t_split + hh
    nd_pois <- ppc_newdata(blocks_test, tt, prev_pois, glm_fits$group_levels)
    out_nb[, hh] <- ppc_predict_nb_marginal_mean(glm_fits, blocks_test, tt)
    out_pois[, hh] <- as.numeric(stats::predict(glm_fits$poisson, newdata = nd_pois, type = "response"))
    prev_pois <- out_pois[, hh]
  }
  list(nb = out_nb, poisson = out_pois)
}

ppc_score_prediction_set <- function(rep_id, task, time, horizon, y_obs, inad_samples_by_fit,
                                     glm_mu, glm_fits, cfg) {
  rows <- list()
  epsilon <- 1 / (cfg$n_sims + 1)
  n <- length(y_obs)

  for (ii in seq_len(n)) {
    sample_vectors <- lapply(inad_samples_by_fit, function(x) as.integer(x[ii, ]))
    support_max <- max(c(y_obs[ii], unlist(sample_vectors, use.names = FALSE)), na.rm = TRUE)
    support <- 0:support_max

    for (fit_name in names(sample_vectors)) {
      sc <- ppc_score_samples(sample_vectors[[fit_name]], y_obs[ii], epsilon = epsilon)
      pred_mean <- unname(sc["mean"])
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, model = paste0("inad_", fit_name), task = task,
        time = time, horizon = horizon, subject_index = ii, observed = y_obs[ii],
        as.data.frame(as.list(sc)),
        squared_error = (pred_mean - y_obs[ii])^2,
        absolute_error = abs(pred_mean - y_obs[ii]),
        row.names = NULL, check.names = FALSE
      )
    }

    if (!inherits(glm_fits$nb, "error") && is.finite(glm_mu$nb[ii])) {
      sc_nb <- ppc_score_nb_support(y_obs[ii], support, glm_mu$nb[ii], glm_fits$nb$theta, epsilon)
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, model = "nb_glm", task = task, time = time,
        horizon = horizon, subject_index = ii, observed = y_obs[ii],
        mean = glm_mu$nb[ii], mc_se = NA_real_, as.data.frame(as.list(sc_nb)),
        cover80 = NA_real_, cover95 = NA_real_,
        squared_error = (glm_mu$nb[ii] - y_obs[ii])^2,
        absolute_error = abs(glm_mu$nb[ii] - y_obs[ii]),
        row.names = NULL, check.names = FALSE
      )
    }

    if (is.finite(glm_mu$poisson[ii])) {
      sc_p <- ppc_score_pois_support(y_obs[ii], support, glm_mu$poisson[ii], epsilon)
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, model = "poisson_glm", task = task, time = time,
        horizon = horizon, subject_index = ii, observed = y_obs[ii],
        mean = glm_mu$poisson[ii], mc_se = NA_real_, as.data.frame(as.list(sc_p)),
        cover80 = NA_real_, cover95 = NA_real_,
        squared_error = (glm_mu$poisson[ii] - y_obs[ii])^2,
        absolute_error = abs(glm_mu$poisson[ii] - y_obs[ii]),
        row.names = NULL, check.names = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

ppc_predict_and_score_rep <- function(rep_id, fits, glm_fits, y_test, blocks_test, cfg) {
  scores <- list()
  fit_names <- names(fits)

  for (tt in 2:cfg$n_time) {
    samples <- list()
    for (ff in fit_names) {
      arr <- simulate_inad_forward(
        fits[[ff]],
        history = y_test[, 1:(tt - 1L), drop = FALSE],
        blocks = blocks_test,
        start_time = tt,
        h = 1L,
        n_sims = cfg$n_sims,
        seed = ppc_rep_seed(cfg, rep_id, 1000L + 100L * match(ff, fit_names) + tt)
      )
      samples[[ff]] <- arr[, 1L, , drop = FALSE][, 1L, ]
    }
    glm_mu <- ppc_predict_glm_means_rolling(glm_fits, y_test, blocks_test, tt)
    scores[[length(scores) + 1L]] <- ppc_score_prediction_set(
      rep_id, "rolling_one_step", tt, 1L, y_test[, tt], samples, glm_mu, glm_fits, cfg
    )
  }

  recursive_samples <- list()
  for (ff in fit_names) {
    recursive_samples[[ff]] <- simulate_inad_forward(
      fits[[ff]],
      history = y_test[, 1:cfg$t_split, drop = FALSE],
      blocks = blocks_test,
      start_time = cfg$t_split + 1L,
      h = cfg$recursive_h,
      n_sims = cfg$n_sims,
      seed = ppc_rep_seed(cfg, rep_id, 2000L + 100L * match(ff, fit_names))
    )
  }
  glm_rec <- ppc_predict_glm_means_recursive(glm_fits, y_test, blocks_test, cfg)
  for (hh in seq_len(cfg$recursive_h)) {
    samples_h <- lapply(recursive_samples, function(x) x[, hh, ])
    glm_mu <- list(nb = glm_rec$nb[, hh], poisson = glm_rec$poisson[, hh])
    tt <- cfg$t_split + hh
    scores[[length(scores) + 1L]] <- ppc_score_prediction_set(
      rep_id, "recursive_multi_step", tt, hh, y_test[, tt], samples_h, glm_mu, glm_fits, cfg
    )
  }

  do.call(rbind, scores)
}

ppc_recovery_rep <- function(rep_id, fits, truth) {
  rows <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    alpha <- as.numeric(fit$alpha)
    alpha_idx <- 2:length(truth$alpha)
    rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_id,
      fit = nm,
      alpha_mae = mean(abs(alpha[alpha_idx] - truth$alpha[alpha_idx]), na.rm = TRUE),
      alpha_mean = mean(alpha[alpha_idx], na.rm = TRUE),
      alpha_min = min(alpha[alpha_idx], na.rm = TRUE),
      alpha_max = max(alpha[alpha_idx], na.rm = TRUE),
      theta_bias_mean = mean(as.numeric(fit$theta) - truth$theta, na.rm = TRUE),
      theta_rmse = sqrt(mean((as.numeric(fit$theta) - truth$theta)^2, na.rm = TRUE)),
      nb_size_bias_mean = mean(as.numeric(fit$nb_inno_size) - truth$nb_inno_size, na.rm = TRUE),
      nb_size_rmse = sqrt(mean((as.numeric(fit$nb_inno_size) - truth$nb_inno_size)^2, na.rm = TRUE)),
      tau2_bias = if (length(fit$tau) >= 2L) as.numeric(fit$tau[2]) - truth$tau[2] else NA_real_
    )
  }
  do.call(rbind, rows)
}

ppc_fit_diagnostics_rep <- function(rep_id, fits) {
  rows <- list()
  for (nm in names(fits)) {
    fit <- fits[[nm]]
    rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_id,
      fit = nm,
      convergence = if (!is.null(fit$convergence)) as.character(fit$convergence) else NA_character_,
      log_l = if (!is.null(fit$log_l)) as.numeric(fit$log_l) else NA_real_,
      alpha_mean = if (!is.null(fit$alpha)) mean(as.numeric(fit$alpha)[-1L], na.rm = TRUE) else NA_real_,
      nb_size_mean = if (!is.null(fit$nb_inno_size)) mean(as.numeric(fit$nb_inno_size), na.rm = TRUE) else NA_real_,
      nb_size_min = if (!is.null(fit$nb_inno_size)) min(as.numeric(fit$nb_inno_size), na.rm = TRUE) else NA_real_,
      nb_size_max = if (!is.null(fit$nb_inno_size)) max(as.numeric(fit$nb_inno_size), na.rm = TRUE) else NA_real_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

ppc_run_inad_rep <- function(rep_id, cfg) {
  set.seed(ppc_rep_seed(cfg, rep_id, 0L))
  sim <- ppc_simulate_inad_dgp(cfg, rep_id)
  split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))

  y_train <- sim$y[split$train, , drop = FALSE]
  y_test <- sim$y[split$test, , drop = FALSE]
  blocks_train <- sim$blocks[split$train]
  blocks_test <- sim$blocks[split$test]

  fits <- ppc_fit_inad_rep(y_train, blocks_train, cfg)
  glm_fits <- ppc_fit_glm_rep(y_train, blocks_train)
  scores <- ppc_predict_and_score_rep(rep_id, fits, glm_fits, y_test, blocks_test, cfg)
  recovery <- ppc_recovery_rep(rep_id, fits, sim$truth)
  diagnostics <- ppc_fit_diagnostics_rep(rep_id, fits)

  list(scores = scores, recovery = recovery, diagnostics = diagnostics)
}

ppc_summarize_scores <- function(scores) {
  metric_cols <- c("log_score", "rps", "squared_error")
  agg <- stats::aggregate(scores[, metric_cols], scores[c("rep", "model", "task", "horizon")], mean, na.rm = TRUE)
  agg$rmse <- sqrt(agg$squared_error)
  nb <- agg[agg$model == "nb_glm", c("rep", "task", "horizon", "log_score", "rps", "rmse")]
  names(nb)[4:6] <- paste0(names(nb)[4:6], "_nb")
  merged <- merge(agg[grepl("^inad_", agg$model), ], nb, by = c("rep", "task", "horizon"), all.x = TRUE)
  merged$diff_log_score_vs_nb <- merged$log_score - merged$log_score_nb
  merged$diff_rps_vs_nb <- merged$rps - merged$rps_nb
  merged$diff_rmse_vs_nb <- merged$rmse - merged$rmse_nb
  merged
}
