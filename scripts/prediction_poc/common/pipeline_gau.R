ppc_simulate_gau_dgp <- function(cfg, rep_id) {
  y <- antedep::simulate_gau(
    n_subjects = cfg$n_subjects,
    n_time = cfg$n_time,
    order = cfg$order,
    mu = cfg$mu,
    phi = cfg$phi,
    sigma = cfg$sigma,
    blocks = cfg$blocks,
    tau = cfg$tau,
    seed = ppc_rep_seed(cfg, rep_id, 1L)
  )
  list(
    y = unclass(y),
    blocks = cfg$blocks,
    truth = list(mu = cfg$mu, phi = cfg$phi, sigma = cfg$sigma, tau = cfg$tau)
  )
}

ppc_fit_gau_lm_baseline <- function(y_train, blocks_train) {
  dat <- ppc_long_lag1(y_train, blocks_train)
  fit <- stats::lm(y ~ group + factor(time) + y_lag1, data = dat)
  sigma <- sqrt(mean(stats::residuals(fit)^2, na.rm = TRUE))
  list(fit = fit, sigma = sigma, group_levels = levels(dat$group))
}

ppc_gau_lm_newdata <- function(blocks, time, y_lag1, group_levels) {
  data.frame(
    group = factor(blocks, levels = group_levels),
    time = time,
    y_lag1 = as.numeric(y_lag1)
  )
}

ppc_predict_gau_lm_rolling <- function(lm_fit, y_test, blocks_test, time) {
  nd <- ppc_gau_lm_newdata(blocks_test, time, y_test[, time - 1L], lm_fit$group_levels)
  as.numeric(stats::predict(lm_fit$fit, newdata = nd))
}

ppc_predict_gau_lm_recursive <- function(lm_fit, y_test, blocks_test, cfg) {
  out <- matrix(NA_real_, nrow = nrow(y_test), ncol = cfg$recursive_h)
  prev <- y_test[, cfg$t_split]
  for (hh in seq_len(cfg$recursive_h)) {
    tt <- cfg$t_split + hh
    nd <- ppc_gau_lm_newdata(blocks_test, tt, prev, lm_fit$group_levels)
    out[, hh] <- as.numeric(stats::predict(lm_fit$fit, newdata = nd))
    prev <- out[, hh]
  }
  out
}

ppc_score_gau_prediction_set <- function(rep_id, model, task, time, horizon,
                                         subject_index, observed, mean, sd) {
  sc <- ppc_score_gau_normal(observed, mean, sd)
  data.frame(
    rep = rep_id,
    model = model,
    task = task,
    time = time,
    horizon = horizon,
    subject_index = subject_index,
    observed = observed,
    sc,
    row.names = NULL
  )
}

ppc_predict_and_score_gau_rep <- function(rep_id, fit, lm_fit, y_test, blocks_test, cfg) {
  rows <- list()
  n_test <- nrow(y_test)

  for (tt in 2:cfg$n_time) {
    mom <- ppc_gau_forward_moments(
      fit,
      history = y_test[, 1:(tt - 1L), drop = FALSE],
      blocks = blocks_test,
      h = 1L
    )
    rows[[length(rows) + 1L]] <- ppc_score_gau_prediction_set(
      rep_id, "gau_fit", "rolling_one_step", tt, 1L, seq_len(n_test),
      y_test[, tt], mom$mean[, 1L], mom$sd[, 1L]
    )

    lm_mean <- ppc_predict_gau_lm_rolling(lm_fit, y_test, blocks_test, tt)
    rows[[length(rows) + 1L]] <- ppc_score_gau_prediction_set(
      rep_id, "lag_lm", "rolling_one_step", tt, 1L, seq_len(n_test),
      y_test[, tt], lm_mean, lm_fit$sigma
    )
  }

  mom_rec <- ppc_gau_forward_moments(
    fit,
    history = y_test[, 1:cfg$t_split, drop = FALSE],
    blocks = blocks_test,
    h = cfg$recursive_h
  )
  lm_rec <- ppc_predict_gau_lm_recursive(lm_fit, y_test, blocks_test, cfg)
  for (hh in seq_len(cfg$recursive_h)) {
    tt <- cfg$t_split + hh
    rows[[length(rows) + 1L]] <- ppc_score_gau_prediction_set(
      rep_id, "gau_fit", "recursive_multi_step", tt, hh, seq_len(n_test),
      y_test[, tt], mom_rec$mean[, hh], mom_rec$sd[, hh]
    )
    rows[[length(rows) + 1L]] <- ppc_score_gau_prediction_set(
      rep_id, "lag_lm", "recursive_multi_step", tt, hh, seq_len(n_test),
      y_test[, tt], lm_rec[, hh], lm_fit$sigma
    )
  }

  do.call(rbind, rows)
}

ppc_validate_gau_moments <- function(fit, y_test, blocks_test, cfg) {
  history <- y_test[, 1:cfg$t_split, drop = FALSE]
  analytic <- ppc_gau_forward_moments(fit, history, blocks_test, h = cfg$recursive_h)
  samples <- ppc_simulate_gau_forward(
    fit, history, blocks_test, h = cfg$recursive_h,
    n_sims = 2000L, seed = ppc_rep_seed(cfg, 1L, 7000L)
  )
  sample_mean <- apply(samples, c(1L, 2L), mean)
  sample_sd <- apply(samples, c(1L, 2L), stats::sd)
  data.frame(
    max_abs_mean_diff = max(abs(sample_mean - analytic$mean), na.rm = TRUE),
    max_abs_sd_diff = max(abs(sample_sd - analytic$sd), na.rm = TRUE)
  )
}

ppc_run_gau_rep <- function(rep_id, cfg) {
  sim <- ppc_simulate_gau_dgp(cfg, rep_id)
  split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))
  y_train <- sim$y[split$train, , drop = FALSE]
  y_test <- sim$y[split$test, , drop = FALSE]
  blocks_train <- sim$blocks[split$train]
  blocks_test <- sim$blocks[split$test]

  fit <- antedep::fit_gau(y_train, order = cfg$order, blocks = blocks_train)
  lm_fit <- ppc_fit_gau_lm_baseline(y_train, blocks_train)
  scores <- ppc_predict_and_score_gau_rep(rep_id, fit, lm_fit, y_test, blocks_test, cfg)

  recovery <- data.frame(
    rep = rep_id,
    mu_rmse = sqrt(mean((as.numeric(fit$mu) - cfg$mu)^2)),
    phi_mae = mean(abs(as.numeric(fit$phi) - cfg$phi)),
    sigma_rmse = sqrt(mean((as.numeric(fit$sigma) - cfg$sigma)^2)),
    tau2_bias = if (length(fit$tau) >= 2L) as.numeric(fit$tau[2L]) - cfg$tau[2L] else NA_real_
  )

  list(scores = scores, recovery = recovery, validation = ppc_validate_gau_moments(fit, y_test, blocks_test, cfg))
}

ppc_summarize_gau_scores <- function(scores) {
  agg <- stats::aggregate(
    cbind(log_score, crps, squared_error) ~ rep + model + task + horizon,
    scores,
    mean,
    na.rm = TRUE
  )
  agg$rmse <- sqrt(agg$squared_error)
  ref <- agg[agg$model == "lag_lm", c("rep", "task", "horizon", "log_score", "crps", "rmse")]
  names(ref)[4:6] <- paste0(names(ref)[4:6], "_lm")
  gau <- agg[agg$model == "gau_fit", ]
  paired <- merge(gau, ref, by = c("rep", "task", "horizon"))
  paired$diff_log_score_vs_lm <- paired$log_score - paired$log_score_lm
  paired$diff_crps_vs_lm <- paired$crps - paired$crps_lm
  paired$diff_rmse_vs_lm <- paired$rmse - paired$rmse_lm

  decision <- stats::aggregate(
    cbind(diff_log_score_vs_lm, diff_crps_vs_lm, diff_rmse_vs_lm) ~ task + horizon,
    paired,
    function(z) c(mean = mean(z), se = stats::sd(z) / sqrt(length(z)))
  )
  list(aggregate = agg, paired = paired, decision = decision)
}

ppc_run_gau_smoke <- function() {
  root <- ppc_find_package_root()
  cfg <- ppc_gau_config("smoke")
  started <- Sys.time()
  reps <- seq_len(cfg$n_reps)
  if (cfg$n_workers > 1L && length(reps) > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, length(reps)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("root", "cfg"), envir = environment())
    parallel::clusterEvalQ(cl, {
      setwd(root)
      source(file.path(root, "scripts", "prediction_poc", "common", "config.R"))
      ppc_load_antedep(root)
      ppc_source_common(c("holdout.R", "sim_gau_forward.R", "pipeline_gau.R"), root)
      NULL
    })
    results <- parallel::parLapply(cl, reps, ppc_run_gau_rep, cfg = cfg)
  } else {
    results <- lapply(reps, ppc_run_gau_rep, cfg = cfg)
  }

  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  recovery <- do.call(rbind, lapply(results, `[[`, "recovery"))
  validation <- do.call(rbind, lapply(results, `[[`, "validation"))
  summaries <- ppc_summarize_gau_scores(scores)
  out_dir <- file.path(ppc_poc_dir(root), "output", "gau", "smoke")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_info <- data.frame(
    mode = "gau_smoke",
    n_reps = cfg$n_reps,
    n_subjects = cfg$n_subjects,
    n_time = cfg$n_time,
    n_workers = cfg$n_workers,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
  result <- list(config = cfg, scores = scores, recovery = recovery,
                 validation = validation, summaries = summaries, run_info = run_info)
  saveRDS(result, file.path(out_dir, "gau_prediction_smoke.rds"))
  utils::write.csv(scores, file.path(out_dir, "scores_gau_smoke.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(out_dir, "recovery_gau_smoke.csv"), row.names = FALSE)
  utils::write.csv(validation, file.path(out_dir, "validation_gau_smoke.csv"), row.names = FALSE)
  utils::write.csv(summaries$paired, file.path(out_dir, "paired_gau_smoke.csv"), row.names = FALSE)
  utils::write.csv(summaries$decision, file.path(out_dir, "decision_gau_smoke.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_gau_smoke.csv"), row.names = FALSE)
  result
}
