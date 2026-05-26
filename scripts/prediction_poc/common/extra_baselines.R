ppc_extra_baseline_packages <- function() {
  c("tscount", "glmmTMB")
}

ppc_extra_package_status <- function(pkgs = ppc_extra_baseline_packages()) {
  data.frame(
    package = pkgs,
    installed = vapply(pkgs, requireNamespace, logical(1L), quietly = TRUE),
    version = vapply(
      pkgs,
      function(pkg) {
        if (!requireNamespace(pkg, quietly = TRUE)) return(NA_character_)
        as.character(utils::packageVersion(pkg))
      },
      character(1L)
    ),
    row.names = NULL
  )
}

ppc_install_extra_baseline_packages <- function(pkgs = ppc_extra_baseline_packages(),
                                                repos = "https://cloud.r-project.org") {
  status <- ppc_extra_package_status(pkgs)
  missing <- status$package[!status$installed]
  if (length(missing)) {
    utils::install.packages(missing, repos = repos)
  }
  ppc_extra_package_status(pkgs)
}

ppc_extra_long_lag1 <- function(y, blocks, cfg, subjects = seq_len(nrow(y))) {
  dat <- ppc_long_lag1(y, blocks, subjects = subjects)
  dat$group <- factor(dat$group, levels = sort(unique(as.integer(blocks))))
  dat$time_fac <- factor(dat$time, levels = 2:cfg$n_time)
  dat$subject <- factor(dat$subject)
  dat[order(dat$subject, dat$time), , drop = FALSE]
}

ppc_extra_newdata <- function(blocks, time, y_lag1, subjects, cfg, group_levels = NULL) {
  if (is.null(group_levels)) group_levels <- sort(unique(as.integer(blocks)))
  data.frame(
    subject = factor(subjects),
    group = factor(blocks, levels = group_levels),
    time = as.integer(time),
    time_fac = factor(as.integer(time), levels = 2:cfg$n_time),
    y_lag1 = as.numeric(y_lag1)
  )
}

ppc_fit_nb_glmm_baseline <- function(y_train, blocks_train, cfg) {
  if (!requireNamespace("glmmTMB", quietly = TRUE)) {
    return(structure(
      list(message = "Package glmmTMB is not installed"),
      class = c("ppc_baseline_error", "error", "condition")
    ))
  }

  train_long <- ppc_extra_long_lag1(y_train, blocks_train, cfg)
  fit <- tryCatch(
    glmmTMB::glmmTMB(
      y ~ group + time_fac + log(y_lag1 + 1) + (1 | subject),
      family = glmmTMB::nbinom2(link = "log"),
      data = train_long
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(fit)

  list(
    model = fit,
    family = "nbinom",
    size = as.numeric(stats::sigma(fit)),
    group_levels = levels(train_long$group)
  )
}

ppc_fit_tscount_baseline <- function(y_train, blocks_train, cfg) {
  if (!requireNamespace("tscount", quietly = TRUE)) {
    return(structure(
      list(message = "Package tscount is not installed"),
      class = c("ppc_baseline_error", "error", "condition")
    ))
  }

  train_long <- ppc_extra_long_lag1(y_train, blocks_train, cfg)
  xreg_formula <- ~ group + time_fac
  xreg <- stats::model.matrix(xreg_formula, data = train_long)
  xreg <- xreg[, setdiff(colnames(xreg), "(Intercept)"), drop = FALSE]

  fit <- tryCatch(
    tscount::tsglm(
      ts = as.integer(train_long$y),
      model = list(past_obs = 1),
      xreg = xreg,
      link = "log",
      distr = "nbinom"
    ),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(fit)

  list(
    model = fit,
    xreg_formula = xreg_formula,
    xreg_names = colnames(xreg),
    group_levels = levels(train_long$group),
    family = if (!is.null(fit$distr)) fit$distr else "poisson",
    size = if (!is.null(fit$distrcoefs)) as.numeric(fit$distrcoefs[["size"]]) else NA_real_,
    note = paste(
      "Exploratory panel surrogate: tscount::tsglm is fit to the pooled",
      "subject-by-time training sequence. Subject-boundary lags are therefore",
      "not reset by tscount and this baseline should be interpreted cautiously."
    )
  )
}

ppc_fit_extra_baselines_rep <- function(y_train, blocks_train, cfg) {
  out <- list()
  if (isTRUE(cfg$extra_baselines$nb_glmm)) {
    out$nb_glmm <- ppc_fit_nb_glmm_baseline(y_train, blocks_train, cfg)
  }
  if (isTRUE(cfg$extra_baselines$tscount)) {
    out$tscount_tsglm <- ppc_fit_tscount_baseline(y_train, blocks_train, cfg)
  }
  out
}

ppc_tscount_xreg <- function(fit, newdata) {
  xreg <- stats::model.matrix(fit$xreg_formula, data = newdata)
  xreg <- xreg[, setdiff(colnames(xreg), "(Intercept)"), drop = FALSE]

  missing <- setdiff(fit$xreg_names, colnames(xreg))
  if (length(missing)) {
    xreg <- cbind(xreg, matrix(0, nrow = nrow(xreg), ncol = length(missing),
                               dimnames = list(NULL, missing)))
  }
  extra <- setdiff(colnames(xreg), fit$xreg_names)
  if (length(extra)) xreg <- xreg[, setdiff(colnames(xreg), extra), drop = FALSE]
  xreg[, fit$xreg_names, drop = FALSE]
}

ppc_predict_tscount_mu <- function(fit, blocks, time, y_lag1, cfg) {
  if (inherits(fit, "error")) return(rep(NA_real_, length(y_lag1)))

  nd <- ppc_extra_newdata(
    blocks = blocks,
    time = time,
    y_lag1 = y_lag1,
    subjects = seq_along(y_lag1),
    cfg = cfg,
    group_levels = fit$group_levels
  )
  xreg <- ppc_tscount_xreg(fit, nd)
  beta <- stats::coef(fit$model)
  if (length(beta) < 2L) return(rep(NA_real_, length(y_lag1)))

  xbeta <- if (ncol(xreg)) {
    as.numeric(xreg %*% beta[seq.int(3L, length.out = ncol(xreg))])
  } else {
    rep(0, length(y_lag1))
  }
  eta <- beta[[1L]] + beta[[2L]] * log(as.numeric(y_lag1) + 1) + xbeta
  as.numeric(exp(pmax(pmin(eta, 30), -30)))
}

ppc_predict_nb_glmm_mu <- function(fit, blocks, time, y_lag1, cfg) {
  if (inherits(fit, "error")) return(rep(NA_real_, length(y_lag1)))

  nd <- ppc_extra_newdata(
    blocks = blocks,
    time = time,
    y_lag1 = y_lag1,
    subjects = seq_along(y_lag1),
    cfg = cfg,
    group_levels = fit$group_levels
  )
  as.numeric(stats::predict(
    fit$model,
    newdata = nd,
    type = "response",
    re.form = NA,
    allow.new.levels = TRUE
  ))
}

ppc_predict_extra_baselines_rolling <- function(extra_fits, y_test, blocks_test, time, cfg) {
  out <- list()
  if (!is.null(extra_fits$nb_glmm)) {
    out$nb_glmm <- list(
      mu = ppc_predict_nb_glmm_mu(extra_fits$nb_glmm, blocks_test, time, y_test[, time - 1L], cfg),
      family = "nbinom",
      size = if (inherits(extra_fits$nb_glmm, "error")) NA_real_ else extra_fits$nb_glmm$size
    )
  }
  if (!is.null(extra_fits$tscount_tsglm)) {
    out$tscount_tsglm <- list(
      mu = ppc_predict_tscount_mu(extra_fits$tscount_tsglm, blocks_test, time, y_test[, time - 1L], cfg),
      family = if (inherits(extra_fits$tscount_tsglm, "error")) "nbinom" else extra_fits$tscount_tsglm$family,
      size = if (inherits(extra_fits$tscount_tsglm, "error")) NA_real_ else extra_fits$tscount_tsglm$size
    )
  }
  out
}

ppc_predict_extra_baselines_recursive <- function(extra_fits, y_test, blocks_test, cfg) {
  out <- list()
  for (nm in names(extra_fits)) {
    mu <- matrix(NA_real_, nrow = nrow(y_test), ncol = cfg$recursive_h)
    prev <- y_test[, cfg$t_split]
    for (hh in seq_len(cfg$recursive_h)) {
      tt <- cfg$t_split + hh
      if (identical(nm, "nb_glmm")) {
        mu[, hh] <- ppc_predict_nb_glmm_mu(extra_fits[[nm]], blocks_test, tt, prev, cfg)
      } else if (identical(nm, "tscount_tsglm")) {
        mu[, hh] <- ppc_predict_tscount_mu(extra_fits[[nm]], blocks_test, tt, prev, cfg)
      }
      prev <- mu[, hh]
    }
    out[[nm]] <- list(
      mu = mu,
      family = if (inherits(extra_fits[[nm]], "error")) "nbinom" else extra_fits[[nm]]$family,
      size = if (inherits(extra_fits[[nm]], "error")) NA_real_ else extra_fits[[nm]]$size
    )
  }
  out
}

ppc_score_exact_on_support <- function(observed, support, family, mean, size, epsilon) {
  if (!is.finite(mean)) return(NULL)
  if (identical(family, "poisson") || !is.finite(size)) {
    out <- ppc_score_pois_support(observed, support, mean, epsilon)
  } else {
    out <- ppc_score_nb_support(observed, support, mean, size, epsilon)
  }
  c(mean = mean, mc_se = NA_real_, out, cover80 = NA_real_, cover95 = NA_real_)
}

ppc_score_prediction_set_extra <- function(rep_id, task, time, horizon, y_obs,
                                           inad_samples_by_fit, glm_mu, glm_fits,
                                           extra_pred, cfg) {
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

    exact_models <- list(
      nb_glm = list(mu = glm_mu$nb[ii], family = "nbinom",
                    size = if (inherits(glm_fits$nb, "error")) NA_real_ else glm_fits$nb$theta),
      poisson_glm = list(mu = glm_mu$poisson[ii], family = "poisson", size = NA_real_)
    )
    for (nm in names(extra_pred)) {
      exact_models[[nm]] <- list(
        mu = extra_pred[[nm]]$mu[ii],
        family = extra_pred[[nm]]$family,
        size = extra_pred[[nm]]$size
      )
    }

    for (model_name in names(exact_models)) {
      pred <- exact_models[[model_name]]
      sc <- ppc_score_exact_on_support(
        observed = y_obs[ii],
        support = support,
        family = pred$family,
        mean = pred$mu,
        size = pred$size,
        epsilon = epsilon
      )
      if (is.null(sc)) next
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, model = model_name, task = task,
        time = time, horizon = horizon, subject_index = ii, observed = y_obs[ii],
        as.data.frame(as.list(sc)),
        squared_error = (pred$mu - y_obs[ii])^2,
        absolute_error = abs(pred$mu - y_obs[ii]),
        row.names = NULL, check.names = FALSE
      )
    }
  }

  do.call(rbind, rows)
}

ppc_predict_and_score_rep_extra <- function(rep_id, fits, glm_fits, extra_fits,
                                            y_test, blocks_test, cfg) {
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
    extra_pred <- ppc_predict_extra_baselines_rolling(extra_fits, y_test, blocks_test, tt, cfg)
    scores[[length(scores) + 1L]] <- ppc_score_prediction_set_extra(
      rep_id, "rolling_one_step", tt, 1L, y_test[, tt],
      samples, glm_mu, glm_fits, extra_pred, cfg
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
  extra_rec <- ppc_predict_extra_baselines_recursive(extra_fits, y_test, blocks_test, cfg)
  for (hh in seq_len(cfg$recursive_h)) {
    samples_h <- lapply(recursive_samples, function(x) x[, hh, ])
    glm_mu <- list(nb = glm_rec$nb[, hh], poisson = glm_rec$poisson[, hh])
    extra_pred <- lapply(extra_rec, function(x) {
      list(mu = x$mu[, hh], family = x$family, size = x$size)
    })
    tt <- cfg$t_split + hh
    scores[[length(scores) + 1L]] <- ppc_score_prediction_set_extra(
      rep_id, "recursive_multi_step", tt, hh, y_test[, tt],
      samples_h, glm_mu, glm_fits, extra_pred, cfg
    )
  }

  do.call(rbind, scores)
}

ppc_extra_fit_diagnostics_rep <- function(rep_id, extra_fits) {
  if (!length(extra_fits)) return(NULL)
  rows <- list()
  for (nm in names(extra_fits)) {
    fit <- extra_fits[[nm]]
    rows[[length(rows) + 1L]] <- data.frame(
      rep = rep_id,
      baseline = nm,
      ok = !inherits(fit, "error"),
      message = if (inherits(fit, "error")) conditionMessage(fit) else NA_character_,
      family = if (!inherits(fit, "error") && !is.null(fit$family)) fit$family else NA_character_,
      size = if (!inherits(fit, "error") && !is.null(fit$size)) fit$size else NA_real_,
      note = if (!inherits(fit, "error") && !is.null(fit$note)) fit$note else NA_character_,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

ppc_run_inad_rep_extra <- function(rep_id, cfg) {
  set.seed(ppc_rep_seed(cfg, rep_id, 0L))
  sim <- ppc_simulate_inad_dgp(cfg, rep_id)
  split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))

  y_train <- sim$y[split$train, , drop = FALSE]
  y_test <- sim$y[split$test, , drop = FALSE]
  blocks_train <- sim$blocks[split$train]
  blocks_test <- sim$blocks[split$test]

  fits <- ppc_fit_inad_rep(y_train, blocks_train, cfg)
  glm_fits <- ppc_fit_glm_rep(y_train, blocks_train)
  extra_fits <- ppc_fit_extra_baselines_rep(y_train, blocks_train, cfg)
  scores <- ppc_predict_and_score_rep_extra(rep_id, fits, glm_fits, extra_fits, y_test, blocks_test, cfg)
  recovery <- ppc_recovery_rep(rep_id, fits, sim$truth)
  diagnostics <- ppc_fit_diagnostics_rep(rep_id, fits)
  extra_diagnostics <- ppc_extra_fit_diagnostics_rep(rep_id, extra_fits)

  list(
    scores = scores,
    recovery = recovery,
    diagnostics = diagnostics,
    extra_diagnostics = extra_diagnostics
  )
}

ppc_summarize_scores_vs_references <- function(scores, reference_models = c("nb_glm", "nb_glmm", "tscount_tsglm")) {
  metric_cols <- c("log_score", "rps", "squared_error")
  agg <- stats::aggregate(scores[, metric_cols], scores[c("rep", "model", "task", "horizon")], mean, na.rm = TRUE)
  agg$rmse <- sqrt(agg$squared_error)
  inad <- agg[grepl("^inad_", agg$model), ]

  out <- list()
  for (ref in reference_models) {
    ref_rows <- agg[agg$model == ref, c("rep", "task", "horizon", "log_score", "rps", "rmse")]
    if (!nrow(ref_rows)) next
    names(ref_rows)[4:6] <- paste0(names(ref_rows)[4:6], "_reference")
    merged <- merge(inad, ref_rows, by = c("rep", "task", "horizon"), all.x = FALSE)
    merged$reference <- ref
    merged$diff_log_score_vs_reference <- merged$log_score - merged$log_score_reference
    merged$diff_rps_vs_reference <- merged$rps - merged$rps_reference
    merged$diff_rmse_vs_reference <- merged$rmse - merged$rmse_reference
    out[[ref]] <- merged
  }
  do.call(rbind, out)
}

ppc_summarize_extra_decision <- function(paired) {
  stats::aggregate(
    cbind(diff_log_score_vs_reference, diff_rps_vs_reference, diff_rmse_vs_reference) ~
      model + reference + task + horizon,
    paired,
    function(z) c(
      mean = mean(z, na.rm = TRUE),
      se = stats::sd(z, na.rm = TRUE) / sqrt(sum(is.finite(z))),
      q025 = unname(stats::quantile(z, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(z, 0.975, na.rm = TRUE))
    )
  )
}
