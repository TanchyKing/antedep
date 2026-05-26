.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(.ppc_script_dir, "11_bolus_prediction.R"))

score_bolus_recursive_smoothing_fold <- function(fold, y, blocks, cfg, epsilons) {
  heldout <- fold$heldout_subject
  train <- setdiff(seq_len(nrow(y)), heldout)
  y_train <- y[train, , drop = FALSE]
  y_test <- y[heldout, , drop = FALSE]
  blocks_train <- blocks[train]
  blocks_test <- blocks[heldout]

  fits <- fit_bolus_inad_models(y_train, blocks_train, cfg)
  glm_fits <- ppc_fit_glm_rep(y_train, blocks_train)

  fit_names <- names(fits)
  recursive_samples <- list()
  for (ff in fit_names) {
    recursive_samples[[ff]] <- simulate_inad_forward(
      fits[[ff]],
      history = y_test[, 1:cfg$t_split, drop = FALSE],
      blocks = blocks_test,
      start_time = cfg$t_split + 1L,
      h = cfg$recursive_h,
      n_sims = cfg$n_sims,
      seed = ppc_rep_seed(cfg, fold$fold[[1L]], 9000L + 100L * match(ff, fit_names))
    )
  }
  glm_rec <- ppc_predict_glm_means_recursive(glm_fits, y_test, blocks_test, cfg)

  rows <- list()
  for (hh in seq_len(cfg$recursive_h)) {
    tt <- cfg$t_split + hh
    y_obs <- y_test[, tt]
    for (ii in seq_along(y_obs)) {
      sample_vectors <- lapply(recursive_samples, function(x) as.integer(x[ii, hh, ]))
      support_max <- max(c(y_obs[ii], unlist(sample_vectors, use.names = FALSE)), na.rm = TRUE)
      support <- 0:support_max

      for (epsilon_name in names(epsilons)) {
        epsilon <- epsilons[[epsilon_name]]
        for (fit_name in names(sample_vectors)) {
          sc <- ppc_score_samples(sample_vectors[[fit_name]], y_obs[ii], epsilon = epsilon)
          pred_mean <- unname(sc["mean"])
          rows[[length(rows) + 1L]] <- data.frame(
            fold = fold$fold[[1L]],
            heldout_subject = heldout[[ii]],
            heldout_group = blocks_test[[ii]],
            model = paste0("inad_", fit_name),
            reference = NA_character_,
            task = "recursive_multi_step",
            time = tt,
            horizon = hh,
            epsilon_label = epsilon_name,
            epsilon = epsilon,
            observed = y_obs[ii],
            mean = pred_mean,
            log_score = unname(sc["log_score"]),
            rps = unname(sc["rps"]),
            squared_error = (pred_mean - y_obs[ii])^2,
            stringsAsFactors = FALSE
          )
        }

        if (!inherits(glm_fits$nb, "error") && is.finite(glm_rec$nb[ii, hh])) {
          sc_nb <- ppc_score_nb_support(y_obs[ii], support, glm_rec$nb[ii, hh], glm_fits$nb$theta, epsilon)
          rows[[length(rows) + 1L]] <- data.frame(
            fold = fold$fold[[1L]],
            heldout_subject = heldout[[ii]],
            heldout_group = blocks_test[[ii]],
            model = "nb_glm",
            reference = NA_character_,
            task = "recursive_multi_step",
            time = tt,
            horizon = hh,
            epsilon_label = epsilon_name,
            epsilon = epsilon,
            observed = y_obs[ii],
            mean = glm_rec$nb[ii, hh],
            log_score = unname(sc_nb["log_score"]),
            rps = unname(sc_nb["rps"]),
            squared_error = (glm_rec$nb[ii, hh] - y_obs[ii])^2,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  do.call(rbind, rows)
}

summarize_smoothing_sensitivity <- function(scores) {
  agg <- stats::aggregate(
    cbind(log_score, rps, squared_error) ~
      fold + epsilon_label + model + task + horizon,
    scores,
    mean,
    na.rm = TRUE
  )
  agg$rmse <- sqrt(agg$squared_error)
  refs <- agg[agg$model == "nb_glm", c("fold", "epsilon_label", "task", "horizon", "log_score", "rps", "rmse")]
  names(refs)[5:7] <- paste0(names(refs)[5:7], "_nb")
  inad <- agg[grepl("^inad_", agg$model), ]
  merged <- merge(inad, refs, by = c("fold", "epsilon_label", "task", "horizon"))
  merged$diff_log_score_vs_nb <- merged$log_score - merged$log_score_nb
  merged$diff_rps_vs_nb <- merged$rps - merged$rps_nb
  merged$diff_rmse_vs_nb <- merged$rmse - merged$rmse_nb
  stats::aggregate(
    cbind(diff_log_score_vs_nb, diff_rps_vs_nb, diff_rmse_vs_nb) ~
      epsilon_label + model + task + horizon,
    merged,
    function(z) c(mean = mean(z), se = stats::sd(z) / sqrt(length(z)))
  )
}

run_bolus_smoothing_sensitivity <- function(fold_limit = 3L, n_sims = 100L,
                                            fit_max_iter = 8L, smoke = TRUE,
                                            fold_design = c("all_pairs", "rotation")) {
  fold_design <- match.arg(fold_design)
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(if (!is.null(bolus_inad$blocks)) bolus_inad$blocks else bolus_inad$bolus)
  cfg <- ppc_inad_config("full")
  cfg$mode <- "bolus_smoothing"
  cfg$n_sims <- as.integer(n_sims)
  cfg$n_time <- ncol(y)
  cfg$n_subjects <- nrow(y)
  cfg$t_split <- 8L
  cfg$recursive_h <- 4L
  cfg$fit_max_iter <- as.integer(fit_max_iter)
  cfg$thinning <- "nbinom"
  cfg$innovation <- "nbinom"
  cfg$order <- 1L

  folds <- if (identical(fold_design, "rotation")) {
    make_bolus_leave_one_per_group_folds(blocks, seed = cfg$seed)
  } else {
    make_bolus_all_pair_folds(blocks)
  }
  folds <- folds[seq_len(min(length(folds), fold_limit))]
  epsilons <- c(
    small = 1 / (10 * cfg$n_sims + 1),
    main = 1 / (cfg$n_sims + 1),
    large = 1 / (max(1, cfg$n_sims / 10) + 1)
  )

  if (cfg$n_workers > 1L && length(folds) > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, length(folds)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c("root"),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      setwd(file.path(root, "scripts", "prediction_poc", "inad"))
      source("15_bolus_smoothing_sensitivity.R")
      NULL
    })
    scores <- do.call(rbind, parallel::parLapply(
      cl, folds, score_bolus_recursive_smoothing_fold,
      y = y, blocks = blocks, cfg = cfg, epsilons = epsilons
    ))
  } else {
    scores <- do.call(rbind, lapply(
      folds, score_bolus_recursive_smoothing_fold,
      y = y, blocks = blocks, cfg = cfg, epsilons = epsilons
    ))
  }
  summary <- summarize_smoothing_sensitivity(scores)

  out_dir <- file.path(
    ppc_poc_dir(root), "output", "inad",
    if (smoke) "bolus_smoothing_sensitivity_smoke" else "bolus_smoothing_sensitivity"
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_info <- data.frame(
    mode = basename(out_dir),
    fold_design = fold_design,
    n_folds = length(folds),
    n_sims = cfg$n_sims,
    fit_max_iter = cfg$fit_max_iter,
    epsilon_small = epsilons[["small"]],
    epsilon_main = epsilons[["main"]],
    epsilon_large = epsilons[["large"]]
  )
  saveRDS(list(config = cfg, scores = scores, summary = summary, run_info = run_info),
          file.path(out_dir, "bolus_smoothing_sensitivity.rds"))
  utils::write.csv(scores, file.path(out_dir, "scores_bolus_smoothing_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(out_dir, "summary_bolus_smoothing_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_bolus_smoothing_sensitivity.csv"), row.names = FALSE)
  list(out_dir = out_dir, run_info = run_info, summary = summary, scores = scores)
}

.ppc_is_main_script <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  length(file_arg) &&
    identical(
      normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = FALSE),
      normalizePath(file.path(.ppc_script_dir, "15_bolus_smoothing_sensitivity.R"), winslash = "/", mustWork = FALSE)
    )
})

if (.ppc_is_main_script && identical(environment(), globalenv()) && !interactive()) {
  result <- run_bolus_smoothing_sensitivity()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
