.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
  "score.R", "pipeline_inad.R", "extra_baselines.R"
), root)

make_bolus_leave_one_per_group_folds <- function(blocks, seed = 20260518L) {
  set.seed(seed)
  g1 <- sample(which(blocks == sort(unique(blocks))[1L]))
  g2 <- sample(which(blocks == sort(unique(blocks))[2L]))
  n_folds <- max(length(g1), length(g2))

  out <- vector("list", n_folds)
  for (ii in seq_len(n_folds)) {
    out[[ii]] <- data.frame(
      fold = ii,
      heldout_subject = c(g1[((ii - 1L) %% length(g1)) + 1L], g2[((ii - 1L) %% length(g2)) + 1L]),
      heldout_group = c(blocks[g1[((ii - 1L) %% length(g1)) + 1L]], blocks[g2[((ii - 1L) %% length(g2)) + 1L]])
    )
  }
  out
}

make_bolus_all_pair_folds <- function(blocks) {
  groups <- sort(unique(blocks))
  if (length(groups) != 2L) stop("Expected exactly two bolus groups")
  g1 <- which(blocks == groups[1L])
  g2 <- which(blocks == groups[2L])
  grid <- expand.grid(group1_subject = g1, group2_subject = g2)

  out <- vector("list", nrow(grid))
  for (ii in seq_len(nrow(grid))) {
    heldout <- c(grid$group1_subject[[ii]], grid$group2_subject[[ii]])
    out[[ii]] <- data.frame(
      fold = ii,
      heldout_subject = heldout,
      heldout_group = blocks[heldout]
    )
  }
  out
}

fit_bolus_inad_models <- function(y_train, blocks_train, cfg) {
  unconstrained <- antedep::fit_inad(
    y = y_train,
    order = cfg$order,
    thinning = cfg$thinning,
    innovation = cfg$innovation,
    blocks = blocks_train,
    max_iter = cfg$fit_max_iter
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

run_bolus_fold <- function(fold, y, blocks, cfg) {
  heldout <- fold$heldout_subject
  train <- setdiff(seq_len(nrow(y)), heldout)

  y_train <- y[train, , drop = FALSE]
  y_test <- y[heldout, , drop = FALSE]
  blocks_train <- blocks[train]
  blocks_test <- blocks[heldout]

  fits <- fit_bolus_inad_models(y_train, blocks_train, cfg)
  glm_fits <- ppc_fit_glm_rep(y_train, blocks_train)
  extra_fits <- ppc_fit_extra_baselines_rep(y_train, blocks_train, cfg)
  scores <- ppc_predict_and_score_rep_extra(
    rep_id = fold$fold[[1L]],
    fits = fits,
    glm_fits = glm_fits,
    extra_fits = extra_fits,
    y_test = y_test,
    blocks_test = blocks_test,
    cfg = cfg
  )
  scores$fold <- fold$fold[[1L]]
  scores$heldout_subject <- heldout[scores$subject_index]
  scores$heldout_group <- blocks_test[scores$subject_index]

  diagnostics <- ppc_fit_diagnostics_rep(fold$fold[[1L]], fits)
  diagnostics$fold <- fold$fold[[1L]]
  diagnostics$heldout_subjects <- paste(heldout, collapse = ",")

  extra_diagnostics <- ppc_extra_fit_diagnostics_rep(fold$fold[[1L]], extra_fits)
  if (!is.null(extra_diagnostics)) {
    extra_diagnostics$fold <- fold$fold[[1L]]
    extra_diagnostics$heldout_subjects <- paste(heldout, collapse = ",")
  }

  list(scores = scores, diagnostics = diagnostics, extra_diagnostics = extra_diagnostics)
}

summarize_bolus_scores <- function(scores) {
  paired <- ppc_summarize_scores_vs_references(
    scores,
    reference_models = c("nb_glm", "poisson_glm", "nb_glmm", "tscount_tsglm")
  )
  decision <- ppc_summarize_extra_decision(paired)

  by_group <- stats::aggregate(
    cbind(log_score, rps, squared_error) ~ heldout_group + model + task + horizon,
    scores,
    mean,
    na.rm = TRUE
  )
  by_group$rmse <- sqrt(by_group$squared_error)

  list(paired = paired, decision_summary = decision, by_group = by_group)
}

run_bolus_prediction <- function(n_sims = 1000L, fit_max_iter = 50L,
                                 include_tscount = TRUE, include_nb_glmm = TRUE) {
  run_bolus_prediction_design(
    fold_design = "rotation",
    n_sims = n_sims,
    fit_max_iter = fit_max_iter,
    include_tscount = include_tscount,
    include_nb_glmm = include_nb_glmm
  )
}

run_bolus_prediction_design <- function(fold_design = c("rotation", "all_pairs"),
                                        n_sims = 1000L, fit_max_iter = 50L,
                                        include_tscount = TRUE, include_nb_glmm = TRUE,
                                        out_label = NULL) {
  fold_design <- match.arg(fold_design)
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- if (!is.null(bolus_inad$blocks)) bolus_inad$blocks else bolus_inad$bolus
  blocks <- as.integer(blocks)

  cfg <- ppc_inad_config("full")
  cfg$mode <- "bolus"
  cfg$n_reps <- NA_integer_
  cfg$n_sims <- as.integer(n_sims)
  cfg$n_time <- ncol(y)
  cfg$n_subjects <- nrow(y)
  cfg$t_split <- 8L
  cfg$recursive_h <- 4L
  cfg$fit_max_iter <- as.integer(fit_max_iter)
  cfg$thinning <- "nbinom"
  cfg$innovation <- "nbinom"
  cfg$order <- 1L
  cfg$extra_baselines <- list(nb_glmm = include_nb_glmm, tscount = include_tscount)

  folds <- if (identical(fold_design, "all_pairs")) {
    make_bolus_all_pair_folds(blocks)
  } else {
    make_bolus_leave_one_per_group_folds(blocks, seed = cfg$seed)
  }
  out_name <- if (identical(fold_design, "all_pairs")) {
    "bolus_leave_one_per_group_all_pairs"
  } else {
    "bolus_leave_one_per_group"
  }
  if (!is.null(out_label) && nzchar(out_label)) {
    out_name <- paste(out_name, out_label, sep = "_")
  }
  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", out_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  started <- Sys.time()
  if (cfg$n_workers > 1L && length(folds) > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, length(folds)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c(
        "root", "cfg", "y", "blocks",
        "fit_bolus_inad_models", "run_bolus_fold"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      setwd(root)
      source(file.path(root, "scripts", "prediction_poc", "common", "config.R"))
      ppc_load_antedep(root)
      ppc_source_common(c(
        "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
        "score.R", "pipeline_inad.R", "extra_baselines.R"
      ), root)
      NULL
    })
    results <- parallel::parLapply(cl, folds, run_bolus_fold, y = y, blocks = blocks, cfg = cfg)
  } else {
    results <- lapply(folds, run_bolus_fold, y = y, blocks = blocks, cfg = cfg)
  }

  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  extra_diagnostics <- do.call(rbind, lapply(results, `[[`, "extra_diagnostics"))
  summaries <- summarize_bolus_scores(scores)

  folds_df <- do.call(rbind, folds)
  run_info <- data.frame(
    mode = out_name,
    fold_design = fold_design,
    n_folds = length(folds),
    n_subjects = nrow(y),
    n_group_1 = sum(blocks == sort(unique(blocks))[1L]),
    n_group_2 = sum(blocks == sort(unique(blocks))[2L]),
    n_sims = cfg$n_sims,
    fit_max_iter = cfg$fit_max_iter,
    include_tscount = include_tscount,
    include_nb_glmm = include_nb_glmm,
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )

  result <- list(
    config = cfg,
    folds = folds_df,
    scores = scores,
    diagnostics = diagnostics,
    extra_diagnostics = extra_diagnostics,
    paired = summaries$paired,
    decision_summary = summaries$decision_summary,
    by_group = summaries$by_group,
    run_info = run_info
  )

  saveRDS(result, file.path(out_dir, "inad_prediction_bolus_leave_one_per_group.rds"))
  utils::write.csv(folds_df, file.path(out_dir, "folds_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(scores, file.path(out_dir, "scores_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(diagnostics, file.path(out_dir, "diagnostics_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(extra_diagnostics, file.path(out_dir, "extra_diagnostics_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(summaries$paired, file.path(out_dir, "paired_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(summaries$decision_summary, file.path(out_dir, "decision_summary_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(summaries$by_group, file.path(out_dir, "by_group_bolus_leave_one_per_group.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_bolus_leave_one_per_group.csv"), row.names = FALSE)

  list(out_dir = out_dir, run_info = run_info, decision_summary = summaries$decision_summary)
}

.ppc_is_main_script <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  length(file_arg) &&
    identical(
      normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = FALSE),
      normalizePath(file.path(.ppc_script_dir, "11_bolus_prediction.R"), winslash = "/", mustWork = FALSE)
    )
})

if (.ppc_is_main_script && identical(environment(), globalenv()) && !interactive()) {
  result <- run_bolus_prediction()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
