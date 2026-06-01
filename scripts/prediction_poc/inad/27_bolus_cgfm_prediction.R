.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c("score.R", "cgfm_prediction.R"), root)

make_bolus_all_pair_folds_local <- function(blocks) {
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

read_existing_bolus_folds <- function(root, blocks,
                                      fold_source = c("rotation", "all_pairs")) {
  fold_source <- match.arg(fold_source)
  candidates <- if (identical(fold_source, "rotation")) {
    c(
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_marginal_nb_glm"),
                "folds_bolus_leave_one_per_group.csv"),
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group"),
                "folds_bolus_leave_one_per_group.csv")
    )
  } else {
    c(
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_all_pairs_marginal_nb_glm"),
                "folds_bolus_leave_one_per_group.csv"),
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_all_pairs"),
                "folds_bolus_leave_one_per_group_all_pairs.csv")
    )
  }
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) return(make_bolus_all_pair_folds_local(blocks))
  split(utils::read.csv(path), utils::read.csv(path)$fold)
}

run_cgfm_bolus_fold <- function(fold, y, blocks, include_time_varying = TRUE,
                                maxit = 3000L, epsilon = 1 / 1001) {
  heldout <- fold$heldout_subject
  train <- setdiff(seq_len(nrow(y)), heldout)
  fit <- ppc_cgfm_fit(
    y = y[train, , drop = FALSE],
    blocks = blocks[train],
    include_time_varying = include_time_varying,
    maxit = maxit
  )
  scores <- ppc_cgfm_score_rolling(
    fit = fit,
    y_test = y[heldout, , drop = FALSE],
    blocks_test = blocks[heldout],
    fold_id = fold$fold[[1L]],
    epsilon = epsilon
  )
  scores$fold <- fold$fold[[1L]]
  scores$heldout_subject <- heldout[scores$subject_index]
  scores$heldout_group <- blocks[heldout][scores$subject_index]

  diagnostics <- ppc_cgfm_diagnostics(fit, fold$fold[[1L]])
  diagnostics$heldout_subjects <- paste(heldout, collapse = ",")

  list(scores = scores, diagnostics = diagnostics)
}

existing_constrained_scores <- function(root, folds_keep = NULL,
                                        fold_source = c("rotation", "all_pairs")) {
  fold_source <- match.arg(fold_source)
  candidates <- if (identical(fold_source, "rotation")) {
    c(
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_marginal_nb_glm"),
                "scores_bolus_leave_one_per_group.csv"),
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group"),
                "scores_bolus_leave_one_per_group.csv")
    )
  } else {
    c(
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_all_pairs_marginal_nb_glm"),
                "scores_bolus_leave_one_per_group.csv"),
      file.path(ppc_output_dir(root, "bolus_leave_one_per_group_all_pairs"),
                "scores_bolus_leave_one_per_group_all_pairs.csv")
    )
  }
  path <- candidates[file.exists(candidates)][1L]
  if (is.na(path)) stop("Could not find existing bolus all-pairs INAD score file")
  scores <- utils::read.csv(path)
  scores <- scores[scores$model == "inad_constrained_alpha" &
                     scores$task == "rolling_one_step" &
                     scores$horizon == 1L, ]
  if (!is.null(folds_keep)) scores <- scores[scores$fold %in% folds_keep, ]
  scores
}

summarize_mean_se <- function(z) {
  z <- z[is.finite(z)]
  c(
    mean = mean(z),
    se = stats::sd(z) / sqrt(length(z)),
    t = mean(z) / (stats::sd(z) / sqrt(length(z))),
    q025 = unname(stats::quantile(z, 0.025)),
    q975 = unname(stats::quantile(z, 0.975))
  )
}

summarize_cgfm_vs_constrained <- function(cgfm_scores, constrained_scores) {
  keys <- c("fold", "heldout_subject", "heldout_group", "task", "horizon", "time")
  keep <- c(keys, "model", "log_score", "rps", "squared_error")
  all_scores <- rbind(
    constrained_scores[, keep],
    cgfm_scores[, keep]
  )

  refs <- setdiff(unique(all_scores$model), "inad_constrained_alpha")
  rows <- list()
  lhs <- all_scores[all_scores$model == "inad_constrained_alpha", ]
  names(lhs)[match(c("log_score", "rps", "squared_error"), names(lhs))] <-
    c("log_score_model", "rps_model", "squared_error_model")
  lhs$model <- NULL

  for (ref in refs) {
    rhs <- all_scores[all_scores$model == ref, ]
    names(rhs)[match(c("log_score", "rps", "squared_error"), names(rhs))] <-
      c("log_score_reference", "rps_reference", "squared_error_reference")
    rhs$model <- NULL
    merged <- merge(lhs, rhs, by = keys)
    merged$model <- "inad_constrained_alpha"
    merged$reference <- ref
    merged$diff_log_score <- merged$log_score_model - merged$log_score_reference
    merged$diff_rps <- merged$rps_model - merged$rps_reference
    merged$diff_squared_error <- merged$squared_error_model - merged$squared_error_reference
    rows[[ref]] <- merged
  }
  obs_deltas <- do.call(rbind, rows)

  patient_scores <- stats::aggregate(
    cbind(log_score, rps, squared_error) ~
      heldout_subject + heldout_group + model + task + horizon,
    all_scores,
    mean,
    na.rm = TRUE
  )
  patient_scores$rmse <- sqrt(patient_scores$squared_error)
  rmse_rows <- list()
  cfit <- patient_scores[patient_scores$model == "inad_constrained_alpha",
                         c("heldout_subject", "heldout_group", "task", "horizon", "rmse")]
  names(cfit)[names(cfit) == "rmse"] <- "rmse_model"
  for (ref in refs) {
    rhs <- patient_scores[patient_scores$model == ref,
                          c("heldout_subject", "heldout_group", "task", "horizon", "rmse")]
    names(rhs)[names(rhs) == "rmse"] <- "rmse_reference"
    merged <- merge(cfit, rhs, by = c("heldout_subject", "heldout_group", "task", "horizon"))
    merged$model <- "inad_constrained_alpha"
    merged$reference <- ref
    merged$diff_rmse <- merged$rmse_model - merged$rmse_reference
    rmse_rows[[ref]] <- merged
  }
  rmse_deltas <- do.call(rbind, rmse_rows)

  patient_log_rps <- stats::aggregate(
    cbind(diff_log_score, diff_rps) ~ model + reference + task + horizon + heldout_subject + heldout_group,
    obs_deltas,
    mean,
    na.rm = TRUE
  )
  summary_log_rps <- stats::aggregate(
    cbind(diff_log_score, diff_rps) ~ model + reference + task + horizon,
    patient_log_rps,
    summarize_mean_se
  )
  summary_rmse <- stats::aggregate(
    diff_rmse ~ model + reference + task + horizon,
    rmse_deltas,
    summarize_mean_se
  )

  list(
    observation_deltas = obs_deltas,
    patient_log_rps = patient_log_rps,
    patient_rmse = rmse_deltas,
    summary_log_rps = summary_log_rps,
    summary_rmse = summary_rmse
  )
}

run_bolus_cgfm_prediction <- function(max_folds = 5L,
                                      fold_source = c("rotation", "all_pairs"),
                                      include_time_varying = TRUE,
                                      maxit = 3000L,
                                      out_label = NULL) {
  fold_source <- match.arg(fold_source)
  if (is.null(out_label)) {
    out_label <- if (identical(fold_source, "rotation")) {
      "bolus_cgfm_leave_one_per_group_rotation"
    } else {
      "bolus_cgfm_leave_one_per_group_all_pairs"
    }
  }
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(bolus_inad$blocks)
  folds <- read_existing_bolus_folds(root, blocks, fold_source = fold_source)
  if (!is.null(max_folds) && is.finite(max_folds) && max_folds > 0L) {
    folds <- folds[seq_len(min(length(folds), as.integer(max_folds)))]
  }

  out_dir <- ppc_output_dir(root, out_label)
  started <- Sys.time()
  message(sprintf("[%s] Running %d CGFM fold(s); time-varying=%s",
                  started, length(folds), include_time_varying))
  results <- lapply(folds, run_cgfm_bolus_fold,
                    y = y, blocks = blocks,
                    include_time_varying = include_time_varying,
                    maxit = maxit)
  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  folds_df <- do.call(rbind, folds)
  constrained <- existing_constrained_scores(
    root,
    folds_keep = unique(folds_df$fold),
    fold_source = fold_source
  )
  summaries <- summarize_cgfm_vs_constrained(scores, constrained)

  run_info <- data.frame(
    mode = out_label,
    fold_source = fold_source,
    n_folds = length(folds),
    include_time_varying = include_time_varying,
    maxit = maxit,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    completed_at = as.character(Sys.time())
  )

  result <- list(
    folds = folds_df,
    scores = scores,
    constrained_scores = constrained,
    diagnostics = diagnostics,
    summaries = summaries,
    run_info = run_info
  )
  saveRDS(result, file.path(out_dir, "cgfm_prediction_bolus_leave_one_per_group_all_pairs.rds"))
  utils::write.csv(folds_df, file.path(out_dir, "folds_cgfm_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(scores, file.path(out_dir, "scores_cgfm_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(diagnostics, file.path(out_dir, "diagnostics_cgfm_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(summaries$observation_deltas, file.path(out_dir, "observation_deltas_constrained_vs_cgfm.csv"), row.names = FALSE)
  utils::write.csv(summaries$patient_log_rps, file.path(out_dir, "patient_log_rps_constrained_vs_cgfm.csv"), row.names = FALSE)
  utils::write.csv(summaries$patient_rmse, file.path(out_dir, "patient_rmse_constrained_vs_cgfm.csv"), row.names = FALSE)
  utils::write.csv(summaries$summary_log_rps, file.path(out_dir, "summary_log_rps_constrained_vs_cgfm.csv"), row.names = FALSE)
  utils::write.csv(summaries$summary_rmse, file.path(out_dir, "summary_rmse_constrained_vs_cgfm.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_cgfm_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)

  list(out_dir = out_dir, run_info = run_info,
       summary_log_rps = summaries$summary_log_rps,
       summary_rmse = summaries$summary_rmse)
}

if (identical(environment(), globalenv()) && !interactive()) {
  max_folds_env <- Sys.getenv("PPC_CGFM_MAX_FOLDS", unset = "5")
  max_folds <- if (tolower(max_folds_env) %in% c("all", "0", "inf")) Inf else as.integer(max_folds_env)
  include_tv <- !tolower(Sys.getenv("PPC_CGFM_INCLUDE_TV", unset = "true")) %in% c("0", "false", "no")
  maxit <- as.integer(Sys.getenv("PPC_CGFM_MAXIT", unset = "3000"))
  fold_source <- match.arg(Sys.getenv("PPC_CGFM_FOLD_SOURCE", unset = "rotation"),
                           c("rotation", "all_pairs"))
  result <- run_bolus_cgfm_prediction(
    max_folds = max_folds,
    fold_source = fold_source,
    include_time_varying = include_tv,
    maxit = maxit
  )
  print(result$run_info)
  print(result$summary_log_rps)
  print(result$summary_rmse)
  cat("Output:", result$out_dir, "\n")
}
