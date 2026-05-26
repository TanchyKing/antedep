.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(.ppc_script_dir, "11_bolus_prediction.R"))

combine_bolus_all_pair_chunks <- function(out_dir, cfg) {
  chunk_dir <- file.path(out_dir, "chunks")
  chunk_paths <- sort(list.files(chunk_dir, pattern = "\\.rds$", full.names = TRUE))
  if (!length(chunk_paths)) stop("No chunk files found in ", chunk_dir)

  chunks <- lapply(chunk_paths, readRDS)
  scores <- do.call(rbind, lapply(chunks, `[[`, "scores"))
  diagnostics <- do.call(rbind, lapply(chunks, `[[`, "diagnostics"))
  extra_diagnostics <- do.call(rbind, lapply(chunks, `[[`, "extra_diagnostics"))
  folds_df <- do.call(rbind, lapply(chunks, `[[`, "folds"))
  summaries <- summarize_bolus_scores(scores)

  list(
    config = cfg,
    folds = folds_df,
    scores = scores,
    diagnostics = diagnostics,
    extra_diagnostics = extra_diagnostics,
    paired = summaries$paired,
    decision_summary = summaries$decision_summary,
    by_group = summaries$by_group
  )
}

run_bolus_prediction_all_pairs_chunked <- function(n_sims = 1000L,
                                                   fit_max_iter = 50L,
                                                   include_tscount = TRUE,
                                                   include_nb_glmm = TRUE,
                                                   chunk_size = 20L) {
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

  folds <- make_bolus_all_pair_folds(blocks)
  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", "bolus_leave_one_per_group_all_pairs")
  chunk_dir <- file.path(out_dir, "chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  started <- Sys.time()
  fold_chunks <- split(folds, ceiling(seq_along(folds) / chunk_size))

  cl <- NULL
  if (cfg$n_workers > 1L && length(folds) > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, chunk_size))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(
      cl,
      c("root", "cfg", "y", "blocks", "fit_bolus_inad_models", "run_bolus_fold"),
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
  }

  for (ii in seq_along(fold_chunks)) {
    folds_i <- fold_chunks[[ii]]
    fold_ids <- vapply(folds_i, function(x) x$fold[[1L]], integer(1L))
    chunk_path <- file.path(
      chunk_dir,
      sprintf("chunk_%03d_folds_%04d_%04d.rds", ii, min(fold_ids), max(fold_ids))
    )
    if (file.exists(chunk_path)) {
      message(sprintf("[%s] Skipping existing chunk %d/%d: folds %d-%d",
                      Sys.time(), ii, length(fold_chunks), min(fold_ids), max(fold_ids)))
      next
    }

    message(sprintf("[%s] Running chunk %d/%d: folds %d-%d",
                    Sys.time(), ii, length(fold_chunks), min(fold_ids), max(fold_ids)))
    flush.console()

    results <- if (!is.null(cl)) {
      parallel::parLapply(cl, folds_i, run_bolus_fold, y = y, blocks = blocks, cfg = cfg)
    } else {
      lapply(folds_i, run_bolus_fold, y = y, blocks = blocks, cfg = cfg)
    }

    chunk <- list(
      folds = do.call(rbind, folds_i),
      scores = do.call(rbind, lapply(results, `[[`, "scores")),
      diagnostics = do.call(rbind, lapply(results, `[[`, "diagnostics")),
      extra_diagnostics = do.call(rbind, lapply(results, `[[`, "extra_diagnostics")),
      completed_at = Sys.time()
    )
    saveRDS(chunk, chunk_path)
    message(sprintf("[%s] Saved %s", Sys.time(), basename(chunk_path)))
    flush.console()
  }

  combined <- combine_bolus_all_pair_chunks(out_dir, cfg)
  run_info <- data.frame(
    mode = "bolus_leave_one_per_group_all_pairs",
    fold_design = "all_pairs",
    n_folds = length(unique(combined$folds$fold)),
    n_subjects = nrow(y),
    n_group_1 = sum(blocks == sort(unique(blocks))[1L]),
    n_group_2 = sum(blocks == sort(unique(blocks))[2L]),
    n_sims = cfg$n_sims,
    fit_max_iter = cfg$fit_max_iter,
    include_tscount = include_tscount,
    include_nb_glmm = include_nb_glmm,
    chunk_size = chunk_size,
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    completed_at = as.character(Sys.time())
  )

  saveRDS(c(combined, list(run_info = run_info)),
          file.path(out_dir, "inad_prediction_bolus_leave_one_per_group_all_pairs.rds"))
  utils::write.csv(combined$folds, file.path(out_dir, "folds_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$scores, file.path(out_dir, "scores_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$diagnostics, file.path(out_dir, "diagnostics_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$extra_diagnostics, file.path(out_dir, "extra_diagnostics_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$paired, file.path(out_dir, "paired_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$decision_summary, file.path(out_dir, "decision_summary_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(combined$by_group, file.path(out_dir, "by_group_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_bolus_leave_one_per_group_all_pairs.csv"), row.names = FALSE)

  list(out_dir = out_dir, run_info = run_info, decision_summary = combined$decision_summary)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_bolus_prediction_all_pairs_chunked()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
