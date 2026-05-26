.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
  "score.R", "pipeline_inad.R"
), root)

summarize_final_decision <- function(paired) {
  stats::aggregate(
    cbind(diff_log_score_vs_nb, diff_rps_vs_nb, diff_rmse_vs_nb) ~ model + task + horizon,
    paired,
    function(z) c(
      mean = mean(z, na.rm = TRUE),
      se = stats::sd(z, na.rm = TRUE) / sqrt(sum(is.finite(z))),
      q025 = unname(stats::quantile(z, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(z, 0.975, na.rm = TRUE))
    )
  )
}

combine_r1000_outputs <- function(out_dir) {
  r200_path <- file.path(
    ppc_poc_dir(root), "output", "inad", "intermediate_r200",
    "inad_prediction_intermediate_r200.rds"
  )
  if (!file.exists(r200_path)) stop("Missing intermediate R=200 output: ", r200_path)

  r200 <- readRDS(r200_path)
  chunk_paths <- sort(list.files(file.path(out_dir, "chunks"), pattern = "\\.rds$", full.names = TRUE))
  if (length(chunk_paths) == 0L) stop("No remaining-run chunk files found")
  chunks <- lapply(chunk_paths, readRDS)

  scores <- do.call(rbind, c(list(r200$scores), lapply(chunks, `[[`, "scores")))
  recovery <- do.call(rbind, c(list(r200$recovery), lapply(chunks, `[[`, "recovery")))
  diagnostics <- do.call(rbind, c(list(r200$diagnostics), lapply(chunks, `[[`, "diagnostics")))
  paired <- ppc_summarize_scores(scores)
  decision_summary <- summarize_final_decision(paired)

  list(
    config = r200$config,
    scores = scores,
    recovery = recovery,
    diagnostics = diagnostics,
    paired = paired,
    decision_summary = decision_summary
  )
}

run_inad_remaining_r800 <- function(chunk_size = 20L) {
  cfg <- ppc_inad_config("full")
  cfg$n_reps <- 1000L
  cfg$n_sims <- 1000L
  cfg$n_per_group <- 50L
  cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
  cfg$n_subjects <- length(cfg$blocks)
  cfg$fit_max_iter <- 50L

  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", "remaining_r800")
  chunk_dir <- file.path(out_dir, "chunks")
  dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

  started <- Sys.time()
  reps <- 201:1000
  chunks <- split(reps, ceiling(seq_along(reps) / chunk_size))

  cl <- NULL
  if (cfg$n_workers > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, chunk_size))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("root", "cfg"), envir = environment())
    parallel::clusterEvalQ(cl, {
      setwd(root)
      source(file.path(root, "scripts", "prediction_poc", "common", "config.R"))
      ppc_load_antedep(root)
      ppc_source_common(c(
        "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
        "score.R", "pipeline_inad.R"
      ), root)
      NULL
    })
  }

  for (ii in seq_along(chunks)) {
    reps_i <- chunks[[ii]]
    chunk_path <- file.path(chunk_dir, sprintf("chunk_%03d_reps_%04d_%04d.rds", ii, min(reps_i), max(reps_i)))
    if (file.exists(chunk_path)) {
      message(sprintf("[%s] Skipping existing chunk %d/%d: reps %d-%d",
                      Sys.time(), ii, length(chunks), min(reps_i), max(reps_i)))
      next
    }

    message(sprintf("[%s] Running chunk %d/%d: reps %d-%d",
                    Sys.time(), ii, length(chunks), min(reps_i), max(reps_i)))
    flush.console()

    results <- if (!is.null(cl)) {
      parallel::parLapply(cl, reps_i, function(rr) ppc_run_inad_rep(rr, cfg))
    } else {
      lapply(reps_i, ppc_run_inad_rep, cfg = cfg)
    }

    chunk <- list(
      reps = reps_i,
      scores = do.call(rbind, lapply(results, `[[`, "scores")),
      recovery = do.call(rbind, lapply(results, `[[`, "recovery")),
      diagnostics = do.call(rbind, lapply(results, `[[`, "diagnostics")),
      completed_at = Sys.time()
    )
    saveRDS(chunk, chunk_path)
    message(sprintf("[%s] Saved %s", Sys.time(), basename(chunk_path)))
    flush.console()
  }

  combined <- combine_r1000_outputs(out_dir)
  run_info <- data.frame(
    mode = "final_r1000_from_r200_plus_r800",
    n_reps = length(unique(combined$scores$rep)),
    n_sims = cfg$n_sims,
    n_per_group = cfg$n_per_group,
    fit_max_iter = cfg$fit_max_iter,
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds_remaining = as.numeric(difftime(Sys.time(), started, units = "secs")),
    completed_at = as.character(Sys.time())
  )

  final_dir <- file.path(ppc_poc_dir(root), "output", "inad", "final_r1000")
  dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(c(combined, list(run_info = run_info)), file.path(final_dir, "inad_prediction_final_r1000.rds"))
  utils::write.csv(combined$recovery, file.path(final_dir, "recovery_final_r1000.csv"), row.names = FALSE)
  utils::write.csv(combined$diagnostics, file.path(final_dir, "diagnostics_final_r1000.csv"), row.names = FALSE)
  utils::write.csv(combined$paired, file.path(final_dir, "paired_final_r1000.csv"), row.names = FALSE)
  utils::write.csv(combined$decision_summary, file.path(final_dir, "decision_summary_final_r1000.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(final_dir, "run_info_final_r1000.csv"), row.names = FALSE)

  list(out_dir = out_dir, final_dir = final_dir, run_info = run_info)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_inad_remaining_r800()
  print(result$run_info)
  cat("Remaining output:", result$out_dir, "\n")
  cat("Final output:", result$final_dir, "\n")
}
