.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "config.R", "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
  "score.R", "pipeline_inad.R"
), root)

run_inad_prediction_poc <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  cfg <- ppc_inad_config(mode)
  out_dir <- ppc_output_dir(root, mode)

  started <- Sys.time()
  reps <- seq_len(cfg$n_reps)

  if (cfg$n_workers > 1L && length(reps) > 1L) {
    cl <- parallel::makeCluster(min(cfg$n_workers, length(reps)))
    on.exit(parallel::stopCluster(cl), add = TRUE)
    parallel::clusterExport(cl, c("root", "mode"), envir = environment())
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
    results <- parallel::parLapply(cl, reps, function(rr) {
      cfg <- ppc_inad_config(mode)
      ppc_run_inad_rep(rr, cfg)
    })
  } else {
    results <- lapply(reps, ppc_run_inad_rep, cfg = cfg)
  }

  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  recovery <- do.call(rbind, lapply(results, `[[`, "recovery"))
  paired <- ppc_summarize_scores(scores)
  elapsed <- difftime(Sys.time(), started, units = "secs")
  run_info <- data.frame(
    mode = mode,
    n_reps = cfg$n_reps,
    n_sims = cfg$n_sims,
    n_per_group = cfg$n_per_group,
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds = as.numeric(elapsed)
  )

  saveRDS(list(config = cfg, scores = scores, recovery = recovery, paired = paired, run_info = run_info),
          file.path(out_dir, paste0("inad_prediction_", mode, ".rds")))
  utils::write.csv(scores, file.path(out_dir, paste0("scores_", mode, ".csv")), row.names = FALSE)
  utils::write.csv(recovery, file.path(out_dir, paste0("recovery_", mode, ".csv")), row.names = FALSE)
  utils::write.csv(paired, file.path(out_dir, paste0("paired_", mode, ".csv")), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, paste0("run_info_", mode, ".csv")), row.names = FALSE)

  list(scores = scores, recovery = recovery, paired = paired, run_info = run_info, out_dir = out_dir)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[[1L]] else "smoke"
  result <- run_inad_prediction_poc(mode)
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
