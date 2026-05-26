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

run_inad_smoke_iter50 <- function() {
  cfg <- ppc_inad_config("smoke")
  cfg$fit_max_iter <- 50L
  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", "smoke_iter50")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
      ppc_source_common(c(
        "holdout.R", "sim_inad_forward.R", "constrained_inad.R",
        "score.R", "pipeline_inad.R"
      ), root)
      NULL
    })
    results <- parallel::parLapply(cl, reps, function(rr) ppc_run_inad_rep(rr, cfg))
  } else {
    results <- lapply(reps, ppc_run_inad_rep, cfg = cfg)
  }

  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  recovery <- do.call(rbind, lapply(results, `[[`, "recovery"))
  paired <- ppc_summarize_scores(scores)
  run_info <- data.frame(
    mode = "smoke_iter50",
    n_reps = cfg$n_reps,
    n_sims = cfg$n_sims,
    n_per_group = cfg$n_per_group,
    fit_max_iter = cfg$fit_max_iter,
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )

  saveRDS(
    list(config = cfg, scores = scores, recovery = recovery, paired = paired, run_info = run_info),
    file.path(out_dir, "inad_prediction_smoke_iter50.rds")
  )
  utils::write.csv(scores, file.path(out_dir, "scores_smoke_iter50.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(out_dir, "recovery_smoke_iter50.csv"), row.names = FALSE)
  utils::write.csv(paired, file.path(out_dir, "paired_smoke_iter50.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_smoke_iter50.csv"), row.names = FALSE)

  list(scores = scores, recovery = recovery, paired = paired, run_info = run_info, out_dir = out_dir)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_inad_smoke_iter50()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
