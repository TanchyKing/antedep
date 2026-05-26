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

run_inad_smoke_n50_iter50 <- function() {
  cfg <- ppc_inad_config("smoke")
  cfg$n_per_group <- 50L
  cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
  cfg$n_subjects <- length(cfg$blocks)
  cfg$fit_max_iter <- 50L

  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", "smoke_n50_iter50")
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
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  paired <- ppc_summarize_scores(scores)

  # One-rep adjunct: omit truth-based initialization and let fit_inad use its defaults.
  sim1 <- ppc_simulate_inad_dgp(cfg, 1L)
  split1 <- ppc_stratified_split(sim1$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, 1L, 2L))
  y_train1 <- sim1$y[split1$train, , drop = FALSE]
  blocks_train1 <- sim1$blocks[split1$train]
  default_fit <- ppc_fit_inad_rep_default_init(y_train1, blocks_train1, cfg)
  default_recovery <- ppc_recovery_rep(1L, list(default_init = default_fit), sim1$truth)
  default_diagnostics <- ppc_fit_diagnostics_rep(1L, list(default_init = default_fit))

  run_info <- data.frame(
    mode = "smoke_n50_iter50",
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
    list(
      config = cfg, scores = scores, recovery = recovery,
      diagnostics = diagnostics, paired = paired,
      default_recovery = default_recovery,
      default_diagnostics = default_diagnostics,
      run_info = run_info
    ),
    file.path(out_dir, "inad_prediction_smoke_n50_iter50.rds")
  )
  utils::write.csv(scores, file.path(out_dir, "scores_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(out_dir, "recovery_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(diagnostics, file.path(out_dir, "diagnostics_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(paired, file.path(out_dir, "paired_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(default_recovery, file.path(out_dir, "default_recovery_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(default_diagnostics, file.path(out_dir, "default_diagnostics_smoke_n50_iter50.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_smoke_n50_iter50.csv"), row.names = FALSE)

  list(
    scores = scores, recovery = recovery, diagnostics = diagnostics,
    paired = paired, default_recovery = default_recovery,
    default_diagnostics = default_diagnostics, run_info = run_info,
    out_dir = out_dir
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_inad_smoke_n50_iter50()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
