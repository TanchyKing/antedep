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

run_inad_extra_baselines_r200 <- function(include_tscount = TRUE, include_nb_glmm = TRUE,
                                          install_missing = FALSE) {
  if (install_missing) {
    ppc_install_extra_baseline_packages()
  }
  pkg_status <- ppc_extra_package_status()
  requested <- c(tscount = include_tscount, glmmTMB = include_nb_glmm)
  missing <- names(requested)[requested & !pkg_status$installed[match(names(requested), pkg_status$package)]]
  if (length(missing)) {
    stop(
      "Missing optional baseline package(s): ", paste(missing, collapse = ", "),
      ". Run scripts/prediction_poc/inad/08_install_extra_baselines.R after the current long job finishes, ",
      "or call run_inad_extra_baselines_r200(install_missing = TRUE)."
    )
  }

  cfg <- ppc_inad_config("full")
  cfg$n_reps <- 200L
  cfg$n_sims <- 1000L
  cfg$n_per_group <- 50L
  cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
  cfg$n_subjects <- length(cfg$blocks)
  cfg$fit_max_iter <- 50L
  cfg$extra_baselines <- list(
    nb_glmm = include_nb_glmm,
    tscount = include_tscount
  )

  out_dir <- file.path(ppc_poc_dir(root), "output", "inad", "extra_baselines_r200")
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
        "score.R", "pipeline_inad.R", "extra_baselines.R"
      ), root)
      NULL
    })
    results <- parallel::parLapply(cl, reps, function(rr) ppc_run_inad_rep_extra(rr, cfg))
  } else {
    results <- lapply(reps, ppc_run_inad_rep_extra, cfg = cfg)
  }

  scores <- do.call(rbind, lapply(results, `[[`, "scores"))
  recovery <- do.call(rbind, lapply(results, `[[`, "recovery"))
  diagnostics <- do.call(rbind, lapply(results, `[[`, "diagnostics"))
  extra_diagnostics <- do.call(rbind, lapply(results, `[[`, "extra_diagnostics"))
  paired <- ppc_summarize_scores_vs_references(scores)
  decision_summary <- ppc_summarize_extra_decision(paired)

  run_info <- data.frame(
    mode = "extra_baselines_r200",
    n_reps = cfg$n_reps,
    n_sims = cfg$n_sims,
    n_per_group = cfg$n_per_group,
    fit_max_iter = cfg$fit_max_iter,
    include_tscount = include_tscount,
    include_nb_glmm = include_nb_glmm,
    tscount_version = pkg_status$version[pkg_status$package == "tscount"],
    glmmTMB_version = pkg_status$version[pkg_status$package == "glmmTMB"],
    n_cores = cfg$n_cores,
    n_workers = cfg$n_workers,
    seed = cfg$seed,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )

  result <- list(
    config = cfg,
    package_status = pkg_status,
    scores = scores,
    recovery = recovery,
    diagnostics = diagnostics,
    extra_diagnostics = extra_diagnostics,
    paired = paired,
    decision_summary = decision_summary,
    run_info = run_info
  )

  saveRDS(result, file.path(out_dir, "inad_prediction_extra_baselines_r200.rds"))
  utils::write.csv(scores, file.path(out_dir, "scores_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(recovery, file.path(out_dir, "recovery_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(diagnostics, file.path(out_dir, "diagnostics_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(extra_diagnostics, file.path(out_dir, "extra_diagnostics_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(paired, file.path(out_dir, "paired_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(decision_summary, file.path(out_dir, "decision_summary_extra_baselines_r200.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_extra_baselines_r200.csv"), row.names = FALSE)

  list(out_dir = out_dir, run_info = run_info, decision_summary = decision_summary)
}

.ppc_is_main_script <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  length(file_arg) &&
    identical(
      normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/", mustWork = FALSE),
      normalizePath(file.path(.ppc_script_dir, "09_extra_baselines_r200.R"), winslash = "/", mustWork = FALSE)
    )
})

if (.ppc_is_main_script && identical(environment(), globalenv()) && !interactive()) {
  result <- run_inad_extra_baselines_r200()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
