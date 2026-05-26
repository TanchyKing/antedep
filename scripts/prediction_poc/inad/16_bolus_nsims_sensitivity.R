.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(.ppc_script_dir, "15_bolus_smoothing_sensitivity.R"))

run_bolus_nsims_sensitivity <- function(fold_limit = 2L, n_sims_values = c(50L, 150L),
                                        fit_max_iter = 8L, smoke = TRUE) {
  out <- list()
  for (ns in n_sims_values) {
    res <- run_bolus_smoothing_sensitivity(
      fold_limit = fold_limit,
      n_sims = ns,
      fit_max_iter = fit_max_iter,
      smoke = TRUE
    )
    scores <- res$scores
    scores$n_sims_setting <- ns
    out[[as.character(ns)]] <- scores
  }
  scores_all <- do.call(rbind, out)

  agg <- stats::aggregate(
    cbind(log_score, rps, squared_error) ~
      n_sims_setting + fold + epsilon_label + model + task + horizon,
    scores_all[scores_all$epsilon_label == "main", ],
    mean,
    na.rm = TRUE
  )
  agg$rmse <- sqrt(agg$squared_error)
  summary <- stats::aggregate(
    cbind(log_score, rps, rmse) ~ n_sims_setting + model + task + horizon,
    agg,
    function(z) c(mean = mean(z), se = stats::sd(z) / sqrt(length(z)))
  )

  out_dir <- file.path(
    ppc_poc_dir(root), "output", "inad",
    if (smoke) "bolus_nsims_sensitivity_smoke" else "bolus_nsims_sensitivity"
  )
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  run_info <- data.frame(
    mode = basename(out_dir),
    n_folds = fold_limit,
    n_sims_values = paste(n_sims_values, collapse = ","),
    fit_max_iter = fit_max_iter
  )
  saveRDS(list(scores = scores_all, summary = summary, run_info = run_info),
          file.path(out_dir, "bolus_nsims_sensitivity.rds"))
  utils::write.csv(scores_all, file.path(out_dir, "scores_bolus_nsims_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(summary, file.path(out_dir, "summary_bolus_nsims_sensitivity.csv"), row.names = FALSE)
  utils::write.csv(run_info, file.path(out_dir, "run_info_bolus_nsims_sensitivity.csv"), row.names = FALSE)
  list(out_dir = out_dir, run_info = run_info, summary = summary)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_bolus_nsims_sensitivity()
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
