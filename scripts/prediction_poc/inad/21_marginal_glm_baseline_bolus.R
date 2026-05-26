.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(.ppc_script_dir, "11_bolus_prediction.R"))

ppc_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(as.integer(default))
  as.integer(value)
}

run_marginal_glm_bolus <- function(fold_design = c("rotation", "all_pairs"),
                                   n_sims = 1000L,
                                   fit_max_iter = 50L) {
  fold_design <- match.arg(fold_design)
  result <- run_bolus_prediction_design(
    fold_design = fold_design,
    n_sims = n_sims,
    fit_max_iter = fit_max_iter,
    out_label = "marginal_nb_glm"
  )
  result$run_info$nb_glm_baseline <- "marginal: y ~ group + factor(time)"
  utils::write.csv(
    result$run_info,
    file.path(result$out_dir, "run_info_bolus_marginal_nb_glm.csv"),
    row.names = FALSE
  )
  result
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  fold_design <- if ("all_pairs" %in% args) "all_pairs" else "rotation"
  result <- run_marginal_glm_bolus(
    fold_design = fold_design,
    n_sims = ppc_env_int("PPC_N_SIMS", 1000L),
    fit_max_iter = ppc_env_int("PPC_FIT_MAX_ITER", 50L)
  )
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
