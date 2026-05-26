.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "config.R", "holdout.R", "constrained_inad.R", "pipeline_inad.R"
), root)

fit_one_inad_prediction_rep <- function(mode = c("smoke", "full"), rep_id = 1L) {
  cfg <- ppc_inad_config(match.arg(mode))
  sim <- ppc_simulate_inad_dgp(cfg, rep_id)
  split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))
  y_train <- sim$y[split$train, , drop = FALSE]
  blocks_train <- sim$blocks[split$train]
  list(
    fits = ppc_fit_inad_rep(y_train, blocks_train, cfg),
    glm = ppc_fit_glm_rep(y_train, blocks_train),
    split = split,
    truth = sim$truth
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[[1L]] else "smoke"
  fit <- fit_one_inad_prediction_rep(mode)
  print(lapply(fit$fits, function(x) x$log_l))
}
