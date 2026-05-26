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

predict_one_inad_prediction_rep <- function(mode = c("smoke", "full"), rep_id = 1L) {
  cfg <- ppc_inad_config(match.arg(mode))
  sim <- ppc_simulate_inad_dgp(cfg, rep_id)
  split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))
  y_train <- sim$y[split$train, , drop = FALSE]
  y_test <- sim$y[split$test, , drop = FALSE]
  blocks_train <- sim$blocks[split$train]
  blocks_test <- sim$blocks[split$test]
  fits <- ppc_fit_inad_rep(y_train, blocks_train, cfg)
  glm_fits <- ppc_fit_glm_rep(y_train, blocks_train)
  ppc_predict_and_score_rep(rep_id, fits, glm_fits, y_test, blocks_test, cfg)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[[1L]] else "smoke"
  scores <- predict_one_inad_prediction_rep(mode)
  print(utils::head(scores))
}
