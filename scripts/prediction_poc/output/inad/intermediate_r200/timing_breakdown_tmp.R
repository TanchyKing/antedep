source("scripts/prediction_poc/common/config.R")
root <- ppc_find_package_root()
ppc_load_antedep(root)
ppc_source_common(c("holdout.R", "sim_inad_forward.R", "constrained_inad.R", "score.R", "pipeline_inad.R"), root)
cfg <- ppc_inad_config("full")
cfg$n_reps <- 1L
cfg$n_sims <- 1000L
cfg$n_per_group <- 50L
cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
cfg$n_subjects <- length(cfg$blocks)
cfg$fit_max_iter <- 50L
rep_id <- 1L
sim <- ppc_simulate_inad_dgp(cfg, rep_id)
split <- ppc_stratified_split(sim$blocks, cfg$train_prop, seed = ppc_rep_seed(cfg, rep_id, 2L))
y_train <- sim$y[split$train, , drop = FALSE]
y_test <- sim$y[split$test, , drop = FALSE]
blocks_train <- sim$blocks[split$train]
blocks_test <- sim$blocks[split$test]

t_fit_inad <- system.time({ fits <- ppc_fit_inad_rep(y_train, blocks_train, cfg) })
t_fit_glm <- system.time({ glm_fits <- ppc_fit_glm_rep(y_train, blocks_train) })

# INAD forward prediction only: rolling + recursive, both INAD fits.
fit_names <- names(fits)
t_inad_forward <- system.time({
  for (tt in 2:cfg$n_time) {
    for (ff in fit_names) {
      simulate_inad_forward(
        fits[[ff]],
        history = y_test[, 1:(tt - 1L), drop = FALSE],
        blocks = blocks_test,
        start_time = tt,
        h = 1L,
        n_sims = cfg$n_sims,
        seed = ppc_rep_seed(cfg, rep_id, 1000L + 100L * match(ff, fit_names) + tt)
      )
    }
  }
  for (ff in fit_names) {
    simulate_inad_forward(
      fits[[ff]],
      history = y_test[, 1:cfg$t_split, drop = FALSE],
      blocks = blocks_test,
      start_time = cfg$t_split + 1L,
      h = cfg$recursive_h,
      n_sims = cfg$n_sims,
      seed = ppc_rep_seed(cfg, rep_id, 2000L + 100L * match(ff, fit_names))
    )
  }
})

# GLM mean prediction only: rolling + recursive, both GLMs.
t_glm_predict <- system.time({
  for (tt in 2:cfg$n_time) {
    ppc_predict_glm_means_rolling(glm_fits, y_test, blocks_test, tt)
  }
  ppc_predict_glm_means_recursive(glm_fits, y_test, blocks_test, cfg)
})

t_all_predict_score <- system.time({
  scores <- ppc_predict_and_score_rep(rep_id, fits, glm_fits, y_test, blocks_test, cfg)
})

out <- data.frame(
  component = c("fit_two_inad", "fit_two_glm", "inad_forward_only", "glm_predict_only", "combined_prediction_and_scoring"),
  elapsed_sec = c(t_fit_inad[["elapsed"]], t_fit_glm[["elapsed"]], t_inad_forward[["elapsed"]], t_glm_predict[["elapsed"]], t_all_predict_score[["elapsed"]]),
  user_sec = c(t_fit_inad[["user.self"]], t_fit_glm[["user.self"]], t_inad_forward[["user.self"]], t_glm_predict[["user.self"]], t_all_predict_score[["user.self"]]),
  system_sec = c(t_fit_inad[["sys.self"]], t_fit_glm[["sys.self"]], t_inad_forward[["sys.self"]], t_glm_predict[["sys.self"]], t_all_predict_score[["sys.self"]])
)
print(out)
