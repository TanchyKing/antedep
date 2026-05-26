.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_source_common("extra_baselines.R", root)

bolus_inad_param_counts <- function(n_time = 12L, n_blocks = 2L) {
  list(
    unconstrained = (n_time - 1L) + n_time + (n_blocks - 1L) + 1L,
    constrained_alpha = 1L + n_time + (n_blocks - 1L) + 1L
  )
}

add_bic_selected_bolus_all_pairs <- function() {
  out_dir <- file.path(
    ppc_poc_dir(root),
    "output", "inad", "bolus_leave_one_per_group_all_pairs"
  )
  rds_path <- file.path(out_dir, "inad_prediction_bolus_leave_one_per_group_all_pairs.rds")
  result <- readRDS(rds_path)

  scores <- result$scores
  diagnostics <- result$diagnostics
  n_subjects_train <- result$run_info$n_subjects[[1L]] - 2L
  counts <- bolus_inad_param_counts(
    n_time = result$config$n_time,
    n_blocks = length(unique(result$folds$heldout_group))
  )

  diag_wide <- reshape(
    diagnostics[, c("fold", "fit", "log_l")],
    idvar = "fold",
    timevar = "fit",
    direction = "wide"
  )
  names(diag_wide) <- sub("^log_l\\.", "log_l_", names(diag_wide))
  diag_wide$bic_unconstrained <- -2 * diag_wide$log_l_unconstrained +
    counts$unconstrained * log(n_subjects_train)
  diag_wide$bic_constrained_alpha <- -2 * diag_wide$log_l_constrained_alpha +
    counts$constrained_alpha * log(n_subjects_train)
  diag_wide$selected_fit <- ifelse(
    diag_wide$bic_constrained_alpha <= diag_wide$bic_unconstrained,
    "constrained_alpha",
    "unconstrained"
  )

  selected_rows <- list()
  for (ii in seq_len(nrow(diag_wide))) {
    fold <- diag_wide$fold[[ii]]
    model_name <- paste0("inad_", diag_wide$selected_fit[[ii]])
    rows <- scores[scores$fold == fold & scores$model == model_name, , drop = FALSE]
    rows$model <- "inad_bic_selected"
    selected_rows[[ii]] <- rows
  }
  selected_scores <- do.call(rbind, selected_rows)
  scores_augmented <- rbind(scores, selected_scores)

  paired_augmented <- ppc_summarize_scores_vs_references(
    scores_augmented,
    reference_models = c("nb_glm", "poisson_glm", "nb_glmm", "tscount_tsglm")
  )
  decision_augmented <- ppc_summarize_extra_decision(paired_augmented)

  selection_rate <- data.frame(
    selected_fit = names(table(diag_wide$selected_fit)),
    n_folds = as.integer(table(diag_wide$selected_fit)),
    proportion = as.numeric(table(diag_wide$selected_fit)) / nrow(diag_wide),
    row.names = NULL
  )

  saveRDS(
    list(
      selection = diag_wide,
      selection_rate = selection_rate,
      selected_scores = selected_scores,
      scores_augmented = scores_augmented,
      paired_augmented = paired_augmented,
      decision_augmented = decision_augmented,
      param_counts = counts,
      n_subjects_train = n_subjects_train
    ),
    file.path(out_dir, "bic_selected_bolus_all_pairs.rds")
  )
  utils::write.csv(diag_wide, file.path(out_dir, "bic_selection_by_fold_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(selection_rate, file.path(out_dir, "bic_selection_rate_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(selected_scores, file.path(out_dir, "scores_bic_selected_only_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(paired_augmented, file.path(out_dir, "paired_with_bic_selected_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(decision_augmented, file.path(out_dir, "decision_with_bic_selected_bolus_all_pairs.csv"), row.names = FALSE)

  list(out_dir = out_dir, selection_rate = selection_rate, decision_augmented = decision_augmented)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- add_bic_selected_bolus_all_pairs()
  print(result$selection_rate)
  cat("Output:", result$out_dir, "\n")
}
