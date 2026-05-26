.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)

summarize_mean_se <- function(z) {
  z <- z[is.finite(z)]
  c(
    mean = mean(z),
    se = stats::sd(z) / sqrt(length(z)),
    t = mean(z) / (stats::sd(z) / sqrt(length(z))),
    q025 = unname(stats::quantile(z, 0.025)),
    q975 = unname(stats::quantile(z, 0.975))
  )
}

make_pairwise_observation_deltas <- function(scores,
                                             inad_models = c("inad_constrained_alpha", "inad_unconstrained"),
                                             reference_models = c("nb_glm", "poisson_glm", "nb_glmm", "tscount_tsglm")) {
  keys <- c("fold", "heldout_subject", "heldout_group", "task", "horizon", "time")
  keep <- c(keys, "model", "log_score", "rps", "squared_error")
  scores <- scores[, keep]
  out <- list()

  for (model_name in inad_models) {
    lhs <- scores[scores$model == model_name, ]
    names(lhs)[match(c("log_score", "rps", "squared_error"), names(lhs))] <-
      c("log_score_model", "rps_model", "squared_error_model")
    lhs$model <- NULL

    for (ref_name in reference_models) {
      rhs <- scores[scores$model == ref_name, ]
      if (!nrow(rhs)) next
      names(rhs)[match(c("log_score", "rps", "squared_error"), names(rhs))] <-
        c("log_score_reference", "rps_reference", "squared_error_reference")
      rhs$model <- NULL

      merged <- merge(lhs, rhs, by = keys)
      merged$model <- model_name
      merged$reference <- ref_name
      merged$diff_log_score <- merged$log_score_model - merged$log_score_reference
      merged$diff_rps <- merged$rps_model - merged$rps_reference
      merged$diff_rmse_proxy <- merged$squared_error_model - merged$squared_error_reference
      out[[paste(model_name, ref_name, sep = "_vs_")]] <- merged
    }
  }

  do.call(rbind, out)
}

summarize_patient_level <- function(deltas) {
  patient <- stats::aggregate(
    cbind(diff_log_score, diff_rps, diff_rmse_proxy) ~
      model + reference + task + horizon + heldout_subject + heldout_group,
    deltas,
    mean,
    na.rm = TRUE
  )

  out <- stats::aggregate(
    cbind(diff_log_score, diff_rps, diff_rmse_proxy) ~ model + reference + task + horizon,
    patient,
    summarize_mean_se
  )
  out
}

summarize_patient_rmse <- function(scores,
                                   inad_models = c("inad_constrained_alpha", "inad_unconstrained"),
                                   reference_models = c("nb_glm", "poisson_glm", "nb_glmm", "tscount_tsglm")) {
  keys <- c("heldout_subject", "heldout_group", "model", "task", "horizon")
  patient <- stats::aggregate(squared_error ~ heldout_subject + heldout_group + model + task + horizon,
                              scores, mean, na.rm = TRUE)
  patient$rmse <- sqrt(patient$squared_error)

  out <- list()
  for (model_name in inad_models) {
    lhs <- patient[patient$model == model_name, c("heldout_subject", "heldout_group", "task", "horizon", "rmse")]
    names(lhs)[names(lhs) == "rmse"] <- "rmse_model"
    for (ref_name in reference_models) {
      rhs <- patient[patient$model == ref_name, c("heldout_subject", "heldout_group", "task", "horizon", "rmse")]
      if (!nrow(rhs)) next
      names(rhs)[names(rhs) == "rmse"] <- "rmse_reference"
      merged <- merge(lhs, rhs, by = c("heldout_subject", "heldout_group", "task", "horizon"))
      merged$model <- model_name
      merged$reference <- ref_name
      merged$diff_rmse <- merged$rmse_model - merged$rmse_reference
      out[[paste(model_name, ref_name, sep = "_vs_")]] <- merged
    }
  }
  diffs <- do.call(rbind, out)
  stats::aggregate(diff_rmse ~ model + reference + task + horizon, diffs, summarize_mean_se)
}

summarize_constrained_vs_unconstrained <- function(scores) {
  keys <- c("heldout_subject", "heldout_group", "model", "task", "horizon")
  patient <- stats::aggregate(cbind(log_score, rps, squared_error) ~
                                heldout_subject + heldout_group + model + task + horizon,
                              scores, mean, na.rm = TRUE)
  patient$rmse <- sqrt(patient$squared_error)

  cfit <- patient[patient$model == "inad_constrained_alpha", c("heldout_subject", "heldout_group", "task", "horizon", "log_score", "rps", "rmse")]
  ufit <- patient[patient$model == "inad_unconstrained", c("heldout_subject", "heldout_group", "task", "horizon", "log_score", "rps", "rmse")]
  names(cfit)[5:7] <- paste0(names(cfit)[5:7], "_constrained")
  names(ufit)[5:7] <- paste0(names(ufit)[5:7], "_unconstrained")
  merged <- merge(cfit, ufit, by = c("heldout_subject", "heldout_group", "task", "horizon"))
  merged$diff_log_score <- merged$log_score_constrained - merged$log_score_unconstrained
  merged$diff_rps <- merged$rps_constrained - merged$rps_unconstrained
  merged$diff_rmse <- merged$rmse_constrained - merged$rmse_unconstrained

  stats::aggregate(cbind(diff_log_score, diff_rps, diff_rmse) ~ task + horizon,
                   merged, summarize_mean_se)
}

analyze_bolus_all_pairs <- function() {
  out_dir <- file.path(
    ppc_poc_dir(root),
    "output", "inad", "bolus_leave_one_per_group_all_pairs"
  )
  rds <- readRDS(file.path(out_dir, "inad_prediction_bolus_leave_one_per_group_all_pairs.rds"))
  scores <- rds$scores

  deltas <- make_pairwise_observation_deltas(scores)
  patient_summary <- summarize_patient_level(deltas)
  patient_rmse <- summarize_patient_rmse(scores)
  constrained_vs_unconstrained <- summarize_constrained_vs_unconstrained(scores)

  utils::write.csv(patient_summary, file.path(out_dir, "patient_level_score_deltas_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(patient_rmse, file.path(out_dir, "patient_level_rmse_deltas_bolus_all_pairs.csv"), row.names = FALSE)
  utils::write.csv(constrained_vs_unconstrained, file.path(out_dir, "patient_level_constrained_vs_unconstrained_bolus_all_pairs.csv"), row.names = FALSE)

  list(
    out_dir = out_dir,
    patient_summary = patient_summary,
    patient_rmse = patient_rmse,
    constrained_vs_unconstrained = constrained_vs_unconstrained
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- analyze_bolus_all_pairs()
  cat("Output:", result$out_dir, "\n")
}
