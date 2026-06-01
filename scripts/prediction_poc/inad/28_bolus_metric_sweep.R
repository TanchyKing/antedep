.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)

fmt_delta_t <- function(mean, se) {
  t_stat <- mean / se
  t_label <- if (abs(abs(t_stat) - 2) < 0.02) sprintf("%.2f", t_stat) else sprintf("%.1f", t_stat)
  sprintf("%.3f (%s)", mean, t_label)
}

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

read_rolling_scores <- function(root) {
  inad_path <- file.path(
    ppc_output_dir(root, "bolus_leave_one_per_group_all_pairs_marginal_nb_glm"),
    "scores_bolus_leave_one_per_group.csv"
  )
  cgfm_path <- file.path(
    ppc_output_dir(root, "bolus_cgfm_leave_one_per_group_all_pairs"),
    "scores_cgfm_bolus_leave_one_per_group_all_pairs.csv"
  )
  if (!file.exists(inad_path)) stop("Missing INAD score file: ", inad_path)
  if (!file.exists(cgfm_path)) stop("Missing CGFM score file: ", cgfm_path)

  inad <- utils::read.csv(inad_path)
  cgfm <- utils::read.csv(cgfm_path)
  keep_models <- c(
    "inad_constrained_alpha", "inad_unconstrained", "poisson_glm", "nb_glm", "nb_glmm",
    "cgfm_independent", "cgfm_shared", "cgfm_time_varying"
  )
  dat <- rbind(
    inad[inad$model %in% keep_models, ],
    cgfm[cgfm$model %in% keep_models, ]
  )
  dat <- dat[dat$task == "rolling_one_step" & dat$horizon == 1L, ]
  dat$error <- dat$mean - dat$observed
  dat
}

patient_metric_values <- function(scores) {
  stats::aggregate(
    cbind(log_score, rps, squared_error, absolute_error, error, pit) ~
      heldout_subject + heldout_group + model,
    scores,
    mean,
    na.rm = TRUE
  )
}

make_metric_sweep <- function(patient) {
  patient$mse <- patient$squared_error
  patient$rmse <- sqrt(patient$squared_error)
  patient$mae <- patient$absolute_error
  patient$abs_bias <- abs(patient$error)
  patient$abs_pit_bias <- abs(patient$pit - 0.5)

  metric_map <- list(
    log_score = "Log score",
    rps = "RPS",
    rmse = "RMSE"
  )
  references <- c(
    "inad_unconstrained", "poisson_glm", "nb_glm", "nb_glmm",
    "cgfm_shared", "cgfm_time_varying"
  )
  labels <- c(
    inad_unconstrained = "Unconstrained INAD",
    poisson_glm = "Poisson GLM",
    nb_glm = "NB GLM / CGFM independent",
    nb_glmm = "NB GLMM",
    cgfm_shared = "CGFM shared",
    cgfm_time_varying = "CGFM time-varying"
  )

  anchor <- patient[patient$model == "inad_constrained_alpha", ]
  out <- list()
  for (ref in references) {
    rhs <- patient[patient$model == ref, ]
    merged <- merge(
      anchor,
      rhs,
      by = c("heldout_subject", "heldout_group"),
      suffixes = c("_inad", "_ref")
    )
    for (metric in names(metric_map)) {
      diff <- merged[[paste0(metric, "_inad")]] - merged[[paste0(metric, "_ref")]]
      s <- summarize_mean_se(diff)
      out[[length(out) + 1L]] <- data.frame(
        reference = ref,
        reference_label = unname(labels[ref]),
        metric = metric,
        metric_label = metric_map[[metric]],
        diff_mean = unname(s["mean"]),
        diff_se = unname(s["se"]),
        diff_t = unname(s["t"]),
        diff_q025 = unname(s["q025"]),
        diff_q975 = unname(s["q975"]),
        n_patients = length(diff),
        lower_is_better = TRUE,
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

wide_metric_table <- function(sweep) {
  order_refs <- c(
    "Unconstrained INAD", "Poisson GLM", "NB GLM / CGFM independent", "NB GLMM",
    "CGFM shared", "CGFM time-varying"
  )
  order_metrics <- c(
    "RPS", "RMSE", "Log score"
  )

  rows <- vector("list", length(order_refs))
  for (ii in seq_along(order_refs)) {
    ref <- order_refs[[ii]]
    sub <- sweep[sweep$reference_label == ref, ]
    vals <- setNames(rep(NA_character_, length(order_metrics)), order_metrics)
    for (metric in order_metrics) {
      z <- sub[sub$metric_label == metric, ]
      if (nrow(z)) vals[[metric]] <- fmt_delta_t(z$diff_mean, z$diff_se)
    }
    rows[[ii]] <- data.frame(Reference = ref, as.list(vals), check.names = FALSE)
  }
  do.call(rbind, rows)
}

run_bolus_metric_sweep <- function() {
  out_dir <- ppc_output_dir(root, "bolus_metric_sweep")
  scores <- read_rolling_scores(root)
  patient <- patient_metric_values(scores)
  sweep <- make_metric_sweep(patient)
  wide <- wide_metric_table(sweep)

  utils::write.csv(scores, file.path(out_dir, "scores_used_bolus_metric_sweep.csv"), row.names = FALSE)
  utils::write.csv(patient, file.path(out_dir, "patient_metric_values_bolus_metric_sweep.csv"), row.names = FALSE)
  utils::write.csv(sweep, file.path(out_dir, "metric_sweep_long.csv"), row.names = FALSE)
  utils::write.csv(wide, file.path(out_dir, "metric_sweep_wide.csv"), row.names = FALSE)
  list(out_dir = out_dir, long = sweep, wide = wide)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_bolus_metric_sweep()
  print(result$wide, row.names = FALSE)
  cat("Output:", result$out_dir, "\n")
}
