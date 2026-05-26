#!/usr/bin/env Rscript

root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)
if (basename(root) != "antedep") {
  stop("Run this script from the antedep package root.", call. = FALSE)
}

out_dir <- file.path(root, "scripts", "prediction_poc", "output", "inad")
dev_dir <- file.path(root, "inst", "dev-notes")

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "", formatC(as.numeric(x), format = "f", digits = digits))
}

fmt_t <- function(mean, se, digits = 3) {
  paste0(fmt(mean, digits), " (", fmt(mean / se, 1), ")")
}

md_table <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  names(x) <- gsub("\\.", " ", names(x))
  rows <- apply(x, 1, function(z) paste0("| ", paste(z, collapse = " | "), " |"))
  c(
    paste0("| ", paste(names(x), collapse = " | "), " |"),
    paste0("| ", paste(rep("---", ncol(x)), collapse = " | "), " |"),
    rows
  )
}

read_csv <- function(...) {
  read.csv(file.path(...), stringsAsFactors = FALSE, check.names = FALSE)
}

metric_table <- function(d, model_id) {
  sub <- d[d$model == model_id & d$task == "rolling_one_step", ]
  sub <- sub[order(sub$reference), ]
  data.frame(
    Reference = c(
      nb_glm = "NB GLM",
      nb_glmm = "NB GLMM",
      poisson_glm = "Poisson GLM",
      tscount_tsglm = "tscount"
    )[sub$reference],
    `Delta RPS (t)` = fmt_t(sub$`diff_rps_vs_reference.mean`, sub$`diff_rps_vs_reference.se`),
    `Delta RMSE (t)` = fmt_t(sub$`diff_rmse_vs_reference.mean`, sub$`diff_rmse_vs_reference.se`),
    `Delta log score (t)` = fmt_t(sub$`diff_log_score_vs_reference.mean`, sub$`diff_log_score_vs_reference.se`),
    check.names = FALSE
  )
}

mard_table <- function() {
  files <- list(
    `G1 mean` = "relative_error_mean_group_1_by_time.csv",
    `G1 variance` = "relative_error_variance_group_1_by_time.csv",
    `G2 mean` = "relative_error_mean_group_2_by_time.csv",
    `G2 variance` = "relative_error_variance_group_2_by_time.csv"
  )
  models <- c(
    "INAD constrained alpha",
    "INAD unconstrained",
    "NB GLM marginal",
    "NB GLMM marginal",
    "tscount tsglm MC",
    "CGFM independent frailty",
    "CGFM shared frailty"
  )
  res <- lapply(names(files), function(cell) {
    z <- read_csv(out_dir, "marginal_moments", files[[cell]])
    vals <- vapply(models, function(m) mean(abs(z[[m]]), na.rm = TRUE), numeric(1))
    data.frame(Cell = cell, Model = names(vals), MARD = vals, check.names = FALSE)
  })
  do.call(rbind, res)
}

rel_table <- function(file, title) {
  z <- read_csv(out_dir, "marginal_moments", file)
  keep <- c(
    "time", "empirical", "INAD constrained alpha", "INAD unconstrained",
    "NB GLM marginal", "NB GLMM marginal", "tscount tsglm MC",
    "CGFM independent frailty", "CGFM shared frailty"
  )
  z <- z[, keep]
  names(z) <- c(
    "Time", "Empirical", "C-alpha INAD", "U-INAD", "NB GLM",
    "NB GLMM", "tscount", "CGFM ind.", "CGFM shared"
  )
  z[-1] <- lapply(z[-1], function(v) fmt(v, 2))
  c(paste0("### ", title), "", md_table(z))
}

bic <- read_csv(out_dir, "validation_revision", "bolus_bic_comparison.csv")
precision <- read_csv(out_dir, "validation_revision", "simulation_inad_precision_r1000.csv")
smoothing <- read_csv(out_dir, "validation_revision", "smoothing_sensitivity_context.csv")
decision <- read_csv(out_dir, "bolus_leave_one_per_group_all_pairs_marginal_nb_glm",
                     "decision_summary_bolus_leave_one_per_group.csv")
recovery <- read_csv(out_dir, "final_r1000", "recovery_final_r1000.csv")

recovery_summary <- aggregate(
  cbind(alpha_mae, theta_rmse, nb_size_rmse, tau2_abs_bias = abs(tau2_bias)) ~ fit,
  data = recovery,
  FUN = mean
)
recovery_summary$fit <- ifelse(
  recovery_summary$fit == "constrained_alpha",
  "INAD constrained-alpha",
  "INAD unconstrained"
)
recovery_summary <- recovery_summary[, c("fit", "alpha_mae", "theta_rmse", "nb_size_rmse", "tau2_abs_bias")]
names(recovery_summary) <- c("Fit", "alpha MAE", "theta RMSE", "NB size RMSE", "|tau2 bias|")
recovery_summary[-1] <- lapply(recovery_summary[-1], fmt, digits = 3)

precision_show <- precision[precision$metric %in% c("log_score", "rps", "rmse", "pit", "cover80", "cover95"), ]
precision_show <- precision_show[, c("metric", "mean", "se", "q025", "q975", "n_reps")]
names(precision_show) <- c("Metric", "Mean", "SE", "2.5%", "97.5%", "Reps")
precision_show[, 2:5] <- lapply(precision_show[, 2:5], fmt, digits = 3)

bic_show <- bic[, c("model", "specification", "BIC", "comparable", "note")]
names(bic_show) <- c("Model", "Specification", "BIC", "Comparable", "Note")
bic_show$BIC <- fmt(bic_show$BIC, 1)
bic_show$Comparable <- ifelse(bic_show$Comparable, "yes", "no")

smooth_h2 <- smoothing[
  smoothing$model == "inad_constrained_alpha" &
    smoothing$task == "recursive_multi_step" &
    smoothing$horizon == 2,
]
eps_label <- c(large = "1/101", main = "1/1001", small = "1/10001")
smooth_show <- data.frame(
  Epsilon = eps_label[smooth_h2$epsilon_label],
  `Delta log score vs NB GLM (t)` = fmt_t(
    smooth_h2$`diff_log_score_vs_nb.mean`,
    smooth_h2$`diff_log_score_vs_nb.se`
  ),
  `Delta RPS vs NB GLM (t)` = fmt_t(
    smooth_h2$`diff_rps_vs_nb.mean`,
    smooth_h2$`diff_rps_vs_nb.se`
  ),
  check.names = FALSE
)

mards <- mard_table()
mard_wide <- reshape(mards, idvar = "Model", timevar = "Cell", direction = "wide")
names(mard_wide) <- sub("^MARD\\.", "", names(mard_wide))
mard_wide[,-1] <- lapply(mard_wide[,-1], fmt, digits = 2)

g1var <- read_csv(out_dir, "marginal_moments", "relative_error_variance_group_1_by_time.csv")
g1var_models <- names(g1var)[!(names(g1var) %in% c("time", "empirical"))]
g1var_sens <- data.frame(
  Model = g1var_models,
  `MARD all 12` = vapply(g1var_models, function(m) mean(abs(g1var[[m]])), numeric(1)),
  `MARD excl. t=11` = vapply(g1var_models, function(m) mean(abs(g1var[g1var$time != 11, m])), numeric(1)),
  check.names = FALSE
)
g1var_sens[, -1] <- lapply(g1var_sens[, -1], fmt, digits = 2)

constrained_roll <- metric_table(decision, "inad_constrained_alpha")
unconstrained_roll <- metric_table(decision, "inad_unconstrained")

lines <- c(
  "# INAD Validation Summary v2 - Prediction and Marginal Moments",
  "",
  "This v2 note implements the revision plan in `inad-validation-revision-todo.md`.",
  "The main change from v1 is scope discipline: the simulation section is",
  "now INAD-only correct-specification validation, and the bolus predictive",
  "section is restricted to rolling one-step forecasting. Recursive multi-step",
  "tables are retained only as smoothing-sensitivity context.",
  "",
  "## 1. Scope and setup",
  "",
  "### 1.1 What this note covers",
  "",
  "- Simulation: absolute prediction precision for constrained-alpha INAD under its own NBT-NBI-INADFE(1) DGP.",
  "- Bolus in-sample fit: BIC comparison across INAD and non-INAD count models.",
  "- Bolus prediction: rolling one-step leave-one-per-group-out prediction on all 30 x 35 = 1050 held-out pairs.",
  "- Bolus marginal moments: Henderson-and-Shimakura-style time-by-time mean and variance comparison.",
  "- Package implications: what already landed in `antedep`, and what remains.",
  "",
  "CGFM is included only in the marginal-moment comparison, using published Henderson-Shimakura estimates. No R implementation is available for refitting CGFM inside each cross-validation fold, so CGFM is not included in the predictive comparison.",
  "",
  "### 1.2 Evaluation framework",
  "",
  "The predictive target is conditional new-subject forecasting. A model is fit on training subjects, early observations of held-out subjects are revealed, and later observations are predicted from the fitted conditional distribution. The v2 headline uses rolling one-step prediction at t = 2, ..., 12, conditioning on the realized history through t - 1.",
  "",
  "All predictive comparisons are plug-in MLE conditional. INAD uses Monte Carlo forward simulation from the fitted process; NB and Poisson GLMs use analytic plug-in PMFs; NB GLMM uses the fitted random-intercept model for prediction; tscount uses a pooled log-linear count time-series fit. Lower RMSE, lower log score, and lower RPS are better.",
  "",
  "### 1.3 Metric primer",
  "",
  "- RMSE measures point prediction error of the predictive mean, on the original count scale.",
  "- Log score is `-log p(Y_obs)`. It rewards sharp probability at the realized count but is sensitive to finite-support smoothing when the PMF is simulation-based.",
  "- RPS compares the predictive CDF with the realized count indicator over the full integer support. It is a proper probabilistic score and is less sensitive than log score to a single tail probability.",
  "",
  "## 2. Models compared - formulations",
  "",
  "Notation: subject i has group g(i), time t = 1, ..., T, and group 1 is the reference.",
  "",
  "### 2.1 NB GLM",
  "",
  "`Y_it ~ NB(mu_it, theta)`, with `log mu_it = beta0 + beta_group[g(i)] + beta_time[t]`. This is the marginal no-lag form used for the BIC and marginal-moment comparison.",
  "",
  "### 2.2 NB GLMM",
  "",
  "Two forms are used. The predictive form is lag-aware: `log mu_it = beta0 + beta_group[g(i)] + beta_time[t] + gamma log(Y_i,t-1 + 1) + b_i`, with `b_i ~ N(0, sigma_b^2)`. The marginal-moment/BIC form drops the lag term so closed-form marginal moments can be obtained by integrating the random intercept.",
  "",
  "### 2.3 Poisson GLM",
  "",
  "The predictive pipeline uses the lag-aware form `Y_it ~ Poisson(mu_it)`, `log mu_it = beta0 + beta_group[g(i)] + beta_time[t] + gamma log(Y_i,t-1 + 1)`. Rows at t = 1 are dropped because no lag is available.",
  "",
  "### 2.4 tscount",
  "",
  "The implemented comparator is an exploratory pooled `tscount::tsglm` fit. Training subjects are concatenated into one artificial sequence, with `model = list(past_obs = 1)` and `xreg = model.matrix(~ group + time_fac, ...)`. This introduces subject-boundary lag artifacts, so the BIC row is labelled non-comparable and the model is treated as a surrogate count-time-series baseline rather than a clean panel likelihood.",
  "",
  "### 2.5 INAD constrained-alpha",
  "",
  "NBT-NBI-INADFE(1): `Y_it = alpha o Y_i,t-1 + epsilon_it`, with constant alpha, negative-binomial thinning, and NB innovations with time-varying mean/size plus group effect tau. This is the BIC-supported primary INAD variant for bolus.",
  "",
  "### 2.6 INAD unconstrained",
  "",
  "Same model class, but alpha is allowed to vary over time. It is reported as a sensitivity because the bolus LRT detects mild non-stationarity even though BIC prefers the constrained model.",
  "",
  "### 2.7 CGFM independent/shared frailty",
  "",
  "CGFM rows use Henderson-Shimakura published full-data estimates only. In the local notation, the gamma frailty has mean 1 and frailty variance psi. Independent frailty has time-specific independent frailties; shared frailty uses a single subject-level frailty across all times. These rows are used for marginal mean and variance only.",
  "",
  "## 3. Simulation - INAD prediction precision",
  "",
  "DGP: NBT-NBI-INADFE(1), constant alpha = 0.35, time-varying NB innovation parameters, tau = c(0, 1.25), 100 subjects, 12 time points, 1000 simulation replicates. Because the DGP is INAD, this section is correct-specification validation rather than a neutral benchmark; cross-model baseline comparisons are not retained here.",
  "",
  "### 3.1 Parameter recovery",
  "",
  md_table(recovery_summary),
  "",
  "The constrained-alpha fit improves alpha and theta recovery, as expected under a constant-alpha DGP. NB innovation size is shared by construction in the current constrained-alpha POC.",
  "",
  "### 3.2 Rolling one-step precision",
  "",
  md_table(precision_show),
  "",
  "PIT is close to 0.5 and the 80%/95% empirical coverages are slightly conservative, which is acceptable for the simulation validation target.",
  "",
  "## 4. Bolus in-sample fit - BIC comparison",
  "",
  "BIC is computed on the full bolus data. The tscount row has the best numeric BIC but is not directly comparable because it is a pooled artificial sequence rather than a proper subject-level panel likelihood. Among comparable rows, constrained-alpha INAD has the lowest BIC, with NB GLMM close behind.",
  "",
  md_table(bic_show),
  "",
  "## 5. Bolus prediction - rolling one-step",
  "",
  "Design: all 30 x 35 cross-group held-out patient pairs, each fold holding out one group-1 and one group-2 patient and fitting on the remaining 63 patients. Standard errors are patient-level clustered, so overlapping folds do not receive naive independent-fold uncertainty.",
  "",
  "### 5.1 Stationarity and primary-model choice",
  "",
  "The full-data LRT for constant alpha rejects (p = 0.000234), but BIC strongly prefers constrained-alpha INAD over unconstrained INAD (4244 vs 4298). The interpretation is mild detectable non-stationarity that is not worth 11 extra alpha parameters under BIC. Constrained-alpha INAD is therefore primary; unconstrained INAD is sensitivity.",
  "",
  "### 5.2 Constrained-alpha INAD vs baselines",
  "",
  "Deltas are INAD minus reference; negative is INAD better.",
  "",
  md_table(constrained_roll),
  "",
  "### 5.3 Unconstrained INAD vs baselines",
  "",
  md_table(unconstrained_roll),
  "",
  "### 5.4 Smoothing sensitivity context",
  "",
  "The v2 headline is rolling one-step, but the existing recursive h = 2 smoothing sweep is retained as methodological context. It shows why log score is treated more cautiously than RPS: changing epsilon changes the single-point tail mass more than it changes integrated scores.",
  "",
  md_table(smooth_show),
  "",
  "### 5.5 Bolus prediction headline",
  "",
  "On rolling one-step prediction, constrained-alpha INAD beats the marginal NB GLM, NB GLMM, Poisson GLM, and tscount on RPS. It also beats all four on RMSE except that the Poisson RMSE margin is small and only borderline by t-statistic. Log score is better than NB GLM, NB GLMM, and Poisson GLM, and essentially tied with tscount. The unconstrained sensitivity tells the same qualitative story.",
  "",
  "## 6. Bolus marginal mean and variance comparison",
  "",
  "This section compares full-data fitted marginal moments with empirical bolus moments. CGFM is included here only, using published Henderson-Shimakura estimates. Relative discrepancy is `100 * (fitted - empirical) / empirical`; MARD is the mean absolute relative discrepancy.",
  "",
  "### 6.1 MARD summary",
  "",
  md_table(mard_wide),
  "",
  "### 6.2 G1 variance sensitivity",
  "",
  "G1 time 11 has a very small empirical variance, so it can dominate relative-error summaries. The ranking is robust when it is removed.",
  "",
  md_table(g1var_sens),
  "",
  rel_table("relative_error_mean_group_1_by_time.csv", "6.3 Group 1 mean relative discrepancy (%)"),
  "",
  rel_table("relative_error_variance_group_1_by_time.csv", "6.4 Group 1 variance relative discrepancy (%)"),
  "",
  rel_table("relative_error_mean_group_2_by_time.csv", "6.5 Group 2 mean relative discrepancy (%)"),
  "",
  rel_table("relative_error_variance_group_2_by_time.csv", "6.6 Group 2 variance relative discrepancy (%)"),
  "",
  "### 6.7 Marginal-moment headline",
  "",
  "The INAD model class gives the closest overall reproduction of the empirical marginal moments. Unconstrained INAD is best for G1 means, G2 means, and G2 variances; constrained-alpha INAD is best for G1 variances. This supports reporting both INAD variants in the marginal-moment table rather than forcing one bolus variant to carry every criterion.",
  "",
  "## 7. Cross-cutting findings",
  "",
  "- Predictions are plug-in MLE conditional, not posterior predictive.",
  "- INAD predictive samples are integer-valued by construction; predictive means are real-valued and are the right point forecasts for RMSE.",
  "- Tau confidence intervals are profile-likelihood intervals. Tau SEs are derived from profile-CI width for internal consistency. A 3+ group unit test remains pending.",
  "- The engineering pipeline completed the R = 1000 simulation, 1050 bolus all-pairs prediction, and full-data marginal-moment drivers without fit failures.",
  "",
  "## 8. Implications for antedep",
  "",
  "Already landed: `.simulate_inad_forward()` and public `predict.inad_fit()` v1 in `R/predict_inad.R`, with `type = c(\"mean\", \"sample\")` and tests in `tests/testthat/test-predict-inad.R`.",
  "",
  "Remaining: a shared cross-family `type = \"distribution\"` API after the gau/cat POCs; a proper `alpha_constraint = \"constant\"` path in `fit_inad()` with joint re-estimation; exact one-step PMFs as future methodological cleanup; and optional scoring helpers only if distribution-valued prediction becomes public.",
  "",
  "Next sequence: start the gau and cat script-first prediction POCs, stop at smoke-test validation, then decide the shared prediction API before moving gau/cat prediction methods into the package.",
  ""
)

v2_path <- file.path(dev_dir, "inad-validation-summary-v2.md")
canonical_path <- file.path(dev_dir, "inad-validation-summary.md")
writeLines(lines, v2_path, useBytes = TRUE)
writeLines(lines, canonical_path, useBytes = TRUE)

message("Wrote: ", v2_path)
message("Updated: ", canonical_path)
