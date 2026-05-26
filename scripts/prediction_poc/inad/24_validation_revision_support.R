.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "holdout.R", "pipeline_inad.R", "extra_baselines.R"
), root)

ppc_loglik_value <- function(x) as.numeric(stats::logLik(x))

ppc_bic_row <- function(model, specification, fit_source, logLik, df, nobs,
                        comparable = TRUE, note = NA_character_) {
  data.frame(
    model = model,
    specification = specification,
    fit_source = fit_source,
    logLik = as.numeric(logLik),
    df = as.integer(df),
    nobs = as.integer(nobs),
    BIC = -2 * as.numeric(logLik) + as.integer(df) * log(as.integer(nobs)),
    comparable = comparable,
    note = note,
    stringsAsFactors = FALSE
  )
}

ppc_bolus_bic_table <- function(root) {
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(bolus_inad$blocks)
  n_subjects <- nrow(y)

  moment_rds <- file.path(
    ppc_output_dir(root, "marginal_moments"),
    "marginal_moments_bolus.rds"
  )
  if (!file.exists(moment_rds)) {
    stop("Missing marginal-moment RDS: ", moment_rds)
  }
  moment <- readRDS(moment_rds)
  fits <- moment$fits

  rows <- list()
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "INAD constrained-alpha",
    specification = "NBT-NBI-INADFE(1), constant alpha",
    fit_source = "local antedep constrained-alpha POC fit",
    logLik = fits$inad_constrained_alpha$log_l,
    df = 15L,
    nobs = n_subjects,
    comparable = TRUE,
    note = "BIC uses n_subjects; constrained POC inherits nb_inno_size from the unconstrained fit."
  )
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "INAD unconstrained",
    specification = "NBT-NBI-INADFE(1), time-varying alpha",
    fit_source = "local antedep fit",
    logLik = fits$inad_unconstrained$log_l,
    df = fits$inad_unconstrained$n_params,
    nobs = n_subjects,
    comparable = TRUE,
    note = "BIC uses n_subjects, matching antedep fit-object convention."
  )

  nb_ll <- stats::logLik(fits$nb_glm)
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "NB GLM",
    specification = "marginal no-lag: y ~ group + factor(time)",
    fit_source = "local MASS::glm.nb full-bolus fit",
    logLik = nb_ll,
    df = attr(nb_ll, "df"),
    nobs = attr(nb_ll, "nobs"),
    comparable = TRUE,
    note = "Marginal no-lag form, aligned with the dissertation BIC comparison."
  )

  glmm_ll <- stats::logLik(fits$nb_glmm)
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "NB GLMM",
    specification = "marginal-moment no-lag: y ~ group + time_fac + (1 | subject)",
    fit_source = "local glmmTMB full-bolus fit",
    logLik = glmm_ll,
    df = attr(glmm_ll, "df"),
    nobs = attr(glmm_ll, "nobs"),
    comparable = TRUE,
    note = "No-lag form used for in-sample marginal fit and closed-form marginal moments."
  )

  train_long <- ppc_long_lag1(y, blocks)
  pois <- stats::glm(
    y ~ group + factor(time) + log(y_lag1 + 1),
    family = stats::poisson(link = "log"),
    data = train_long
  )
  pois_ll <- stats::logLik(pois)
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "Poisson GLM",
    specification = "lag-aware predictive: y ~ group + factor(time) + log(y_lag1 + 1)",
    fit_source = "local stats::glm full-bolus lag-aware fit",
    logLik = pois_ll,
    df = attr(pois_ll, "df"),
    nobs = attr(pois_ll, "nobs"),
    comparable = TRUE,
    note = "Lag-aware form used in the predictive pipeline; asymmetric with marginal NB GLM."
  )

  if (!requireNamespace("tscount", quietly = TRUE)) {
    rows[[length(rows) + 1L]] <- data.frame(
      model = "tscount tsglm",
      specification = "exploratory pooled tsglm",
      fit_source = "local tscount unavailable",
      logLik = NA_real_, df = NA_integer_, nobs = NA_integer_,
      BIC = NA_real_, comparable = FALSE,
      note = "tscount is not installed.",
      stringsAsFactors = FALSE
    )
  } else {
    tsglm_fit <- fits$tscount$model
    tsglm_ll <- stats::logLik(tsglm_fit)
    rows[[length(rows) + 1L]] <- ppc_bic_row(
      model = "tscount tsglm",
      specification = "exploratory pooled log-linear NB tsglm",
      fit_source = "local tscount pooled subject-by-time fit",
      logLik = tsglm_ll,
      df = attr(tsglm_ll, "df"),
      nobs = attr(tsglm_ll, "nobs"),
      comparable = FALSE,
      note = paste(
        "Pooled artificial sequence; subject-boundary lag artifacts;",
        "not directly comparable to a proper panel likelihood."
      )
    )
  }

  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "CGFM independent frailty",
    specification = "published Henderson-Shimakura independent frailty",
    fit_source = "published Table 4 log-likelihood",
    logLik = -2191.8,
    df = 14L,
    nobs = n_subjects,
    comparable = TRUE,
    note = "Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance."
  )
  rows[[length(rows) + 1L]] <- ppc_bic_row(
    model = "CGFM shared frailty",
    specification = "published Henderson-Shimakura shared frailty",
    fit_source = "published Table 4 log-likelihood",
    logLik = -2249.2,
    df = 14L,
    nobs = n_subjects,
    comparable = TRUE,
    note = "Reconstructed from published full log-likelihood and 12 time effects + group + frailty variance."
  )

  out <- do.call(rbind, rows)
  out[order(out$BIC, na.last = TRUE), ]
}

ppc_simulation_precision <- function(root) {
  final_rds <- file.path(
    ppc_output_dir(root, "final_r1000"),
    "inad_prediction_final_r1000.rds"
  )
  if (!file.exists(final_rds)) stop("Missing final R=1000 RDS: ", final_rds)
  x <- readRDS(final_rds)
  scores <- x$scores
  keep <- scores$model == "inad_constrained_alpha" &
    scores$task == "rolling_one_step" &
    scores$horizon == 1L
  dat <- scores[keep, ]
  rep_summary <- stats::aggregate(
    cbind(log_score, rps, squared_error, absolute_error, pit, cover80, cover95) ~ rep,
    dat,
    mean,
    na.rm = TRUE
  )
  rep_summary$rmse <- sqrt(rep_summary$squared_error)
  metrics <- c("log_score", "rps", "rmse", "absolute_error", "pit", "cover80", "cover95")
  out <- lapply(metrics, function(metric) {
    z <- rep_summary[[metric]]
    data.frame(
      model = "INAD constrained-alpha",
      task = "rolling_one_step",
      horizon = 1L,
      metric = metric,
      mean = mean(z, na.rm = TRUE),
      se = stats::sd(z, na.rm = TRUE) / sqrt(sum(is.finite(z))),
      q025 = unname(stats::quantile(z, 0.025, na.rm = TRUE)),
      q975 = unname(stats::quantile(z, 0.975, na.rm = TRUE)),
      n_reps = sum(is.finite(z)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

ppc_smoothing_context <- function(root) {
  src <- file.path(
    ppc_output_dir(root, "bolus_smoothing_sensitivity"),
    "summary_bolus_smoothing_sensitivity.csv"
  )
  if (!file.exists(src)) stop("Missing smoothing summary: ", src)
  x <- utils::read.csv(src, check.names = FALSE)
  x$context <- paste(
    "Existing sensitivity run: 35-fold balanced rotation, recursive h-specific",
    "NB GLM comparison. Retained as methodological context because v2 headline",
    "uses rolling one-step and rolling is less exposed to recursive tail support."
  )
  x
}

run_validation_revision_support <- function() {
  out_dir <- ppc_output_dir(root, "validation_revision")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  bic <- ppc_bolus_bic_table(root)
  precision <- ppc_simulation_precision(root)
  smoothing <- ppc_smoothing_context(root)

  utils::write.csv(bic, file.path(out_dir, "bolus_bic_comparison.csv"), row.names = FALSE)
  utils::write.csv(precision, file.path(out_dir, "simulation_inad_precision_r1000.csv"), row.names = FALSE)
  utils::write.csv(smoothing, file.path(out_dir, "smoothing_sensitivity_context.csv"), row.names = FALSE)

  run_info <- data.frame(
    mode = "validation_revision_support",
    created = as.character(Sys.time()),
    source_final_r1000 = file.path(ppc_output_dir(root, "final_r1000"), "inad_prediction_final_r1000.rds"),
    source_marginal_moments = file.path(ppc_output_dir(root, "marginal_moments"), "marginal_moments_bolus.rds"),
    source_smoothing = file.path(ppc_output_dir(root, "bolus_smoothing_sensitivity"), "summary_bolus_smoothing_sensitivity.csv"),
    henderson_shimakura_pdf = "D:/UI/papers/Henderson&Shimakura.pdf",
    stringsAsFactors = FALSE
  )
  utils::write.csv(run_info, file.path(out_dir, "run_info_validation_revision.csv"), row.names = FALSE)

  list(out_dir = out_dir, bic = bic, precision = precision, smoothing = smoothing)
}

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_validation_revision_support()
  print(result$bic)
  print(result$precision)
  cat("Output:", result$out_dir, "\n")
}
