#!/usr/bin/env Rscript

# Example: Missing-data workflow for AD models in antedep

load_antedep <- function() {
  if (requireNamespace("antedep", quietly = TRUE)) {
    suppressPackageStartupMessages(library(antedep))
    return(invisible(TRUE))
  }

  if (!file.exists("DESCRIPTION")) {
    stop("Run this script from the package root, or install 'antedep' first.")
  }
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("Install 'pkgload' to run from source: install.packages('pkgload').")
  }

  pkgload::load_all(
    ".",
    export_all = FALSE,
    helpers = FALSE,
    attach_testthat = FALSE,
    quiet = TRUE
  )
  invisible(TRUE)
}

is_monotone_dropout <- function(row_values) {
  missing_idx <- which(is.na(row_values))
  if (length(missing_idx) == 0) {
    return(FALSE)
  }
  identical(missing_idx, seq.int(min(missing_idx), length(row_values)))
}

summarize_missing <- function(y) {
  has_row_missing <- rowSums(is.na(y)) > 0
  n_incomplete <- sum(has_row_missing)
  n_dropout <- 0L
  if (n_incomplete > 0) {
    n_dropout <- sum(apply(y[has_row_missing, , drop = FALSE], 1, is_monotone_dropout))
  }

  list(
    n_subjects = nrow(y),
    n_time = ncol(y),
    n_missing = sum(is.na(y)),
    pct_missing = 100 * mean(is.na(y)),
    n_complete = sum(!has_row_missing),
    n_dropout = n_dropout,
    n_intermittent = n_incomplete - n_dropout
  )
}

print_missing <- function(info) {
  cat("Data dimensions:", info$n_subjects, "x", info$n_time, "\n")
  cat("Missing values:", info$n_missing, "\n")
  cat("Percent missing:", sprintf("%.1f%%", info$pct_missing), "\n")
  cat("Complete subjects:", info$n_complete, "\n")
  cat("Dropout subjects:", info$n_dropout, "\n")
  cat("Intermittent subjects:", info$n_intermittent, "\n\n")
}

print_fit <- function(name, fit_obj) {
  cat(name, "\n", sep = "")
  cat("  log-likelihood:", round(fit_obj$log_l, 2), "\n")
  cat("  AIC:", round(fit_obj$aic, 2), "\n")
  cat("  BIC:", round(fit_obj$bic, 2), "\n")
  if (!is.null(fit_obj$em_converged)) {
    cat("  EM converged:", fit_obj$em_converged, "\n")
    cat("  EM iterations:", fit_obj$em_iterations, "\n")
  }
  cat("\n")
}

load_antedep()
set.seed(123)

cat("\n=== Example 1: Monotone dropout ===\n")
n_subjects <- 80
n_time <- 6
y_complete <- simulate_gau(
  n_subjects = n_subjects,
  n_time = n_time,
  order = 1,
  mu = 10,
  phi = 0.55,
  sigma = 1.4
)

y_dropout <- y_complete
dropout_subjects <- sample(seq_len(n_subjects), size = 24)
for (s in dropout_subjects) {
  n_drop <- sample(1:2, 1)
  y_dropout[s, (n_time - n_drop + 1):n_time] <- NA
}
print_missing(summarize_missing(y_dropout))

fit_em <- fit_gau(
  y_dropout,
  order = 1,
  na_action = "em",
  em_max_iter = 100,
  em_verbose = FALSE
)
fit_cc <- fit_gau(y_dropout, order = 1, na_action = "complete")
print_fit("EM fit:", fit_em)
print_fit("Complete-case fit:", fit_cc)

cat("Subjects used by complete-case fit:", fit_cc$settings$n_subjects, "/", n_subjects, "\n")
cat("mu (EM):", paste(round(fit_em$mu, 2), collapse = ", "), "\n")
cat("mu (complete-case):", paste(round(fit_cc$mu, 2), collapse = ", "), "\n")
cat("mu (true sample mean):", paste(round(colMeans(y_complete), 2), collapse = ", "), "\n\n")

cat("=== Example 2: Intermittent (MCAR) missingness ===\n")
y_mcar <- y_complete
mcar_idx <- sample(length(y_mcar), size = round(0.15 * length(y_mcar)))
y_mcar[mcar_idx] <- NA
print_missing(summarize_missing(y_mcar))

fit_mcar_em <- fit_gau(
  y_mcar,
  order = 1,
  na_action = "em",
  em_max_iter = 100,
  em_verbose = FALSE
)
print_fit("EM fit (MCAR):", fit_mcar_em)

cat("=== Example 3: Block effects with missing data ===\n")
n_per_group <- 40
blocks <- c(rep(1L, n_per_group), rep(2L, n_per_group))
y_blocks <- simulate_gau(
  n_subjects = 2 * n_per_group,
  n_time = n_time,
  order = 1,
  mu = 10,
  phi = 0.5,
  sigma = 1.5,
  blocks = blocks,
  tau = c(0, 2)
)

y_blocks_missing <- y_blocks
blocks_missing_idx <- sample(length(y_blocks_missing), size = round(0.10 * length(y_blocks_missing)))
y_blocks_missing[blocks_missing_idx] <- NA
print_missing(summarize_missing(y_blocks_missing))

fit_blocks <- fit_gau(
  y_blocks_missing,
  order = 1,
  blocks = blocks,
  na_action = "em",
  em_max_iter = 100,
  em_verbose = FALSE
)
print_fit("EM fit with blocks:", fit_blocks)
cat("Estimated tau:", paste(round(fit_blocks$tau, 2), collapse = ", "), "\n\n")

cat("=== Example 4: logL_gau with different missing-data actions ===\n")
ll_marginalize <- logL_gau(
  y_dropout,
  order = 1,
  mu = fit_em$mu,
  phi = fit_em$phi,
  sigma = fit_em$sigma,
  na_action = "marginalize"
)
ll_complete <- logL_gau(
  y_dropout,
  order = 1,
  mu = fit_em$mu,
  phi = fit_em$phi,
  sigma = fit_em$sigma,
  na_action = "complete"
)

cat("logL (marginalize):", round(ll_marginalize, 2), "\n")
cat("logL (complete-case):", round(ll_complete, 2), "\n\n")

cat("Done.\n")
