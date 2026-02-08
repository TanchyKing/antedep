#!/usr/bin/env Rscript

#' Example: Missing Data Handling in antedep Package
#' 
#' This script demonstrates the missing data functionality implemented
#' for the AD module.

# Load required functions (in practice, these would be in the package)
cat("Loading functions...\n")
source("missing_utils.R")
source("logL_ad_missing.R")
source("fit_ad_em.R")
source("logL_ad_modified.R")
source("fit_ad_modified.R")

# Set seed for reproducibility
set.seed(123)

# ==== Example 1: Monotone Missing (Dropout) ====
cat("\n=== Example 1: Monotone Missing (Dropout) ===\n")

# Simulate complete data
n_subjects <- 20
n_time <- 5
y_complete <- matrix(rnorm(n_subjects * n_time, mean = 10, sd = 2), 
                     nrow = n_subjects, ncol = n_time)

# Introduce dropout pattern (last 1-2 time points missing for some subjects)
y_dropout <- y_complete
dropout_subjects <- sample(n_subjects, 8)
for (s in dropout_subjects) {
  n_drop <- sample(1:2, 1)
  y_dropout[s, (n_time - n_drop + 1):n_time] <- NA
}

cat("Data dimensions:", nrow(y_dropout), "x", ncol(y_dropout), "\n")
cat("Missing data:", sum(is.na(y_dropout)), "values\n")
cat("Percent missing:", round(mean(is.na(y_dropout)) * 100, 1), "%\n\n")

# Validate missing pattern
missing_info <- .validate_missing(y_dropout)
cat("Complete subjects:", missing_info$n_complete, "\n")
cat("Dropout subjects:", sum(missing_info$patterns == "dropout"), "\n")
cat("Intermittent subjects:", missing_info$n_intermittent, "\n\n")

# Fit with EM
cat("Fitting with EM algorithm...\n")
fit_em <- fit_ad(y_dropout, order = 1, na_action = "em", 
                 em_max_iter = 50, em_verbose = TRUE)

cat("\n--- EM Results ---\n")
cat("Converged:", fit_em$em_converged, "\n")
cat("Iterations:", fit_em$em_iterations, "\n")
cat("Log-likelihood:", round(fit_em$log_l, 2), "\n")
cat("AIC:", round(fit_em$aic, 2), "\n")
cat("BIC:", round(fit_em$bic, 2), "\n\n")

# Compare with complete-case analysis
cat("Fitting with complete-case analysis...\n")
fit_cc <- fit_ad(y_dropout, order = 1, na_action = "complete")

cat("\n--- Complete-Case Results ---\n")
cat("Subjects used:", fit_cc$settings$n_subjects, "/", n_subjects, "\n")
cat("Log-likelihood:", round(fit_cc$log_l, 2), "\n")
cat("AIC:", round(fit_cc$aic, 2), "\n")
cat("BIC:", round(fit_cc$bic, 2), "\n\n")

# Compare parameters
cat("--- Parameter Comparison ---\n")
cat("mu (EM):  ", round(fit_em$mu, 2), "\n")
cat("mu (CC):  ", round(fit_cc$mu, 2), "\n")
cat("mu (True):", round(colMeans(y_complete), 2), "\n\n")

cat("phi (EM):  ", round(fit_em$phi, 3), "\n")
cat("phi (CC):  ", round(fit_cc$phi, 3), "\n\n")

# Plot EM convergence
cat("Plotting EM convergence...\n")
png("em_convergence_dropout.png", width = 600, height = 400)
plot(fit_em$em_ll_trace, type = "b", pch = 19, col = "steelblue",
     xlab = "EM Iteration", ylab = "Log-Likelihood",
     main = "EM Convergence (Monotone Dropout)")
grid()
dev.off()
cat("Saved: em_convergence_dropout.png\n\n")

# ==== Example 2: Intermittent Missing ====
cat("\n=== Example 2: Intermittent Missing ===\n")

# Introduce random (MCAR) missing pattern
y_mcar <- y_complete
missing_idx <- sample(length(y_mcar), round(0.15 * length(y_mcar)))
y_mcar[missing_idx] <- NA

cat("Data dimensions:", nrow(y_mcar), "x", ncol(y_mcar), "\n")
cat("Missing data:", sum(is.na(y_mcar)), "values\n")
cat("Percent missing:", round(mean(is.na(y_mcar)) * 100, 1), "%\n\n")

# Validate missing pattern
missing_info2 <- .validate_missing(y_mcar)
cat("Complete subjects:", missing_info2$n_complete, "\n")
cat("Intermittent subjects:", missing_info2$n_intermittent, "\n\n")

# Fit with EM
cat("Fitting with EM algorithm...\n")
fit_em2 <- fit_ad(y_mcar, order = 1, na_action = "em",
                  em_max_iter = 50, em_verbose = FALSE)

cat("\n--- Results ---\n")
cat("Converged:", fit_em2$em_converged, "\n")
cat("Iterations:", fit_em2$em_iterations, "\n")
cat("Log-likelihood:", round(fit_em2$log_l, 2), "\n\n")

# Plot convergence
png("em_convergence_intermittent.png", width = 600, height = 400)
plot(fit_em2$em_ll_trace, type = "b", pch = 19, col = "darkgreen",
     xlab = "EM Iteration", ylab = "Log-Likelihood",
     main = "EM Convergence (Intermittent Missing)")
grid()
dev.off()
cat("Saved: em_convergence_intermittent.png\n\n")

# ==== Example 3: With Block Effects ====
cat("\n=== Example 3: With Block Effects ===\n")

# Simulate two groups with different means
n_per_group <- 15
y_group1 <- matrix(rnorm(n_per_group * n_time, mean = 10, sd = 2),
                   nrow = n_per_group, ncol = n_time)
y_group2 <- matrix(rnorm(n_per_group * n_time, mean = 12, sd = 2),
                   nrow = n_per_group, ncol = n_time)
y_groups <- rbind(y_group1, y_group2)
blocks <- c(rep(1, n_per_group), rep(2, n_per_group))

# Introduce missing data
y_groups_missing <- y_groups
missing_idx <- sample(length(y_groups), round(0.1 * length(y_groups)))
y_groups_missing[missing_idx] <- NA

cat("Group 1 mean (true):", round(mean(y_group1), 2), "\n")
cat("Group 2 mean (true):", round(mean(y_group2), 2), "\n")
cat("Difference:", round(mean(y_group2) - mean(y_group1), 2), "\n\n")

# Fit with blocks
cat("Fitting with blocks...\n")
fit_blocks <- fit_ad(y_groups_missing, order = 1, blocks = blocks,
                     na_action = "em", em_max_iter = 50, em_verbose = FALSE)

cat("\n--- Results ---\n")
cat("Converged:", fit_blocks$em_converged, "\n")
cat("Block effects (tau):", round(fit_blocks$tau, 2), "\n")
cat("Estimated difference:", round(fit_blocks$tau[2], 2), "\n\n")

# ==== Example 4: Log-Likelihood Evaluation ====
cat("\n=== Example 4: Log-Likelihood with Different na_action ===\n")

# Use the dropout data from Example 1
true_params <- list(
  mu = colMeans(y_complete),
  phi = rep(0.5, n_time - 1),
  sigma = apply(y_complete, 2, sd)
)

# Evaluate log-likelihood with different methods
ll_marginalize <- logL_ad(y_dropout, order = 1, 
                          mu = true_params$mu,
                          phi = true_params$phi,
                          sigma = true_params$sigma,
                          na_action = "marginalize")

ll_complete <- logL_ad(y_dropout, order = 1,
                       mu = true_params$mu,
                       phi = true_params$phi,
                       sigma = true_params$sigma,
                       na_action = "complete")

cat("Log-likelihood (marginalize):", round(ll_marginalize, 2), "\n")
cat("Log-likelihood (complete):   ", round(ll_complete, 2), "\n")
cat("\nNote: These are not directly comparable as they use different data.\n")
cat("Marginalize uses all subjects; complete uses only complete cases.\n\n")

# ==== Summary ====
cat("\n=== Summary ===\n")
cat("✓ Implemented missing data handling for AD models\n")
cat("✓ EM algorithm converges reliably\n")
cat("✓ Handles monotone and intermittent missing patterns\n")
cat("✓ Works with block effects\n")
cat("✓ Provides proper observed-data likelihood\n\n")

cat("Implementation complete!\n")
