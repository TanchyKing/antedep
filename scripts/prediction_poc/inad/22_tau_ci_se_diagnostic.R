.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common("constrained_inad.R", root)

tau_ci_se_table <- function(ci) {
  tau <- ci$tau
  if (is.null(tau) || !nrow(tau)) return(data.frame())
  z <- stats::qnorm(1 - (1 - ci$level) / 2)
  tau$profile_midpoint <- (tau$lower + tau$upper) / 2
  tau$profile_width_se <- (tau$upper - tau$lower) / (2 * z)
  tau$normal_lower_from_se <- tau$est - z * tau$se
  tau$normal_upper_from_se <- tau$est + z * tau$se
  tau$ci_asymmetry <- (tau$upper - tau$est) - (tau$est - tau$lower)
  tau[, c(
    "param", "est", "se", "lower", "upper", "width",
    "profile_midpoint", "profile_width_se",
    "normal_lower_from_se", "normal_upper_from_se", "ci_asymmetry"
  )]
}

run_tau_ci_se_diagnostic <- function(profile_maxeval = 2500L,
                                     fit_max_iter = 50L) {
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(bolus_inad$blocks)

  fit_unconstrained <- antedep::fit_inad(
    y = y,
    order = 1L,
    thinning = "nbinom",
    innovation = "nbinom",
    blocks = blocks,
    max_iter = fit_max_iter
  )
  fit_constrained <- fit_inad_alpha_constant_poc(
    y = y,
    blocks = blocks,
    fit_unconstrained = fit_unconstrained,
    thinning = "nbinom",
    innovation = "nbinom",
    max_iter = fit_max_iter
  )

  ci <- antedep::ci_inad(
    y = y,
    fit = fit_constrained,
    blocks = blocks,
    profile_maxeval = profile_maxeval
  )
  tau_table <- tau_ci_se_table(ci)
  fit_table <- data.frame(
    fit = c("unconstrained", "constrained_alpha"),
    log_l = c(fit_unconstrained$log_l, fit_constrained$log_l),
    convergence = c(fit_unconstrained$convergence, fit_constrained$convergence),
    tau2 = c(fit_unconstrained$tau[2L], fit_constrained$tau[2L]),
    stringsAsFactors = FALSE
  )

  out_dir <- file.path(ppc_output_dir(root, "bolus_tau_ci_se_diagnostic"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      fit_unconstrained = fit_unconstrained,
      fit_constrained = fit_constrained,
      ci = ci,
      tau_table = tau_table,
      fit_table = fit_table
    ),
    file.path(out_dir, "tau_ci_se_diagnostic.rds")
  )
  utils::write.csv(tau_table, file.path(out_dir, "tau_ci_se_table.csv"), row.names = FALSE)
  utils::write.csv(fit_table, file.path(out_dir, "tau_ci_fit_table.csv"), row.names = FALSE)
  list(out_dir = out_dir, tau_table = tau_table, fit_table = fit_table)
}

if (identical(environment(), globalenv()) && !interactive()) {
  profile_maxeval <- as.integer(Sys.getenv("PPC_PROFILE_MAXEVAL", unset = "2500"))
  result <- run_tau_ci_se_diagnostic(profile_maxeval = profile_maxeval)
  print(result$fit_table)
  print(result$tau_table)
  cat("Output:", result$out_dir, "\n")
}
