.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c(
  "holdout.R", "constrained_inad.R", "pipeline_inad.R",
  "extra_baselines.R", "marginal_moments.R", "cgfm_marginal.R"
), root)

ppc_plot_marginal_panel <- function(moments, metric, out_file) {
  groups <- sort(unique(moments$group))
  models <- setdiff(unique(moments$model), "Empirical bolus")
  colors <- grDevices::hcl.colors(max(length(models), 3L), "Dark 3")[seq_along(models)]
  names(colors) <- models

  grDevices::png(out_file, width = 1400, height = 700, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = c(1, length(groups)), mar = c(4, 4, 3, 1))
  for (group in groups) {
    dat <- moments[moments$group == group, ]
    ylim <- range(dat[[metric]], finite = TRUE)
    emp <- dat[dat$model == "Empirical bolus", ]
    graphics::plot(
      emp$time, emp[[metric]],
      type = "b", pch = 16, lwd = 2, col = "black",
      xlab = "Time", ylab = metric,
      ylim = ylim, main = paste("Group", group)
    )
    for (model in models) {
      cur <- dat[dat$model == model, ]
      graphics::lines(cur$time, cur[[metric]], col = colors[[model]], lwd = 1.6)
    }
  }
  graphics::legend(
    "topright",
    inset = c(-0.02, -0.08),
    legend = c("Empirical bolus", models),
    col = c("black", colors),
    lwd = c(2, rep(1.6, length(models))),
    pch = c(16, rep(NA, length(models))),
    cex = 0.65,
    xpd = TRUE
  )
}

run_marginal_moments_bolus <- function(fit_max_iter = 50L,
                                       tscount_sims = 2000L) {
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(bolus_inad$blocks)
  group_levels <- sort(unique(blocks))

  cfg <- ppc_inad_config("full")
  cfg$n_time <- ncol(y)
  cfg$n_subjects <- nrow(y)
  cfg$t_split <- 8L
  cfg$recursive_h <- 4L

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
  fit_nb_glm <- ppc_fit_marginal_nb_glm_full(y, blocks)
  fit_nb_glmm <- ppc_fit_marginal_nb_glmm_full(y, blocks)
  fit_tscount <- ppc_fit_tscount_baseline(y, blocks, cfg)

  moment_parts <- list(
    empirical = ppc_empirical_bolus_moments(y, blocks),
    inad_constrained = ppc_marginal_moments_inad_order1(
      fit_constrained, group_levels, "INAD constrained alpha"
    ),
    inad_unconstrained = ppc_marginal_moments_inad_order1(
      fit_unconstrained, group_levels, "INAD unconstrained"
    ),
    nb_glm = ppc_marginal_moments_nb_glm(fit_nb_glm, group_levels, ncol(y)),
    nb_glmm = ppc_marginal_moments_nb_glmm(fit_nb_glmm, group_levels, ncol(y)),
    tscount = ppc_marginal_moments_tscount_mc(
      fit_tscount, y, blocks, cfg, n_sims = tscount_sims
    )
  )
  cgfm_specs <- cgfm_bolus_published_specs()
  for (nm in names(cgfm_specs)) {
    moment_parts[[nm]] <- ppc_marginal_moments_cgfm_bolus(cgfm_specs[[nm]], group_levels)
  }
  moment_parts <- Filter(Negate(is.null), moment_parts)
  moments <- do.call(rbind, moment_parts)
  distances <- ppc_moment_l1_distance(moments)
  empirical_table3 <- ppc_empirical_table3_style(moments)
  by_time_tables <- list()
  relative_error_tables <- list()
  for (group in group_levels) {
    by_time_tables[[paste0("mean_group_", group)]] <- ppc_moment_by_time_table(
      moments, group = group, metric = "mean"
    )
    by_time_tables[[paste0("variance_group_", group)]] <- ppc_moment_by_time_table(
      moments, group = group, metric = "variance"
    )
    relative_error_tables[[paste0("mean_group_", group)]] <- ppc_moment_relative_error_table(
      moments, group = group, metric = "mean"
    )
    relative_error_tables[[paste0("variance_group_", group)]] <- ppc_moment_relative_error_table(
      moments, group = group, metric = "variance"
    )
  }

  fit_summary <- data.frame(
    model = c(
      "INAD unconstrained", "INAD constrained alpha",
      "NB GLM marginal", "NB GLMM marginal", "tscount tsglm MC",
      "CGFM independent frailty", "CGFM shared frailty"
    ),
    fit_source = c(
      "local antedep fit", "local antedep constrained-alpha fit",
      "local MASS::glm.nb", "local glmmTMB", "local tscount",
      "published Henderson-Shimakura Table 4",
      "published Henderson-Shimakura Table 4"
    ),
    note = c(
      NA, "nb_inno_size inherited by the POC constrained fit",
      "y ~ group + factor(time)", "random intercept; no lag",
      "panel-recursive MC from pooled tsglm fit",
      "independent frailty marginal moments",
      "shared frailty marginal moments"
    ),
    stringsAsFactors = FALSE
  )

  out_dir <- file.path(ppc_output_dir(root, "marginal_moments"))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(
    list(
      moments = moments,
      distances = distances,
      empirical_table3 = empirical_table3,
      by_time_tables = by_time_tables,
      relative_error_tables = relative_error_tables,
      fit_summary = fit_summary,
      fits = list(
        inad_unconstrained = fit_unconstrained,
        inad_constrained_alpha = fit_constrained,
        nb_glm = fit_nb_glm,
        nb_glmm = fit_nb_glmm,
        tscount = fit_tscount
      )
    ),
    file.path(out_dir, "marginal_moments_bolus.rds")
  )
  utils::write.csv(moments, file.path(out_dir, "marginal_moments_bolus.csv"), row.names = FALSE)
  utils::write.csv(distances, file.path(out_dir, "marginal_moments_l1.csv"), row.names = FALSE)
  utils::write.csv(
    empirical_table3,
    file.path(out_dir, "empirical_bolus_table3_style.csv"),
    row.names = FALSE
  )
  for (nm in names(by_time_tables)) {
    utils::write.csv(
      by_time_tables[[nm]],
      file.path(out_dir, paste0("marginal_", nm, "_by_time.csv")),
      row.names = FALSE
    )
  }
  for (nm in names(relative_error_tables)) {
    utils::write.csv(
      relative_error_tables[[nm]],
      file.path(out_dir, paste0("relative_error_", nm, "_by_time.csv")),
      row.names = FALSE
    )
  }
  utils::write.csv(fit_summary, file.path(out_dir, "marginal_moments_fit_sources.csv"), row.names = FALSE)
  ppc_plot_marginal_panel(moments, "mean", file.path(out_dir, "marginal_mean_bolus.png"))
  ppc_plot_marginal_panel(moments, "variance", file.path(out_dir, "marginal_variance_bolus.png"))
  list(
    out_dir = out_dir,
    empirical_table3 = empirical_table3,
    by_time_tables = by_time_tables,
    relative_error_tables = relative_error_tables,
    distances = distances,
    fit_summary = fit_summary
  )
}

if (identical(environment(), globalenv()) && !interactive()) {
  tscount_sims <- as.integer(Sys.getenv("PPC_TSCOUNT_MOMENT_SIMS", unset = "2000"))
  result <- run_marginal_moments_bolus(tscount_sims = tscount_sims)
  print(result$fit_summary)
  print(result$empirical_table3)
  print(result$distances)
  cat("Output:", result$out_dir, "\n")
}
