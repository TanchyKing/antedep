# Published Henderson and Shimakura bolus fits from Biometrika (2003), Table 4.
# Counts have Poisson means multiplied by gamma frailty with mean one and
# variance phi, so the marginal variance is mu + phi * mu^2.
cgfm_bolus_published_specs <- function() {
  list(
    cgfm_independent = list(
      label = "CGFM independent frailty",
      alpha = c(2.10, 1.61, 1.71, 1.78, 1.96, 1.66, 1.47, 1.37, 1.38, 1.55, 1.40, 1.35),
      beta_group2 = 0.34,
      phi = 0.49
    ),
    cgfm_shared = list(
      label = "CGFM shared frailty",
      alpha = c(2.09, 1.61, 1.72, 1.80, 1.97, 1.67, 1.49, 1.38, 1.37, 1.55, 1.42, 1.37),
      beta_group2 = 0.33,
      phi = 0.24
    )
  )
}

ppc_marginal_moments_cgfm_bolus <- function(spec, group_levels = c(1L, 2L)) {
  out <- vector("list", length(group_levels))
  for (ii in seq_along(group_levels)) {
    group_shift <- if (ii == 1L) 0 else spec$beta_group2
    mu <- exp(spec$alpha + group_shift)
    out[[ii]] <- ppc_moment_frame(
      model = spec$label,
      group = group_levels[[ii]],
      time = seq_along(mu),
      mean = mu,
      variance = mu + spec$phi * mu^2,
      source = "published Henderson-Shimakura Table 4"
    )
  }
  do.call(rbind, out)
}
