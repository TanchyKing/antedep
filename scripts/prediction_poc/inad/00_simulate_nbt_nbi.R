.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)

simulate_one_nbt_nbi_inadfe <- function(mode = c("smoke", "full"), rep_id = 1L) {
  mode <- match.arg(mode)
  cfg <- ppc_inad_config(mode)
  ppc_source_common(c("pipeline_inad.R"), root)
  ppc_simulate_inad_dgp(cfg, rep_id)
}

if (identical(environment(), globalenv()) && !interactive()) {
  args <- commandArgs(trailingOnly = TRUE)
  mode <- if (length(args)) args[[1L]] else "smoke"
  sim <- simulate_one_nbt_nbi_inadfe(mode)
  str(sim)
}
