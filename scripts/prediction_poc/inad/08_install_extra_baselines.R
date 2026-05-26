.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_source_common("extra_baselines.R", root)

install_extra_prediction_baselines <- function() {
  before <- ppc_extra_package_status()
  print(before)
  after <- ppc_install_extra_baseline_packages()
  print(after)
  invisible(after)
}

if (identical(environment(), globalenv()) && !interactive()) {
  install_extra_prediction_baselines()
}
