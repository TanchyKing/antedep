.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(.ppc_script_dir, "09_extra_baselines_r200.R"))

if (identical(environment(), globalenv()) && !interactive()) {
  result <- run_inad_extra_baselines_r200(install_missing = TRUE)
  print(result$run_info)
  cat("Output:", result$out_dir, "\n")
}
