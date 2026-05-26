.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})

source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common(c("holdout.R", "sim_gau_forward.R", "pipeline_gau.R"), root)

result <- ppc_run_gau_smoke()
print(result$run_info)
print(result$summaries$decision)
print(summary(result$recovery))
print(result$validation)
