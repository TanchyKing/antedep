ppc_find_package_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    desc <- file.path(path, "DESCRIPTION")
    if (file.exists(desc)) {
      first <- readLines(desc, n = 1L, warn = FALSE)
      if (identical(first, "Package: antedep")) return(path)
    }
    parent <- dirname(path)
    if (identical(parent, path)) stop("Could not find antedep package root")
    path <- parent
  }
}

ppc_poc_dir <- function(root = ppc_find_package_root()) {
  file.path(root, "scripts", "prediction_poc")
}

ppc_output_dir <- function(root = ppc_find_package_root(), mode = "smoke") {
  out <- file.path(ppc_poc_dir(root), "output", "inad", mode)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  out
}

ppc_load_antedep <- function(root = ppc_find_package_root()) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(root, quiet = TRUE)
  } else {
    library(antedep)
  }
  invisible(TRUE)
}

ppc_source_common <- function(files, root = ppc_find_package_root()) {
  for (file in files) {
    source(file.path(ppc_poc_dir(root), "common", file), local = FALSE)
  }
  invisible(TRUE)
}

ppc_inad_config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  n_cores <- parallel::detectCores(logical = TRUE)
  n_workers <- max(1L, floor(2 * n_cores / 3))

  cfg <- list(
    mode = mode,
    seed = 20260518L,
    n_reps = if (mode == "smoke") 5L else 1000L,
    n_sims = if (mode == "smoke") 100L else 1000L,
    n_per_group = if (mode == "smoke") 20L else 50L,
    n_time = 12L,
    train_prop = 0.70,
    t_split = 8L,
    recursive_h = 4L,
    thinning = "nbinom",
    innovation = "nbinom",
    order = 1L,
    alpha = 0.35,
    theta = c(2.2, 2.4, 2.8, 3.2, 3.6, 3.3, 3.0, 2.7, 2.4, 2.1, 1.9, 1.7),
    nb_inno_size = c(1.4, 1.4, 1.6, 1.8, 2.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.5),
    tau = c(0, 1.25),
    fit_max_iter = if (mode == "smoke") 8L else 50L,
    n_cores = n_cores,
    n_workers = n_workers
  )
  cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
  cfg$n_subjects <- length(cfg$blocks)
  cfg
}

ppc_gau_config <- function(mode = c("smoke", "full")) {
  mode <- match.arg(mode)
  n_cores <- parallel::detectCores(logical = TRUE)
  n_workers <- max(1L, floor(2 * n_cores / 3))

  cfg <- list(
    mode = mode,
    seed = 20260520L,
    n_reps = if (mode == "smoke") 5L else 500L,
    n_per_group = if (mode == "smoke") 20L else 60L,
    n_time = 8L,
    train_prop = 0.70,
    t_split = 5L,
    recursive_h = 3L,
    order = 1L,
    mu = c(0.0, 0.4, 0.8, 1.2, 1.6, 1.5, 1.3, 1.1),
    phi = c(0, rep(0.55, 7L)),
    sigma = c(0.9, 0.9, 1.0, 1.0, 1.1, 1.1, 1.0, 1.0),
    tau = c(0, 1.0),
    n_cores = n_cores,
    n_workers = n_workers
  )
  cfg$blocks <- rep(seq_along(cfg$tau), each = cfg$n_per_group)
  cfg$n_subjects <- length(cfg$blocks)
  cfg
}

ppc_rep_seed <- function(cfg, rep_id, stream = 0L) {
  as.integer(cfg$seed + 10000L * as.integer(rep_id) + as.integer(stream))
}
