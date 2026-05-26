.ppc_script_dir <- local({
  file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]), winslash = "/")) else getwd()
})
source(file.path(dirname(.ppc_script_dir), "common", "config.R"))
root <- ppc_find_package_root(.ppc_script_dir)
ppc_load_antedep(root)
ppc_source_common("constrained_inad.R", root)

get_ns <- function(name) getFromNamespace(name, "antedep")

make_tau_profile_refitter <- function(y, fit, blocks,
                                      maxeval = 2500L,
                                      xtol_rel = 1e-6) {
  if (!requireNamespace("nloptr", quietly = TRUE)) {
    stop("nloptr is required for tau profile refits.", call. = FALSE)
  }

  if (is.null(fit$tau) || length(fit$tau) < 2L) {
    stop("fit$tau must exist with at least two blocks.", call. = FALSE)
  }

  n <- nrow(y)
  n_time <- ncol(y)
  n_blocks <- max(blocks)
  ord <- fit$settings$order
  thinning <- fit$settings$thinning
  innovation <- fit$settings$innovation

  eps_pos <- 1e-10
  eps_alpha <- 1e-10
  alpha_ub <- 1 - eps_alpha
  alpha_share_tol <- 1e-6

  alpha_profile_shared <- FALSE
  if (ord == 1L && n_time >= 2L && length(fit$alpha) == n_time) {
    a <- as.numeric(fit$alpha)
    if (all(is.finite(a[2:n_time])) && max(abs(a[2:n_time] - a[2])) <= alpha_share_tol) {
      alpha_profile_shared <- TRUE
    }
  }

  bell_mean_from_theta <- get_ns(".bell_mean_from_theta")
  bell_theta_from_mean <- get_ns(".bell_theta_from_mean")

  lambda_ok <- function(theta_vec, tau_vec) {
    mu <- if (innovation == "bell") {
      outer(rep(1, n), bell_mean_from_theta(theta_vec)) + matrix(tau_vec[blocks], n, n_time)
    } else {
      outer(rep(1, n), theta_vec) + matrix(tau_vec[blocks], n, n_time)
    }
    all(is.finite(mu)) && all(mu > eps_pos)
  }

  pack_par <- function(alpha, theta, tau, ord, n_blocks, n_time, tau_fix_idx) {
    out <- numeric(0)
    tau_free_idx <- 2:n_blocks
    tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
    if (length(tau_free_idx) > 0L) out <- c(out, tau[tau_free_idx])

    if (ord == 1L && n_time >= 2L) {
      if (alpha_profile_shared) {
        out <- c(out, alpha[2])
      } else {
        out <- c(out, alpha[2:n_time])
      }
    }
    if (ord == 2L) {
      if (n_time >= 2L) out <- c(out, alpha[2:n_time, 1])
      if (n_time >= 3L) out <- c(out, alpha[3:n_time, 2])
    }

    c(out, theta)
  }

  unpack_par <- function(par, alpha0, theta0, tau0, nb0,
                         ord, n_blocks, n_time, tau_fix_idx, tau_fix_val) {
    k <- 0L

    tau <- tau0
    tau[1] <- 0
    tau_free_idx <- 2:n_blocks
    tau_free_idx <- setdiff(tau_free_idx, tau_fix_idx)
    if (length(tau_free_idx) > 0L) {
      tau[tau_free_idx] <- par[(k + 1L):(k + length(tau_free_idx))]
      k <- k + length(tau_free_idx)
    }
    tau[tau_fix_idx] <- tau_fix_val
    tau[1] <- 0

    alpha <- alpha0
    if (ord == 1L && n_time >= 2L) {
      if (alpha_profile_shared) {
        alpha[2:n_time] <- par[k + 1L]
        k <- k + 1L
      } else {
        alpha[2:n_time] <- par[(k + 1L):(k + n_time - 1L)]
        k <- k + n_time - 1L
      }
      alpha[1] <- 0
      alpha <- pmin(pmax(alpha, 0), alpha_ub)
      alpha[1] <- 0
    }
    if (ord == 2L) {
      if (n_time >= 2L) {
        alpha[2:n_time, 1] <- par[(k + 1L):(k + n_time - 1L)]
        k <- k + n_time - 1L
      }
      if (n_time >= 3L) {
        alpha[3:n_time, 2] <- par[(k + 1L):(k + n_time - 2L)]
        k <- k + n_time - 2L
      }
      alpha <- pmin(pmax(alpha, 0), alpha_ub)
      alpha[1, ] <- 0
      if (n_time >= 2L) alpha[2, 2] <- 0
    }

    theta <- pmax(par[(k + 1L):(k + n_time)], eps_pos)
    list(alpha = alpha, theta = theta, tau = tau, nb_inno_size = nb0)
  }

  refit_tau_fixed <- function(tau_fix_idx, tau_fix_val) {
    alpha0 <- fit$alpha
    theta0 <- pmax(as.numeric(fit$theta), eps_pos)
    tau0 <- as.numeric(fit$tau)
    tau0[1] <- 0

    nb0 <- fit$nb_inno_size
    if (innovation == "nbinom") {
      if (is.null(nb0)) nb0 <- rep(1, n_time)
      nb0 <- as.numeric(nb0)
      if (length(nb0) == 1L) nb0 <- rep(nb0, n_time)
      if (length(nb0) != n_time) stop("fit$nb_inno_size must be length 1 or ncol(y).")
      nb0 <- pmax(nb0, eps_pos)
    } else {
      nb0 <- NULL
    }

    par0 <- pack_par(alpha0, theta0, tau0, ord, n_blocks, n_time, tau_fix_idx)

    obj <- function(par) {
      cur <- unpack_par(
        par = par,
        alpha0 = alpha0,
        theta0 = theta0,
        tau0 = tau0,
        nb0 = nb0,
        ord = ord,
        n_blocks = n_blocks,
        n_time = n_time,
        tau_fix_idx = tau_fix_idx,
        tau_fix_val = tau_fix_val
      )
      if (!lambda_ok(cur$theta, cur$tau)) return(1e8)

      ll <- tryCatch(
        antedep::logL_inad(
          y = y,
          order = ord,
          thinning = thinning,
          innovation = innovation,
          alpha = cur$alpha,
          theta = cur$theta,
          nb_inno_size = cur$nb_inno_size,
          blocks = blocks,
          tau = cur$tau
        ),
        error = function(e) NA_real_
      )
      if (!is.finite(ll)) return(1e8)
      -ll
    }

    opt <- nloptr::nloptr(
      x0 = par0,
      eval_f = obj,
      opts = list(
        algorithm = "NLOPT_LN_BOBYQA",
        xtol_rel = xtol_rel,
        maxeval = maxeval
      )
    )

    list(
      log_l = -opt$objective,
      convergence = opt$status,
      message = opt$message
    )
  }

  refit_tau_fixed
}

estimate_tau_profile_curvature <- function(y, fit, blocks,
                                           tau_index = 2L,
                                           deltas = c(-0.20, -0.15, -0.10, -0.075,
                                                      -0.05, 0, 0.05, 0.075,
                                                      0.10, 0.15, 0.20),
                                           maxeval = 2500L,
                                           xtol_rel = 1e-6,
                                           level = 0.95) {
  refit <- make_tau_profile_refitter(y, fit, blocks, maxeval = maxeval, xtol_rel = xtol_rel)
  tau_hat <- as.numeric(fit$tau[tau_index])
  ll_hat <- as.numeric(fit$log_l)

  prof <- do.call(rbind, lapply(deltas, function(d) {
    val <- tau_hat + d
    rr <- refit(tau_fix_idx = tau_index, tau_fix_val = val)
    data.frame(
      tau_index = tau_index,
      tau_value = val,
      delta = d,
      log_l_profile = rr$log_l,
      profile_deviance = 2 * (ll_hat - rr$log_l),
      convergence = rr$convergence,
      message = rr$message,
      stringsAsFactors = FALSE
    )
  }))

  ok <- is.finite(prof$log_l_profile) & is.finite(prof$delta)
  local_prof <- prof[ok, , drop = FALSE]
  if (nrow(local_prof) < 4L) {
    stop("Too few finite profile points to estimate curvature.", call. = FALSE)
  }

  # Local profile curvature: if l_p(tau) is locally quadratic,
  # l_p(tau_hat + d) ~= a + b d + c d^2, with c < 0 and
  # Var(tau_hat) ~= -1 / (2c). The linear term is kept because the
  # constrained-alpha POC object is not a fully joint fit, so the stored
  # tau may be slightly off the local profiled maximum after nuisance refits.
  fit_quad_logl <- stats::lm(log_l_profile ~ delta + I(delta^2), data = local_prof)
  coef_quad <- stats::coef(fit_quad_logl)
  slope <- unname(coef_quad[["delta"]])
  curve <- unname(coef_quad[["I(delta^2)"]])
  if (!is.finite(curve) || curve >= 0) {
    stop("Profile log-likelihood quadratic curvature is not negative.", call. = FALSE)
  }
  vertex_delta <- -slope / (2 * curve)
  est_curvature <- tau_hat + vertex_delta
  se_curvature <- sqrt(-1 / (2 * curve))
  z <- stats::qnorm(1 - (1 - level) / 2)

  ci <- antedep::ci_inad(
    y = y,
    fit = fit,
    blocks = blocks,
    level = level,
    profile_maxeval = maxeval,
    profile_xtol_rel = xtol_rel
  )
  tau_ci <- ci$tau[ci$tau$param == paste0("tau[", tau_index, "]"), , drop = FALSE]

  summary <- data.frame(
    param = paste0("tau[", tau_index, "]"),
    est_fit = tau_hat,
    est_profile_curvature = est_curvature,
    se_profile_curvature = se_curvature,
    lower_profile_curvature = est_curvature - z * se_curvature,
    upper_profile_curvature = est_curvature + z * se_curvature,
    se_profile_ci_width = tau_ci$se,
    lower_profile_ci = tau_ci$lower,
    upper_profile_ci = tau_ci$upper,
    width_profile_ci = tau_ci$width,
    level = level,
    n_profile_points = nrow(local_prof),
    curvature_delta_from_fit = vertex_delta,
    stringsAsFactors = FALSE
  )

  list(summary = summary, profile_points = prof)
}

run_tau_profile_curvature <- function(maxeval = 2500L,
                                      fit_max_iter = 50L) {
  data("bolus_inad", package = "antedep")
  y <- bolus_inad$y
  blocks <- as.integer(bolus_inad$blocks)

  diag_path <- file.path(
    ppc_output_dir(root, "bolus_tau_ci_se_diagnostic"),
    "tau_ci_se_diagnostic.rds"
  )

  if (file.exists(diag_path)) {
    diag <- readRDS(diag_path)
    fit_constrained <- diag$fit_constrained
  } else {
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
  }

  res <- estimate_tau_profile_curvature(
    y = y,
    fit = fit_constrained,
    blocks = blocks,
    tau_index = 2L,
    maxeval = maxeval
  )

  out_dir <- ppc_output_dir(root, "bolus_tau_profile_curvature")
  utils::write.csv(res$summary, file.path(out_dir, "tau_profile_curvature_summary.csv"), row.names = FALSE)
  utils::write.csv(res$profile_points, file.path(out_dir, "tau_profile_curvature_points.csv"), row.names = FALSE)
  saveRDS(res, file.path(out_dir, "tau_profile_curvature.rds"))

  list(out_dir = out_dir, summary = res$summary, profile_points = res$profile_points)
}

if (identical(environment(), globalenv()) && !interactive()) {
  maxeval <- as.integer(Sys.getenv("PPC_PROFILE_MAXEVAL", unset = "2500"))
  result <- run_tau_profile_curvature(maxeval = maxeval)
  print(result$summary)
  print(result$profile_points)
  cat("Output:", result$out_dir, "\n")
}
