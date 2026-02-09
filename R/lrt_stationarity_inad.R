#' Likelihood Ratio Test for INAD Stationarity
#'
#' Tests whether time-varying parameters can be constrained to constants.
#'
#' @param y Integer matrix with n_subjects rows and n_time columns.
#'
#' @importFrom stats pchisq optim
#' @param order Model order (1 or 2).
#' @param thinning Thinning operator: "binom", "pois", or "nbinom".
#' @param innovation Innovation distribution: "pois", "bell", or "nbinom".
#' @param blocks Optional integer vector for block effects.
#' @param constrain Which parameters to constrain: "alpha", "theta", "both" for order 1;
#'   "alpha1", "alpha2", "alpha", "theta", "all" for order 2.
#' @param fit_unconstrained Optional pre-computed unconstrained fit.
#' @param verbose Logical; if TRUE, print progress.
#' @param ... Additional arguments.
#'
#' @return A list with class "lrt_stationarity_inad".
#' @export
lrt_stationarity_inad <- function(y, order = 1, thinning = "binom", innovation = "pois",
                                  blocks = NULL, constrain = "both",
                                  fit_unconstrained = NULL, verbose = FALSE, ...) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (anyNA(y)) {
    stop(
      "lrt_stationarity_inad currently supports complete data only. Missing-data INAD likelihood-ratio tests are not implemented yet.",
      call. = FALSE
    )
  }
  if (any(y < 0) || any(y != floor(y))) stop("y must contain non-negative integers")
  n_subjects <- nrow(y); n_time <- ncol(y)
  if (!order %in% c(1, 2)) stop("order must be 1 or 2 for stationarity testing")
  if (order >= n_time) stop("order must be less than n_time")
  thinning <- match.arg(thinning, c("binom", "pois", "nbinom"))
  innovation <- match.arg(innovation, c("pois", "bell", "nbinom"))
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  n_blocks <- length(unique(blocks))
  
  constraint_info <- .parse_constraint(constrain, order, n_time)
  
  if (is.null(fit_unconstrained)) {
    if (verbose) cat("Fitting unconstrained model...\n")
    fit_unconstrained <- fit_inad(y, order = order, thinning = thinning, innovation = innovation, blocks = blocks, ...)
  }
  if (verbose) cat("Fitting constrained model...\n")
  fit_constrained <- .fit_constrained(y, order, thinning, innovation, blocks, constraint_info, fit_unconstrained, ...)
  
  logL_uncon <- fit_unconstrained$log_l; logL_con <- fit_constrained$log_l
  lrt_stat <- 2 * (logL_uncon - logL_con)
  if (lrt_stat < 0) { warning("LRT negative, setting to 0"); lrt_stat <- 0 }
  df <- constraint_info$df
  p_value <- 1 - pchisq(lrt_stat, df = df)
  
  bic_uncon <- bic_inad(fit_unconstrained, n_subjects)
  n_params_con <- .count_params_inad(order, n_time, n_blocks, thinning, innovation) - df
  bic_con <- -2 * logL_con + n_params_con * log(n_subjects)
  bic_selected <- if (bic_uncon <= bic_con) "unconstrained" else "constrained"
  
  table <- data.frame(model = c("Unconstrained", paste("Constrained:", constraint_info$description)),
                      logLik = c(logL_uncon, logL_con),
                      n_params = c(.count_params_inad(order, n_time, n_blocks, thinning, innovation), n_params_con),
                      BIC = c(bic_uncon, bic_con), stringsAsFactors = FALSE)
  
  result <- list(fit_unconstrained = fit_unconstrained, fit_constrained = fit_constrained,
                 constraint = constraint_info$description, lrt_stat = lrt_stat, df = df, p_value = p_value,
                 bic_unconstrained = bic_uncon, bic_constrained = bic_con, bic_selected = bic_selected,
                 estimates = list(unconstrained = list(alpha = fit_unconstrained$alpha, theta = fit_unconstrained$theta, tau = fit_unconstrained$tau),
                                  constrained = list(alpha = fit_constrained$alpha, theta = fit_constrained$theta, tau = fit_constrained$tau)),
                 table = table, settings = list(order = order, thinning = thinning, innovation = innovation,
                                                n_subjects = n_subjects, n_time = n_time, n_blocks = n_blocks, constrain = constrain))
  class(result) <- "lrt_stationarity_inad"
  result
}

#' @export
print.lrt_stationarity_inad <- function(x, digits = 4, ...) {
  cat("\nLikelihood Ratio Test for INAD Stationarity\n============================================\n\n")
  cat("Constraint: H0:", x$constraint, "\n\n")
  print(x$table, row.names = FALSE, digits = digits)
  cat("\nLRT:", round(x$lrt_stat, digits), " df:", x$df, " p-value:", format.pval(x$p_value, digits = digits))
  cat("\nBIC selects:", x$bic_selected, "\n")
  invisible(x)
}

#' @keywords internal
.parse_constraint <- function(constrain, order, n_time) {
  constrain <- tolower(constrain)
  if (order == 1) {
    if (length(constrain) == 1) {
      if (constrain %in% c("both", "all")) { alpha_const <- TRUE; theta_const <- TRUE }
      else if (constrain == "alpha") { alpha_const <- TRUE; theta_const <- FALSE }
      else if (constrain == "theta") { alpha_const <- FALSE; theta_const <- TRUE }
      else stop("Invalid constraint for order 1. Use: 'alpha', 'theta', or 'both'")
    } else { alpha_const <- "alpha" %in% constrain; theta_const <- "theta" %in% constrain }
    df <- (if (alpha_const) n_time - 2 else 0) + (if (theta_const) n_time - 1 else 0)
    desc <- paste(c(if(alpha_const) "alpha constant", if(theta_const) "theta constant"), collapse = ", ")
    list(order = 1, alpha_const = alpha_const, theta_const = theta_const, df = df, description = desc)
  } else {
    alpha1_const <- alpha2_const <- theta_const <- FALSE
    if (length(constrain) == 1) {
      if (constrain %in% c("both", "all")) { alpha1_const <- alpha2_const <- theta_const <- TRUE }
      else if (constrain == "alpha") { alpha1_const <- alpha2_const <- TRUE }
      else if (constrain == "alpha1") alpha1_const <- TRUE
      else if (constrain == "alpha2") alpha2_const <- TRUE
      else if (constrain == "theta") theta_const <- TRUE
      else stop("Invalid constraint for order 2")
    } else {
      alpha1_const <- "alpha1" %in% constrain; alpha2_const <- "alpha2" %in% constrain
      if ("alpha" %in% constrain) { alpha1_const <- alpha2_const <- TRUE }
      theta_const <- "theta" %in% constrain
    }
    df <- (if (alpha1_const) n_time - 3 else 0) + (if (alpha2_const) n_time - 3 else 0) + (if (theta_const) n_time - 1 else 0)
    desc <- paste(c(if(alpha1_const) "alpha1 constant", if(alpha2_const) "alpha2 constant", if(theta_const) "theta constant"), collapse = ", ")
    list(order = 2, alpha1_const = alpha1_const, alpha2_const = alpha2_const, theta_const = theta_const, df = df, description = desc)
  }
}

#' @keywords internal
.fit_constrained <- function(y, order, thinning, innovation, blocks, cinfo, init_from, ...) {
  n <- nrow(y); N <- ncol(y); B <- length(unique(blocks))
  
  if (order == 1) {
    obj_fn <- function(par) {
      idx <- 1
      if (cinfo$alpha_const) { alpha <- c(0, rep(par[idx], N - 1)); idx <- idx + 1 }
      else { alpha <- c(0, par[idx:(idx + N - 2)]); idx <- idx + N - 1 }
      if (cinfo$theta_const) { theta <- rep(par[idx], N); idx <- idx + 1 }
      else { theta <- par[idx:(idx + N - 1)]; idx <- idx + N }
      tau <- if (B > 1) c(0, par[idx:(idx + B - 2)]) else 0
      alpha <- pmax(alpha, 0); if (thinning == "binom") alpha <- pmin(alpha, 0.999)
      theta <- pmax(theta, 1e-8)
      lam_mat <- outer(rep(1, n), theta) + matrix(tau[blocks], n, N)
      if (any(lam_mat <= 0)) return(1e20)
      -logL_inad(y = y, order = 1, thinning = thinning, innovation = innovation,
                 alpha = alpha, theta = theta, blocks = blocks, tau = tau)
    }
    alpha_init <- if (!is.null(init_from$alpha)) init_from$alpha else c(0, rep(0.3, N - 1))
    theta_init <- if (!is.null(init_from$theta)) init_from$theta else colMeans(y)
    tau_init <- if (!is.null(init_from$tau)) init_from$tau else rep(0, B)
    par0 <- c(); lower <- c(); upper <- c()
    if (cinfo$alpha_const) { par0 <- c(par0, mean(alpha_init[-1])); lower <- c(lower, 0); upper <- c(upper, if(thinning == "binom") 0.999 else 10) }
    else { par0 <- c(par0, alpha_init[-1]); lower <- c(lower, rep(0, N-1)); upper <- c(upper, rep(if(thinning == "binom") 0.999 else 10, N-1)) }
    if (cinfo$theta_const) { par0 <- c(par0, mean(theta_init)); lower <- c(lower, 1e-8); upper <- c(upper, 1e6) }
    else { par0 <- c(par0, theta_init); lower <- c(lower, rep(1e-8, N)); upper <- c(upper, rep(1e6, N)) }
    if (B > 1) { par0 <- c(par0, tau_init[-1]); lower <- c(lower, rep(-1e3, B-1)); upper <- c(upper, rep(1e3, B-1)) }
  } else {
    obj_fn <- function(par) {
      idx <- 1
      if (cinfo$alpha1_const) { alpha1 <- c(0, 0, rep(par[idx], N - 2)); idx <- idx + 1 }
      else { alpha1 <- c(0, 0, par[idx:(idx + N - 3)]); idx <- idx + N - 2 }
      if (cinfo$alpha2_const) { alpha2 <- c(0, 0, rep(par[idx], N - 2)); idx <- idx + 1 }
      else { alpha2 <- c(0, 0, par[idx:(idx + N - 3)]); idx <- idx + N - 2 }
      if (cinfo$theta_const) { theta <- rep(par[idx], N); idx <- idx + 1 }
      else { theta <- par[idx:(idx + N - 1)]; idx <- idx + N }
      tau <- if (B > 1) c(0, par[idx:(idx + B - 2)]) else 0
      alpha1 <- pmax(alpha1, 0); alpha2 <- pmax(alpha2, 0)
      if (thinning == "binom") { alpha1 <- pmin(alpha1, 0.999); alpha2 <- pmin(alpha2, 0.999) }
      theta <- pmax(theta, 1e-8)
      lam_mat <- outer(rep(1, n), theta) + matrix(tau[blocks], n, N)
      if (any(lam_mat <= 0)) return(1e20)
      -logL_inad(y = y, order = 2, thinning = thinning, innovation = innovation,
                 alpha = cbind(alpha1, alpha2), theta = theta, blocks = blocks, tau = tau)
    }
    if (!is.null(init_from$alpha)) { au <- .unpack_alpha(2, init_from$alpha, N); alpha1_init <- au$a1; alpha2_init <- au$a2 }
    else { alpha1_init <- c(0, 0, rep(0.3, N-2)); alpha2_init <- c(0, 0, rep(0.1, N-2)) }
    theta_init <- if (!is.null(init_from$theta)) init_from$theta else colMeans(y)
    tau_init <- if (!is.null(init_from$tau)) init_from$tau else rep(0, B)
    par0 <- c(); lower <- c(); upper <- c()
    if (cinfo$alpha1_const) { par0 <- c(par0, mean(alpha1_init[-(1:2)])); lower <- c(lower, 0); upper <- c(upper, if(thinning == "binom") 0.999 else 10) }
    else { par0 <- c(par0, alpha1_init[-(1:2)]); lower <- c(lower, rep(0, N-2)); upper <- c(upper, rep(if(thinning == "binom") 0.999 else 10, N-2)) }
    if (cinfo$alpha2_const) { par0 <- c(par0, mean(alpha2_init[-(1:2)])); lower <- c(lower, 0); upper <- c(upper, if(thinning == "binom") 0.999 else 10) }
    else { par0 <- c(par0, alpha2_init[-(1:2)]); lower <- c(lower, rep(0, N-2)); upper <- c(upper, rep(if(thinning == "binom") 0.999 else 10, N-2)) }
    if (cinfo$theta_const) { par0 <- c(par0, mean(theta_init)); lower <- c(lower, 1e-8); upper <- c(upper, 1e6) }
    else { par0 <- c(par0, theta_init); lower <- c(lower, rep(1e-8, N)); upper <- c(upper, rep(1e6, N)) }
    if (B > 1) { par0 <- c(par0, tau_init[-1]); lower <- c(lower, rep(-1e3, B-1)); upper <- c(upper, rep(1e3, B-1)) }
  }
  
  opt <- tryCatch(nloptr::nloptr(x0 = par0, eval_f = obj_fn, lb = lower, ub = upper,
                                  opts = list(algorithm = "NLOPT_LN_BOBYQA", xtol_rel = 1e-7, maxeval = 4000)),
                  error = function(e) optim(par0, obj_fn, method = "L-BFGS-B", lower = lower, upper = upper))
  sol <- if (inherits(opt, "nloptr")) opt$solution else opt$par
  log_l <- if (inherits(opt, "nloptr")) -opt$objective else -opt$value
  
  if (order == 1) {
    idx <- 1
    if (cinfo$alpha_const) { alpha <- c(0, rep(sol[idx], N - 1)); idx <- idx + 1 }
    else { alpha <- c(0, sol[idx:(idx + N - 2)]); idx <- idx + N - 1 }
    if (cinfo$theta_const) { theta <- rep(sol[idx], N); idx <- idx + 1 }
    else { theta <- sol[idx:(idx + N - 1)]; idx <- idx + N }
    tau <- if (B > 1) c(0, sol[idx:(idx + B - 2)]) else 0
    list(alpha = pmax(alpha, 0), theta = pmax(theta, 1e-8), tau = tau, log_l = log_l)
  } else {
    idx <- 1
    if (cinfo$alpha1_const) { alpha1 <- c(0, 0, rep(sol[idx], N - 2)); idx <- idx + 1 }
    else { alpha1 <- c(0, 0, sol[idx:(idx + N - 3)]); idx <- idx + N - 2 }
    if (cinfo$alpha2_const) { alpha2 <- c(0, 0, rep(sol[idx], N - 2)); idx <- idx + 1 }
    else { alpha2 <- c(0, 0, sol[idx:(idx + N - 3)]); idx <- idx + N - 2 }
    if (cinfo$theta_const) { theta <- rep(sol[idx], N); idx <- idx + 1 }
    else { theta <- sol[idx:(idx + N - 1)]; idx <- idx + N }
    tau <- if (B > 1) c(0, sol[idx:(idx + B - 2)]) else 0
    list(alpha = cbind(pmax(alpha1, 0), pmax(alpha2, 0)), theta = pmax(theta, 1e-8), tau = tau, log_l = log_l)
  }
}

#' Run All Stationarity Tests
#' @param y Integer matrix.
#' @param order Model order (1 or 2).
#' @param thinning Thinning operator.
#' @param innovation Innovation distribution.
#' @param blocks Optional block assignments.
#' @param verbose Logical.
#' @param ... Additional arguments.
#' @return A list with class "stationarity_tests_inad".
#' @export
run_stationarity_tests_inad <- function(y, order = 1, thinning = "binom", innovation = "pois",
                                        blocks = NULL, verbose = TRUE, ...) {
  if (!is.matrix(y)) y <- as.matrix(y)
  if (anyNA(y)) {
    stop(
      "run_stationarity_tests_inad currently supports complete data only. Missing-data INAD likelihood-ratio tests are not implemented yet.",
      call. = FALSE
    )
  }
  n_subjects <- nrow(y)
  if (is.null(blocks)) blocks <- rep(1L, n_subjects)
  if (verbose) cat("Fitting unconstrained model...\n")
  fit_uncon <- fit_inad(y, order = order, thinning = thinning, innovation = innovation, blocks = blocks, ...)
  tests <- if (order == 1) c("alpha", "theta", "both") else c("alpha1", "alpha2", "alpha", "theta", "all")
  results <- list()
  for (test in tests) {
    if (verbose) cat("Testing:", test, "constant...\n")
    results[[test]] <- lrt_stationarity_inad(y, order, thinning, innovation, blocks, test, fit_uncon, FALSE, ...)
  }
  summary_table <- data.frame(constraint = tests, df = sapply(results, function(r) r$df),
                              LRT = sapply(results, function(r) r$lrt_stat),
                              p_value = sapply(results, function(r) r$p_value),
                              BIC_selected = sapply(results, function(r) r$bic_selected), stringsAsFactors = FALSE)
  output <- list(fit_unconstrained = fit_uncon, tests = results, summary = summary_table,
                 settings = list(order = order, thinning = thinning, innovation = innovation,
                                 n_subjects = n_subjects, n_time = ncol(y)))
  class(output) <- "stationarity_tests_inad"
  output
}

#' @export
print.stationarity_tests_inad <- function(x, digits = 4, ...) {
  cat("\nStationarity Tests for INAD Model\n==================================\n\n")
  print(x$summary, row.names = FALSE, digits = digits)
  invisible(x)
}
