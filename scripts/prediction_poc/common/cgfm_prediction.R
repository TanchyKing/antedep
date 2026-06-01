# Henderson-Shimakura common gamma frailty model (CGFM) prediction helpers.
#
# These helpers are intentionally script-side, not package-side.  They support
# bolus leave-two-out comparisons against INAD without making CGFM part of the
# antedep public model scope.

ppc_cgfm_compile_cpp <- local({
  compiled <- FALSE
  function() {
    if (compiled || !requireNamespace("Rcpp", quietly = TRUE) || !nzchar(Sys.which("make"))) {
      return(compiled)
    }
    Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static double ppc_log_pois_gamma_pair_recur(
    int n1, int n2, double u1, double u2, double theta, double r) {
  if (r <= 1e-12) {
    double size = 1.0 / theta;
    double lp1 = R::lgammafn(n1 + size) - R::lgammafn(size) - R::lgammafn(n1 + 1.0) +
      size * std::log(size) + n1 * std::log(u1) - (n1 + size) * std::log(size + u1);
    double lp2 = R::lgammafn(n2 + size) - R::lgammafn(size) - R::lgammafn(n2 + 1.0) +
      size * std::log(size) + n2 * std::log(u2) - (n2 + size) * std::log(size + u2);
    return lp1 + lp2;
  }
  if (r >= 1.0 - 1e-12) {
    int n = n1 + n2;
    double u = u1 + u2;
    double size = 1.0 / theta;
    return n1 * std::log(u1) + n2 * std::log(u2) -
      R::lgammafn(n1 + 1.0) - R::lgammafn(n2 + 1.0) +
      R::lgammafn(n + size) - R::lgammafn(size) + size * std::log(size) -
      (n + size) * std::log(size + u);
  }

  double a = 1.0 / theta;
  double A = 1.0 + theta * u1 + theta * u2 + theta * theta * u1 * u2 * (1.0 - r);
  double B = theta * u1 * (1.0 + theta * u2 * (1.0 - r));
  double C = theta * u2 * (1.0 + theta * u1 * (1.0 - r));
  double E = theta * theta * u1 * u2 * (1.0 - r);
  int nr = n1 + 1;
  int nc = n2 + 1;
  std::vector<double> g(nr * nc, 0.0);
  auto at = [&](int i, int j) -> double& { return g[i * nc + j]; };

  at(0, 0) = std::exp(-a * std::log(A));
  for (int n = 1; n <= n1; ++n) {
    at(n, 0) = at(n - 1, 0) * B * (n - 1 + a) / (A * n);
  }
  for (int m = 1; m <= n2; ++m) {
    at(0, m) = at(0, m - 1) * C * (m - 1 + a) / (A * m);
  }
  for (int n = 1; n <= n1; ++n) {
    for (int m = 1; m <= n2; ++m) {
      at(n, m) = (B * (n - 1 + a) * at(n - 1, m) +
                  C * n * at(n, m - 1) -
                  E * (n - 1 + a) * at(n - 1, m - 1)) / (A * n);
    }
  }
  double val = at(n1, n2);
  if (!R_FINITE(val) || val <= 0.0) return R_NegInf;
  return std::log(val);
}

// [[Rcpp::export]]
double ppc_cgfm_composite_loglik_recur_cpp(NumericMatrix y, NumericVector x,
                                           NumericVector alpha, double beta,
                                           double theta, double rho) {
  int s = y.nrow();
  int p = y.ncol();
  double total = 0.0;
  for (int j = 0; j < p - 1; ++j) {
    for (int k = j + 1; k < p; ++k) {
      double lag_r = std::pow(rho, k - j);
      for (int i = 0; i < s; ++i) {
        double u1 = std::exp(alpha[j] + beta * x[i]);
        double u2 = std::exp(alpha[k] + beta * x[i]);
        double lp = ppc_log_pois_gamma_pair_recur((int)y(i, j), (int)y(i, k),
                                                  u1, u2, theta, lag_r);
        if (!R_FINITE(lp)) return R_NegInf;
        total += lp;
      }
    }
  }
  return total / (p - 1);
}
',
      env = .GlobalEnv,
      rebuild = TRUE
    )
    compiled <<- exists("ppc_cgfm_composite_loglik_recur_cpp", mode = "function", envir = .GlobalEnv)
    compiled
  }
})

ppc_cgfm_x <- function(blocks) {
  as.integer(blocks == sort(unique(blocks))[2L])
}

ppc_cgfm_nb_logpmf <- function(n, mu, theta) {
  size <- 1 / theta
  stats::dnbinom(n, size = size, mu = mu, log = TRUE)
}

ppc_cgfm_biv_logpmf_table <- function(max1, max2, u1, u2, theta, r) {
  if (r <= 1e-12) {
    lp1 <- ppc_cgfm_nb_logpmf(0:max1, u1, theta)
    lp2 <- ppc_cgfm_nb_logpmf(0:max2, u2, theta)
    return(outer(lp1, lp2, `+`))
  }
  if (r >= 1 - 1e-12) {
    out <- matrix(NA_real_, max1 + 1L, max2 + 1L)
    size <- 1 / theta
    u <- u1 + u2
    for (n1 in 0:max1) {
      for (n2 in 0:max2) {
        n <- n1 + n2
        out[n1 + 1L, n2 + 1L] <-
          n1 * log(u1) + n2 * log(u2) - lgamma(n1 + 1) - lgamma(n2 + 1) +
          lgamma(n + size) - lgamma(size) + size * log(size) -
          (n + size) * log(size + u)
      }
    }
    return(out)
  }

  a <- 1 / theta
  A <- 1 + theta * u1 + theta * u2 + theta^2 * u1 * u2 * (1 - r)
  B <- theta * u1 * (1 + theta * u2 * (1 - r))
  C <- theta * u2 * (1 + theta * u1 * (1 - r))
  E <- theta^2 * u1 * u2 * (1 - r)

  g <- matrix(0, nrow = max1 + 1L, ncol = max2 + 1L)
  g[1L, 1L] <- exp(-a * log(A))
  if (max1 > 0L) {
    for (n in 1:max1) g[n + 1L, 1L] <- g[n, 1L] * B * (n - 1 + a) / (A * n)
  }
  if (max2 > 0L) {
    for (m in 1:max2) g[1L, m + 1L] <- g[1L, m] * C * (m - 1 + a) / (A * m)
  }
  if (max1 > 0L && max2 > 0L) {
    for (n in 1:max1) {
      for (m in 1:max2) {
        g[n + 1L, m + 1L] <-
          (B * (n - 1 + a) * g[n, m + 1L] +
             C * n * g[n + 1L, m] -
             E * (n - 1 + a) * g[n, m]) / (A * n)
      }
    }
  }
  if (any(!is.finite(g)) || any(g < 0)) return(NULL)
  log(pmax(g, .Machine$double.xmin))
}

ppc_cgfm_pair_logpmf <- function(n1, n2, u1, u2, theta, r) {
  tab <- ppc_cgfm_biv_logpmf_table(n1, n2, u1, u2, theta, r)
  if (is.null(tab)) return(-Inf)
  tab[n1 + 1L, n2 + 1L]
}

ppc_cgfm_fit <- function(y, blocks, include_time_varying = TRUE,
                         maxit = 3000L, compile_cpp = TRUE) {
  y <- as.matrix(y)
  storage.mode(y) <- "numeric"
  x <- ppc_cgfm_x(blocks)
  p <- ncol(y)

  unpack_independent <- function(par) {
    list(alpha = par[seq_len(p)], beta = par[p + 1L], theta = exp(par[p + 2L]))
  }
  unpack_correlated <- function(par) {
    z <- unpack_independent(par)
    z$rho <- stats::plogis(par[p + 3L])
    z
  }

  start_from_glm <- function() {
    dat <- data.frame(
      y = as.vector(t(y)),
      time = factor(rep(seq_len(p), times = nrow(y))),
      x = rep(x, each = p)
    )
    fit <- stats::glm(y ~ 0 + time + x, family = stats::poisson(), data = dat)
    co <- stats::coef(fit)
    c(unname(co[grep("^time", names(co))]), unname(co["x"]), log(0.45))
  }

  loglik_independent <- function(par) {
    pp <- unpack_independent(par)
    mu <- exp(outer(rep(1, nrow(y)), pp$alpha) + outer(x * pp$beta, rep(1, p)))
    sum(ppc_cgfm_nb_logpmf(y, mu, pp$theta))
  }

  loglik_shared <- function(par) {
    pp <- unpack_independent(par)
    mu <- exp(outer(rep(1, nrow(y)), pp$alpha) + outer(x * pp$beta, rep(1, p)))
    size <- 1 / pp$theta
    total <- 0
    for (ii in seq_len(nrow(y))) {
      nsum <- sum(y[ii, ])
      usum <- sum(mu[ii, ])
      total <- total +
        sum(y[ii, ] * log(mu[ii, ]) - lgamma(y[ii, ] + 1)) +
        lgamma(nsum + size) - lgamma(size) + size * log(size) -
        (nsum + size) * log(size + usum)
    }
    total
  }

  use_cpp <- isTRUE(compile_cpp) && ppc_cgfm_compile_cpp()
  if (isTRUE(include_time_varying) && !use_cpp) {
    stop(
      "Fold-wise time-varying CGFM requires Rcpp compilation/Rtools for practical runtime. ",
      "Install Rtools or run with include_time_varying = FALSE."
    )
  }
  loglik_timevarying <- function(par) {
    pp <- unpack_correlated(par)
    if (use_cpp) {
      return(ppc_cgfm_composite_loglik_recur_cpp(y, x, pp$alpha, pp$beta, pp$theta, pp$rho))
    }
    total <- 0
    for (j in 1:(p - 1L)) {
      for (k in (j + 1L):p) {
        lag_r <- pp$rho^(k - j)
        u1 <- exp(pp$alpha[j] + pp$beta * x)
        u2 <- exp(pp$alpha[k] + pp$beta * x)
        for (ii in seq_len(nrow(y))) {
          total <- total + ppc_cgfm_pair_logpmf(y[ii, j], y[ii, k], u1[ii], u2[ii], pp$theta, lag_r)
        }
      }
    }
    total / (p - 1L)
  }

  fit_model <- function(loglik_fun, start) {
    lower <- rep(-Inf, length(start))
    upper <- rep(Inf, length(start))
    lower[p + 2L] <- log(1e-3)
    upper[p + 2L] <- log(5)
    if (length(start) >= p + 3L) {
      lower[p + 3L] <- stats::qlogis(1e-4)
      upper[p + 3L] <- stats::qlogis(1 - 1e-4)
    }
    nll <- function(par) {
      below <- pmax(lower - par, 0)
      above <- pmax(par - upper, 0)
      if (any(below > 0 | above > 0)) return(1e100 + 1e20 * sum(below^2 + above^2))
      val <- tryCatch(loglik_fun(par), error = function(e) -Inf)
      if (!is.finite(val)) return(1e100)
      -val
    }
    stats::optim(start, nll, method = "BFGS",
                 control = list(maxit = maxit, reltol = 1e-10))
  }

  fit_with_starts <- function(loglik_fun, starts, model_name) {
    fits <- lapply(starts, function(st) {
      tryCatch(fit_model(loglik_fun, st), error = function(e) NULL)
    })
    fits <- Filter(Negate(is.null), fits)
    if (!length(fits)) stop("All starts failed for ", model_name)
    vals <- vapply(fits, function(z) -z$value, numeric(1))
    fits[[which.max(vals)]]
  }

  base <- start_from_glm()
  ind_starts <- list(base, replace(base, p + 2L, log(0.25)), replace(base, p + 2L, log(0.75)))
  tv_starts <- lapply(c(0.5, 0.75, 0.85, 0.95), function(rho) c(base, stats::qlogis(rho)))
  tv_starts <- c(tv_starts, lapply(c(0.25, 0.75), function(theta) {
    st <- c(base, stats::qlogis(0.85))
    st[p + 2L] <- log(theta)
    st
  }))

  ind <- fit_with_starts(loglik_independent, ind_starts, "cgfm_independent")
  shared <- fit_with_starts(loglik_shared, ind_starts, "cgfm_shared")
  tv <- if (include_time_varying) {
    fit_with_starts(loglik_timevarying, tv_starts, "cgfm_time_varying")
  } else {
    NULL
  }

  structure(
    list(
      independent = c(unpack_independent(ind$par), list(logLik = -ind$value, convergence = ind$convergence)),
      shared = c(unpack_independent(shared$par), list(logLik = -shared$value, convergence = shared$convergence)),
      time_varying = if (!is.null(tv)) c(unpack_correlated(tv$par), list(logLik = -tv$value, convergence = tv$convergence)) else NULL
    ),
    class = "ppc_cgfm_fit"
  )
}

ppc_cgfm_score_distribution <- function(model, observed, support, prob, epsilon) {
  prob <- pmax(as.numeric(prob), 0)
  prob <- prob / sum(prob)
  sc <- ppc_score_pmf(observed, support, prob, epsilon = epsilon)
  mean <- sum(support * prob)
  c(mean = mean, mc_se = NA_real_, sc,
    cover80 = NA_real_, cover95 = NA_real_,
    squared_error = (mean - observed)^2,
    absolute_error = abs(mean - observed))
}

ppc_cgfm_nb_support <- function(observed, mean, size, lower_support = 0L,
                                prob = 0.99999, min_upper = 60L) {
  upper <- max(observed, min_upper, stats::qnbinom(prob, size = size, mu = mean), na.rm = TRUE)
  lower_support:as.integer(ceiling(upper))
}

ppc_cgfm_score_nb <- function(observed, mean, size, epsilon) {
  support <- ppc_cgfm_nb_support(observed, mean, size)
  prob <- stats::dnbinom(support, size = size, mu = mean)
  ppc_cgfm_score_distribution("nb", observed, support, prob, epsilon)
}

ppc_cgfm_score_time_varying <- function(observed, y_prev, mu_prev, mu_now,
                                        theta, rho, epsilon,
                                        min_upper = 80L, tail_tol = 1e-7) {
  upper <- max(observed, min_upper, stats::qnbinom(0.99999, size = 1 / theta, mu = mu_now), na.rm = TRUE)
  repeat {
    support <- 0:as.integer(ceiling(upper))
    tab <- ppc_cgfm_biv_logpmf_table(y_prev, max(support), mu_prev, mu_now, theta, rho)
    if (is.null(tab)) stop("Failed to evaluate time-varying CGFM predictive PMF")
    log_joint <- tab[y_prev + 1L, support + 1L]
    log_denom <- ppc_cgfm_nb_logpmf(y_prev, mu_prev, theta)
    prob <- exp(log_joint - log_denom)
    mass <- sum(prob)
    if (mass > 1 - tail_tol || upper > 500L) break
    upper <- ceiling(upper * 1.5)
  }
  prob <- prob / sum(prob)
  ppc_cgfm_score_distribution("tv", observed, support, prob, epsilon)
}

ppc_cgfm_score_rolling <- function(fit, y_test, blocks_test, fold_id,
                                   epsilon = 1 / 1001) {
  x_test <- ppc_cgfm_x(blocks_test)
  rows <- list()
  n_time <- ncol(y_test)
  for (tt in 2:n_time) {
    observed <- y_test[, tt]
    for (ii in seq_len(nrow(y_test))) {
      for (model_name in c("cgfm_independent", "cgfm_shared", "cgfm_time_varying")) {
        if (identical(model_name, "cgfm_independent")) {
          pp <- fit$independent
          mean <- exp(pp$alpha[tt] + pp$beta * x_test[ii])
          sc <- ppc_cgfm_score_nb(observed[ii], mean, 1 / pp$theta, epsilon)
        } else if (identical(model_name, "cgfm_shared")) {
          pp <- fit$shared
          mu <- exp(pp$alpha + pp$beta * x_test[ii])
          shape <- 1 / pp$theta + sum(y_test[ii, 1:(tt - 1L)])
          rate <- 1 / pp$theta + sum(mu[1:(tt - 1L)])
          mean <- mu[tt] * shape / rate
          sc <- ppc_cgfm_score_nb(observed[ii], mean, shape, epsilon)
        } else {
          if (is.null(fit$time_varying)) next
          pp <- fit$time_varying
          mu <- exp(pp$alpha + pp$beta * x_test[ii])
          sc <- ppc_cgfm_score_time_varying(
            observed = observed[ii],
            y_prev = y_test[ii, tt - 1L],
            mu_prev = mu[tt - 1L],
            mu_now = mu[tt],
            theta = pp$theta,
            rho = pp$rho,
            epsilon = epsilon
          )
        }
        rows[[length(rows) + 1L]] <- data.frame(
          rep = fold_id, model = model_name, task = "rolling_one_step",
          time = tt, horizon = 1L, subject_index = ii, observed = observed[ii],
          as.data.frame(as.list(sc)),
          row.names = NULL, check.names = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

ppc_cgfm_diagnostics <- function(fit, fold_id) {
  parts <- list(
    data.frame(
      fold = fold_id, model = "cgfm_independent",
      logLik = fit$independent$logLik,
      theta = fit$independent$theta,
      rho = 0,
      convergence = fit$independent$convergence
    ),
    data.frame(
      fold = fold_id, model = "cgfm_shared",
      logLik = fit$shared$logLik,
      theta = fit$shared$theta,
      rho = 1,
      convergence = fit$shared$convergence
    )
  )
  if (!is.null(fit$time_varying)) {
    parts[[length(parts) + 1L]] <- data.frame(
      fold = fold_id, model = "cgfm_time_varying",
      logLik = fit$time_varying$logLik,
      theta = fit$time_varying$theta,
      rho = fit$time_varying$rho,
      convergence = fit$time_varying$convergence
    )
  }
  do.call(rbind, parts)
}
