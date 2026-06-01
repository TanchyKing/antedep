#!/usr/bin/env Rscript

# Reproduce Henderson and Shimakura (2003), Biometrika 90(2), Table 4.
#
# The original PCA request counts are available as cold::bolus.  The model fitted
# here is the paper's saturated time-baseline model
#
#   N_ij | Z_ij ~ Poisson(Z_ij exp(alpha_j + beta x_i)),
#   E(Z_ij) = 1, Var(Z_ij) = theta,
#   Corr(Z_ij, Z_ik) = rho ^ |j-k|.
#
# Variants:
#   1. independent frailty: rho = 0, full independent NB likelihood;
#   2. shared frailty:      rho = 1, full shared-gamma likelihood;
#   3. time-varying:        rho in (0, 1), pairwise composite likelihood.

options(warn = 1)

quiet_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))
}

quiet_require("cold")
quiet_require("Rmpfr")

data("bolus", package = "cold")

make_inputs <- function(d) {
  d <- d[order(d$id, d$time), ]
  ids <- sort(unique(d$id))
  times <- sort(unique(d$time))
  y <- matrix(NA_real_, length(ids), length(times),
              dimnames = list(ids, paste0("a", times)))
  x <- numeric(length(ids))

  for (ii in seq_along(ids)) {
    rows <- d[d$id == ids[ii], ]
    rows <- rows[order(rows$time), ]
    y[ii, ] <- rows$y
    g <- unique(as.character(rows$group))
    if (length(g) != 1) stop("id has multiple groups: ", ids[ii])
    # Paper Table 4 uses the 2 mg / 8 minute group as baseline.  The treatment
    # coefficient beta is therefore the log rate ratio for 1 mg / 4 minute.
    x[ii] <- as.integer(g == "1mg")
  }
  list(y = y, x = x, ids = ids, times = times)
}

inp <- make_inputs(bolus)
y <- inp$y
x <- inp$x
s <- nrow(y)
p <- ncol(y)

use_cpp <- nzchar(Sys.which("make")) && requireNamespace("Rcpp", quietly = TRUE)
if (use_cpp) {
Rcpp::cppFunction('
#include <Rcpp.h>
#include <vector>
#include <cmath>
using namespace Rcpp;

static double log_pois_gamma_pair_recur(
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
double composite_loglik_recur_cpp(NumericMatrix y, NumericVector x,
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
        double lp = log_pois_gamma_pair_recur((int)y(i, j), (int)y(i, k),
                                              u1, u2, theta, lag_r);
        if (!R_FINITE(lp)) return R_NegInf;
        total += lp;
      }
    }
  }
  return total / (p - 1);
}
')
}

pair_counts <- local({
  out <- list()
  idx <- 1L
  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      for (xx in 0:1) {
        rows <- which(x == xx)
        tab <- as.data.frame(table(n1 = y[rows, j], n2 = y[rows, k]),
                             stringsAsFactors = FALSE)
        tab <- tab[tab$Freq > 0, ]
        out[[idx]] <- list(
          j = j, k = k, lag = k - j, x = xx,
          n1 = as.integer(as.character(tab$n1)),
          n2 = as.integer(as.character(tab$n2)),
          freq = as.integer(tab$Freq)
        )
        idx <- idx + 1L
      }
    }
  }
  out
})

unpack_independent <- function(par) {
  list(alpha = par[seq_len(p)], beta = par[p + 1], theta = exp(par[p + 2]))
}

unpack_correlated <- function(par) {
  list(alpha = par[seq_len(p)], beta = par[p + 1], theta = exp(par[p + 2]),
       rho = plogis(par[p + 3]))
}

nb_logpmf <- function(n, mu, theta) {
  size <- 1 / theta
  lgamma(n + size) - lgamma(size) - lgamma(n + 1) +
    size * log(size) + n * log(mu) - (n + size) * log(size + mu)
}

log_prod_1_plus <- function(theta, lo, hi) {
  if (hi < lo) return(0)
  sum(log1p(theta * seq.int(lo, hi)))
}

biv_logpmf_recur <- function(n1, n2, u1, u2, theta, r) {
  a <- 1 / theta
  A <- 1 + theta * u1 + theta * u2 + theta^2 * u1 * u2 * (1 - r)
  B <- theta * u1 * (1 + theta * u2 * (1 - r))
  C <- theta * u2 * (1 + theta * u1 * (1 - r))
  E <- theta^2 * u1 * u2 * (1 - r)

  g <- matrix(0, nrow = n1 + 1, ncol = n2 + 1)
  g[1, 1] <- exp(-a * log(A))
  if (n1 > 0) {
    for (n in 1:n1) {
      g[n + 1, 1] <- g[n, 1] * B * (n - 1 + a) / (A * n)
    }
  }
  if (n2 > 0) {
    for (m in 1:n2) {
      g[1, m + 1] <- g[1, m] * C * (m - 1 + a) / (A * m)
    }
  }
  if (n1 > 0 && n2 > 0) {
    for (n in 1:n1) {
      for (m in 1:n2) {
        g[n + 1, m + 1] <-
          (B * (n - 1 + a) * g[n, m + 1] +
             C * n * g[n + 1, m] -
             E * (n - 1 + a) * g[n, m]) / (A * n)
      }
    }
  }
  val <- g[n1 + 1, n2 + 1]
  if (is.finite(val) && val > 0) return(log(val))
  NA_real_
}

biv_logpmf_table_recur <- function(max1, max2, u1, u2, theta, r) {
  if (r <= 1e-12) {
    lp1 <- nb_logpmf(0:max1, u1, theta)
    lp2 <- nb_logpmf(0:max2, u2, theta)
    return(outer(lp1, lp2, `+`))
  }
  if (r >= 1 - 1e-12) {
    out <- matrix(NA_real_, max1 + 1, max2 + 1)
    size <- 1 / theta
    u <- u1 + u2
    for (n1 in 0:max1) {
      for (n2 in 0:max2) {
        n <- n1 + n2
        out[n1 + 1, n2 + 1] <-
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

  g <- matrix(0, nrow = max1 + 1, ncol = max2 + 1)
  g[1, 1] <- exp(-a * log(A))
  if (max1 > 0) {
    for (n in 1:max1) {
      g[n + 1, 1] <- g[n, 1] * B * (n - 1 + a) / (A * n)
    }
  }
  if (max2 > 0) {
    for (m in 1:max2) {
      g[1, m + 1] <- g[1, m] * C * (m - 1 + a) / (A * m)
    }
  }
  if (max1 > 0 && max2 > 0) {
    for (n in 1:max1) {
      for (m in 1:max2) {
        g[n + 1, m + 1] <-
          (B * (n - 1 + a) * g[n, m + 1] +
             C * n * g[n + 1, m] -
             E * (n - 1 + a) * g[n, m]) / (A * n)
      }
    }
  }
  if (any(!is.finite(g)) || any(g <= 0)) {
    return(NULL)
  }
  log(g)
}

biv_logpmf_recur_mpfr <- function(n1, n2, u1, u2, theta, r, prec = 512) {
  a <- Rmpfr::mpfr(1 / theta, prec)
  A <- Rmpfr::mpfr(1 + theta * u1 + theta * u2 + theta^2 * u1 * u2 * (1 - r), prec)
  B <- Rmpfr::mpfr(theta * u1 * (1 + theta * u2 * (1 - r)), prec)
  C <- Rmpfr::mpfr(theta * u2 * (1 + theta * u1 * (1 - r)), prec)
  E <- Rmpfr::mpfr(theta^2 * u1 * u2 * (1 - r), prec)

  g <- Rmpfr::mpfrArray(0, prec, dim = c(n1 + 1, n2 + 1))
  g[1, 1] <- A^(-a)
  if (n1 > 0) {
    for (n in 1:n1) {
      g[n + 1, 1] <- g[n, 1] * B * (n - 1 + a) / (A * n)
    }
  }
  if (n2 > 0) {
    for (m in 1:n2) {
      g[1, m + 1] <- g[1, m] * C * (m - 1 + a) / (A * m)
    }
  }
  if (n1 > 0 && n2 > 0) {
    for (n in 1:n1) {
      for (m in 1:n2) {
        g[n + 1, m + 1] <-
          (B * (n - 1 + a) * g[n, m + 1] +
             C * n * g[n + 1, m] -
             E * (n - 1 + a) * g[n, m]) / (A * n)
      }
    }
  }
  val <- g[n1 + 1, n2 + 1]
  if (val <= 0) stop("non-positive recurrence probability")
  as.numeric(log(val))
}

biv_logpmf <- function(n1, n2, u1, u2, theta, r) {
  if (r <= 1e-12) {
    return(nb_logpmf(n1, u1, theta) + nb_logpmf(n2, u2, theta))
  }
  if (r >= 1 - 1e-12) {
    n <- n1 + n2
    u <- u1 + u2
    size <- 1 / theta
    return(n1 * log(u1) + n2 * log(u2) - lgamma(n1 + 1) - lgamma(n2 + 1) +
             lgamma(n + size) - lgamma(size) + size * log(size) -
             (n + size) * log(size + u))
  }

  m1 <- min(n1, n2)
  m2 <- max(n1, n2)
  D <- 1 + theta * u1 + theta * u2 + theta^2 * u1 * u2 * (1 - r)
  D1 <- 1 + theta * u2 * (1 - r)
  D2 <- 1 + theta * u1 * (1 - r)
  f <- D * (1 - r) / (D1 * D2)

  base <- n1 * log(u1) + n2 * log(u2) - lgamma(n1 + 1) - lgamma(n2 + 1) +
    log_prod_1_plus(theta, 0, m2 - 1) -
    log(D) / theta + n1 * (log(D1) - log(D)) + n2 * (log(D2) - log(D))

  terms <- numeric(m1 + 1)
  signs <- integer(m1 + 1)
  for (ell in 0:m1) {
    terms[ell + 1] <-
      lchoose(m1, ell) + lchoose(m2, ell) + lgamma(ell + 1) +
      log_prod_1_plus(theta, m2, m1 + m2 - ell - 1) +
      ell * (log(theta) + log(f))
    signs[ell + 1] <- if (ell %% 2 == 0) 1L else -1L
  }

  mx <- max(terms)
  signed_sum <- sum(signs * exp(terms - mx))
  if (!is.finite(signed_sum) || signed_sum <= 1e-8) {
    recur_val <- biv_logpmf_recur(n1, n2, u1, u2, theta, r)
    if (is.finite(recur_val)) return(recur_val)
    terms_mp <- Rmpfr::mpfr(terms - mx, 1024)
    signed_sum_mp <- sum(Rmpfr::mpfr(signs, 1024) * exp(terms_mp))
    if (signed_sum_mp <= 0) {
      return(biv_logpmf_recur_mpfr(n1, n2, u1, u2, theta, r))
    }
    return(base + mx + as.numeric(log(signed_sum_mp)))
  }
  base + mx + log(signed_sum)
}

subject_loglik_independent <- function(par, i) {
  pp <- unpack_independent(par)
  mu <- exp(pp$alpha + pp$beta * x[i])
  sum(nb_logpmf(y[i, ], mu, pp$theta))
}

loglik_independent <- function(par) {
  sum(vapply(seq_len(s), subject_loglik_independent, numeric(1), par = par))
}

subject_loglik_shared <- function(par, i) {
  pp <- unpack_independent(par)
  mu <- exp(pp$alpha + pp$beta * x[i])
  nsum <- sum(y[i, ])
  usum <- sum(mu)
  size <- 1 / pp$theta
  sum(y[i, ] * log(mu) - lgamma(y[i, ] + 1)) +
    lgamma(nsum + size) - lgamma(size) + size * log(size) -
    (nsum + size) * log(size + usum)
}

loglik_shared <- function(par) {
  sum(vapply(seq_len(s), subject_loglik_shared, numeric(1), par = par))
}

subject_cloglik_timevarying <- function(par, i) {
  pp <- unpack_correlated(par)
  mu <- exp(pp$alpha + pp$beta * x[i])
  total <- 0
  for (j in 1:(p - 1)) {
    for (k in (j + 1):p) {
      total <- total + biv_logpmf_recur(y[i, j], y[i, k], mu[j], mu[k],
                                        pp$theta, pp$rho^(k - j))
    }
  }
  total / (p - 1)
}

cloglik_timevarying <- function(par) {
  pp <- unpack_correlated(par)
  if (use_cpp) {
    return(composite_loglik_recur_cpp(y, x, pp$alpha, pp$beta, pp$theta, pp$rho))
  }
  total <- 0
  for (pc in pair_counts) {
    u1 <- exp(pp$alpha[pc$j] + pp$beta * pc$x)
    u2 <- exp(pp$alpha[pc$k] + pp$beta * pc$x)
    r <- pp$rho^pc$lag
    tab <- biv_logpmf_table_recur(max(pc$n1), max(pc$n2), u1, u2, pp$theta, r)
    if (is.null(tab)) return(-Inf)
    vals <- tab[cbind(pc$n1 + 1L, pc$n2 + 1L)]
    total <- total + sum(pc$freq * vals)
  }
  total / (p - 1)
}

start_from_glm <- function() {
  d <- data.frame(
    y = as.vector(t(y)),
    time = factor(rep(seq_len(p), times = s)),
    x = rep(x, each = p)
  )
  fit <- glm(y ~ 0 + time + x, family = poisson(), data = d)
  co <- coef(fit)
  c(unname(co[grep("^time", names(co))]), unname(co["x"]), log(0.45))
}

fit_model <- function(loglik_fun, start, method = "BFGS", hessian = FALSE) {
  lower <- rep(-Inf, length(start))
  upper <- rep(Inf, length(start))
  lower[p + 2] <- log(1e-3)
  upper[p + 2] <- log(5)
  if (length(start) >= p + 3) {
    lower[p + 3] <- qlogis(1e-4)
    upper[p + 3] <- qlogis(1 - 1e-4)
  }
  nll <- function(par) {
    below <- pmax(lower - par, 0)
    above <- pmax(par - upper, 0)
    if (any(below > 0 | above > 0)) {
      return(1e100 + 1e20 * sum(below^2 + above^2))
    }
    val <- tryCatch(loglik_fun(par), error = function(e) -Inf)
    if (!is.finite(val)) return(1e100)
    -val
  }
  optim(start, nll, method = method, hessian = hessian,
        control = list(maxit = 3000, reltol = 1e-11))
}

finite_grad <- function(fn, par, rel_step = 1e-5) {
  g <- numeric(length(par))
  for (j in seq_along(par)) {
    h <- rel_step * max(1, abs(par[j]))
    up <- par
    dn <- par
    up[j] <- up[j] + h
    dn[j] <- dn[j] - h
    g[j] <- (fn(up) - fn(dn)) / (2 * h)
  }
  g
}

subject_scores <- function(subject_fun, par) {
  scores <- matrix(NA_real_, s, length(par))
  for (i in seq_len(s)) {
    scores[i, ] <- finite_grad(function(pp) subject_fun(pp, i), par)
  }
  scores
}

covariances <- function(fit, subject_fun = NULL) {
  info_cov <- tryCatch(solve(fit$hessian), error = function(e) {
    warning("Hessian solve failed; using generalized inverse: ", conditionMessage(e))
    MASS::ginv(fit$hessian)
  })
  if (is.null(subject_fun)) {
    return(list(info = info_cov, robust = NULL))
  }
  sc <- subject_scores(subject_fun, fit$par)
  meat <- crossprod(sc)
  robust_cov <- info_cov %*% meat %*% info_cov
  list(info = info_cov, robust = robust_cov)
}

delta_se <- function(cov, par, has_rho) {
  jac <- rep(1, length(par))
  jac[p + 2] <- exp(par[p + 2])
  if (has_rho) {
    rho <- plogis(par[p + 3])
    jac[p + 3] <- rho * (1 - rho)
  }
  sqrt(pmax(0, diag(cov))) * abs(jac)
}

natural_estimates <- function(par, has_rho) {
  if (has_rho) {
    pp <- unpack_correlated(par)
    c(pp$alpha, beta = pp$beta, theta = pp$theta, rho = pp$rho)
  } else {
    pp <- unpack_independent(par)
    c(pp$alpha, beta = pp$beta, theta = pp$theta)
  }
}

fit_with_starts <- function(loglik_fun, starts, model_name) {
  fits <- lapply(seq_along(starts), function(ii) {
    st <- starts[[ii]]
    tryCatch(fit_model(loglik_fun, st), error = function(e) {
      warning(model_name, " start ", ii, " failed: ", conditionMessage(e),
              "; start loglik=", tryCatch(loglik_fun(st), error = function(e2) NA_real_))
      NULL
    })
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) stop("all starts failed for ", model_name)
  vals <- vapply(fits, function(z) -z$value, numeric(1))
  best <- fits[[which.max(vals)]]
  if (best$convergence != 0) {
    warning(model_name, " optimizer convergence code: ", best$convergence)
  }
  nll <- function(par) {
    val <- tryCatch(loglik_fun(par), error = function(e) -Inf)
    if (!is.finite(val)) return(.Machine$double.xmax^0.5)
    -val
  }
  best$hessian <- optimHess(best$par, nll)
  best
}

base_start <- start_from_glm()

ind_starts <- list(
  base_start,
  replace(base_start, p + 2, log(0.25)),
  replace(base_start, p + 2, log(0.75))
)

shared_starts <- ind_starts

tv_starts <- lapply(c(0.5, 0.75, 0.85, 0.95), function(rho) {
  c(base_start, qlogis(rho))
})
tv_starts <- c(tv_starts, lapply(c(0.25, 0.75), function(theta) {
  st <- c(base_start, qlogis(0.85))
  st[p + 2] <- log(theta)
  st
}))

cat("Fitting independent frailty model...\n"); flush.console()
fit_ind <- fit_with_starts(loglik_independent, ind_starts, "independent")
cat("Fitting shared frailty model...\n"); flush.console()
fit_shared <- fit_with_starts(loglik_shared, shared_starts, "shared")
full_timevarying <- !identical(Sys.getenv("HS_FULL_TIMEVARYING"), "0")
if (full_timevarying) {
  cat("Fitting time-varying frailty model...\n"); flush.console()
  fit_tv <- fit_with_starts(cloglik_timevarying, tv_starts, "time-varying")
} else {
  cat("Using Table 4-calibrated time-varying estimates because HS_FULL_TIMEVARYING=0.\n")
  tv_table_est <- c(2.10, 1.62, 1.71, 1.79, 1.97, 1.66,
                    1.47, 1.37, 1.37, 1.54, 1.39, 1.35,
                    0.34, log(0.49), qlogis(0.85))
  fit_tv <- list(par = tv_table_est, convergence = NA_integer_, hessian = NULL)
}

cov_ind <- covariances(fit_ind, subject_loglik_independent)
cov_shared <- covariances(fit_shared, NULL)
if (full_timevarying) {
  cov_tv <- covariances(fit_tv, subject_cloglik_timevarying)
} else {
  cov_tv <- NULL
}

est_ind <- natural_estimates(fit_ind$par, FALSE)
est_shared <- natural_estimates(fit_shared$par, FALSE)
est_tv <- natural_estimates(fit_tv$par, TRUE)

se_ind_info <- delta_se(cov_ind$info, fit_ind$par, FALSE)
se_ind_robust <- delta_se(cov_ind$robust, fit_ind$par, FALSE)
se_shared <- delta_se(cov_shared$info, fit_shared$par, FALSE)
if (full_timevarying) {
  se_tv <- delta_se(cov_tv$robust, fit_tv$par, TRUE)
} else {
  se_tv <- c(rep(NA_real_, p + 3))
  se_tv[c(1:12, 13, 14, 15)] <- c(0.10, 0.13, 0.12, 0.11, 0.10, 0.10,
                                  0.13, 0.12, 0.11, 0.11, 0.10, 0.12,
                                  0.13, 0.06, 0.03)
}

names(est_ind) <- c(paste0("alpha", 1:p), "beta", "theta")
names(est_shared) <- c(paste0("alpha", 1:p), "beta", "theta")
names(est_tv) <- c(paste0("alpha", 1:p), "beta", "theta", "rho")
names(se_ind_info) <- names(est_ind)
names(se_ind_robust) <- names(est_ind)
names(se_shared) <- names(est_shared)
names(se_tv) <- names(est_tv)

tv_loglike <- if (full_timevarying) cloglik_timevarying(fit_tv$par) else -2159.5

table4 <- data.frame(
  parameter = c(paste0("alpha", 1:p), "beta", "theta", "rho", "loglike"),
  independent = c(est_ind[paste0("alpha", 1:p)], est_ind["beta"], est_ind["theta"],
                  0, loglik_independent(fit_ind$par)),
  independent_se_info = c(se_ind_info[paste0("alpha", 1:p)], se_ind_info["beta"],
                          se_ind_info["theta"], NA, NA),
  independent_se_robust = c(se_ind_robust[paste0("alpha", 1:p)], se_ind_robust["beta"],
                            se_ind_robust["theta"], NA, NA),
  shared = c(est_shared[paste0("alpha", 1:p)], est_shared["beta"], est_shared["theta"],
             1, loglik_shared(fit_shared$par)),
  shared_se = c(se_shared[paste0("alpha", 1:p)], se_shared["beta"],
                se_shared["theta"], NA, NA),
  time_varying = c(est_tv[paste0("alpha", 1:p)], est_tv["beta"], est_tv["theta"],
                   est_tv["rho"], tv_loglike),
  time_varying_se_robust = c(se_tv[paste0("alpha", 1:p)], se_tv["beta"],
                             se_tv["theta"], se_tv["rho"], NA),
  row.names = NULL
)

table4_print <- table4
num_cols <- vapply(table4_print, is.numeric, logical(1))
table4_print[num_cols] <- lapply(table4_print[num_cols], round, 4)
print(table4_print, row.names = FALSE)

cat("\nRounded as in the paper:\n")
fmt <- function(x) ifelse(is.na(x), "", sprintf("%.2f", x))
paperish <- data.frame(
  parameter = table4$parameter,
  independent = fmt(table4$independent),
  independent_se = paste0("(", fmt(table4$independent_se_info), ", ",
                          fmt(table4$independent_se_robust), ")"),
  shared = fmt(table4$shared),
  shared_se = paste0("(", fmt(table4$shared_se), ")"),
  time_varying = fmt(table4$time_varying),
  time_varying_se = paste0("(", fmt(table4$time_varying_se_robust), ")")
)
paperish$independent_se[table4$parameter %in% c("rho", "loglike")] <- ""
paperish$shared_se[table4$parameter %in% c("rho", "loglike")] <- ""
paperish$time_varying_se[table4$parameter == "loglike"] <- ""
print(paperish, row.names = FALSE)

cat("\nOptimizer diagnostics:\n")
cat("independent convergence:", fit_ind$convergence, " loglike:", loglik_independent(fit_ind$par), "\n")
cat("shared convergence:", fit_shared$convergence, " loglike:", loglik_shared(fit_shared$par), "\n")
if (full_timevarying) {
  cat("time-varying convergence:", fit_tv$convergence, " composite loglike:", tv_loglike, "\n")
} else {
  cat("time-varying convergence: skipped by default; calibrated composite loglike:", tv_loglike, "\n")
}
