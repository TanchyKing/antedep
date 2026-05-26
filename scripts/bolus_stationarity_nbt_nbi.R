## Reproduce stationarity tests for the bolus data under
## the NBT-NBI-INADFE(1) model (nbinom thinning, nbinom innovation, order 1).
##
## Unconstrained model parameters (order 1, n_time = 12):
##   alpha : 11  (alpha_1 = 0, then alpha_2,...,alpha_12)
##   theta : 12  (innovation mean at each time)
##   kappa : 12  (innovation nbinom size at each time)
##   total : 35
##
## Three stationarity hypotheses are tested via LRT:
##   H0a : alpha constant over t
##   H0b : theta AND kappa both constant over t
##   H0c : all three (alpha, theta, kappa) constant over t

suppressPackageStartupMessages({
  library(antedep)
  library(nloptr)
})

data("bolus_inad")
y      <- bolus_inad$y        ## 65 x 12
blocks <- bolus_inad$blocks   ## 1 = 2mg, 2 = 1mg
N      <- ncol(y)
n_subj <- nrow(y)
B      <- length(unique(blocks))
stopifnot(N == 12, B == 2)

ORDER <- 1
THINNING <- "nbinom"
INNOVATION <- "nbinom"

## ---------------------------------------------------------------
## 1. Unconstrained fit
## ---------------------------------------------------------------
set.seed(1)
fit_u <- fit_inad(
  y,
  order = ORDER,
  thinning = THINNING,
  innovation = INNOVATION,
  blocks = blocks,
  max_iter = 200,
  tol = 1e-8
)

logL_u <- fit_u$log_l
k_u    <- 36                          ## 11 + 12 + 12 + 1 (tau)
bic_u  <- -2 * logL_u + k_u * log(n_subj)

cat("=== Unconstrained NBT-NBI-INADFE(1) ===\n")
cat(sprintf("logL = %.4f   k = %d   BIC = %.4f\n", logL_u, k_u, bic_u))
cat("alpha:\n"); print(round(fit_u$alpha, 4))
cat("theta:\n"); print(round(fit_u$theta, 4))
cat("kappa:\n"); print(round(fit_u$nb_inno_size, 4))
cat("tau:  \n"); print(round(fit_u$tau,   4))

## ---------------------------------------------------------------
## 2. Helper for constrained fits
## ---------------------------------------------------------------
## par layout:
##   alpha part: 1 if alpha const else (N-1) (alpha_1 is fixed at 0)
##   theta part: 1 if theta const else N
##   kappa part: 1 if kappa const else N
##
fit_constrained <- function(alpha_const, theta_const, kappa_const,
                            init_alpha = fit_u$alpha,
                            init_theta = fit_u$theta,
                            init_kappa = fit_u$nb_inno_size,
                            init_tau   = fit_u$tau) {

  ## Build initial par and bounds
  par0 <- numeric(); lo <- numeric(); up <- numeric()

  if (alpha_const) {
    par0 <- c(par0, mean(init_alpha[-1]))
    lo   <- c(lo, 1e-8); up <- c(up, 50)
  } else {
    par0 <- c(par0, pmax(init_alpha[-1], 1e-4))
    lo   <- c(lo, rep(1e-8, N - 1)); up <- c(up, rep(50, N - 1))
  }

  if (theta_const) {
    par0 <- c(par0, mean(init_theta))
    lo   <- c(lo, 1e-8); up <- c(up, 1e6)
  } else {
    par0 <- c(par0, pmax(init_theta, 1e-4))
    lo   <- c(lo, rep(1e-8, N)); up <- c(up, rep(1e6, N))
  }

  if (kappa_const) {
    par0 <- c(par0, mean(init_kappa))
    lo   <- c(lo, 1e-6); up <- c(up, 1e6)
  } else {
    par0 <- c(par0, pmax(init_kappa, 1e-4))
    lo   <- c(lo, rep(1e-6, N)); up <- c(up, rep(1e6, N))
  }

  ## tau: first block fixed at 0, free params are tau[2..B]
  par0 <- c(par0, init_tau[-1])
  lo   <- c(lo, rep(-1e3, B - 1)); up <- c(up, rep(1e3, B - 1))

  unpack <- function(par) {
    idx <- 1
    if (alpha_const) { alpha <- c(0, rep(par[idx], N - 1)); idx <- idx + 1 }
    else             { alpha <- c(0, par[idx:(idx + N - 2)]); idx <- idx + N - 1 }
    if (theta_const) { theta <- rep(par[idx], N); idx <- idx + 1 }
    else             { theta <- par[idx:(idx + N - 1)]; idx <- idx + N }
    if (kappa_const) { kappa <- rep(par[idx], N); idx <- idx + 1 }
    else             { kappa <- par[idx:(idx + N - 1)]; idx <- idx + N }
    tau <- c(0, par[idx:(idx + B - 2)])
    list(alpha = alpha, theta = theta, kappa = kappa, tau = tau)
  }

  obj <- function(par) {
    p <- unpack(par)
    ll <- tryCatch(
      logL_inad(y, order = ORDER, thinning = THINNING, innovation = INNOVATION,
                alpha = p$alpha, theta = p$theta, nb_inno_size = p$kappa,
                blocks = blocks, tau = p$tau),
      error = function(e) -Inf
    )
    if (!is.finite(ll)) return(1e20)
    -ll
  }

  opt <- nloptr::nloptr(
    x0 = par0, eval_f = obj, lb = lo, ub = up,
    opts = list(algorithm = "NLOPT_LN_BOBYQA",
                xtol_rel = 1e-8, maxeval = 20000)
  )
  est <- unpack(opt$solution)
  list(alpha = est$alpha, theta = est$theta, kappa = est$kappa, tau = est$tau,
       log_l = -opt$objective, convergence = opt$status)
}

## Free-parameter counts for each constrained model (+ (B-1) tau params)
n_params <- function(alpha_const, theta_const, kappa_const) {
  (if (alpha_const) 1 else (N - 1)) +
  (if (theta_const) 1 else N) +
  (if (kappa_const) 1 else N) +
  (B - 1)
}

## ---------------------------------------------------------------
## 3. Fit the three constrained models
## ---------------------------------------------------------------
cases <- list(
  "H0a: alpha constant"              = list(a = TRUE,  t = FALSE, k = FALSE),
  "H0b: theta and kappa constant"    = list(a = FALSE, t = TRUE,  k = TRUE),
  "H0c: alpha, theta, kappa constant"= list(a = TRUE,  t = TRUE,  k = TRUE)
)

results <- list()
for (nm in names(cases)) {
  cat("\n--- Fitting ", nm, " ---\n", sep = "")
  cs <- cases[[nm]]
  fit_c <- fit_constrained(cs$a, cs$t, cs$k)
  k_c   <- n_params(cs$a, cs$t, cs$k)
  bic_c <- -2 * fit_c$log_l + k_c * log(n_subj)
  lrt   <- 2 * (logL_u - fit_c$log_l)
  df    <- k_u - k_c
  pval  <- 1 - pchisq(lrt, df = df)
  results[[nm]] <- list(
    logL = fit_c$log_l, k = k_c, BIC = bic_c,
    LRT = lrt, df = df, p = pval, fit = fit_c
  )
  cat(sprintf("  logL = %.4f   k = %d   BIC = %.4f\n",
              fit_c$log_l, k_c, bic_c))
  cat(sprintf("  LRT = %.4f   df = %d   p-value = %.4g\n",
              lrt, df, pval))
}

## ---------------------------------------------------------------
## 4. Summary table
## ---------------------------------------------------------------
summary_df <- data.frame(
  Model   = c("Unconstrained (NBT-NBI-INADFE(1))", names(results)),
  logLik  = c(logL_u, vapply(results, `[[`, numeric(1), "logL")),
  k       = c(k_u,    vapply(results, `[[`, numeric(1), "k")),
  BIC     = c(bic_u,  vapply(results, `[[`, numeric(1), "BIC")),
  LRT     = c(NA,     vapply(results, `[[`, numeric(1), "LRT")),
  df      = c(NA,     vapply(results, `[[`, numeric(1), "df")),
  p_value = c(NA,     vapply(results, `[[`, numeric(1), "p")),
  row.names = NULL, check.names = FALSE
)

cat("\n================ SUMMARY ================\n")
print(summary_df, digits = 5, row.names = FALSE)

invisible(list(unconstrained = fit_u, constrained = results,
               summary = summary_df))
