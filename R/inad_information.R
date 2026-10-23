# Observed information for INAD models.
#
# Louis' identity gives the observed information exactly: for every time point
# the latent thinned counts have a finite conditional distribution, so the
# conditional moments of the complete-data score and Hessian are finite sums.
# The time-point contributions are mapped onto the free parameters of the fit,
# which may share parameters across time (time-invariance constraints) and
# include block effects. A numerical Hessian with fixed absolute steps is kept
# as an independent check.

#' Free parameters of a fitted INAD model
#'
#' @return A list with \code{tab}, a data frame with one row per free
#'   parameter (columns \code{component}, \code{group}, \code{param},
#'   \code{est}, \code{index}), and the unpacked parameter vectors.
#' @keywords internal
.inad_free_params <- function(fit, N, B) {
    ord <- fit$settings$order
    innovation <- fit$settings$innovation
    cinfo <- fit$settings$constraint_info
    if (is.null(cinfo)) {
        cinfo <- list(alpha1_const = FALSE, alpha2_const = FALSE,
                      theta_const = FALSE, kappa_const = FALSE)
    }
    au <- if (ord >= 1) .unpack_alpha(ord, fit$alpha, N) else list(a1 = NULL, a2 = NULL)
    theta <- rep_len(as.numeric(fit$theta), N)
    kappa <- if (innovation == "nbinom") rep_len(as.numeric(fit$nb_inno_size), N) else NULL
    tau <- if (B > 1L) as.numeric(fit$tau)[seq_len(B)] else 0

    rows <- list()
    add <- function(component, group, param, est, index) {
        rows[[length(rows) + 1L]] <<- data.frame(
            component = component, group = group, param = param, est = est,
            index = index, stringsAsFactors = FALSE
        )
    }
    lab1 <- if (ord == 2) "alpha1" else "alpha"
    if (ord >= 1) {
        if (cinfo$alpha1_const) add("alpha", "a1", lab1, au$a1[2], NA_integer_)
        else add("alpha", "a1", paste0(lab1, "[", 2:N, "]"), au$a1[2:N], 2:N)
    }
    if (ord == 2) {
        if (cinfo$alpha2_const) add("alpha", "a2", "alpha2", au$a2[3], NA_integer_)
        else add("alpha", "a2", paste0("alpha2[", 3:N, "]"), au$a2[3:N], 3:N)
    }
    if (cinfo$theta_const) add("theta", "theta", "theta", theta[1], NA_integer_)
    else add("theta", "theta", paste0("theta[", 1:N, "]"), theta, 1:N)
    if (!is.null(kappa)) {
        if (cinfo$kappa_const) add("nb_inno_size", "kappa", "kappa", kappa[1], NA_integer_)
        else add("nb_inno_size", "kappa", paste0("kappa[", 1:N, "]"), kappa, 1:N)
    }
    if (B > 1L) add("tau", "tau", paste0("tau[", 2:B, "]"), tau[2:B], 2:B)
    tab <- do.call(rbind, rows)
    # Estimates on the boundary of the parameter space (a thinning coefficient
    # at zero, or a negative binomial size at its upper bound) have no Wald
    # interval; they are held fixed in information calculations.
    kappa_ub <- if (!is.null(fit$settings$nb_inno_size_ub)) fit$settings$nb_inno_size_ub else 50
    tab$boundary <- (tab$group %in% c("a1", "a2") & tab$est <= 1e-8) |
        (tab$group == "kappa" & tab$est >= kappa_ub * (1 - 1e-6))
    list(tab = tab, a1 = au$a1, a2 = au$a2, theta = theta,
         kappa = kappa, tau = tau, cinfo = cinfo)
}

#' Row of the free-parameter table for a local (time or block) parameter
#' @keywords internal
.inad_free_index <- function(tab, group, index) {
    k <- which(tab$group == group & (is.na(tab$index) | tab$index == index))
    if (length(k) != 1L) stop("internal error: cannot map ", group, "[", index, "]")
    k
}

#' Write free-parameter values back into a fit-like list
#' @keywords internal
.inad_set_free <- function(fit, fp, values, N) {
    tab <- fp$tab
    a1 <- fp$a1; a2 <- fp$a2; theta <- fp$theta; kappa <- fp$kappa; tau <- fp$tau
    for (k in seq_len(nrow(tab))) {
        idx <- tab$index[k]
        v <- values[k]
        switch(tab$group[k],
               a1 = if (is.na(idx)) a1[2:N] <- v else a1[idx] <- v,
               a2 = if (is.na(idx)) a2[3:N] <- v else a2[idx] <- v,
               theta = if (is.na(idx)) theta[] <- v else theta[idx] <- v,
               kappa = if (is.na(idx)) kappa[] <- v else kappa[idx] <- v,
               tau = tau[idx] <- v)
    }
    ord <- fit$settings$order
    fit$alpha <- if (ord == 2) cbind(a1, a2) else a1
    fit$theta <- theta
    if (!is.null(kappa)) fit$nb_inno_size <- kappa
    if (length(tau) > 1L) fit$tau <- tau
    fit
}

#' Bell parameter for a given mean, to machine precision
#' @keywords internal
.bell_theta_newton <- function(mu, start = NULL) {
    th <- if (is.null(start) || !is.finite(start)) .bell_theta_from_mean(mu) else start
    for (k in seq_len(50)) {
        step <- (th * exp(th) - mu) / (exp(th) * (1 + th))
        th <- th - step
        if (abs(step) <= 1e-15 * max(1, abs(th))) break
    }
    th
}

#' Complete-data score and negative Hessian of one thinned count
#'
#' @return A list with the score \code{s} and the negative second derivative
#'   \code{b} of \code{log P(alpha odot yprev = z)} with respect to alpha.
#' @keywords internal
.louis_thin_terms <- function(z, yprev, alpha, thinning) {
    a <- max(alpha, 1e-10)
    if (thinning == "binom") {
        a <- min(a, 1 - 1e-10)
        list(s = z / a - (yprev - z) / (1 - a), b = z / a^2 + (yprev - z) / (1 - a)^2)
    } else if (thinning == "pois") {
        list(s = z / a - yprev, b = z / a^2)
    } else {
        list(s = z / a - (yprev + z) / (1 + a), b = z / a^2 - (yprev + z) / (1 + a)^2)
    }
}

#' Louis score and information for the parameters entering time i
#'
#' Local parameters are \code{a1} (i >= 2), \code{a2} (order 2, i >= 3),
#' \code{theta}, \code{kappa} (negative binomial innovations), and the block
#' effects \code{tau2}, ..., \code{tauB}. Returns the observed score (the
#' conditional mean of the complete-data score) and the observed information
#' \eqn{E(B \mid y) - \mathrm{Cov}(S \mid y)}, both summed over subjects.
#'
#' @keywords internal
.louis_time_terms <- function(y, i, ord, thinning, innovation, a1, a2, theta_i,
                              kappa_i, tau, blocks) {
    n <- nrow(y)
    B <- if (is.null(blocks)) 1L else length(tau)
    has_a1 <- ord >= 1 && i >= 2
    has_a2 <- ord == 2 && i >= 3
    is_nb <- innovation == "nbinom"
    tau_nm <- if (B > 1L) paste0("tau", 2:B) else character(0)
    nm <- c(if (has_a1) "a1", if (has_a2) "a2", "theta", if (is_nb) "kappa", tau_nm)
    P <- length(nm)
    info <- matrix(0, P, P, dimnames = list(nm, nm))
    score <- stats::setNames(numeric(P), nm)
    theta_bell <- if (innovation == "bell") theta_i else NA_real_

    for (s in seq_len(n)) {
        yi <- y[s, i]
        b <- if (is.null(blocks)) 1L else blocks[s]
        tb <- if (b > 1L) tau[b] else 0
        tnm <- if (b > 1L) paste0("tau", b) else NULL

        if (has_a2) {
            g <- expand.grid(z1 = 0:yi, z2 = 0:yi)
            g <- g[g$z1 + g$z2 <= yi, , drop = FALSE]
            z1 <- g$z1; z2 <- g$z2
        } else if (has_a1) {
            z1 <- 0:yi; z2 <- 0 * z1
        } else {
            z1 <- 0; z2 <- 0
        }
        w <- yi - z1 - z2

        if (innovation == "bell") {
            v <- if (tb == 0) theta_bell else
                .bell_theta_newton(.bell_mean_from_theta(theta_bell) + tb)
            par_inno <- v
        } else {
            m <- theta_i + tb
            par_inno <- m
        }

        lp <- log(.innov_vec(w, par_inno, innovation, if (is_nb) kappa_i else NA_real_))
        if (has_a1) lp <- lp + log(.thin_vec(z1, y[s, i - 1], a1[i], thinning))
        if (has_a2) lp <- lp + log(.thin_vec(z2, y[s, i - 2], a2[i], thinning))
        lZ <- .log_sum_exp(lp)
        if (!is.finite(lZ)) stop("Zero likelihood contribution at time ", i, ", subject ", s)
        post <- exp(lp - lZ)

        S <- matrix(0, length(w), P, dimnames = list(NULL, nm))
        Bm <- matrix(0, P, P, dimnames = list(nm, nm))
        if (has_a1) {
            tt <- .louis_thin_terms(z1, y[s, i - 1], a1[i], thinning)
            S[, "a1"] <- tt$s; Bm["a1", "a1"] <- sum(post * tt$b)
        }
        if (has_a2) {
            tt <- .louis_thin_terms(z2, y[s, i - 2], a2[i], thinning)
            S[, "a2"] <- tt$s; Bm["a2", "a2"] <- sum(post * tt$b)
        }

        # Innovation terms. For Poisson and negative binomial innovations the
        # mean m = theta_i + tau_b, so d m / d theta_i = d m / d tau_b = 1.
        mean_nm <- c("theta", tnm)
        if (innovation == "pois") {
            mm <- max(m, 1e-12)
            s_m <- w / mm - 1
            b_mm <- sum(post * w / mm^2)
            for (a in mean_nm) S[, a] <- s_m
            Bm[mean_nm, mean_nm] <- b_mm
        } else if (innovation == "nbinom") {
            mm <- max(m, 1e-12); k <- max(kappa_i, 1e-12)
            s_m <- w / mm - (w + k) / (k + mm)
            for (a in mean_nm) S[, a] <- s_m
            S[, "kappa"] <- digamma(w + k) - digamma(k) + log(k / (k + mm)) + (mm - w) / (k + mm)
            Bm[mean_nm, mean_nm] <- sum(post * (w / mm^2 - (w + k) / (k + mm)^2))
            Bm["kappa", "kappa"] <- -sum(post * (trigamma(w + k) - trigamma(k) + 1 / k -
                                                    1 / (k + mm) - (mm - w) / (k + mm)^2))
            b_mk <- sum(post * (mm - w) / (k + mm)^2)
            Bm["kappa", mean_nm] <- b_mk
            Bm[mean_nm, "kappa"] <- b_mk
        } else {
            # Bell: the Bell parameter v solves v exp(v) = theta exp(theta) + tau_b.
            th <- theta_bell
            q <- exp(v) * (1 + v); q1 <- exp(v) * (2 + v)
            c_th <- exp(th) * (1 + th) / q
            c_thth <- (exp(th) * (2 + th) - c_th^2 * q1) / q
            dl <- w / v - exp(v)
            d2l <- -w / v^2 - exp(v)
            S[, "theta"] <- dl * c_th
            Bm["theta", "theta"] <- -sum(post * (d2l * c_th^2 + dl * c_thth))
            if (!is.null(tnm)) {
                c_t <- 1 / q
                c_tt <- -q1 * c_t^3
                c_tht <- -q1 * c_th * c_t^2
                S[, tnm] <- dl * c_t
                Bm[tnm, tnm] <- -sum(post * (d2l * c_t^2 + dl * c_tt))
                Bm["theta", tnm] <- Bm[tnm, "theta"] <- -sum(post * (d2l * c_th * c_t + dl * c_tht))
            }
        }
        ES <- colSums(post * S)
        score <- score + ES
        info <- info + Bm - (crossprod(S, post * S) - tcrossprod(ES))
    }
    list(score = score, info = info)
}

#' Observed score and information for all free parameters, by Louis' identity
#'
#' @keywords internal
.louis_full_info <- function(y, fit, blocks = NULL) {
    N <- ncol(y)
    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    fp <- .inad_free_params(fit, N, B)
    tab <- fp$tab
    P <- nrow(tab)
    info <- matrix(0, P, P, dimnames = list(tab$param, tab$param))
    score <- stats::setNames(numeric(P), tab$param)
    for (i in seq_len(N)) {
        tt <- .louis_time_terms(y, i, ord, thinning, innovation, fp$a1, fp$a2, fp$theta[i],
                                if (!is.null(fp$kappa)) fp$kappa[i] else NULL, fp$tau, blocks)
        loc <- names(tt$score)
        map <- vapply(loc, function(p) {
            if (p %in% c("a1", "a2", "theta", "kappa")) .inad_free_index(tab, p, i)
            else .inad_free_index(tab, "tau", as.integer(sub("tau", "", p, fixed = TRUE)))
        }, integer(1))
        score[map] <- score[map] + tt$score
        info[map, map] <- info[map, map] + tt$info
    }
    list(score = score, info = info, fp = fp, B = B, blocks = blocks)
}

#' Wald intervals from an observed information matrix
#'
#' @keywords internal
.inad_ci_from_info <- function(tab, info, level = 0.95, log_scale = TRUE, ridge = 0) {
    inner <- if (is.null(tab$boundary)) rep(TRUE, nrow(tab)) else !tab$boundary
    if (any(!inner)) {
        warning("Estimates on the boundary of the parameter space (",
                paste(tab$param[!inner], collapse = ", "),
                ") are held fixed and have no Wald interval.", call. = FALSE)
    }
    I_in <- info[inner, inner, drop = FALSE]
    if (ridge > 0) I_in <- I_in + diag(ridge, nrow(I_in))
    I_in <- (I_in + t(I_in)) / 2
    V_in <- tryCatch(chol2inv(chol(I_in)), error = function(e) NULL)
    if (is.null(V_in)) stop("Observed information is not positive definite at the fitted values.")
    V <- matrix(NA_real_, nrow(tab), nrow(tab), dimnames = list(tab$param, tab$param))
    V[inner, inner] <- V_in
    info <- I_in
    z <- stats::qnorm(1 - (1 - level) / 2)
    se <- sqrt(pmax(diag(V), 0))
    positive <- tab$group != "tau"
    use_log <- log_scale & positive & tab$est > 0
    lower <- ifelse(use_log, tab$est * exp(-z * se / tab$est), tab$est - z * se)
    upper <- ifelse(use_log, tab$est * exp(z * se / tab$est), tab$est + z * se)
    lower[positive & !use_log] <- pmax(lower[positive & !use_log], 0)
    part <- function(component) {
        k <- tab$component == component
        if (!any(k)) return(NULL)
        data.frame(param = tab$param[k], est = tab$est[k], se = se[k], lower = lower[k],
                   upper = upper[k], width = upper[k] - lower[k], level = level,
                   scale = ifelse(use_log[k], "log", "natural"),
                   row.names = NULL, stringsAsFactors = FALSE)
    }
    list(alpha = part("alpha"), theta = part("theta"), nb_inno_size = part("nb_inno_size"),
         tau = part("tau"), vcov = V,
         info_eigen_range = range(eigen(info, symmetric = TRUE, only.values = TRUE)$values))
}

#' Numerical observed-data Hessian with fixed absolute steps
#'
#' Positive parameters are differentiated on the log scale and block effects on
#' their natural scale. Steps are absolute (\code{step}, refined by Richardson
#' extrapolation), because steps proportional to the parameter value become too
#' small for log-scale parameters near zero.
#'
#' @return The observed information for the free parameters on their natural
#'   scale.
#' @keywords internal
.inad_numeric_info <- function(y, fit, blocks = NULL, step = 1e-3) {
    if (!requireNamespace("numDeriv", quietly = TRUE)) stop("numDeriv is required for method = 'hessian'")
    N <- ncol(y)
    ord <- fit$settings$order
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    fp <- .inad_free_params(fit, N, B)
    tab <- fp$tab
    use_log <- tab$group != "tau" & tab$est > 0
    eta_hat <- ifelse(use_log, log(pmax(tab$est, 1e-300)), tab$est)
    inner <- which(!tab$boundary)  # boundary estimates are held fixed
    negll <- function(delta_inner) {
        delta <- numeric(nrow(tab)); delta[inner] <- delta_inner
        eta <- eta_hat + delta
        f <- .inad_set_free(fit, fp, ifelse(use_log, exp(eta), eta), N)
        ll <- tryCatch(
            logL_inad(y = y, order = ord, thinning = fit$settings$thinning,
                      innovation = fit$settings$innovation, alpha = f$alpha, theta = f$theta,
                      nb_inno_size = f$nb_inno_size, blocks = blocks,
                      tau = if (B > 1L) f$tau else 0),
            error = function(e) -Inf
        )
        if (!is.finite(ll)) return(1e20)
        -ll
    }
    H <- numDeriv::hessian(negll, rep(0, length(inner)), method = "Richardson",
                           method.args = list(eps = step, d = step, r = 4, v = 2))
    H <- (H + t(H)) / 2
    # Back to the natural scale (the score vanishes at the MLE, so only the
    # Jacobian enters). Rows and columns of boundary estimates are NA.
    jac <- ifelse(use_log, tab$est, 1)[inner]
    info <- matrix(NA_real_, nrow(tab), nrow(tab), dimnames = list(tab$param, tab$param))
    info[inner, inner] <- H / outer(jac, jac)
    list(info = info, fp = fp)
}

#' Newton polish of an INAD fit using the exact score and information
#'
#' Takes Newton-Raphson steps on the free parameters, using the Louis score and
#' observed information, with step halving so that the log-likelihood never
#' decreases and parameters stay in their parameter space. Used to finish
#' derivative-free fits to full precision.
#'
#' @return A list with the polished \code{fit} (parameters only) and its
#'   \code{log_l}.
#' @keywords internal
.inad_newton_polish <- function(y, fit, blocks = NULL, nb_inno_size_ub = 50,
                                max_iter = 25, tol = 1e-10) {
    N <- ncol(y)
    ord <- fit$settings$order
    thinning <- fit$settings$thinning
    innovation <- fit$settings$innovation
    B <- if (!is.null(blocks) && !is.null(fit$tau) && length(fit$tau) > 1L) length(fit$tau) else 1L
    if (B == 1L) blocks <- NULL
    loglik <- function(f) {
        tryCatch(
            logL_inad(y = y, order = ord, thinning = thinning, innovation = innovation,
                      alpha = f$alpha, theta = f$theta, nb_inno_size = f$nb_inno_size,
                      blocks = blocks, tau = if (B > 1L) f$tau else 0),
            error = function(e) -Inf
        )
    }
    in_space <- function(tab, v) {
        a <- v[tab$group %in% c("a1", "a2")]
        th <- v[tab$group == "theta"]
        kp <- v[tab$group == "kappa"]
        all(a >= 0) && (thinning != "binom" || all(a < 1)) && all(a < 1 - 1e-10) &&
            all(th > 0) && all(kp > 0) && all(kp <= nb_inno_size_ub)
    }
    ll <- loglik(fit)
    for (it in seq_len(max_iter)) {
        li <- tryCatch(.louis_full_info(y, fit, blocks), error = function(e) NULL)
        if (is.null(li)) break
        tab <- li$fp$tab
        est <- tab$est
        # Active set: parameters on the boundary whose score points out of the
        # parameter space stay fixed (Karush-Kuhn-Tucker conditions).
        at_lower <- tab$group %in% c("a1", "a2", "kappa") & est <= 1e-8 & li$score < 0
        at_upper <- tab$group == "kappa" & est >= nb_inno_size_ub * (1 - 1e-8) & li$score > 0
        free <- !(at_lower | at_upper)
        if (!any(free)) break
        step <- numeric(length(est))
        sol <- tryCatch(solve(li$info[free, free, drop = FALSE], li$score[free]),
                        error = function(e) NULL)
        if (is.null(sol) || any(!is.finite(sol))) break
        step[free] <- sol
        lam <- 1
        accepted <- FALSE
        while (lam >= 1 / 1024) {
            v <- est + lam * step
            if (in_space(li$fp$tab, v)) {
                cand <- .inad_set_free(fit, li$fp, v, N)
                ll_c <- loglik(cand)
                if (is.finite(ll_c) && ll_c >= ll - 1e-9) {
                    fit <- cand; ll <- max(ll, ll_c); accepted <- TRUE
                    break
                }
            }
            lam <- lam / 2
        }
        if (!accepted || max(abs(lam * step)) < tol) break
    }
    list(fit = fit, log_l = ll)
}
