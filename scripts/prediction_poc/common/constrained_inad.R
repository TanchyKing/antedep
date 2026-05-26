fit_inad_alpha_constant_poc <- function(y,
                                        blocks,
                                        fit_unconstrained = NULL,
                                        thinning = "nbinom",
                                        innovation = "nbinom",
                                        max_iter = 50L,
                                        ...) {
  if (is.null(fit_unconstrained)) {
    fit_unconstrained <- antedep::fit_inad(
      y = y,
      order = 1L,
      thinning = thinning,
      innovation = innovation,
      blocks = blocks,
      max_iter = max_iter,
      ...
    )
  }

  st <- antedep::test_stationarity_inad(
    y = y,
    order = 1L,
    thinning = thinning,
    innovation = innovation,
    blocks = blocks,
    constrain = "alpha",
    fit_unconstrained = fit_unconstrained,
    max_iter = max_iter,
    ...
  )

  fit <- st$fit_constrained
  n_time <- ncol(y)
  blocks_int <- as.integer(factor(blocks))

  fit$settings <- modifyList(
    list(
      order = 1L,
      thinning = thinning,
      innovation = innovation,
      n_subjects = nrow(y),
      n_time = n_time,
      blocks = blocks_int,
      n_blocks = length(unique(blocks_int)),
      alpha_constraint = "constant",
      nb_inno_size_source = "unconstrained_fit"
    ),
    fit$settings %||% list()
  )
  fit$n_obs <- sum(!is.na(y))
  fit$n_missing <- sum(is.na(y))
  fit$convergence <- fit$convergence %||% NA_integer_
  fit$convergence_info <- fit$convergence_info %||% list(source = "test_stationarity_inad")
  class(fit) <- c("inad_fit", "ad_fit")
  attr(fit, "stationarity_test") <- st
  fit
}

`%||%` <- function(x, y) if (is.null(x)) y else x
