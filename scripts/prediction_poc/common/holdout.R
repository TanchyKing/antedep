ppc_stratified_split <- function(blocks, train_prop = 0.70, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  blocks <- as.integer(blocks)
  train <- logical(length(blocks))

  for (b in sort(unique(blocks))) {
    idx <- which(blocks == b)
    n_train <- max(1L, floor(length(idx) * train_prop))
    n_train <- min(n_train, length(idx) - 1L)
    train[sample(idx, n_train)] <- TRUE
  }

  list(train = which(train), test = which(!train))
}

ppc_long_lag1 <- function(y, blocks, subjects = seq_len(nrow(y))) {
  y <- as.matrix(y)
  out <- vector("list", max(ncol(y) - 1L, 0L))
  pos <- 0L
  for (tt in 2:ncol(y)) {
    pos <- pos + 1L
    out[[pos]] <- data.frame(
      subject = subjects,
      group = factor(blocks),
      time = tt,
      y = as.numeric(y[, tt]),
      y_lag1 = as.numeric(y[, tt - 1L])
    )
  }
  do.call(rbind, out)
}
