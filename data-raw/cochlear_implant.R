path <- file.path("data-raw", "external", "speech_recognition_data.txt")
raw <- utils::read.table(
  path,
  header = FALSE,
  stringsAsFactors = FALSE,
  na.strings = "NA"
)

if (ncol(raw) != 5L) {
  stop("speech_recognition_data.txt is expected to have 5 columns")
}

group <- raw[[1L]]
y <- as.matrix(raw[, -1L, drop = FALSE])
storage.mode(y) <- "double"

is_A <- group == "A"
is_B <- group == "B"

cochlear_implant <- list(
  y = y,
  y_A = y[is_A, , drop = FALSE],
  y_B = y[is_B, , drop = FALSE],
  blocks = ifelse(is_A, 1L, 2L),
  group = group,
  time = seq_len(ncol(y))
)

save(cochlear_implant, file = file.path("data", "cochlear_implant.rda"), compress = "bzip2")
