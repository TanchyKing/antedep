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
y_pct <- as.matrix(raw[, -1L, drop = FALSE])
storage.mode(y_pct) <- "double"

# Map percentages to ordered categories 1..11 using 10-point bins.
to_cat <- function(x) {
  ifelse(is.na(x), NA_integer_, as.integer(cut(
    x,
    breaks = c(seq(0, 100, by = 10), Inf),
    right = FALSE,
    include.lowest = TRUE,
    labels = FALSE
  )))
}

y_cat <- apply(y_pct, 2, to_cat)
storage.mode(y_cat) <- "integer"

is_A <- group == "A"
is_B <- group == "B"

cochlear_implant_cat <- list(
  y = y_cat,
  y_A = y_cat[is_A, , drop = FALSE],
  y_B = y_cat[is_B, , drop = FALSE],
  blocks = ifelse(is_A, 1L, 2L),
  group = group,
  n_categories = 11L,
  category_breaks = c(seq(0, 100, by = 10), Inf)
)

save(cochlear_implant_cat, file = file.path("data", "cochlear_implant_cat.rda"), compress = "bzip2")
