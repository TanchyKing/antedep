path <- file.path('data-raw', 'external', 'labor_force_table1.csv')
raw <- utils::read.csv(path, stringsAsFactors = FALSE)

required_cols <- c('Y1', 'Y2', 'Y3', 'Y4', 'Y5', 'Count')
if (!all(required_cols %in% names(raw))) {
  stop('labor_force_table1.csv must contain columns: Y1,Y2,Y3,Y4,Y5,Count')
}

for (nm in required_cols) {
  raw[[nm]] <- as.integer(raw[[nm]])
}

if (anyNA(raw)) {
  stop('labor_force_table1.csv contains NA values; expected complete table counts')
}
if (any(raw$Count < 0L)) {
  stop('Count must be nonnegative')
}

seq_cols <- c('Y1', 'Y2', 'Y3', 'Y4', 'Y5')
if (any(raw[seq_cols] < 1L | raw[seq_cols] > 2L)) {
  stop('Y1..Y5 must be coded as 1/2 categories')
}

row_index <- rep.int(seq_len(nrow(raw)), raw$Count)
y <- as.matrix(raw[row_index, seq_cols, drop = FALSE])
storage.mode(y) <- 'integer'

labor_force_cat <- list(
  y = y,
  counts = raw,
  n_categories = 2L,
  time = 1967:1971,
  status_labels = c('employed', 'unemployed')
)

save(labor_force_cat, file = file.path('data', 'labor_force_cat.rda'), compress = 'bzip2')
