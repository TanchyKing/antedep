path <- file.path("data-raw", "external", "100km_race_combined_extracted.csv")
raw <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("Subject", "Age", paste0("Section_", 1:10))
if (!all(required_cols %in% names(raw))) {
  stop("100km_race_combined_extracted.csv must contain Subject, Age, Section_1..Section_10")
}

raw$Subject <- as.integer(raw$Subject)
raw$Age <- as.numeric(raw$Age)

section_cols <- paste0("Section_", 1:10)
y <- as.matrix(raw[, section_cols, drop = FALSE])
storage.mode(y) <- "double"

race_100km <- list(
  y = y,
  age = raw$Age,
  subject = raw$Subject,
  time = seq_len(ncol(y)),
  section_labels = colnames(y)
)

save(race_100km, file = file.path("data", "race_100km.rda"), compress = "bzip2")
