read_cattle_block <- function(path, n_time = 11L) {
  txt <- readLines(path, warn = FALSE)
  txt <- trimws(txt)
  txt <- txt[nzchar(txt)]

  rows <- strsplit(txt, "\\s+")
  mat <- t(vapply(rows, function(x) as.numeric(x), numeric(n_time)))
  storage.mode(mat) <- "double"
  mat
}

path_a <- file.path("data-raw", "external", "cattle_growth_data_Treatment_A.txt")
path_b <- file.path("data-raw", "external", "cattle_growth_data_Treatment_B.txt")

y_A <- read_cattle_block(path_a, n_time = 11L)
y_B <- read_cattle_block(path_b, n_time = 11L)

cattle_growth <- list(
  y = rbind(y_A, y_B),
  y_A = y_A,
  y_B = y_B,
  blocks = c(rep.int(1L, nrow(y_A)), rep.int(2L, nrow(y_B))),
  time = seq_len(ncol(y_A))
)

save(cattle_growth, file = file.path("data", "cattle_growth.rda"), compress = "bzip2")
