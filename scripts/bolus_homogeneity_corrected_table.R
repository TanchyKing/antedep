## Regenerate Table 4.2 of the bolus homogeneity analysis with the fixed
## .count_params_homogeneity() (which now counts the NB innovation size
## parameter honestly from the fit).

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
})

data("bolus_inad")
y      <- bolus_inad$y
blocks <- bolus_inad$blocks
n_subj <- nrow(y); N <- ncol(y)

families <- list(
  list(name = "NBT-NBI", thin = "nbinom", inno = "nbinom"),
  list(name = "PT-NBI",  thin = "pois",   inno = "nbinom"),
  list(name = "BT-NBI",  thin = "binom",  inno = "nbinom"),
  list(name = "NBT-BI",  thin = "nbinom", inno = "bell"),
  list(name = "PT-BI",   thin = "pois",   inno = "bell"),
  list(name = "BT-BI",   thin = "binom",  inno = "bell"),
  list(name = "NBT-PI",  thin = "nbinom", inno = "pois"),
  list(name = "PT-PI",   thin = "pois",   inno = "pois"),
  list(name = "BT-PI",   thin = "binom",  inno = "pois")
)

rows <- list()
tests_rows <- list()
for (f in families) {
  cat(">>> ", f$name, "\n", sep = "")
  set.seed(1)
  res <- run_homogeneity_tests_inad(
    y, blocks, order = 1,
    thinning = f$thin, innovation = f$inno,
    max_iter = 200, tol = 1e-8
  )
  mt <- res$model_table
  rows[[f$name]] <- data.frame(
    Family = f$name,
    M1     = mt$BIC[mt$model == "M1 (Pooled)"],
    M2     = mt$BIC[mt$model == "M2 (INADFE)"],
    M3     = mt$BIC[mt$model == "M3 (Heterogeneous)"],
    stringsAsFactors = FALSE
  )
  tm <- res$test_mean
  td <- res$test_dependence
  tests_rows[[f$name]] <- data.frame(
    Family       = f$name,
    LRT_M1_M2    = tm$lrt_stat,
    df_M1_M2     = tm$df,
    p_M1_M2      = tm$p_value,
    LRT_M2_M3    = td$lrt_stat,
    df_M2_M3     = td$df,
    p_M2_M3      = td$p_value,
    stringsAsFactors = FALSE
  )
}

bic_table <- do.call(rbind, rows)
bic_table$MinBIC <- pmin(bic_table$M1, bic_table$M2, bic_table$M3)
bic_table <- bic_table[order(bic_table$MinBIC), ]

test_table <- do.call(rbind, tests_rows)
test_table <- test_table[match(bic_table$Family, test_table$Family), ]

cat("\n==== Corrected Table 4.2 (ordered by smallest BIC) ====\n")
print(bic_table, row.names = FALSE, digits = 7)

cat("\n==== Homogeneity LRTs (M1 vs M2, M2 vs M3) ====\n")
print(test_table, row.names = FALSE, digits = 5)

saveRDS(list(bic_table = bic_table, test_table = test_table),
        "scripts/bolus_homogeneity_corrected.rds")
