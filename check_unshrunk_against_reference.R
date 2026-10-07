# check_unshrunk_against_reference.R
# Compares the unshrunk tables (scripts 01-04) and the unshrunk pooled result
# (script 08) with the reference files. Run from the repository root.
tol <- 1e-6
ds <- list(
  GSE167459 = "ENSEMBL", GSE253476 = "ENSEMBL", GSE276647 = "Gene", GSE283560 = "ENSEMBL")
num_cols <- c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")
all_ok <- TRUE

cat("##### Per-dataset unshrunk tables\n")
for (nm in names(ds)) {
  nf <- file.path("results/deg_tables_unshrunk", paste0(nm, "_7d_unshrunk.csv"))
  rf <- file.path("reference", paste0(nm, "_7d_unshrunk_reference.csv"))
  cat("\n=== ", nm, " ===\n", sep = "")
  if (!file.exists(nf)) { cat("SKIPPED: result file not found\n"); next }
  new <- read.csv(nf, check.names = FALSE); ref <- read.csv(rf, check.names = FALSE)
  key <- ds[[nm]]
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")
  same <- setequal(new[[key]], ref[[key]]); cat("Same genes:", same, "\n")
  if (!same) { all_ok <- FALSE; next }
  ref <- ref[match(new[[key]], ref[[key]]), ]
  ok <- TRUE
  for (col in num_cols) {
    d <- suppressWarnings(max(abs(new[[col]] - ref[[col]]), na.rm = TRUE))
    na_same <- identical(is.na(new[[col]]), is.na(ref[[col]]))
    cat(sprintf("  %-15s max difference = %.3g   NA pattern identical: %s\n", col, d, na_same))
    if ((is.finite(d) && d > tol) || !na_same) ok <- FALSE
  }
  sym_same <- identical(as.character(new$SYMBOL), as.character(ref$SYMBOL))
  cat("  SYMBOL identical:", sym_same, "\n"); if (!sym_same) ok <- FALSE
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n")
  if (!ok) all_ok <- FALSE
}

cat("\n##### Unshrunk pooled result\n")
mcols <- c("k", "meta_logFC", "meta_SE", "meta_z", "meta_p", "ci_lb", "ci_ub",
           "I2", "tau2", "meta_padj", "meta_stat")
for (nm in c("Meta7dpi_all_genes_k2", "Meta7dpi_metaDEGs_FDR0.05")) {
  nf <- file.path("results/meta_unshrunk", paste0(nm, "_unshrunk.csv"))
  rf <- file.path("reference", paste0(nm, "_unshrunk_reference.csv"))
  cat("\n=== ", nm, " ===\n", sep = "")
  if (!file.exists(nf)) { cat("SKIPPED: result file not found\n"); next }
  new <- read.csv(nf); ref <- read.csv(rf)
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")
  same <- setequal(new$SYMBOL, ref$SYMBOL); cat("Same genes:", same, "\n")
  if (!same) { all_ok <- FALSE; next }
  ref <- ref[match(new$SYMBOL, ref$SYMBOL), ]
  ok <- TRUE
  for (col in mcols) {
    d <- suppressWarnings(max(abs(new[[col]] - ref[[col]]), na.rm = TRUE))
    na_same <- identical(is.na(new[[col]]), is.na(ref[[col]]))
    cat(sprintf("  %-11s max difference = %.3g   NA pattern identical: %s\n", col, d, na_same))
    if ((is.finite(d) && d > tol) || !na_same) ok <- FALSE
  }
  m_same <- identical(as.character(new$method), as.character(ref$method))
  cat("  method identical:", m_same, "\n"); if (!m_same) ok <- FALSE
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n")
  if (!ok) all_ok <- FALSE
}
cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one file failed\n")
