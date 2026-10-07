# check_meta_against_reference.R
# Compares the pooled results of 05_meta_analysis.R with the reference files.
# Run from the repository root after running scripts 01 to 05.
tol <- 1e-6
pairs <- list(
  all_genes = c("results/meta/Meta7dpi_all_genes_k2.csv",
                "reference/Meta7dpi_all_genes_k2_reference.csv"),
  meta_degs = c("results/meta/Meta7dpi_metaDEGs_FDR0.05.csv",
                "reference/Meta7dpi_metaDEGs_FDR0.05_reference.csv"))
num_cols <- c("k", "meta_logFC", "meta_SE", "meta_z", "meta_p", "ci_lb",
              "ci_ub", "I2", "tau2", "meta_padj", "meta_stat")
all_ok <- TRUE
for (nm in names(pairs)) {
  cat("\n=== ", nm, " ===\n", sep = "")
  if (!file.exists(pairs[[nm]][1])) { cat("SKIPPED: result file not found\n"); next }
  new <- read.csv(pairs[[nm]][1]); ref <- read.csv(pairs[[nm]][2])
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")
  same <- setequal(new$SYMBOL, ref$SYMBOL); cat("Same genes:", same, "\n")
  if (!same) { all_ok <- FALSE; next }
  ref <- ref[match(new$SYMBOL, ref$SYMBOL), ]
  ok <- TRUE
  for (col in num_cols) {
    d <- suppressWarnings(max(abs(new[[col]] - ref[[col]]), na.rm = TRUE))
    na_same <- identical(is.na(new[[col]]), is.na(ref[[col]]))
    cat(sprintf("  %-11s max difference = %.3g   NA pattern identical: %s\n", col, d, na_same))
    if ((is.finite(d) && d > tol) || !na_same) ok <- FALSE
  }
  m_same <- identical(as.character(new$method), as.character(ref$method))
  e_same <- identical(as.character(new$ENSEMBL), as.character(ref$ENSEMBL))
  cat("  method identical:", m_same, "| ENSEMBL identical:", e_same, "\n")
  if (!m_same || !e_same) ok <- FALSE
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n")
  if (!ok) all_ok <- FALSE
}
cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one file failed\n")
