# =============================================================================
# check_against_reference.R
# Compares each per-dataset result produced by the deposited scripts with a
# reference result file. Run from the repository root after running the scripts.
#
#   Rscript check_against_reference.R
#
# A dataset passes if it has the same genes and every numeric column differs by
# less than 1e-6. Datasets whose result file is not present yet are skipped.
# =============================================================================

tol <- 1e-6

# dataset = c(new file, reference file, key column; "rownames" or a column name)
checks <- list(
  GSE167459 = list(new = "results/deg_tables/GSE167459_7d.csv",
                   ref = "reference/GSE167459_7d_reference.csv",
                   key = "rownames"),
  GSE253476 = list(new = "results/deg_tables/GSE253476_7d.csv",
                   ref = "reference/GSE253476_7d_reference.csv",
                   key = "ENSEMBL"),
  GSE276647 = list(new = "results/deg_tables/GSE276647_7d.csv",
                   ref = "reference/GSE276647_7d_reference.csv",
                   key = "Gene"),
  GSE283560 = list(new = "results/deg_tables/GSE283560_7d.csv",
                   ref = "reference/GSE283560_7d_reference.csv",
                   key = "ENSEMBL")
)

read_keyed <- function(path, key) {
  if (key == "rownames") {
    df <- read.csv(path, row.names = 1, check.names = FALSE)
  } else {
    df <- read.csv(path, check.names = FALSE)
    rownames(df) <- df[[key]]
  }
  df
}

num_cols <- c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
all_ok <- TRUE

for (ds in names(checks)) {
  cf <- checks[[ds]]
  cat("\n=== ", ds, " ===\n", sep = "")
  if (!file.exists(cf$new)) { cat("SKIPPED: result file not found yet\n"); next }
  if (!file.exists(cf$ref)) { cat("SKIPPED: reference file not found\n");  next }

  new <- read_keyed(cf$new, cf$key)
  ref <- read_keyed(cf$ref, cf$key)
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")

  same_genes <- setequal(rownames(new), rownames(ref))
  cat("Same genes:", same_genes, "\n")
  if (!same_genes) { all_ok <- FALSE; next }

  ref <- ref[rownames(new), ]
  ok <- TRUE
  for (col in intersect(num_cols, colnames(new))) {
    d <- max(abs(new[[col]] - ref[[col]]), na.rm = TRUE)
    na_same <- identical(is.na(new[[col]]), is.na(ref[[col]]))
    cat(sprintf("  %-15s max difference = %.3g   NA pattern identical: %s\n",
                col, d, na_same))
    if (!is.finite(d) || d > tol || !na_same) ok <- FALSE
  }
  sym_same <- identical(as.character(new$SYMBOL), as.character(ref$SYMBOL))
  cat("  SYMBOL identical:", sym_same, "\n")
  if (!sym_same) ok <- FALSE

  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n")
  if (!ok) all_ok <- FALSE
}

cat("\nOverall:", if (all_ok) "all checked datasets pass\n" else "at least one dataset failed\n")
