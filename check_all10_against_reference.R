# check_all10_against_reference.R
# Compares the all-ten-samples sensitivity analysis of GSE167459 (script 20)
# with the reference files.
# Run from the repository root after running script 20.
#
# Notes
#   - Gene counts, Cybb, Fcer1g, Btk and the other genes are compared with
#     tolerance 1e-6.
#   - The Hallmark and panel GSEA values (NES, padj) come from fgsea with
#     permutations (set.seed(1)); they are compared with relative tolerance 1e-3.
#   - Only the columns that exist in both files are compared.
tol_deg  <- 1e-6
tol_gsea <- 1e-3
all_ok <- TRUE

cat("##### Sample table (10 microglia samples)\n")
st <- "results/sensitivity_all10/GSE167459_all10_sample_table.csv"
if (file.exists(st)) {
  s <- read.csv(st, stringsAsFactors = FALSE)
  ok <- nrow(s) == 10 && all(table(s$group) == 5) && sum(s$batch == "J-A") == 4 &&
        all(table(s$group, s$batch)["Sham", ] == table(s$group, s$batch)["FPI", ])
  cat("  samples:", nrow(s), "| per group:", paste(names(table(s$group)), table(s$group), collapse = ", "),
      "| batch J-A:", sum(s$batch == "J-A"), "\n")
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
} else { cat("SKIPPED: run script 20 first\n") }

cat("\n##### Summary: primary versus all ten samples\n")
nf <- "results/sensitivity_all10/Summary_primary_vs_all10.csv"
rf <- "reference/sensitivity_all10/Summary_primary_vs_all10_reference.csv"
if (file.exists(nf) && file.exists(rf)) {
  new <- read.csv(nf, stringsAsFactors = FALSE)
  ref <- read.csv(rf, stringsAsFactors = FALSE)
  cols <- intersect(names(new), names(ref))
  cols <- setdiff(cols, "analysis")
  ok <- identical(new$analysis, ref$analysis)
  for (col in cols) {
    is_gsea <- grepl("NES|GSEA|IFNG_padj|IFNA_padj", col)
    d <- if (is_gsea) max(abs(new[[col]] - ref[[col]]) / pmax(abs(ref[[col]]), 1e-12), na.rm = TRUE)
         else max(abs(new[[col]] - ref[[col]]), na.rm = TRUE)
    lim <- if (is_gsea) tol_gsea else tol_deg
    if (!is.finite(d) || d > lim) { ok <- FALSE; cat(sprintf("  DIFFERENT: %-28s %.3g\n", col, d)) }
  }
  cat("  columns compared:", length(cols), "\n")
  print(new[, c("analysis", "meta_DEGs", "overlap_with_primary", "Cybb_padj")], row.names = FALSE)
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 20 first, or reference file missing\n")

cat("\n##### Output files\n")
for (f in c("results/sensitivity_all10/Figure_S9_GSE167459_all10_PCA.tiff",
            "results/sensitivity_all10/Table_S18_GSE167459_all10_sensitivity.csv")) {
  ok <- file.exists(f) && file.size(f) > 0
  cat(if (ok) "  found   " else "  MISSING ", f, "\n")
  if (!ok) all_ok <- FALSE
}

cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one difference, see the output above\n")
