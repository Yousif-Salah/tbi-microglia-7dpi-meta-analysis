# check_qc_against_reference.R
# Compares the sample QC output of script 04b with the reference files.
# Run from the repository root after running script 04b.
#
# Notes
#   - The sign of a principal component can flip between computers, so PC1 and
#     PC2 are compared as absolute values.
#   - Tolerance is 1e-4 for the PCA and correlation values (they come from
#     vst and prcomp, which can differ in the last digits between R versions).
tol <- 1e-4
all_ok <- TRUE

cat("##### QC sample table\n")
nf <- "results/qc/QC_sample_table.csv"
rf <- "reference/qc/QC_sample_table_reference.csv"
if (file.exists(nf) && file.exists(rf)) {
  new <- read.csv(nf, stringsAsFactors = FALSE)
  ref <- read.csv(rf, stringsAsFactors = FALSE)
  ok <- identical(names(new), names(ref)) && identical(dim(new), dim(ref)) &&
        identical(paste(new$dataset, new$sample), paste(ref$dataset, ref$sample))
  if (ok) {
    ok <- identical(new$condition, ref$condition) &&
          identical(new$flagged_in_script16, ref$flagged_in_script16) &&
          all(new$library_size == ref$library_size)
    d1 <- max(abs(abs(new$PC1) - abs(ref$PC1)))
    d2 <- max(abs(abs(new$PC2) - abs(ref$PC2)))
    d3 <- max(abs(new$mean_within_group_spearman - ref$mean_within_group_spearman))
    cat(sprintf("  max difference |PC1| = %.3g, |PC2| = %.3g, within-group correlation = %.3g\n",
                d1, d2, d3))
    ok <- ok && all(c(d1, d2, d3) < tol)
  }
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 04b first, or reference file missing\n")

cat("\n##### QC summary (PC1 and PC2 percent variance, flagged sample)\n")
nf <- "results/qc/QC_summary.csv"
rf <- "reference/qc/QC_summary_reference.csv"
if (file.exists(nf) && file.exists(rf)) {
  new <- read.csv(nf, stringsAsFactors = FALSE)
  ref <- read.csv(rf, stringsAsFactors = FALSE)
  print(new[, c("dataset", "n_samples", "n_genes", "PC1_percent", "PC2_percent",
                "flagged_sample", "lowest_within_group_correlation")], row.names = FALSE)
  ok <- isTRUE(all.equal(new, ref, check.attributes = FALSE))
  cat(if (ok) "RESULT: PASS (identical to the reference)\n" else "RESULT: FAIL\n")
  if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 04b first, or reference file missing\n")

cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one difference, see the output above\n")
