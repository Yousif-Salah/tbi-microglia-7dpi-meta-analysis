# check_k4_against_reference.R
# Compares the k = 4 sensitivity results (script 09) with the reference files.
# Run from the repository root after running script 09.
tol <- 1e-9
nf <- "results/sensitivity_k4/k4_sensitivity_analysis.csv"
bf <- "results/sensitivity_k4/k_breakdown_summary.csv"
if (!file.exists(nf)) stop("Run script 09 first.")
new <- read.csv(nf); ref <- read.csv("reference/k4_sensitivity_analysis_reference.csv")
cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")
same <- setequal(new$SYMBOL, ref$SYMBOL); cat("Same genes:", same, "\n")
ok <- same
if (same) {
  ref <- ref[match(new$SYMBOL, ref$SYMBOL), ]
  for (col in c("k", "meta_logFC", "meta_SE", "meta_z", "meta_p", "ci_lb", "ci_ub",
                "I2", "tau2", "meta_padj", "meta_stat", "meta_padj_k4")) {
    d <- suppressWarnings(max(abs(new[[col]] - ref[[col]]), na.rm = TRUE))
    na_same <- identical(is.na(new[[col]]), is.na(ref[[col]]))
    cat(sprintf("  %-13s max difference = %.3g   NA pattern identical: %s\n", col, d, na_same))
    if ((is.finite(d) && d > tol) || !na_same) ok <- FALSE
  }
}
bnew <- read.csv(bf); bref <- read.csv("reference/k_breakdown_summary_reference.csv")
b_ok <- isTRUE(all.equal(bnew, bref))
cat("k breakdown table identical:", b_ok, "\n")
cat("\nOverall:", if (ok && b_ok) "PASS\n" else "FAIL, see the differences printed above\n")
