# check_lodo_against_reference.R
# Compares the LODO results (scripts 10 and 11) with the reference files.
# Reference for script 11 = the harmonized run (same input rules as script 05).
# Run from the repository root.
tol <- 1e-6
all_ok <- TRUE

cat("##### Cybb and Fcer1g LODO\n")
nf <- "results/lodo/LODO_Cybb_Fcer1g.csv"
if (file.exists(nf)) {
  new <- read.csv(nf); ref <- read.csv("reference/LODO_Cybb_Fcer1g_reference.csv")
  ref <- ref[match(paste(new$dropped_dataset, new$gene), paste(ref$dropped_dataset, ref$gene)), ]
  ok <- nrow(new) == 8 && !anyNA(ref$gene)
  for (col in c("n_datasets", "pooled_log2FC", "pval", "I2", "tau2")) {
    d <- max(abs(new[[col]] - ref[[col]]), na.rm = TRUE)
    cat(sprintf("  %-14s max difference = %.3g\n", col, d))
    if (!is.finite(d) || d > tol) ok <- FALSE
  }
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 10 first\n")

cat("\n##### Genome-wide LODO\n")
sf <- "results/lodo/LODO_GenomeWide_summary.csv"
if (file.exists(sf)) {
  new <- read.csv(sf); ref <- read.csv("reference/LODO_GenomeWide_summary_reference.csv")
  cmp <- merge(new, ref, by = "dropped_dataset", suffixes = c("_new", "_reference"))
  print(cmp[, c("dropped_dataset", "n_genes_tested_new", "n_genes_tested_reference",
                "n_meta_DEGs_new", "n_meta_DEGs_reference", "n_up_new", "n_up_reference",
                "n_down_new", "n_down_reference")], row.names = FALSE)
  same <- isTRUE(all.equal(new[order(new$dropped_dataset), ], ref[order(ref$dropped_dataset), ],
                           check.attributes = FALSE))
  cat("Summary identical to the reference:", same, "\n")
  if (!same) all_ok <- FALSE
} else cat("SKIPPED: run script 11 first\n")

cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one difference, see the output above\n")
