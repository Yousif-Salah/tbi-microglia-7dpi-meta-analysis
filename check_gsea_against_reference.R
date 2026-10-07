# check_gsea_against_reference.R
# Compares the GSEA results (scripts 13 and 14) with the reference files.
# Run from the repository root.
tol <- 1e-6
all_ok <- TRUE

cmp_tab <- function(new_file, ref_file, key, num_cols, label) {
  cat("\n=== ", label, " ===\n", sep = "")
  if (!file.exists(new_file)) { cat("SKIPPED: file not found, run the script first\n"); return(NA) }
  new <- read.csv(new_file, stringsAsFactors = FALSE)
  ref <- read.csv(ref_file, stringsAsFactors = FALSE)
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "\n")
  same_keys <- setequal(new[[key]], ref[[key]])
  cat("Same terms:", same_keys, "\n")
  if (!same_keys) return(FALSE)
  ref <- ref[match(new[[key]], ref[[key]]), ]
  ok <- TRUE
  for (col in num_cols) {
    a <- new[[col]]; b <- ref[[col]]
    d <- max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE)   # relative difference
    na_same <- identical(is.na(a), is.na(b))
    cat(sprintf("  %-8s max relative difference = %.3g   NA pattern identical: %s\n", col, d, na_same))
    if (!is.finite(d) || d > tol || !na_same) ok <- FALSE
  }
  if ("leadingEdge" %in% names(new)) {
    le <- identical(as.character(new$leadingEdge), as.character(ref$leadingEdge))
    cat("  leadingEdge identical:", le, "\n"); if (!le) ok <- FALSE
  }
  cat("RESULT:", if (ok) "PASS" else "FAIL, see the differences printed above", "\n")
  ok
}

res <- list(
  cmp_tab("results/gsea/GSEA_GOBP.csv",     "reference/gsea/GSEA_GOBP_reference.csv",
          "pathway", c("pval", "padj", "ES", "NES", "size"), "GO Biological Process"),
  cmp_tab("results/gsea/GSEA_Hallmark.csv", "reference/gsea/GSEA_Hallmark_reference.csv",
          "pathway", c("pval", "padj", "ES", "NES", "size"), "Hallmark"),
  cmp_tab("results/gsea/GSEA_Reactome.csv", "reference/gsea/GSEA_Reactome_reference.csv",
          "pathway", c("pval", "padj", "ES", "NES", "size"), "Reactome"),
  cmp_tab("results/gsea_panel/BTK_NOX2_panel_per_gene.csv", "reference/gsea/BTK_NOX2_panel_per_gene_reference.csv",
          "SYMBOL", c("meta_logFC", "meta_p", "meta_padj", "rank_position"), "Panel per gene"),
  cmp_tab("results/gsea_panel/BTK_NOX2_panel_fgsea.csv", "reference/gsea/BTK_NOX2_panel_fgsea_reference.csv",
          "pathway", c("size", "ES", "NES", "pval"), "Panel fgsea"),
  cmp_tab("results/gsea_panel/BTK_NOX2_module_t_tests.csv", "reference/gsea/BTK_NOX2_module_t_tests_reference.csv",
          "module", c("n_genes", "mean_log2FC", "t", "df", "p", "ci_low", "ci_high"), "Module t tests")
)
res <- unlist(res)
cat("\nOverall:", if (all(res, na.rm = TRUE)) "all checked files pass\n" else "at least one difference, see the output above\n")
