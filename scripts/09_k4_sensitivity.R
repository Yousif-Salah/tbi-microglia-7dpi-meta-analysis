# =============================================================================
# 09_k4_sensitivity.R
# Dataset-representation sensitivity analysis: genes measured in all 4 datasets
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 9. Tests whether the result depends on genes that were present
#           in only 2 or 3 of the 4 datasets.
#
# Method
#   - Start from the pooled result of script 05 (genes present in at least 2
#     datasets). k is the number of datasets in which a gene was measured.
#   - Count genes by k (2, 3, 4).
#   - Keep only genes with k = 4 and repeat the Benjamini-Hochberg correction
#     within this subset (meta_padj_k4). The pooled estimates and P values are
#     not recalculated.
#   - Report the number of meta-DEGs in the k = 4 subset (meta_padj_k4 < 0.05)
#     and the values of the BTK/NOX2 panel genes Cybb, Fcer1g and Btk.
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv           (script 05)
#
# Output (results/sensitivity_k4/)
#   k_breakdown_summary.csv         number and percent of genes for each k
#   k4_sensitivity_analysis.csv     all k = 4 genes with meta_padj_k4
#   results/session_info/09_k4_sensitivity_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/09_k4_sensitivity.R
# Packages: base R only
# =============================================================================


# 0. Paths and input ----------------------------------------------------------
meta_file <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
out_dir   <- file.path("results", "sensitivity_k4")
info_dir  <- file.path("results", "session_info")

for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(meta_file))
  stop("Pooled result not found: ", meta_file,
       ". Run script 05 first. Current folder: ", getwd())

meta_res <- read.csv(meta_file)
stopifnot(all(c("SYMBOL", "k", "meta_p", "meta_padj", "meta_logFC") %in% names(meta_res)))
cat("Analyzable genes:", nrow(meta_res), "\n")


# 1. Genes by number of datasets ----------------------------------------------
k_breakdown <- as.data.frame(table(meta_res$k))
colnames(k_breakdown) <- c("k", "n_genes")
k_breakdown$pct <- round(100 * k_breakdown$n_genes / nrow(meta_res), 1)
cat("\nGenes by k:\n")
print(k_breakdown, row.names = FALSE)


# 2. Genes present in all 4 datasets, BH correction within this subset --------
k4_only <- subset(meta_res, k == 4)
cat("\nk = 4 subset size:", nrow(k4_only), "\n")
stopifnot(nrow(k4_only) > 0)

k4_only$meta_padj_k4 <- p.adjust(k4_only$meta_p, method = "BH")

cat("\nCybb, Fcer1g and Btk with BH correction inside the k = 4 subset:\n")
print(k4_only[k4_only$SYMBOL %in% c("Cybb", "Fcer1g", "Btk"),
              c("SYMBOL", "meta_logFC", "meta_p", "meta_padj", "meta_padj_k4", "I2")],
      row.names = FALSE)

cat("\nMeta-DEGs (BH < 0.05) in the k = 4 subset:", sum(k4_only$meta_padj_k4 < 0.05), "\n")


# 3. Save ---------------------------------------------------------------------
write.csv(k_breakdown, file.path(out_dir, "k_breakdown_summary.csv"), row.names = FALSE)
write.csv(k4_only,     file.path(out_dir, "k4_sensitivity_analysis.csv"), row.names = FALSE)
cat("\nSaved to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "09_k4_sensitivity_sessionInfo.txt"))
