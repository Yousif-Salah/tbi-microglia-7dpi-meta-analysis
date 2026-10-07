# =============================================================================
# 08_Sensitivity_unshrunk.R
# Sensitivity analysis: pooled result with unshrunk instead of shrunk fold changes
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 8. Repeats the pooling of script 05 with the unshrunk (maximum
#           likelihood) fold changes of all four datasets and compares the
#           meta-DEG lists of the two analyses.
#
# Method
#   - Same samples, count filter (rowSums >= 10), DESeq2 model and pooling
#     (random effects, REML, k >= 2, BH) as the main analysis. Only the fold
#     change differs: not shrunk with apeglm.
#   - The unshrunk tables are written by scripts 01 to 04
#     (results/deg_tables_unshrunk/). Pooling is done by script 05 in its
#     "unshrunk" mode, which this script switches on.
#   - Comparison: number of analyzable genes and meta-DEGs, overlap of the two
#     meta-DEG lists, agreement of direction, and the BTK/NOX2 panel genes.
#
# Input
#   results/deg_tables_unshrunk/*_7d_unshrunk.csv     (scripts 01-04)
#   results/meta/Meta7dpi_all_genes_k2.csv            (script 05, main analysis)
#   results/meta/Meta7dpi_metaDEGs_FDR0.05.csv        (script 05, main analysis)
#
# Output
#   results/meta_unshrunk/Meta7dpi_all_genes_k2_unshrunk.csv
#   results/meta_unshrunk/Meta7dpi_metaDEGs_FDR0.05_unshrunk.csv
#   results/meta_unshrunk/shrinkage_comparison_summary.txt
#   results/meta_unshrunk/metaDEGs_only_in_main.csv
#   results/meta_unshrunk/metaDEGs_only_in_unshrunk.csv
#   results/session_info/08_Sensitivity_unshrunk_sessionInfo.txt
#
# Usage   : run from the repository root, after scripts 01 to 05, for example
#           Rscript scripts/08_Sensitivity_unshrunk.R
# Packages: dplyr, purrr, tibble, metafor (used by script 05)
# =============================================================================


# 0. Paths and checks ---------------------------------------------------------
main_all  <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
main_deg  <- file.path("results", "meta", "Meta7dpi_metaDEGs_FDR0.05.csv")
unsh_dir  <- file.path("results", "deg_tables_unshrunk")
out_dir   <- file.path("results", "meta_unshrunk")
info_dir  <- file.path("results", "session_info")

for (f in c(main_all, main_deg))
  if (!file.exists(f))
    stop("File not found: ", f, ". Run scripts 01 to 05 first. Current folder: ", getwd())
if (!dir.exists(unsh_dir) || length(list.files(unsh_dir, pattern = "_7d_unshrunk\\.csv$")) != 4)
  stop("Expected 4 unshrunk tables in ", unsh_dir,
       ". Run scripts 01 to 04 first. Current folder: ", getwd())


# 1. Pool the unshrunk results (script 05 in "unshrunk" mode) ------------------
meta_arm <- "unshrunk"
source(file.path("scripts", "05_meta_analysis.R"))
rm(meta_arm)   # so that a later run of script 05 uses the main analysis again


# 2. Compare the two analyses ---------------------------------------------------
suppressPackageStartupMessages(library(dplyr))

main_res <- read.csv(main_all)
main_sig <- read.csv(main_deg)
unsh_res <- read.csv(file.path(out_dir, "Meta7dpi_all_genes_k2_unshrunk.csv"))
unsh_sig <- read.csv(file.path(out_dir, "Meta7dpi_metaDEGs_FDR0.05_unshrunk.csv"))

both     <- intersect(main_sig$SYMBOL, unsh_sig$SYMBOL)
only_m   <- setdiff(main_sig$SYMBOL, unsh_sig$SYMBOL)
only_u   <- setdiff(unsh_sig$SYMBOL, main_sig$SYMBOL)

dir_main <- sign(main_sig$meta_logFC[match(both, main_sig$SYMBOL)])
dir_unsh <- sign(unsh_sig$meta_logFC[match(both, unsh_sig$SYMBOL)])
same_dir <- sum(dir_main == dir_unsh)

summary_lines <- c(
  sprintf("Analyzable genes (k >= 2)     main: %d | unshrunk: %d", nrow(main_res), nrow(unsh_res)),
  sprintf("Meta-DEGs                     main: %d (up %d, down %d) | unshrunk: %d (up %d, down %d)",
          nrow(main_sig), sum(main_sig$meta_logFC > 0), sum(main_sig$meta_logFC < 0),
          nrow(unsh_sig), sum(unsh_sig$meta_logFC > 0), sum(unsh_sig$meta_logFC < 0)),
  sprintf("Meta-DEGs in both analyses    %d of %d main meta-DEGs", length(both), nrow(main_sig)),
  sprintf("Same direction in both        %d of %d", same_dir, length(both)),
  sprintf("Only in main                  %d", length(only_m)),
  sprintf("Only in unshrunk              %d", length(only_u))
)
cat("\n", paste(summary_lines, collapse = "\n"), "\n", sep = "")
writeLines(summary_lines, file.path(out_dir, "shrinkage_comparison_summary.txt"))

write.csv(main_sig[main_sig$SYMBOL %in% only_m, ],
          file.path(out_dir, "metaDEGs_only_in_main.csv"), row.names = FALSE)
write.csv(unsh_sig[unsh_sig$SYMBOL %in% only_u, ],
          file.path(out_dir, "metaDEGs_only_in_unshrunk.csv"), row.names = FALSE)

cat("\nBTK/NOX2 panel genes, main versus unshrunk analysis:\n")
panel <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1",
           "Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
panel_tab <- merge(
  main_res[main_res$SYMBOL %in% panel, c("SYMBOL", "k", "meta_logFC", "meta_padj")],
  unsh_res[unsh_res$SYMBOL %in% panel, c("SYMBOL", "meta_logFC", "meta_padj")],
  by = "SYMBOL", all = TRUE, suffixes = c("_main", "_unshrunk"))
print(panel_tab, row.names = FALSE)

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "08_Sensitivity_unshrunk_sessionInfo.txt"))
