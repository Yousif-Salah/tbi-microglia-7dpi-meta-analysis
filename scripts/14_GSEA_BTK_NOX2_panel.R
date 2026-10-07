# =============================================================================
# 14_GSEA_BTK_NOX2_panel.R
# Tests of the 14-gene BTK/NOX2 panel on the pooled 7-dpi analysis
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 14. Table 2 of the manuscript.
#
# Tests
#   (1) fgsea of the 14-gene panel on the ranked list of all analyzable genes
#   (2) one-sample t test of the pooled log2 fold changes, BTK module (8 genes)
#   (3) one-sample t test of the pooled log2 fold changes, NOX2 module (6 genes)
#
# Settings for the panel, and how they differ from script 13
#   - The ranking statistic is the same: sign(pooled log2 fold change) *
#     -log10(pooled P value), over the same 16,636 genes.
#   - Script 13 maps symbols to Entrez IDs and tests collections of gene sets
#     with minSize = 15. Here the ranking is kept by gene SYMBOL (if a symbol
#     repeated, the entry with the largest absolute statistic would be kept)
#     and the panel is one fixed set of 14 genes, tested with
#     minSize = 5. The panel is smaller than the minimum size of 15 used for
#     the genome-wide collections, so the larger minimum cannot apply to it.
#     A minimum of 15 would leave the 14-gene panel untested.
#   - Other fgsea settings are identical: maxSize = 500, eps = 0,
#     nPermSimple = 10000, SerialParam, set.seed(1).
#   - The panel is a single set, so there is no multiple-testing correction
#     within this test; the nominal P value is reported.
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv                      (script 05)
#
# Output (results/gsea_panel/)
#   BTK_NOX2_panel_per_gene.csv, BTK_NOX2_panel_fgsea.csv,
#   BTK_NOX2_module_t_tests.csv
#   results/session_info/14_GSEA_BTK_NOX2_panel_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/14_GSEA_BTK_NOX2_panel.R
# Packages: fgsea, dplyr, BiocParallel
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(fgsea)
  library(dplyr)
})

set.seed(1)

in_file  <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
out_dir  <- file.path("results", "gsea_panel")
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file))
  stop("Cannot find ", in_file, ". Run script 05 first. Current folder: ", getwd())


# 1. Load the pooled analysis -------------------------------------------------
meta_res <- read.csv(in_file, stringsAsFactors = FALSE)
stopifnot(nrow(meta_res) == 16636)


# 2. Panel definition (same in all scripts) -------------------------------------
btk_genes  <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1")
nox2_genes <- c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
panel_all  <- c(btk_genes, nox2_genes)

panel <- meta_res[meta_res$SYMBOL %in% panel_all, ]
cat("Panel genes found:", nrow(panel), "of 14\n")
stopifnot(nrow(panel) == 14)


# 3. Ranked vector (same rule as script 13) -----------------------------------
ranked <- meta_res %>%
  filter(!is.na(SYMBOL), !is.na(meta_logFC), !is.na(meta_p)) %>%
  mutate(rank_stat = sign(meta_logFC) * -log10(meta_p)) %>%
  group_by(SYMBOL) %>%
  slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
  ungroup()
ranked_vec <- setNames(ranked$rank_stat, ranked$SYMBOL)
ranked_vec <- sort(ranked_vec, decreasing = TRUE)
cat("Genes in ranked list:", length(ranked_vec), "\n")


# 4. Per-gene table (rank position out of the full list) ------------------------
panel$rank_position <- match(panel$SYMBOL, names(ranked_vec))
panel$module <- ifelse(panel$SYMBOL %in% btk_genes, "BTK", "NOX2")
cols <- intersect(c("SYMBOL", "module", "meta_logFC", "meta_p", "meta_padj", "rank_position"),
                  colnames(panel))
cat("\n--- Per-gene results (rank out of", length(ranked_vec), ") ---\n")
print(panel[order(panel$module, -panel$meta_logFC), cols], row.names = FALSE)


# 5. fgsea on the 14-gene panel -------------------------------------------------
set.seed(1)
gsea_panel <- fgsea(
  pathways    = list(BTK_NOX2_panel = panel_all),
  stats       = ranked_vec,
  minSize     = 5,
  maxSize     = 500,
  eps         = 0,
  nPermSimple = 10000,
  BPPARAM     = BiocParallel::SerialParam()
)
cat("\n--- fgsea, 14-gene panel ---\n")
print(as.data.frame(gsea_panel[, c("pathway", "size", "ES", "NES", "pval")]), row.names = FALSE)
cat("Leading edge:", paste(unlist(gsea_panel$leadingEdge), collapse = ", "), "\n")


# 6. One-sample t tests on pooled log2 fold change (H0: mean = 0) ---------------
run_t <- function(genes, label) {
  x <- panel$meta_logFC[panel$SYMBOL %in% genes]
  tt <- t.test(x, mu = 0)
  cat(sprintf("\n%s module: n = %d genes, mean log2FC = %.3f, t(%d) = %.2f, p = %.3f\n",
              label, length(x), mean(x), tt$parameter, tt$statistic, tt$p.value))
  cat(sprintf("  95%% CI of the mean: %.3f to %.3f\n", tt$conf.int[1], tt$conf.int[2]))
  data.frame(module = label, n_genes = length(x), mean_log2FC = mean(x),
             t = unname(tt$statistic), df = unname(tt$parameter), p = tt$p.value,
             ci_low = tt$conf.int[1], ci_high = tt$conf.int[2])
}
tests <- rbind(run_t(btk_genes, "BTK"), run_t(nox2_genes, "NOX2"))


# 7. Save ---------------------------------------------------------------------
write.csv(panel[order(panel$module, -panel$meta_logFC), cols],
          file.path(out_dir, "BTK_NOX2_panel_per_gene.csv"), row.names = FALSE)
write.csv(tests, file.path(out_dir, "BTK_NOX2_module_t_tests.csv"), row.names = FALSE)
gp <- as.data.frame(gsea_panel[, c("pathway", "size", "ES", "NES", "pval")])
gp$leadingEdge <- paste(unlist(gsea_panel$leadingEdge), collapse = ";")
write.csv(gp, file.path(out_dir, "BTK_NOX2_panel_fgsea.csv"), row.names = FALSE)
cat("\nSaved to", out_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "14_GSEA_BTK_NOX2_panel_sessionInfo.txt"))
