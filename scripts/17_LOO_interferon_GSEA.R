# =============================================================================
# 17_LOO_interferon_GSEA.R
# Interferon signatures and homeostatic markers after removing the flagged samples
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 17. After the flagged sample of each dataset is removed (the
#           four leave-one-sample-out runs of script 16), this script checks
#             (1) that the Hallmark Interferon Gamma Response and Interferon
#                 Alpha Response sets stay enriched, and
#             (2) that Cx3cr1, Sall1 and Tjp1 keep their direction.
#           The full main analysis is run too, as the reference.
#
# Settings: the same as script 13. Ranking = sign(log2 fold change) *
#           -log10(P value); fgsea with minSize = 15, maxSize = 500, eps = 0,
#           nPermSimple = 10000, one worker, set.seed(1); the full Hallmark
#           collection is tested, so adjusted P values are comparable with the
#           main analysis.
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv                       (script 05)
#   results/loo/LOO_pooled_<dataset>.csv  (four files)           (script 16)
#
# Output (results/loo/)
#   LOO_interferon_GSEA.csv        NES and adjusted P of the two interferon sets
#   LOO_Cx3cr1_Sall1_Tjp1.csv      pooled estimate of the three genes
#   results/session_info/17_LOO_interferon_GSEA_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/17_LOO_interferon_GSEA.R
# Packages: fgsea, msigdbr (25.1.1), clusterProfiler, org.Mm.eg.db, dplyr
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
})

set.seed(1)  # fixed seed for reproducible NES and P values

out_dir  <- file.path("results", "loo")
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

run_files <- c(
  baseline  = file.path("results", "meta", "Meta7dpi_all_genes_k2.csv"),
  GSE167459 = file.path(out_dir, "LOO_pooled_GSE167459.csv"),
  GSE253476 = file.path(out_dir, "LOO_pooled_GSE253476.csv"),
  GSE276647 = file.path(out_dir, "LOO_pooled_GSE276647.csv"),
  GSE283560 = file.path(out_dir, "LOO_pooled_GSE283560.csv")
)
missing_files <- run_files[!file.exists(run_files)]
if (length(missing_files) > 0)
  stop("Missing input file(s):\n", paste(missing_files, collapse = "\n"),
       "\nRun scripts 05 and 16 first. Current folder: ", getwd())

genes_to_check <- c("Cx3cr1", "Sall1", "Tjp1")
ifn_sets <- c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE")


# 1. Hallmark gene sets (mouse, same call as script 13) --------------------------
gs_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
hallmark_list <- split(as.character(gs_hallmark$ncbi_gene), gs_hallmark$gs_name)
stopifnot(all(ifn_sets %in% names(hallmark_list)))


# 2. Ranking helper (same steps as script 13) --------------------------------------
build_ranks <- function(meta_res) {
  stopifnot(all(c("SYMBOL", "meta_logFC", "meta_p") %in% names(meta_res)))
  meta_res <- meta_res[!is.na(meta_res$meta_logFC) & !is.na(meta_res$meta_p), ]
  meta_res$rank_stat <- sign(meta_res$meta_logFC) * -log10(meta_res$meta_p)

  n_inf <- sum(is.infinite(meta_res$rank_stat))
  if (n_inf > 0) {
    smallest_nonzero_p <- min(meta_res$meta_p[meta_res$meta_p > 0], na.rm = TRUE)
    inf_idx <- is.infinite(meta_res$rank_stat)
    meta_res$rank_stat[inf_idx] <- sign(meta_res$meta_logFC[inf_idx]) * -log10(smallest_nonzero_p)
  }
  stopifnot(!any(is.na(meta_res$rank_stat)))

  id_map <- bitr(meta_res$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  mapped <- merge(meta_res, id_map, by = "SYMBOL") %>%
    group_by(ENTREZID) %>%
    slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup()

  ranked_vec <- setNames(mapped$rank_stat, mapped$ENTREZID)
  ranked_vec <- sort(ranked_vec, decreasing = TRUE)
  stopifnot(!any(duplicated(names(ranked_vec))))
  list(ranks = ranked_vec, n_inf = n_inf)
}


# 3. Run the main analysis and the four leave-one-sample-out analyses ---------------
gsea_rows <- list()
gene_rows <- list()

for (run in names(run_files)) {
  cat("\n=== Run:", run, "===\n")
  df <- read.csv(run_files[[run]], stringsAsFactors = FALSE)
  cat("Genes in file:", nrow(df), "\n")

  if (run == "baseline") {
    stopifnot(nrow(df) == 16636)
    stopifnot(abs(df$meta_logFC[df$SYMBOL == "Cybb"] - 1.66) < 0.01)
    cat("Baseline verified: 16,636 genes, Cybb log2FC = 1.66.\n")
  }
  if (!"meta_padj" %in% names(df)) df$meta_padj <- p.adjust(df$meta_p, method = "BH")
  if (!"k" %in% names(df)) df$k <- NA_integer_

  cat("Global meta-DEG count (BH < 0.05):", sum(df$meta_padj < 0.05, na.rm = TRUE), "\n")

  rk <- build_ranks(df)
  cat("Ranked vector length:", length(rk$ranks), "| infinite ranking statistics:", rk$n_inf, "\n")

  res <- fgsea(pathways = hallmark_list, stats = rk$ranks,
               minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000,
               BPPARAM = BiocParallel::SerialParam())
  cat("Hallmark sets tested:", nrow(res), "\n")

  ifn <- as.data.frame(res[res$pathway %in% ifn_sets, c("pathway", "size", "NES", "pval", "padj")])
  stopifnot(nrow(ifn) == 2)
  ifn$run <- run
  gsea_rows[[run]] <- ifn
  print(ifn, row.names = FALSE)

  g <- df[df$SYMBOL %in% genes_to_check, c("SYMBOL", "meta_logFC", "meta_p", "meta_padj", "k")]
  g$run <- run
  gene_rows[[run]] <- g
  if (nrow(g) < length(genes_to_check))
    cat("WARNING: not all of", paste(genes_to_check, collapse = ", "), "are in this run.\n")
  print(g, row.names = FALSE)
}

gsea_tab <- bind_rows(gsea_rows)
gene_tab <- bind_rows(gene_rows)
gsea_tab$run <- factor(gsea_tab$run, levels = names(run_files))
gene_tab$run <- factor(gene_tab$run, levels = names(run_files))

write.csv(gsea_tab, file.path(out_dir, "LOO_interferon_GSEA.csv"), row.names = FALSE)
write.csv(gene_tab, file.path(out_dir, "LOO_Cx3cr1_Sall1_Tjp1.csv"), row.names = FALSE)


# 4. Summary (four leave-one-sample-out runs only) ----------------------------------
loo_gsea <- gsea_tab[gsea_tab$run != "baseline", ]
loo_gene <- gene_tab[gene_tab$run != "baseline", ]

cat("\n=============== SUMMARY (four leave-one-sample-out runs) ===============\n")
for (s in ifn_sets) {
  d <- loo_gsea[loo_gsea$pathway == s, ]
  cat(s, ":\n")
  cat("  NES range:", round(min(d$NES), 2), "to", round(max(d$NES), 2), "\n")
  cat("  padj range:", signif(min(d$padj), 3), "to", signif(max(d$padj), 3), "\n")
  cat("  significant (padj < 0.05) in", sum(d$padj < 0.05), "of", nrow(d), "runs\n")
}

cat("\nDirection of the three genes in the baseline and in each run:\n")
dir_tab <- gene_tab %>%
  mutate(direction = ifelse(meta_logFC > 0, "up", "down")) %>%
  select(SYMBOL, run, meta_logFC, meta_padj, direction) %>%
  arrange(SYMBOL, run)
print(as.data.frame(dir_tab), row.names = FALSE)

for (g in genes_to_check) {
  base_dir <- sign(gene_tab$meta_logFC[gene_tab$SYMBOL == g & gene_tab$run == "baseline"])
  d <- loo_gene[loo_gene$SYMBOL == g, ]
  cat(g, ": same direction as baseline in", sum(sign(d$meta_logFC) == base_dir),
      "of", nrow(d), "runs | nominal p < 0.05 in", sum(d$meta_p < 0.05),
      "| BH < 0.05 in", sum(d$meta_padj < 0.05), "runs\n")
}

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "17_LOO_interferon_GSEA_sessionInfo.txt"))
