# =============================================================================
# 13_GSEA.R
# Gene set enrichment analysis (GSEA) of the pooled 7-dpi ranking, and Figure 7
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 13. Pre-ranked GSEA of all analyzable genes against
#           GO Biological Process, Hallmark and Reactome gene sets.
#
# Method
#   - Ranking statistic: sign(pooled log2 fold change) * -log10(pooled P value)
#     for every gene in the analyzable universe. A P value of exactly 0 is
#     replaced by the smallest non-zero P value, so no gene is dropped.
#   - Symbols are mapped to Entrez IDs (org.Mm.eg.db). If several symbols share
#     one Entrez ID, the one with the largest absolute ranking statistic is kept.
#   - Gene sets: msigdbr (mouse orthologs of the human MSigDB collections),
#     GO:BP (C5), Hallmark (H), Reactome (C2:CP:REACTOME).
#   - fgsea settings: minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000,
#     one worker (SerialParam), set.seed(1) once at the start.
#   - Figure 7: top significant terms per collection, split by direction.
#     Virus-named Reactome terms are left out of the figure because they are
#     driven by host ribosomal proteins (they stay in the CSV). GO terms with
#     "viral" in the name are left out of panel A because they repeat
#     "defense response to virus" (they stay in the CSV). Removed Reactome
#     terms are printed to the console.
#
# Note: the 14-gene BTK/NOX2 panel is tested separately in script 14 with its
#       own settings (see the header of that script).
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv                      (script 05)
#
# Output
#   results/gsea/GSEA_GOBP.csv, GSEA_Hallmark.csv, GSEA_Reactome.csv
#   figures/Figure7_GSEA.tiff                       (TIFF, 600 dpi, 180 mm wide)
#   results/session_info/13_GSEA_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/13_GSEA.R
# Packages: fgsea, msigdbr (25.1.1), clusterProfiler, org.Mm.eg.db, dplyr,
#           ggplot2, patchwork, BiocParallel
# Note    : the Arial font must be installed for the figure text.
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(1)  # fixed seed for reproducible NES and P values

in_file  <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
out_dir  <- file.path("results", "gsea")
fig_dir  <- "figures"
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, fig_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_file))
  stop("Cannot find ", in_file, ". Run script 05 first. Current folder: ", getwd())

cat("fgsea version:  ", as.character(packageVersion("fgsea")), "\n")
cat("msigdbr version:", as.character(packageVersion("msigdbr")), "\n")


# 1. Ranked gene list ---------------------------------------------------------
meta_res <- read.csv(in_file, stringsAsFactors = FALSE)

stopifnot(all(c("SYMBOL", "meta_logFC", "meta_p") %in% names(meta_res)))
stopifnot(nrow(meta_res) == 16636)
stopifnot(abs(meta_res$meta_logFC[meta_res$SYMBOL == "Cybb"]   - 1.66) < 0.01)
stopifnot(abs(meta_res$meta_logFC[meta_res$SYMBOL == "Fcer1g"] - 0.30) < 0.01)
cat("Input verified: 16,636 genes, Cybb and Fcer1g match expected values.\n")

meta_res$rank_stat <- sign(meta_res$meta_logFC) * -log10(meta_res$meta_p)

# A P value of exactly 0 gives an infinite ranking statistic. Replace it with
# the smallest non-zero P value in the data instead of dropping the gene.
n_inf <- sum(is.infinite(meta_res$rank_stat))
cat("Genes with infinite ranking statistic (P = 0 exactly):", n_inf, "\n")
if (n_inf > 0) {
  smallest_nonzero_p <- min(meta_res$meta_p[meta_res$meta_p > 0], na.rm = TRUE)
  inf_idx <- is.infinite(meta_res$rank_stat)
  meta_res$rank_stat[inf_idx] <- sign(meta_res$meta_logFC[inf_idx]) * -log10(smallest_nonzero_p)
}
stopifnot(!any(is.na(meta_res$rank_stat)))

id_map <- bitr(meta_res$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("SYMBOL to ENTREZ rows returned by bitr:", nrow(id_map), "\n")
cat("Unique symbols mapped:", length(unique(id_map$SYMBOL)), "of", nrow(meta_res),
    "(", round(100 * length(unique(id_map$SYMBOL)) / nrow(meta_res), 1), "%)\n")

meta_res_mapped <- merge(meta_res, id_map, by = "SYMBOL")

# If several symbols map to the same Entrez ID, keep the one with the largest
# absolute ranking statistic.
meta_res_mapped <- meta_res_mapped %>%
  group_by(ENTREZID) %>%
  slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
  ungroup()

ranked_vec <- setNames(meta_res_mapped$rank_stat, meta_res_mapped$ENTREZID)
ranked_vec <- sort(ranked_vec, decreasing = TRUE)
cat("Final ranked vector length:", length(ranked_vec), "\n")
stopifnot(!any(duplicated(names(ranked_vec))))


# 2. Run GSEA (GO BP, Hallmark, Reactome; minSize = 15) -----------------------
gs_gobp     <- msigdbr(species = "Mus musculus", collection = "C5", subcollection = "GO:BP")
gs_hallmark <- msigdbr(species = "Mus musculus", collection = "H")
gs_reactome <- msigdbr(species = "Mus musculus", collection = "C2", subcollection = "CP:REACTOME")
# msigdbr 25.1.1 uses the human MSigDB collections mapped to mouse orthologs
# (default db_species = "HS").

to_fgsea_list <- function(gs_df) split(as.character(gs_df$ncbi_gene), gs_df$gs_name)

run_fgsea <- function(pathways, ranks, label) {
  res <- fgsea(pathways = pathways, stats = ranks,
               minSize = 15, maxSize = 500, eps = 0, nPermSimple = 10000,
               BPPARAM = BiocParallel::SerialParam())  # one core: repeatable results
  res <- res[order(res$padj), ]
  cat("\n---", label, "---\n")
  cat("Total pathways tested:", nrow(res), "\n")
  cat("Significant (padj < 0.05):", sum(res$padj < 0.05, na.rm = TRUE), "\n")
  cat("  positive NES:", sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE),
      "| negative NES:", sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE), "\n")
  res
}

# The three calls must stay in this order: the random seed is set once.
gsea_gobp     <- run_fgsea(to_fgsea_list(gs_gobp),     ranked_vec, "GO Biological Process")
gsea_hallmark <- run_fgsea(to_fgsea_list(gs_hallmark), ranked_vec, "Hallmark")
gsea_reactome <- run_fgsea(to_fgsea_list(gs_reactome), ranked_vec, "Reactome")

flatten_leadingEdge <- function(res) {
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = "/")
  as.data.frame(res)
}

write.csv(flatten_leadingEdge(gsea_gobp),     file.path(out_dir, "GSEA_GOBP.csv"),     row.names = FALSE)
write.csv(flatten_leadingEdge(gsea_hallmark), file.path(out_dir, "GSEA_Hallmark.csv"), row.names = FALSE)
write.csv(flatten_leadingEdge(gsea_reactome), file.path(out_dir, "GSEA_Reactome.csv"), row.names = FALSE)


# 3. Figure 7 -----------------------------------------------------------------
clean_name <- function(x, prefix) {
  x <- gsub(prefix, "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  x <- paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
  # Restore capital letters for acronyms and gene names
  fixes <- c("\\bbmp\\b" = "BMP", "\\bnotch\\b" = "Notch", "\\bsrp\\b" = "SRP",
             "\\bfoxo\\b" = "FOXO", "\\bnmd\\b" = "NMD", "\\beif2ak4\\b" = "EIF2AK4",
             "\\bgcn2\\b" = "GCN2", "\\brrna\\b" = "rRNA", "\\bmrna\\b" = "mRNA", "\\bn linked\\b" = "N linked")
  for (i in seq_along(fixes)) x <- gsub(names(fixes)[i], fixes[i], x, ignore.case = TRUE)
  x
}

# Virus-named Reactome terms are driven by ribosomal proteins (checked from their
# leading edge genes). Every significant term this pattern removes is printed.
artifact_pattern <- "virus|viral|influenza|hepatitis|herpes|coronavirus|papillomavirus|immunodeficiency|measles|carcinogenesis|sars_cov"

gobp_clean     <- as.data.frame(gsea_gobp)
hallmark_clean <- as.data.frame(gsea_hallmark)
reactome_all   <- as.data.frame(gsea_reactome)

is_artifact <- grepl(artifact_pattern, reactome_all$pathway, ignore.case = TRUE)
removed <- reactome_all[is_artifact & !is.na(reactome_all$padj) & reactome_all$padj < 0.05,
                        c("pathway", "NES", "padj")]
cat("\nSignificant Reactome terms removed from Figure 7 by the viral filter:", nrow(removed), "\n")
print(removed)

reactome_clean <- reactome_all[!is_artifact, ]
cat("Reactome terms after filter:", nrow(reactome_clean),
    "(significant:", sum(reactome_clean$padj < 0.05, na.rm = TRUE), ")\n")

# Top significant terms per collection, split by direction, so that
# downregulated terms are shown even when upregulated terms are more numerous.
select_top <- function(df, prefix, n_up = 8, n_down = 5) {
  df <- df[!is.na(df$padj) & df$padj < 0.05, ]
  up   <- df[df$NES > 0, ]
  down <- df[df$NES < 0, ]
  up   <- up[order(up$padj), ][seq_len(min(n_up, nrow(up))), ]
  down <- down[order(down$padj), ][seq_len(min(n_down, nrow(down))), ]
  combo <- rbind(up, down)
  combo$label <- clean_name(combo$pathway, prefix)
  combo$direction <- ifelse(combo$NES > 0, "Up (TBI vs sham)", "Down (TBI vs sham)")
  combo[order(combo$NES), ]
}

# Panel A only: GO terms with "viral" in the name are left out (see header).
gobp_top     <- select_top(gobp_clean[!grepl("VIRAL", gobp_clean$pathway), ], "GOBP_")
hallmark_top <- select_top(hallmark_clean, "HALLMARK_")
reactome_top <- select_top(reactome_clean, "REACTOME_")

cat("\nGO BP:", sum(gobp_top$NES > 0), "up /", sum(gobp_top$NES < 0), "down shown\n")
cat("Hallmark:", sum(hallmark_top$NES > 0), "up /", sum(hallmark_top$NES < 0), "down shown\n")
cat("Reactome:", sum(reactome_top$NES > 0), "up /", sum(reactome_top$NES < 0), "down shown\n")

UP_COLOR   <- "#C0392B"
DOWN_COLOR <- "#2471A3"

fig_theme <- theme_minimal(base_family = "Arial", base_size = 9) +
  theme(axis.text.y = element_text(size = 9, color = "black"),
        axis.text.x = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 9), axis.title.y = element_blank(),
        plot.title = element_text(size = 10, face = "bold", hjust = 0),
        legend.position = "none",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())

make_panel <- function(df, title) {
  ggplot(df, aes(x = reorder(label, NES), y = NES, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.3) +
    coord_flip() +
    scale_fill_manual(values = c("Up (TBI vs sham)" = UP_COLOR, "Down (TBI vs sham)" = DOWN_COLOR)) +
    labs(title = title, x = NULL, y = "Normalized Enrichment Score (NES)") +
    fig_theme
}

panel_A <- make_panel(gobp_top,     "A. GO Biological Process")
panel_B <- make_panel(hallmark_top, "B. Hallmark")
panel_C <- make_panel(reactome_top, "C. Reactome")

fig7 <- panel_A / panel_B / panel_C + plot_layout(heights = c(1, 0.8, 1))

ggsave(file.path(fig_dir, "Figure7_GSEA.tiff"), fig7, width = 180, height = 250, units = "mm",
       dpi = 600, compression = "lzw")
cat("\nFigure 7 saved.\n")


# 4. Print all significant terms (for the Results text) -----------------------
for (nm in c("gobp_clean", "hallmark_clean", "reactome_all")) {
  d <- get(nm)
  d <- d[!is.na(d$padj) & d$padj < 0.05, c("pathway", "size", "NES", "padj")]
  cat("\n---", nm, ": significant terms ---\n")
  print(d[order(-d$NES), ], row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(info_dir, "13_GSEA_sessionInfo.txt"))
