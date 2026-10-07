# =============================================================================
# 12b_Figures_enrichment.R
# Dot plots of the enrichment results (overall and per MCODE cluster)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 12b. Draws the enrichment figures from the tables of script 12a.
#
# Method
#   - Only significant terms (p.adjust < 0.05) are drawn; the saved tables keep
#     all terms.
#   - Overall figure (four panels): GO BP up, GO BP down, Reactome up, KEGG up.
#     The 10 terms with the lowest p.adjust are shown per panel.
#     KEGG up: viral-disease and other annotation-artifact terms are removed
#     (see artifact_terms below). These pathways appear because of shared
#     immune genes in the KEGG annotation and not because of the biology.
#   - Cluster figure (three panels): simplified GO BP of MCODE clusters 1 to 3.
#
# Input   results/enrichment/*.csv                      (script 12a)
# Output  results/figures/Enrichment_overall.tiff       180 x 350 mm, 300 dpi
#         results/figures/Enrichment_clusters.tiff      180 x 260 mm, 300 dpi
#         results/session_info/12b_Figures_enrichment_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/12b_Figures_enrichment.R
# Packages: ggplot2, stringr, patchwork
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
  library(patchwork)
})

in_dir   <- file.path("results", "enrichment")
fig_dir  <- file.path("results", "figures")
info_dir <- file.path("results", "session_info")
for (d in c(fig_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

need <- c("GO_BP_UP.csv", "GO_BP_DOWN.csv", "KEGG_UP.csv", "Reactome_UP.csv")
if (!all(file.exists(file.path(in_dir, need))))
  stop("Enrichment tables not found in ", in_dir, ". Run script 12 first. ",
       "Current folder: ", getwd())


# 1. Dot plot function --------------------------------------------------------
make_dotplot <- function(df, title, n_show = 10, legend_order = "color_first") {
  if (is.null(df) || nrow(df) == 0) {
    return(ggplot() + labs(title = paste0(title, " (no significant terms)")) +
             theme_void() +
             theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)))
  }
  df <- df[order(df$p.adjust), ][1:min(n_show, nrow(df)), ]
  df$Description <- str_wrap(df$Description, width = 55)
  df$Description <- factor(df$Description, levels = rev(unique(df$Description)))

  p <- ggplot(df, aes(x = GeneRatio, y = Description, size = Count, color = p.adjust)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(title = title, x = "GeneRatio", y = NULL) +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5,
                                    margin = margin(b = 6, t = 15)),
          axis.text.y = element_text(size = 7, lineheight = 0.85))

  if (legend_order == "color_first") {
    p + guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))
  } else {
    p + guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2))
  }
}


# 2. Overall figure -----------------------------------------------------------
go_up   <- read.csv(file.path(in_dir, "GO_BP_UP.csv"),    stringsAsFactors = FALSE)
go_down <- read.csv(file.path(in_dir, "GO_BP_DOWN.csv"),  stringsAsFactors = FALSE)
kegg_up <- read.csv(file.path(in_dir, "KEGG_UP.csv"),     stringsAsFactors = FALSE)
rt_up   <- read.csv(file.path(in_dir, "Reactome_UP.csv"), stringsAsFactors = FALSE)

go_up   <- go_up[go_up$p.adjust < 0.05, ]
go_down <- go_down[go_down$p.adjust < 0.05, ]
rt_up   <- rt_up[rt_up$p.adjust < 0.05, ]

artifact_terms <- c("Coronavirus", "Herpes", "Epstein", "Influenza", "Measles",
                    "Hepatitis", "papillomavirus", "carcinogenesis", "[Vv]irus", "[Vv]iral",
                    "Toxoplasma", "Graft-versus-host", "diabetes mellitus",
                    "Autoimmune thyroid", "Allograft rejection", "HIV")
kegg_up <- kegg_up[kegg_up$p.adjust < 0.05, ]
kegg_up <- kegg_up[!grepl(paste(artifact_terms, collapse = "|"),
                          kegg_up$Description, ignore.case = TRUE), ]
cat("Significant terms drawn  GO up:", nrow(go_up), "| GO down:", nrow(go_down),
    "| Reactome up:", nrow(rt_up), "| KEGG up (after filter):", nrow(kegg_up), "\n")

p_A <- make_dotplot(go_up,   "GO BP: Upregulated")
p_B <- make_dotplot(go_down, "GO BP: Downregulated", legend_order = "size_first")
p_C <- make_dotplot(rt_up,   "Reactome: Upregulated")
p_D <- make_dotplot(kegg_up, "KEGG: Upregulated")

combined_ABCD <- (p_A / p_B / p_C / p_D) +
  plot_layout(heights = c(1, 1, 1, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(text = element_text(family = "sans"))

ggsave(file.path(fig_dir, "Enrichment_overall.tiff"), combined_ABCD,
       width = 180, height = 350, units = "mm", dpi = 300, compression = "lzw")
cat("Saved: Enrichment_overall.tiff\n")


# 3. Cluster figure -----------------------------------------------------------
cluster_titles <- c("Cluster 1: Interferon/ISG",
                    "Cluster 2: MHC Class I",
                    "Cluster 3: Ribosomal fragment")
cluster_keys <- c("Cluster1_ISG", "Cluster2_MHC1", "Cluster3_Ribo")

go_clusters <- lapply(cluster_keys, function(k) {
  df <- read.csv(file.path(in_dir, paste0("MCODE_", k, "_GO_BP_simplified.csv")),
                 stringsAsFactors = FALSE)
  df[df$p.adjust < 0.05, ]
})
for (i in seq_along(cluster_keys))
  cat(cluster_titles[i], "- significant simplified GO BP terms:", nrow(go_clusters[[i]]), "\n")

panels <- lapply(seq_along(cluster_keys),
                 function(i) make_dotplot(go_clusters[[i]], cluster_titles[i]))

combined_3 <- wrap_plots(panels, ncol = 1) +
  plot_annotation(tag_levels = "A") &
  theme(text = element_text(family = "sans"))

ggsave(file.path(fig_dir, "Enrichment_clusters.tiff"), combined_3,
       width = 180, height = 260, units = "mm", dpi = 300, compression = "lzw")
cat("Saved: Enrichment_clusters.tiff\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "12b_Figures_enrichment_sessionInfo.txt"))
