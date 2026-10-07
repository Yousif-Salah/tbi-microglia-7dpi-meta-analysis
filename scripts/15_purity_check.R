# =============================================================================
# 15_purity_check.R
# Microglial purity check across the four datasets
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 15. The isolation methods differ among the datasets, so this
#           check shows whether each dataset behaves like a microglial
#           sample, using two microglial markers (P2ry12, Cx3cr1), one pan-
#           leukocyte marker (Ptprc) and one T-cell marker (Cd3e).
#
# Method
#   - For each dataset, the fitted DESeq2 object saved by scripts 01 to 04 is
#     read (same samples, same low-count filter, same design as the main
#     analysis).
#   - Expression is the variance-stabilizing transformation (VST,
#     blind = FALSE). The values of the four markers are taken for every
#     sample.
#   - Ensembl IDs of the markers come from org.Mm.eg.db (GSE276647 uses gene
#     symbols as row names).
#   - A marker that was removed by the low-count filter is reported as absent
#     (NA); the raw counts of every marker are printed so that detection can
#     be checked.
#
# Input
#   results/dds/GSE167459_dds.rds ... GSE283560_dds.rds        (scripts 01-04)
#
# Output
#   results/purity/purity_check_markers.csv        (one row per sample)
#   figures/FigureS6_purity_check.tiff             (TIFF, 300 dpi, 180 mm wide)
#   results/session_info/15_purity_check_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/15_purity_check.R
# Packages: DESeq2, org.Mm.eg.db, AnnotationDbi, dplyr, tidyr, ggplot2
# Note    : the Arial font must be installed for the figure text.
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

dds_dir  <- file.path("results", "dds")
out_dir  <- file.path("results", "purity")
fig_dir  <- "figures"
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, fig_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE276647", "GSE253476", "GSE167459", "GSE283560")
dds_files <- file.path(dds_dir, paste0(datasets, "_dds.rds"))
if (!all(file.exists(dds_files)))
  stop("Missing fitted DESeq2 objects in ", dds_dir,
       ". Run scripts 01 to 04 first. Current folder: ", getwd())

# The condition column of each dataset, and whether its rows are Ensembl IDs
condition_column <- c(GSE276647 = "condition", GSE253476 = "condition",
                      GSE167459 = "group",     GSE283560 = "condition")
rows_are_symbols <- c(GSE276647 = TRUE, GSE253476 = FALSE,
                      GSE167459 = FALSE, GSE283560 = FALSE)


# 1. Marker genes --------------------------------------------------------------
marker_genes <- c("P2ry12", "Cx3cr1", "Ptprc", "Cd3e")
ens_map <- mapIds(org.Mm.eg.db, keys = marker_genes, keytype = "SYMBOL", column = "ENSEMBL")
cat("Symbol to Ensembl map:\n"); print(ens_map)
stopifnot(!anyNA(ens_map))


# 2. VST values of the markers, per dataset --------------------------------------
purity_for_dataset <- function(ds) {
  dds <- readRDS(file.path(dds_dir, paste0(ds, "_dds.rds")))
  stopifnot(is(dds, "DESeqDataSet"))
  cat("\n", ds, ": genes", nrow(dds), "| samples", ncol(dds), "\n", sep = "")

  ids <- if (rows_are_symbols[[ds]]) setNames(marker_genes, marker_genes) else ens_map
  raw <- counts(dds, normalized = FALSE)
  for (sym in names(ids)) {
    if (ids[[sym]] %in% rownames(raw)) {
      cat(ds, "-", sym, ": raw counts =", paste(raw[ids[[sym]], ], collapse = ","), "\n")
    } else {
      cat(ds, "-", sym, ": not in the object (removed by the low-count filter or absent)\n")
    }
  }

  mat <- assay(vst(dds, blind = FALSE))
  present <- ids[ids %in% rownames(mat)]
  sub_mat <- mat[present, , drop = FALSE]
  rownames(sub_mat) <- names(ids)[match(present, ids)]
  df <- as.data.frame(t(sub_mat))
  for (m in marker_genes) if (!(m %in% colnames(df))) df[[m]] <- NA_real_
  df <- df[, marker_genes, drop = FALSE]
  df$sample    <- rownames(df)
  df$condition <- as.character(colData(dds)[[condition_column[[ds]]]])
  df$dataset   <- ds
  df
}

purity_all <- bind_rows(lapply(datasets, purity_for_dataset))

stopifnot(nrow(purity_all) == 31,
          all(table(purity_all$dataset)[datasets] == c(8, 6, 6, 11)))
cat("\nSamples per dataset:\n"); print(table(purity_all$dataset))
print(purity_all)

write.csv(purity_all, file.path(out_dir, "purity_check_markers.csv"), row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "purity_check_markers.csv"), "\n")


# 3. Figure S6: four-panel purity plot --------------------------------------------
plot_df <- purity_all %>%
  pivot_longer(cols = all_of(marker_genes), names_to = "gene", values_to = "vst_expr") %>%
  mutate(
    gene = factor(gene, levels = marker_genes),
    condition_group = ifelse(condition %in% c("Sham", "Control"), "Sham/Control", "TBI"),
    condition_group = factor(condition_group, levels = c("Sham/Control", "TBI")),
    dataset = factor(dataset, levels = c("GSE276647", "GSE283560", "GSE253476", "GSE167459"))
  )

p <- ggplot(plot_df, aes(x = dataset, y = vst_expr, color = condition_group)) +
  geom_jitter(width = 0.15, size = 1.6, alpha = 0.85) +
  stat_summary(fun = median, geom = "crossbar", width = 0.4, fatten = 1,
               color = "grey20", linewidth = 0.3) +
  facet_wrap(~ gene, nrow = 1, scales = "free_y") +
  scale_color_manual(values = c("Sham/Control" = "grey50", "TBI" = "#C0392B")) +
  labs(x = NULL, y = "VST-normalized expression", color = NULL,
       title = "Microglial purity check across four 7-dpi datasets") +
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title = element_text(size = 9, face = "bold")
  )

ggsave(file.path(fig_dir, "FigureS6_purity_check.tiff"), plot = p,
       width = 180, height = 65, units = "mm", dpi = 300, compression = "lzw")
cat("Figure saved:", file.path(fig_dir, "FigureS6_purity_check.tiff"), "\n")

writeLines(capture.output(sessionInfo()), file.path(info_dir, "15_purity_check_sessionInfo.txt"))
