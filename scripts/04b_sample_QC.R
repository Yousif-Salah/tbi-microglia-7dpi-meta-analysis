# =============================================================================
# 04b_sample_QC.R
# Sample quality control of the four datasets (library size, log2 CPM boxplot,
# PCA, Spearman correlation heatmap)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 4b. Quality control of the samples used in the primary analysis
#           (Methods 2.2.3; Supplementary Figures S1 to S4: S1 = GSE253476,
#           S2 = GSE283560, S3 = GSE276647, S4 = GSE167459). It shows which
#           sample of each dataset separates most from the others. Script 16
#           tests those samples by leaving them out.
#
# Method (the same for all four datasets)
#   1. The fitted DESeq2 object of the dataset is read (scripts 01 to 04). It
#      holds the samples and genes of the primary analysis (low-count filter
#      rowSums >= 10 already applied).
#   2. Counts are transformed with the variance-stabilising transformation
#      (vst, blind = TRUE, so the group labels do not influence it).
#   3. Boxplot of log2 CPM (log2 of counts per million + 1) for every sample.
#   4. PCA on the 500 genes with the highest variance after vst (scaled).
#   5. Spearman correlation of all sample pairs (vst values), shown as a
#      heatmap with hierarchical clustering.
#   6. A sample table is written: library size, PC1, PC2, mean Spearman
#      correlation with the other samples of the same group, and the sample
#      that script 16 leaves out. The script prints whether that sample has
#      the lowest within-group correlation of its dataset.
#
# Input
#   results/dds/<dataset>_dds.rds                 (scripts 01-04)
#
# Output (results/qc/)
#   QC_<dataset>_a_boxplot.png, QC_<dataset>_b_PCA.png,
#   QC_<dataset>_c_correlation_heatmap.png        (300 dpi)
#   QC_sample_table.csv       one row per sample
#   QC_summary.csv            one row per dataset (PC1 and PC2 % variance, ...)
#   results/session_info/04b_sample_QC_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/04b_sample_QC.R
# Packages: DESeq2, ggplot2, ggrepel, pheatmap
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(pheatmap)
})

dds_dir  <- file.path("results", "dds")
out_dir  <- file.path("results", "qc")
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

datasets <- c("GSE167459", "GSE253476", "GSE276647", "GSE283560")

# Column of the sample table that holds the Sham/Control versus injury group
condition_column <- c(GSE167459 = "group",     GSE253476 = "condition",
                      GSE276647 = "condition", GSE283560 = "condition")

# Sample that script 16 leaves out in each dataset (flagged in this QC)
flagged_sample <- c(GSE167459 = "1304M_S77", GSE253476 = "Sham_Mg_R3",
                    GSE276647 = "TBI_Iso_1", GSE283560 = "GSM8666331")

dds_files <- file.path(dds_dir, paste0(datasets, "_dds.rds"))
if (!all(file.exists(dds_files)))
  stop("Missing input files. Run scripts 01 to 04 first. Current folder: ", getwd())


# 1. QC of one dataset -------------------------------------------------------
qc_one <- function(ds) {
  dds <- readRDS(file.path(dds_dir, paste0(ds, "_dds.rds")))
  cond <- as.character(colData(dds)[[condition_column[[ds]]]])
  n    <- ncol(dds)
  cat("\n====", ds, ":", n, "samples,", nrow(dds), "genes ====\n")
  print(table(cond))

  # Library size and log2 CPM
  lib    <- colSums(counts(dds))
  logcpm <- log2(sweep(counts(dds), 2, lib / 1e6, "/") + 1)

  # Variance-stabilising transformation, blind to the group labels
  vst_mat <- assay(vst(dds, blind = TRUE))

  # PCA on the 500 most variable genes
  top_genes <- names(sort(apply(vst_mat, 1, var), decreasing = TRUE))[1:500]
  pca     <- prcomp(t(vst_mat[top_genes, ]), scale. = TRUE)
  pct_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2))[1:2], 1)
  pca_df  <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                        condition = cond, sample = colnames(dds))

  # Spearman correlation of all sample pairs
  cors <- cor(vst_mat, method = "spearman")
  ann  <- data.frame(condition = cond, row.names = colnames(dds))

  # Mean correlation with the other samples of the same group
  within_cor <- sapply(seq_len(n), function(i) {
    same <- which(cond == cond[i] & seq_len(n) != i)
    mean(cors[i, same])
  })

  tag <- paste0(ds, ", 7dpi (n=", n, ")")

  # (a) boxplot
  png(file.path(out_dir, paste0("QC_", ds, "_a_boxplot.png")),
      width = 6, height = 5, units = "in", res = 300)
  boxplot(logcpm, las = 2, main = paste0("Boxplot: log2 CPM (", tag, ")"),
          ylab = "log2 CPM")
  dev.off()

  # (b) PCA
  p <- ggplot(pca_df, aes(PC1, PC2, color = condition, label = sample)) +
    geom_point(size = 3) +
    geom_text_repel(size = 3, max.overlaps = 30) +
    theme_minimal() +
    ggtitle(paste0("PCA (VST, top 500 var genes): ", tag)) +
    xlab(paste0("PC1: ", pct_var[1], "% variance")) +
    ylab(paste0("PC2: ", pct_var[2], "% variance"))
  ggsave(file.path(out_dir, paste0("QC_", ds, "_b_PCA.png")), p,
         width = 6, height = 5, dpi = 300)

  # (c) correlation heatmap
  pheatmap(cors, annotation_col = ann, border_color = NA,
           main = paste0("Spearman correlation (VST): ", tag),
           filename = file.path(out_dir, paste0("QC_", ds, "_c_correlation_heatmap.png")),
           width = 6, height = 5)

  cat("Library size (sorted):\n"); print(sort(lib))
  cat("PC1:", pct_var[1], "% | PC2:", pct_var[2], "%\n")

  tab <- data.frame(dataset = ds, sample = colnames(dds), condition = cond,
                    library_size = as.numeric(lib),
                    PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                    mean_within_group_spearman = within_cor,
                    flagged_in_script16 = colnames(dds) == flagged_sample[[ds]],
                    stringsAsFactors = FALSE)

  # The flagged sample must exist; report whether it has the lowest
  # within-group correlation (a check for the reader, not a filter).
  stopifnot(flagged_sample[[ds]] %in% colnames(dds))
  lowest <- colnames(dds)[which.min(within_cor)]
  cat("Sample flagged for script 16:", flagged_sample[[ds]],
      "| lowest within-group correlation:", lowest,
      if (lowest == flagged_sample[[ds]]) "(same)\n" else "(DIFFERENT)\n")

  summ <- data.frame(dataset = ds, n_samples = n, n_genes = nrow(dds),
                     PC1_percent = pct_var[1], PC2_percent = pct_var[2],
                     flagged_sample = flagged_sample[[ds]],
                     lowest_within_group_correlation = lowest,
                     stringsAsFactors = FALSE)
  list(tab = tab, summ = summ)
}

res <- lapply(datasets, qc_one)

write.csv(do.call(rbind, lapply(res, `[[`, "tab")),
          file.path(out_dir, "QC_sample_table.csv"), row.names = FALSE)
write.csv(do.call(rbind, lapply(res, `[[`, "summ")),
          file.path(out_dir, "QC_summary.csv"), row.names = FALSE)
cat("\nSaved: results/qc/QC_sample_table.csv and QC_summary.csv\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "04b_sample_QC_sessionInfo.txt"))
