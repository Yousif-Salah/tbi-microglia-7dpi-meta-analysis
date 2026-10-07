# =============================================================================
# 18_BTK_NOX2_forest_plots.R
# Figure 7a, Figure 7b, Figure S6 and Table S17 (BTK/NOX2 panel)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 18.
#
# Outputs
#   Table S17  per-gene table of the 14-gene panel (pooled log2FC, 95% CI,
#              P, BH-adjusted P, I2, k)
#   Figure 7a  pooled log2FC with 95% CI for the 14 panel genes, by module
#   Figure 7b  per-dataset forest plot of Cybb (primary analysis)
#   Figure S6  per-dataset forest plots of Fcer1g, primary (shrunk) analysis
#              and unshrunk sensitivity analysis
#
# Inputs
#   results/meta/Meta7dpi_all_genes_k2.csv                (script 05)
#   results/deg_tables/<dataset>_7d.csv                   (scripts 01 to 04)
#   results/meta_unshrunk/Meta7dpi_all_genes_k2_unshrunk.csv  (script 08)
#   results/deg_tables_unshrunk/<dataset>_7d_unshrunk.csv     (scripts 01 to 04)
#
# Outputs
#   results/panel/Table_S7_panel_per_gene.csv
#   results/panel/Table_S7_formatted.csv
#   figures/Figure8a_panel_forest.tiff, Figure8b_Cybb_forest.tiff,
#   figures/FigureS7_Fcer1g_forest.tiff           (300 dpi, 180 mm wide)
#   The file names keep the earlier numbers: Table_S7 = Table S17,
#   Figure8a = Figure 7a, Figure8b = Figure 7b, FigureS7 = Figure S6.
#   results/session_info/18_BTK_NOX2_forest_plots_sessionInfo.txt
#
# Notes
#   - Per-dataset rows: one row per gene and dataset (the first row of the
#     table written by scripts 01 to 04, which is sorted by adjusted P value),
#     and lfcSE > 0.
#   - The weight of each dataset is re-computed here with the same random-
#     effects model (REML) and the pooled estimate is checked against the
#     pooled table. The script stops if they differ by more than 0.02.
#   - Every number printed on the figures comes from the data.
#   - Text uses the Arial font. If Arial is not installed, change "Arial" in
#     the FONT line below.
#
# Usage   : run from the repository root, for example
#           Rscript scripts/18_BTK_NOX2_forest_plots.R
# Packages: ggplot2, dplyr, metafor
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(metafor)
})

FONT <- "Arial"

meta_file   <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
deg_dir     <- file.path("results", "deg_tables")
meta_unsh   <- file.path("results", "meta_unshrunk", "Meta7dpi_all_genes_k2_unshrunk.csv")
deg_dir_uns <- file.path("results", "deg_tables_unshrunk")
out_dir     <- file.path("results", "panel")
fig_dir     <- "figures"
info_dir    <- file.path("results", "session_info")
for (d in c(out_dir, fig_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

study_order <- c("GSE167459", "GSE253476", "GSE276647", "GSE283560")
files_main  <- setNames(paste0(study_order, "_7d.csv"), study_order)
files_unsh  <- setNames(paste0(study_order, "_7d_unshrunk.csv"), study_order)

need_files <- c(meta_file, meta_unsh,
                file.path(deg_dir, files_main), file.path(deg_dir_uns, files_unsh))
missing <- need_files[!file.exists(need_files)]
if (length(missing) > 0)
  stop("Missing input file(s):\n", paste(missing, collapse = "\n"),
       "\nRun scripts 01 to 05 and 08 first. Current folder: ", getwd())

btk_genes   <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1")
nox2_genes  <- c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
panel_genes <- c(btk_genes, nox2_genes)
stopifnot(length(panel_genes) == 14, !any(duplicated(panel_genes)))

COL_BTK  <- "#2166AC"
COL_NOX2 <- "#B2182B"


# 1. Pooled table and Table S7 -------------------------------------------------
read_pooled <- function(path) {
  m <- read.csv(path, stringsAsFactors = FALSE)
  m <- m[!is.na(m$SYMBOL) & m$SYMBOL != "", ]
  need <- c("SYMBOL", "k", "meta_logFC", "meta_p", "meta_padj",
            "ci_lb", "ci_ub", "I2")
  if (!all(need %in% names(m)))
    stop(path, " is missing: ", paste(setdiff(need, names(m)), collapse = ", "))
  m
}

meta <- read_pooled(meta_file)
stopifnot(nrow(meta) == 16636)

panel <- meta[meta$SYMBOL %in% panel_genes,
              c("SYMBOL", "meta_logFC", "ci_lb", "ci_ub", "meta_p", "meta_padj", "I2", "k")]
stopifnot(nrow(panel) == 14, !any(duplicated(panel$SYMBOL)), !any(is.na(panel)))
panel$module <- ifelse(panel$SYMBOL %in% btk_genes, "BTK", "NOX2")
panel <- panel[order(panel$module, panel$meta_padj),
               c("SYMBOL", "module", "meta_logFC", "ci_lb", "ci_ub",
                 "meta_p", "meta_padj", "I2", "k")]
stopifnot(sum(panel$module == "BTK") == 8, sum(panel$module == "NOX2") == 6)
write.csv(panel, file.path(out_dir, "Table_S7_panel_per_gene.csv"), row.names = FALSE)

tbl_pub <- panel %>%
  arrange(module, meta_padj) %>%
  transmute(Gene = SYMBOL, Module = module,
            `Pooled log2FC` = round(meta_logFC, 3),
            `95% CI lower`  = round(ci_lb, 3),
            `95% CI upper`  = round(ci_ub, 3),
            `p value`       = signif(meta_p, 3),
            `BH-adjusted p value` = signif(meta_padj, 3),
            `I2 (%)`        = round(I2, 1),
            `Datasets (k)`  = k)
write.csv(tbl_pub, file.path(out_dir, "Table_S7_formatted.csv"), row.names = FALSE)
cat("\nPanel genes with BH-adjusted P < 0.05:",
    paste(panel$SYMBOL[panel$meta_padj < 0.05], collapse = ", "), "\n")


# 2. Figure 8a: pooled log2FC of the 14 panel genes ----------------------------
d8a <- panel %>%
  mutate(stars = case_when(meta_padj < 0.001 ~ "***",
                           meta_padj < 0.01  ~ "**",
                           meta_padj < 0.05  ~ "*",
                           TRUE              ~ "ns"),
         module_lab = factor(ifelse(module == "BTK", "BTK module (8 genes)",
                                    "NOX2 module (6 genes)"),
                             levels = c("BTK module (8 genes)", "NOX2 module (6 genes)"))) %>%
  arrange(module, meta_logFC)
d8a$gene <- factor(d8a$SYMBOL, levels = d8a$SYMBOL)
cat("\nSignificance marks (Figure 8a):\n")
print(d8a[, c("SYMBOL", "module", "meta_padj", "stars")], row.names = FALSE)

fig8a <- ggplot(d8a, aes(x = meta_logFC, y = gene, color = module)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), orientation = "y",
                width = 0, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(x = ci_ub, label = stars), hjust = -0.4, vjust = 0.3,
            size = 9 / ggplot2::.pt, family = FONT, color = "black") +
  scale_color_manual(values = c(BTK = COL_BTK, NOX2 = COL_NOX2), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.03, 0.10))) +
  facet_grid(module_lab ~ ., scales = "free_y", space = "free_y") +
  labs(title = "BTK/NOX2 candidate panel: pooled log2 fold change (14 genes)",
       x = "Pooled log2 fold change (TBI vs sham/control), 95% CI", y = NULL) +
  theme_minimal(base_family = FONT, base_size = 9) +
  theme(axis.text   = element_text(size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black", face = "italic"),
        axis.title  = element_text(size = 9, color = "black"),
        plot.title  = element_text(size = 9, face = "bold", color = "black", hjust = 0),
        strip.background = element_rect(fill = "grey85", color = NA),
        strip.text = element_text(size = 9, face = "bold", color = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        plot.margin = margin(t = 8, r = 12, b = 5, l = 5))

ggsave(file.path(fig_dir, "Figure8a_panel_forest.tiff"), fig8a,
       width = 180, height = 110, units = "mm", dpi = 300, compression = "lzw")


# 3. Forest-plot helpers (Figure 8b and Figure S7) ----------------------------
# One gene, one analysis: four dataset rows (with random-effects weights)
# and one pooled row.
build_forest_data <- function(gene, dir, files, pooled_tbl, label) {
  st <- bind_rows(lapply(names(files), function(nm) {
    df <- read.csv(file.path(dir, files[[nm]]), stringsAsFactors = FALSE)
    need <- c("SYMBOL", "log2FoldChange", "lfcSE")
    if (!all(need %in% names(df)))
      stop(nm, " is missing: ", paste(setdiff(need, names(df)), collapse = ", "))
    row <- df[!is.na(df$SYMBOL) & df$SYMBOL == gene &
                !is.na(df$log2FoldChange) & !is.na(df$lfcSE) & df$lfcSE > 0, need]
    if (nrow(row) < 1) stop(nm, ": no usable row for ", gene)
    row <- row[1, ]                       # first row, same rule as script 05
    row$dataset <- nm
    row
  }))
  stopifnot(nrow(st) == length(files))

  fit <- rma(yi = st$log2FoldChange, sei = st$lfcSE, method = "REML")
  st$weight <- as.numeric(weights(fit))   # percent, sums to 100

  pr <- pooled_tbl[pooled_tbl$SYMBOL == gene, ]
  stopifnot(nrow(pr) == 1)
  cat(sprintf("\n%s | %s\n  re-fit pooled log2FC = %.3f | table pooled log2FC = %.3f | table I2 = %.1f%%\n",
              gene, label, as.numeric(fit$beta), pr$meta_logFC, pr$I2))
  if (abs(as.numeric(fit$beta) - pr$meta_logFC) > 0.02)
    stop("Re-fit pooled estimate differs from the table for ", gene, " (", label, ").")
  print(st[, c("dataset", "log2FoldChange", "lfcSE", "weight")], row.names = FALSE)

  studies <- data.frame(dataset = st$dataset, log2FC = st$log2FoldChange,
                        ci_lb = st$log2FoldChange - 1.96 * st$lfcSE,
                        ci_ub = st$log2FoldChange + 1.96 * st$lfcSE,
                        weight = st$weight, row_type = "study",
                        stringsAsFactors = FALSE)
  pooled <- data.frame(dataset = "Pooled (random-effects)", log2FC = pr$meta_logFC,
                       ci_lb = pr$ci_lb, ci_ub = pr$ci_ub,
                       weight = NA_real_, row_type = "pooled",
                       stringsAsFactors = FALSE)
  out <- rbind(studies, pooled)
  out$dataset <- factor(out$dataset, levels = rev(c(study_order, "Pooled (random-effects)")))
  out$est_lab <- sprintf("%.2f [%.2f, %.2f]", out$log2FC, out$ci_lb, out$ci_ub)
  out$panel <- sprintf("%s\nPooled log2FC = %.2f [%.2f, %.2f]; p = %s; I² = %.1f%%; τ² = %.2f; k = %d",
                       label, pr$meta_logFC, pr$ci_lb, pr$ci_ub,
                       format(signif(pr$meta_p, 3), scientific = TRUE),
                       pr$I2, fit$tau2, nrow(st))
  out
}

draw_forest <- function(dat, col) {
  x_rng <- range(c(dat$ci_lb, dat$ci_ub, 0))
  x_txt <- x_rng[2] + 0.05 * diff(x_rng)
  ggplot(dat, aes(x = log2FC, y = dataset)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbar(aes(xmin = ci_lb, xmax = ci_ub), orientation = "y",
                  width = 0, linewidth = 0.8, color = col) +
    geom_point(data = dat[dat$row_type == "study", ], aes(size = weight),
               shape = 15, color = col) +
    geom_point(data = dat[dat$row_type == "pooled", ], shape = 18, size = 5, color = col) +
    geom_text(aes(x = x_txt, label = est_lab), hjust = 0,
              size = 9 / ggplot2::.pt, family = FONT) +
    scale_size_area(max_size = 5, guide = "none") +
    scale_x_continuous(breaks = pretty(x_rng), expand = expansion(mult = c(0.03, 0.60))) +
    facet_wrap(~ panel, ncol = 1, scales = "free_y", strip.position = "top") +
    labs(x = "log2 fold change (TBI vs sham/control), 95% CI", y = NULL) +
    theme_minimal(base_family = FONT, base_size = 9) +
    theme(axis.text  = element_text(size = 9, color = "black"),
          axis.title = element_text(size = 9, color = "black"),
          strip.background = element_rect(fill = "grey85", color = NA),
          strip.text = element_text(size = 9, face = "plain", color = "black", hjust = 0),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          plot.margin = margin(t = 8, r = 12, b = 5, l = 5))
}


# 4. Figure 8b: Cybb, primary analysis ------------------------------------------
cybb_dat <- build_forest_data("Cybb", deg_dir, files_main, meta, "Cybb")
fig8b <- draw_forest(cybb_dat, col = COL_NOX2)
ggsave(file.path(fig_dir, "Figure8b_Cybb_forest.tiff"), fig8b,
       width = 180, height = 75, units = "mm", dpi = 300, compression = "lzw")


# 5. Figure S7: Fcer1g, primary vs unshrunk -------------------------------------
meta_u <- read_pooled(meta_unsh)
fcer_h <- build_forest_data("Fcer1g", deg_dir,     files_main, meta,   "Fcer1g, primary (shrunk)")
fcer_u <- build_forest_data("Fcer1g", deg_dir_uns, files_unsh, meta_u, "Fcer1g, unshrunk")
fcer_dat <- rbind(fcer_h, fcer_u)
fcer_dat$panel <- factor(fcer_dat$panel, levels = unique(fcer_dat$panel))
figS7 <- draw_forest(fcer_dat, col = COL_BTK)
ggsave(file.path(fig_dir, "FigureS7_Fcer1g_forest.tiff"), figS7,
       width = 180, height = 120, units = "mm", dpi = 300, compression = "lzw")

cat("\nDone. Saved Table S7 to", out_dir, "and three figures to", fig_dir, "\n")
writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "18_BTK_NOX2_forest_plots_sessionInfo.txt"))
