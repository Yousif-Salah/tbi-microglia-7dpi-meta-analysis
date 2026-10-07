# =============================================================================
# 07_Figure2_volcano_forest.R
# Figure 2: genome-wide volcano plot (panel A) and forest plots of four
# representative genes (panels B to E)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 7. Draws Figure 2 from the pooled result (script 05) and the
#           four per-dataset result tables (scripts 01 to 04).
#
# Method
#   - Panel A : volcano plot of all pooled genes. x = pooled log2 fold change,
#               y = -log10 of the BH-adjusted P value. Meta-DEGs (meta_padj <
#               0.05) are coloured by direction. Labels: the 10 most significant
#               up genes and the 5 most significant down genes.
#   - Panels B to E : forest plots for Lpl, Spp1 (up) and Sall1, Tjp1 (down).
#               Each plot shows the per-dataset log2 fold change with a 95%
#               confidence interval (estimate +/- 1.96 x lfcSE) and the pooled
#               random-effects estimate (REML), as in script 05.
#   - Check   : the pooled estimate drawn in each forest plot must equal
#               meta_logFC of that gene in the pooled table.
#   - Colours : up #C0392B, down #2471A3.
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv           (script 05)
#   results/deg_tables/GSE167459_7d.csv  ... GSE283560_7d.csv   (scripts 01-04)
#
# Output (results/figures/)
#   Figure2A_volcano.tiff                  180 x 130 mm, 300 dpi
#   Fig2_Forest_<gene>.tiff                one per gene, 180 x 60 mm, 300 dpi
#   Fig2_Forest_Combined_BE.tiff           four forest plots, 180 x 220 mm
#   results/session_info/07_Figure2_volcano_forest_sessionInfo.txt
#
# Note    : text uses the Arial font. If Arial is not installed on your system,
#           R falls back to a default font and the plots look slightly different.
#
# Usage   : run from the repository root, for example
#           Rscript scripts/07_Figure2_volcano_forest.R
# Packages: ggplot2, ggrepel, dplyr, tibble, metafor, patchwork
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tibble)
  library(metafor)
  library(patchwork)
})

meta_file   <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
input_dir   <- file.path("results", "deg_tables")
out_fig_dir <- file.path("results", "figures")
info_dir    <- file.path("results", "session_info")

for (d in c(out_fig_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(meta_file))
  stop("Pooled result not found: ", meta_file,
       ". Run scripts 01 to 05 first. Current folder: ", getwd())


# 1. Panel A: volcano plot ----------------------------------------------------
meta_res <- read.csv(meta_file, stringsAsFactors = FALSE)
cat("Analyzable genes:", nrow(meta_res), "\n")

meta_res$direction <- "NS"
meta_res$direction[meta_res$meta_padj < 0.05 & meta_res$meta_logFC > 0] <- "UP"
meta_res$direction[meta_res$meta_padj < 0.05 & meta_res$meta_logFC < 0] <- "DOWN"
meta_res$direction <- factor(meta_res$direction, levels = c("UP", "DOWN", "NS"))

n_up   <- sum(meta_res$direction == "UP")
n_down <- sum(meta_res$direction == "DOWN")
cat("Meta-DEGs: up =", n_up, "| down =", n_down, "\n")
stopifnot(n_up >= 10, n_down >= 5)

top_up   <- meta_res[meta_res$direction == "UP", ]
top_up   <- top_up[order(top_up$meta_padj), ][1:10, ]
top_down <- meta_res[meta_res$direction == "DOWN", ]
top_down <- top_down[order(top_down$meta_padj), ][1:5, ]
label_genes <- rbind(top_up, top_down)

p_volcano <- ggplot(meta_res, aes(x = meta_logFC, y = -log10(meta_padj),
                                  color = direction)) +
  geom_point(size = 1.2, alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "black",
             linewidth = 0.4) +
  geom_text_repel(data = label_genes, aes(label = SYMBOL), size = 2.6,
                  color = "black", max.overlaps = 20, segment.size = 0.2,
                  family = "Arial") +
  scale_color_manual(values = c(UP = "#C0392B", DOWN = "#2471A3", NS = "grey70"),
                     labels = c(paste0("Up (n=", n_up, ")"),
                                paste0("Down (n=", n_down, ")"), "NS"),
                     name = "Direction") +
  labs(title = "Meta-analysis Volcano Plot (7 dpi)",
       x = expression(paste("Meta log"[2], " Fold Change")),
       y = expression(paste(-log[10], "(BH-adj ", italic(p), ")"))) +
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold", hjust = 0.5),
        legend.position = "right")

ggsave(file.path(out_fig_dir, "Figure2A_volcano.tiff"), p_volcano,
       width = 180, height = 130, units = "mm", dpi = 300, compression = "lzw")
cat("Saved:", file.path(out_fig_dir, "Figure2A_volcano.tiff"), "\n")


# 2. Panels B to E: forest plots ---------------------------------------------
# Read one dataset table and standardize the column names, because the four
# tables come from different scripts and do not share identical names.
read_deg_clean <- function(fpath) {
  df <- read.csv(fpath, check.names = FALSE)
  colnames(df) <- trimws(colnames(df))
  sym_col  <- c("SYMBOL", "symbol", "Gene", "gene")[c("SYMBOL", "symbol", "Gene", "gene") %in% colnames(df)][1]
  lfc_col  <- c("log2FoldChange", "logFC", "LOG2FOLDCHANGE")[c("log2FoldChange", "logFC", "LOG2FOLDCHANGE") %in% colnames(df)][1]
  se_col   <- c("lfcSE", "SE", "se")[c("lfcSE", "SE", "se") %in% colnames(df)][1]
  stat_col <- c("stat", "WaldStat", "STAT")[c("stat", "WaldStat", "STAT") %in% colnames(df)][1]
  if (is.na(sym_col) | is.na(lfc_col)) { warning("Missing SYMBOL or log2FC; skipping."); return(NULL) }
  if (is.na(se_col)) {
    if (!is.na(stat_col)) { df$lfcSE <- abs(df[[lfc_col]] / df[[stat_col]]); se_col <- "lfcSE" }
    else { warning("No lfcSE/SE and no stat; skipping."); return(NULL) }
  }
  df2 <- df[, c(sym_col, lfc_col, se_col)]
  colnames(df2) <- c("SYMBOL", "log2FoldChange", "lfcSE")
  df2
}

deg_files_vec <- list.files(input_dir, pattern = "_7d\\.csv$", full.names = TRUE)
names(deg_files_vec) <- gsub("_7d\\.csv$", "", basename(deg_files_vec))
deg_list <- lapply(deg_files_vec, read_deg_clean)
names(deg_list) <- names(deg_files_vec)
deg_list <- deg_list[!vapply(deg_list, is.null, logical(1))]
cat("Datasets read for forest plots:", paste(names(deg_list), collapse = ", "), "\n")
stopifnot(length(deg_list) == 4)

# Per-gene random-effects pooling (REML), the same model as script 05
run_gene_meta_core <- function(gene_symbol, gene_data) {
  n_studies <- nrow(gene_data)
  if (n_studies < 2) {
    return(tibble(gene = gene_symbol,
                  pooled_log2FC = if (n_studies == 1) gene_data$log2FC else NA_real_,
                  pooled_SE = if (n_studies == 1) gene_data$SE else NA_real_,
                  ci_lower = if (n_studies == 1) gene_data$log2FC - 1.96 * gene_data$SE else NA_real_,
                  ci_upper = if (n_studies == 1) gene_data$log2FC + 1.96 * gene_data$SE else NA_real_,
                  pval = NA_real_, I2 = NA_real_, tau2 = NA_real_, n_datasets = n_studies,
                  datasets_used = paste(gene_data$dataset, collapse = ", "),
                  note = if (n_studies == 1) "Single dataset" else "Not detected"))
  }
  fit <- tryCatch(metafor::rma(yi = gene_data$log2FC, sei = gene_data$SE, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) return(tibble(gene = gene_symbol, pooled_log2FC = NA_real_, pooled_SE = NA_real_,
                                  ci_lower = NA_real_, ci_upper = NA_real_, pval = NA_real_, I2 = NA_real_, tau2 = NA_real_,
                                  n_datasets = n_studies, datasets_used = paste(gene_data$dataset, collapse = ", "), note = "rma() failed"))
  tibble(gene = gene_symbol, pooled_log2FC = as.numeric(fit$beta), pooled_SE = fit$se,
         ci_lower = fit$ci.lb, ci_upper = fit$ci.ub, pval = fit$pval,
         I2 = round(fit$I2, 1), tau2 = round(fit$tau2, 4), n_datasets = n_studies,
         datasets_used = paste(gene_data$dataset, collapse = ", "), note = "OK")
}

explain_missing_datasets <- function(gene_symbol, deg_list_in, datasets_used) {
  missing_ds <- setdiff(names(deg_list_in), datasets_used)
  if (length(missing_ds) == 0) return(NULL)
  reasons <- sapply(missing_ds, function(ds) {
    df <- deg_list_in[[ds]]
    if (is.null(df) || !(gene_symbol %in% df$SYMBOL)) "gene not present in this dataset's output"
    else "present but excluded (failed QC)"
  })
  paste0(missing_ds, " (", reasons, ")", collapse = "; ")
}

# Forest plot for one gene: per-dataset estimates (squares) and the pooled
# random-effects estimate (diamond), coloured by the pooled direction.
make_forest_plot <- function(gene_symbol, deg_list_in) {
  indiv <- bind_rows(lapply(names(deg_list_in), function(ds) {
    df <- deg_list_in[[ds]]
    row <- df[!is.na(df$SYMBOL) & df$SYMBOL == gene_symbol & !is.na(df$log2FoldChange) & !is.na(df$lfcSE), ]
    if (nrow(row) == 0) return(NULL)
    data.frame(dataset = ds, log2FC = row$log2FoldChange[1], SE = row$lfcSE[1])
  }))
  if (is.null(indiv) || nrow(indiv) == 0) { message("Gene ", gene_symbol, " not found."); return(NULL) }

  pooled_row <- run_gene_meta_core(gene_symbol, indiv)

  col_up   <- "#C0392B"
  col_down <- "#2471A3"
  point_color <- if (!is.na(pooled_row$pooled_log2FC) && pooled_row$pooled_log2FC > 0) col_up else col_down

  indiv_plot <- indiv %>%
    mutate(ci_lo = log2FC - 1.96 * SE, ci_hi = log2FC + 1.96 * SE, weight = 1 / SE^2,
           label_txt = sprintf("%.2f [%.2f, %.2f]", log2FC, ci_lo, ci_hi)) %>%
    mutate(y_pos = rev(row_number()))

  w_range <- range(indiv_plot$weight)
  indiv_plot$pt_size <- if (diff(w_range) == 0) 4 else
    2.5 + (indiv_plot$weight - w_range[1]) / diff(w_range) * 3.5

  has_pool <- !is.na(pooled_row$pooled_log2FC) && pooled_row$n_datasets >= 2
  y_labels <- c(indiv_plot$dataset, if (has_pool) "Random-effects model")
  y_breaks <- c(indiv_plot$y_pos, if (has_pool) 0)

  missing_note <- explain_missing_datasets(gene_symbol, deg_list_in, indiv$dataset)
  if (!is.na(pooled_row$pval)) {
    pv_label <- if (pooled_row$pval < 0.001) "p < 0.001" else sprintf("p = %.3f", pooled_row$pval)
    sub_text <- sprintf("Pooled log2FC = %.2f [%.2f, %.2f]   %s   I² = %.0f%%   τ² = %.3f   k = %d datasets",
                        pooled_row$pooled_log2FC, pooled_row$ci_lower, pooled_row$ci_upper,
                        pv_label, pooled_row$I2, pooled_row$tau2, pooled_row$n_datasets)
  } else sub_text <- sprintf("Single dataset (k = %d)", pooled_row$n_datasets)
  if (!is.null(missing_note)) sub_text <- paste0(sub_text, "\nNot detected in: ", missing_note)

  x_max <- max(indiv_plot$ci_hi, if (has_pool) pooled_row$ci_upper else -Inf)
  x_min <- min(indiv_plot$ci_lo, if (has_pool) pooled_row$ci_lower else Inf)
  label_x <- x_max + 0.12 * (x_max - x_min)

  p <- ggplot() +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.4) +
    geom_errorbarh(data = indiv_plot, aes(xmin = ci_lo, xmax = ci_hi, y = y_pos),
                   height = 0.15, color = point_color, linewidth = 0.5) +
    geom_point(data = indiv_plot, aes(x = log2FC, y = y_pos, size = pt_size),
               shape = 15, color = point_color) +
    geom_text(data = indiv_plot, aes(x = label_x, y = y_pos, label = label_txt),
              hjust = 0, size = 2.5, family = "Arial", color = "black") +
    scale_size_identity()

  if (has_pool) {
    diamond_df <- data.frame(
      x = c(pooled_row$ci_lower, pooled_row$pooled_log2FC, pooled_row$ci_upper, pooled_row$pooled_log2FC),
      y = c(0, 0.18, 0, -0.18))
    pooled_label <- sprintf("%.2f [%.2f, %.2f]", pooled_row$pooled_log2FC, pooled_row$ci_lower, pooled_row$ci_upper)
    p <- p +
      geom_polygon(data = diamond_df, aes(x = x, y = y), fill = point_color, color = point_color, alpha = 0.85) +
      geom_text(aes(x = label_x, y = 0, label = pooled_label), hjust = 0, size = 2.5,
                family = "Arial", color = "black", fontface = "bold")
  }

  p <- p +
    scale_y_continuous(breaks = y_breaks, labels = y_labels,
                       expand = expansion(mult = c(0.18, 0.18))) +
    coord_cartesian(xlim = c(x_min - 0.05 * (x_max - x_min), label_x + 0.35 * (x_max - x_min))) +
    labs(title = paste0("Meta-analysis Forest Plot – ", gene_symbol),
         subtitle = sub_text, x = "Log2 Fold Change (TBI vs Sham)", y = NULL) +
    theme_classic(base_size = 9, base_family = "Arial") +
    theme(plot.title = element_text(face = "bold", size = 10, hjust = 0),
          plot.subtitle = element_text(size = 7, color = "black"),
          axis.text.y = element_text(size = 8, color = "black"),
          axis.line.y = element_blank(), axis.ticks.y = element_blank())

  attr(p, "pooled_log2FC") <- pooled_row$pooled_log2FC   # used for the check below
  p
}


# 3. Build and save the four panels ------------------------------------------
# Lpl and Spp1 are up-regulated; Sall1 and Tjp1 are down-regulated meta-DEGs.
figure2_genes <- c("Lpl", "Spp1", "Sall1", "Tjp1")
figure2_forest_plots <- lapply(figure2_genes, function(g) make_forest_plot(g, deg_list))
names(figure2_forest_plots) <- figure2_genes
stopifnot(all(!sapply(figure2_forest_plots, is.null)))

# The pooled estimate in each forest plot must equal the pooled table
for (g in figure2_genes) {
  tab_val  <- meta_res$meta_logFC[meta_res$SYMBOL == g]
  plot_val <- attr(figure2_forest_plots[[g]], "pooled_log2FC")
  cat(sprintf("%-6s pooled log2FC: plot %.6f | table %.6f\n", g, plot_val, tab_val))
  stopifnot(length(tab_val) == 1, abs(plot_val - tab_val) < 1e-6)
}

for (g in figure2_genes) {
  ggsave(file.path(out_fig_dir, paste0("Fig2_Forest_", g, ".tiff")),
         figure2_forest_plots[[g]], width = 180, height = 60, units = "mm",
         dpi = 300, compression = "lzw")
}

combined_fig2 <- wrap_plots(figure2_forest_plots, ncol = 1) +
  plot_annotation(title = "Figure 2B–E: Representative Meta-DEG Forest Plots (7 dpi)") &
  theme(text = element_text(family = "Arial"))

ggsave(file.path(out_fig_dir, "Fig2_Forest_Combined_BE.tiff"), combined_fig2,
       width = 180, height = 220, units = "mm", dpi = 300, compression = "lzw")

cat("Saved 1 volcano, 4 forest and 1 combined TIFF to", out_fig_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "07_Figure2_volcano_forest_sessionInfo.txt"))
