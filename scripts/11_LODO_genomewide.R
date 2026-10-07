# =============================================================================
# 11_LODO_genomewide.R
# Genome-wide leave-one-dataset-out (LODO) analysis
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 11. Tests whether the number of meta-DEGs depends on any single
#           dataset.
#
# Method
#   - For each of the 4 datasets in turn, the dataset is left out and the whole
#     random-effects meta-analysis (REML, fixed-effect fallback, genes in at
#     least 2 of the remaining 3 datasets) is repeated with Benjamini-Hochberg
#     correction. The number of meta-DEGs (adjusted P < 0.05) is counted.
#   - Same input rules as script 05: a row needs a symbol, a fold change and
#     lfcSE > 0, and each gene is used once per dataset (the first row, which
#     is the one with the lowest adjusted P value, is kept if a symbol repeats).
#   - Check: the same pooling with all 4 datasets must give exactly the result
#     of script 05 (same genes, same meta-DEG count). The script stops if not.
#
# Input
#   results/deg_tables/*_7d.csv                      (scripts 01-04)
#   results/meta/Meta7dpi_all_genes_k2.csv           (script 05, for the check)
#
# Output (results/lodo/)
#   LODO_GenomeWide_summary.csv        one row per left-out dataset: genes
#                                      tested, meta-DEGs, up, down
#   LODO_GenomeWide_drop_<dataset>.csv pooled result of each iteration
#   results/session_info/11_LODO_genomewide_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/11_LODO_genomewide.R      (takes several minutes)
# Packages: dplyr, metafor, tibble, purrr
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(tibble)
  library(purrr)
})

input_dir <- file.path("results", "deg_tables")
meta_file <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
out_dir   <- file.path("results", "lodo")
info_dir  <- file.path("results", "session_info")
for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

deg_files <- list.files(input_dir, pattern = "_7d\\.csv$", full.names = TRUE)
if (length(deg_files) != 4)
  stop("Expected 4 per-dataset tables in ", input_dir, " but found ", length(deg_files),
       ". Run scripts 01 to 04 first. Current folder: ", getwd())
names(deg_files) <- gsub("_7d\\.csv$", "", basename(deg_files))


# 1. Read the per-dataset tables ----------------------------------------------
read_deg_clean <- function(fpath) {
  df <- read.csv(fpath, check.names = FALSE)
  colnames(df) <- trimws(colnames(df))
  sym_col <- c("SYMBOL","symbol","Gene","gene")[c("SYMBOL","symbol","Gene","gene") %in% colnames(df)][1]
  lfc_col <- c("log2FoldChange","logFC","LOG2FOLDCHANGE")[c("log2FoldChange","logFC","LOG2FOLDCHANGE") %in% colnames(df)][1]
  se_col  <- c("lfcSE","SE","se")[c("lfcSE","SE","se") %in% colnames(df)][1]
  df2 <- df[, c(sym_col, lfc_col, se_col)]
  colnames(df2) <- c("SYMBOL","log2FoldChange","lfcSE")
  df2 <- df2[!is.na(df2$SYMBOL) & df2$SYMBOL != "", ]
  df2
}

deg_list <- lapply(deg_files, read_deg_clean)
datasets <- names(deg_list)
cat("Datasets read:", paste(datasets, collapse = ", "), "\n")

# Long table: dataset, SYMBOL, log2FC, SE
full_long <- bind_rows(lapply(datasets, function(ds) {
  deg_list[[ds]] %>% transmute(dataset = ds, SYMBOL, log2FC = log2FoldChange, SE = lfcSE)
}))

full_long <- full_long %>%
  filter(!is.na(log2FC), !is.na(SE), SE > 0) %>%
  group_by(dataset, SYMBOL) %>%
  slice(1) %>%
  ungroup()
cat("Gene-dataset rows used:", nrow(full_long), "\n")


# 2. Pooling (random effects, REML, fixed-effect fallback) ----------------------
pool_one_gene <- function(gene_data) {
  n_studies <- nrow(gene_data)
  if (n_studies < 2) return(NULL)
  fit <- tryCatch(metafor::rma(yi = gene_data$log2FC, sei = gene_data$SE, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit))
    fit <- tryCatch(metafor::rma(yi = gene_data$log2FC, sei = gene_data$SE, method = "FE"),
                    error = function(e) NULL)
  if (is.null(fit))
    return(tibble(meta_logFC = NA_real_, meta_p = NA_real_, I2 = NA_real_,
                  tau2 = NA_real_, k = n_studies, method = "failed"))
  tibble(meta_logFC = as.numeric(fit$beta), meta_p = fit$pval,
         I2 = round(fit$I2, 1), tau2 = round(fit$tau2, 4), k = n_studies,
         method = if (inherits(fit, "rma.uni") && !is.null(fit$method) && fit$method == "FE") "FE" else "REML")
}

run_genome_wide_meta <- function(long_table) {
  long_table %>%
    group_by(SYMBOL) %>%
    filter(n() >= 2) %>%
    group_split() %>%
    map_dfr(function(gd) {
      res <- pool_one_gene(gd)
      if (is.null(res)) return(NULL)
      res$SYMBOL <- gd$SYMBOL[1]
      res
    })
}


# 3. Check: all 4 datasets must reproduce the result of script 05 ---------------
cat("\n=== All 4 datasets (check against script 05) ===\n")
baseline <- run_genome_wide_meta(full_long)
baseline$meta_padj <- p.adjust(baseline$meta_p, method = "BH")
cat("Analyzable genes:", nrow(baseline),
    "| meta-DEGs:", sum(baseline$meta_padj < 0.05, na.rm = TRUE), "\n")

if (file.exists(meta_file)) {
  main <- read.csv(meta_file)
  cmp  <- merge(baseline, main, by = "SYMBOL", suffixes = c("_lodo", "_main"))
  stopifnot(nrow(baseline) == nrow(main),
            nrow(cmp) == nrow(main),
            sum(baseline$meta_padj < 0.05, na.rm = TRUE) == sum(main$meta_padj < 0.05, na.rm = TRUE),
            max(abs(cmp$meta_logFC_lodo - cmp$meta_logFC_main), na.rm = TRUE) < 1e-6)
  cat("Check passed: identical to script 05.\n")
}


# 4. Leave one dataset out ------------------------------------------------------
lodo_summary <- list()

for (drop_ds in datasets) {
  cat("\n=== Leaving out", drop_ds, "===\n")
  res <- run_genome_wide_meta(full_long %>% filter(dataset != drop_ds))
  res$meta_padj <- p.adjust(res$meta_p, method = "BH")

  n_up   <- sum(res$meta_padj < 0.05 & res$meta_logFC > 0, na.rm = TRUE)
  n_down <- sum(res$meta_padj < 0.05 & res$meta_logFC < 0, na.rm = TRUE)
  cat("  Genes tested:", nrow(res), "| meta-DEGs:", n_up + n_down,
      "(", n_up, "up /", n_down, "down )\n")

  lodo_summary[[drop_ds]] <- tibble(dropped_dataset = drop_ds,
                                    n_genes_tested = nrow(res),
                                    n_meta_DEGs = n_up + n_down,
                                    n_up = n_up, n_down = n_down)
  write.csv(res, file.path(out_dir, paste0("LODO_GenomeWide_drop_", drop_ds, ".csv")),
            row.names = FALSE)
}

lodo_summary_df <- bind_rows(lodo_summary)
cat("\nGenome-wide LODO summary:\n")
print(lodo_summary_df)

write.csv(lodo_summary_df, file.path(out_dir, "LODO_GenomeWide_summary.csv"), row.names = FALSE)
cat("\nSaved to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "11_LODO_genomewide_sessionInfo.txt"))
