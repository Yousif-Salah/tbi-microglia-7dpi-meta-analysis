# =============================================================================
# 10_LODO_Cybb_Fcer1g.R
# Leave-one-dataset-out (LODO) analysis for Cybb and Fcer1g
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 10a. Tests whether the pooled result for the two candidate
#           genes depends on any single dataset.
#
# Method
#   - For each of the 4 datasets in turn, the dataset is left out and the
#     pooled estimate of Cybb and Fcer1g is calculated again from the other 3
#     (random-effects model, REML, metafor::rma), as in script 05.
#   - Output per gene and left-out dataset: pooled log2 fold change, P value,
#     I2 and tau2.
#   - These P values are not corrected for multiple testing; they show how
#     stable the estimate is.
#
# Input
#   results/deg_tables/GSE167459_7d.csv ... GSE283560_7d.csv   (scripts 01-04)
#
# Output
#   results/lodo/LODO_Cybb_Fcer1g.csv
#   results/session_info/10_LODO_Cybb_Fcer1g_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/10_LODO_Cybb_Fcer1g.R
# Packages: dplyr, metafor, tibble
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(metafor)
  library(tibble)
})

input_dir <- file.path("results", "deg_tables")
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
  df2
}

deg_list <- lapply(deg_files, read_deg_clean)
datasets <- names(deg_list)
cat("Datasets read:", paste(datasets, collapse = ", "), "\n")


# 2. Pooling function (random effects, REML) ----------------------------------
run_gene_meta_core <- function(gene_symbol, gene_data) {
  n_studies <- nrow(gene_data)
  if (n_studies < 2) {
    return(tibble(gene = gene_symbol, pooled_log2FC = NA_real_, pooled_SE = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_, pval = NA_real_,
                  I2 = NA_real_, tau2 = NA_real_, n_datasets = n_studies,
                  datasets_used = paste(gene_data$dataset, collapse = ", "),
                  note = "Single dataset, not pooled"))
  }
  fit <- tryCatch(
    metafor::rma(yi = gene_data$log2FC, sei = gene_data$SE, method = "REML"),
    error = function(e) { message("  rma() failed for ", gene_symbol, ": ", e$message); NULL }
  )
  if (is.null(fit)) {
    return(tibble(gene = gene_symbol, pooled_log2FC = NA_real_, pooled_SE = NA_real_,
                  ci_lower = NA_real_, ci_upper = NA_real_, pval = NA_real_,
                  I2 = NA_real_, tau2 = NA_real_, n_datasets = n_studies,
                  datasets_used = paste(gene_data$dataset, collapse = ", "), note = "rma() failed"))
  }
  tibble(gene = gene_symbol, pooled_log2FC = as.numeric(fit$beta), pooled_SE = fit$se,
         ci_lower = fit$ci.lb, ci_upper = fit$ci.ub, pval = fit$pval,
         I2 = round(fit$I2, 1), tau2 = round(fit$tau2, 4), n_datasets = n_studies,
         datasets_used = paste(gene_data$dataset, collapse = ", "), note = "OK")
}


# 3. Leave one dataset out -----------------------------------------------------
lodo_genes <- c("Cybb", "Fcer1g")

gene_dataset_table <- bind_rows(lapply(datasets, function(ds) {
  deg_list[[ds]] %>%
    filter(SYMBOL %in% lodo_genes) %>%
    transmute(dataset = ds, SYMBOL, log2FC = log2FoldChange, SE = lfcSE)
}))

cat("\nPer-dataset values used:\n")
print(gene_dataset_table %>% arrange(SYMBOL, dataset))

# Each gene must appear once per dataset
stopifnot(all(table(gene_dataset_table$SYMBOL, gene_dataset_table$dataset) == 1))

lodo_results <- list()
for (drop_ds in datasets) {
  remaining <- gene_dataset_table %>% filter(dataset != drop_ds)
  for (g in lodo_genes) {
    res <- run_gene_meta_core(g, remaining %>% filter(SYMBOL == g))
    res$dropped_dataset <- drop_ds
    lodo_results[[paste(drop_ds, g)]] <- res
  }
}

lodo_summary <- bind_rows(lodo_results) %>%
  select(dropped_dataset, gene, n_datasets, pooled_log2FC, pval, I2, tau2, note)

cat("\nLeave-one-dataset-out summary:\n")
print(lodo_summary, n = Inf)


# 4. Save ---------------------------------------------------------------------
write.csv(lodo_summary, file.path(out_dir, "LODO_Cybb_Fcer1g.csv"), row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "LODO_Cybb_Fcer1g.csv"), "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "10_LODO_Cybb_Fcer1g_sessionInfo.txt"))
