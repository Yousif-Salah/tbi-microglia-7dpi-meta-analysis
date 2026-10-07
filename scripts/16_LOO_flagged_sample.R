# =============================================================================
# 16_LOO_flagged_sample.R
# Leave-one-sample-out (LOO) analysis of the flagged sample in each dataset
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 16. Tests whether the pooled result depends on the one sample
#           per dataset that separated from its group in the quality control
#           (PCA and sample correlation; see the QC scripts).
#
# Method (the same for all four datasets)
#   1. The fitted DESeq2 object of the dataset (scripts 01 to 04) is read. Its
#      raw counts and sample table are used to build a new object without the
#      flagged sample. The genes are the ones kept by the low-count filter of
#      the full dataset (the filter is not repeated).
#   2. DESeq2 is run again on the remaining samples (same design; GSE283560
#      with minReplicatesForReplace = Inf, as in script 04).
#   3. Log2 fold changes are shrunk with apeglm. If apeglm gives more than five
#      convergence warnings, shrinkage falls back to type = "normal". The same
#      rule applies to every dataset.
#   4. The new table of this dataset is pooled again with the unchanged tables
#      of the other three datasets (random effects, REML, genes in at least two
#      datasets, Benjamini-Hochberg correction).
#   5. The number of meta-DEGs and the adjusted P values of Cybb, Fcer1g and Btk
#      are recorded.
#
# Pooling rule for repeated gene symbols (the same rule as script 05): rows with
# a missing symbol, fold change or standard error, or with a standard error of 0,
# are removed. The table is sorted by adjusted P value and the first row is kept
# for each gene symbol.
#
# Input
#   results/dds/<dataset>_dds.rds                              (scripts 01-04)
#   results/deg_tables/<dataset>_7d.csv                        (scripts 01-04)
#
# Output (results/loo/)
#   LOO_pooled_<dataset>.csv     full pooled table of each run
#                                (read by script 17)
#   LOO_global_summary.csv       one row per left-out sample
#   results/session_info/16_LOO_flagged_sample_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/16_LOO_flagged_sample.R   (takes several minutes)
# Packages: DESeq2, apeglm, metafor, dplyr
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
  library(metafor)
  library(dplyr)
})

dds_dir  <- file.path("results", "dds")
deg_dir  <- file.path("results", "deg_tables")
out_dir  <- file.path("results", "loo")
info_dir <- file.path("results", "session_info")
for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

all_datasets <- c("GSE167459", "GSE253476", "GSE276647", "GSE283560")

dds_files <- file.path(dds_dir, paste0(all_datasets, "_dds.rds"))
deg_files <- file.path(deg_dir, paste0(all_datasets, "_7d.csv"))
if (!all(file.exists(c(dds_files, deg_files))))
  stop("Missing input files. Run scripts 01 to 04 first. Current folder: ", getwd())

# Flagged sample of each dataset, from the quality control of each dataset
flagged_sample <- c(GSE167459 = "1304M_S77", GSE253476 = "Sham_Mg_R3",
                    GSE276647 = "TBI_Iso_1", GSE283560 = "GSM8666331")

# Model coefficient of the TBI-versus-control contrast of each dataset
coef_name <- c(GSE167459 = "group_FPI_vs_Sham",
               GSE253476 = "condition_TBI_7d_vs_Sham",
               GSE276647 = "condition_TBI_vs_Sham",
               GSE283560 = "condition_TBI_vs_Control")
# Datasets whose rows are Ensembl IDs (GSE276647 rows are gene symbols)
rows_are_ensembl <- c(GSE167459 = TRUE, GSE253476 = TRUE, GSE276647 = FALSE, GSE283560 = TRUE)
min_replicates_for_replace <- c(GSE167459 = 7, GSE253476 = 7, GSE276647 = 7, GSE283560 = Inf)

main_baseline_degs <- 160   # meta-DEGs of the main analysis (script 05)


# 1. Helpers ------------------------------------------------------------------
shrink_with_fallback <- function(dds_loo, coef) {
  warn_count <- 0
  res <- withCallingHandlers(
    lfcShrink(dds_loo, coef = coef, type = "apeglm"),
    warning = function(w) {
      if (grepl("line search", conditionMessage(w))) warn_count <<- warn_count + 1
      invokeRestart("muffleWarning")
    }
  )
  if (warn_count > 5) {
    message(sprintf("apeglm failed to converge (%d warnings); using type = 'normal'.", warn_count))
    res <- lfcShrink(dds_loo, coef = coef, type = "normal")
    method_used <- "normal (apeglm failed to converge)"
  } else {
    method_used <- "apeglm"
  }
  list(res = res, method = method_used)
}

# Same rule as script 05: remove unusable rows, sort by adjusted P value (rows
# without an adjusted P value come last) and keep the first row of each symbol.
dedup_first_by_padj <- function(df) {
  df <- df[!is.na(df$SYMBOL) & df$SYMBOL != "" &
           !is.na(df$log2FoldChange) & !is.na(df$lfcSE) & df$lfcSE > 0, ]
  df <- df[order(df$padj), ]
  df[!duplicated(df$SYMBOL), ]
}

# New DESeq2 object without one sample. The size factors of the full dataset
# are removed so that DESeq2 estimates them again from the remaining samples.
drop_sample <- function(dds_full, sample_id) {
  keep <- colnames(dds_full) != sample_id
  cd <- as.data.frame(colData(dds_full))[keep, , drop = FALSE]
  cd$sizeFactor <- NULL
  for (v in names(cd)) if (is.factor(cd[[v]])) cd[[v]] <- droplevels(cd[[v]])
  DESeqDataSetFromMatrix(countData = counts(dds_full)[, keep, drop = FALSE],
                         colData = cd, design = design(dds_full))
}

pool_tables <- function(deg_list) {
  long <- bind_rows(lapply(names(deg_list), function(ds) {
    d <- deg_list[[ds]]
    d <- d[!is.na(d$SYMBOL) & d$SYMBOL != "" & !is.na(d$log2FoldChange), ]
    data.frame(dataset = ds, SYMBOL = d$SYMBOL, log2FC = d$log2FoldChange, SE = d$lfcSE,
               stringsAsFactors = FALSE)
  }))
  by_gene <- split(long, long$SYMBOL)
  by_gene <- by_gene[vapply(by_gene, nrow, integer(1)) >= 2]
  pooled <- bind_rows(lapply(names(by_gene), function(g) {
    rows <- by_gene[[g]]
    fit <- tryCatch(suppressWarnings(rma(yi = log2FC, sei = SE, data = rows, method = "REML")),
                    error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    data.frame(SYMBOL = g, meta_logFC = fit$b[1], meta_p = fit$pval, k = nrow(rows),
               stringsAsFactors = FALSE)
  }))
  pooled$meta_padj <- p.adjust(pooled$meta_p, method = "BH")
  pooled
}

# Tables of the main analysis, one row per symbol (same rule as script 05)
deg_main <- setNames(lapply(deg_files, function(f) dedup_first_by_padj(read.csv(f, stringsAsFactors = FALSE))),
                     all_datasets)

# Safety check: the four unchanged tables must give the result of script 05
# (160 meta-DEGs). If not, the pooling rule of this script differs from script 05.
main_check <- pool_tables(deg_main)
n_main <- sum(main_check$meta_padj < 0.05, na.rm = TRUE)
cat("Check: all four tables, no sample removed:", nrow(main_check), "genes,",
    n_main, "meta-DEGs (script 05: 16636 genes, 160 meta-DEGs)\n")
stopifnot(n_main == main_baseline_degs)
rm(main_check)


# 2. One run per dataset ---------------------------------------------------------
summary_rows <- list()

for (ds in all_datasets) {
  cat("\n=== Leaving out", flagged_sample[[ds]], "from", ds, "===\n")
  dds_full <- readRDS(file.path(dds_dir, paste0(ds, "_dds.rds")))
  flagged  <- flagged_sample[[ds]]
  stopifnot(flagged %in% colnames(dds_full), coef_name[[ds]] %in% resultsNames(dds_full))
  cat("Genes:", nrow(dds_full), "| samples:", ncol(dds_full), "| after drop:", ncol(dds_full) - 1, "\n")

  dds_loo <- drop_sample(dds_full, flagged)
  dds_loo <- DESeq(dds_loo, minReplicatesForReplace = min_replicates_for_replace[[ds]])
  stopifnot(coef_name[[ds]] %in% resultsNames(dds_loo))

  shrink_out <- shrink_with_fallback(dds_loo, coef_name[[ds]])
  cat("Shrinkage method used:", shrink_out$method, "\n")
  res_df <- as.data.frame(shrink_out$res)

  # Symbols: Ensembl IDs are mapped with the table of the main analysis
  if (rows_are_ensembl[[ds]]) {
    canon <- read.csv(file.path(deg_dir, paste0(ds, "_7d.csv")), stringsAsFactors = FALSE)
    res_df$SYMBOL <- setNames(canon$SYMBOL, canon$ENSEMBL)[rownames(res_df)]
    cat("Unmapped IDs:", sum(is.na(res_df$SYMBOL)), "of", nrow(res_df), "\n")
  } else {
    res_df$SYMBOL <- rownames(res_df)
  }
  res_df <- res_df[!is.na(res_df$SYMBOL) & res_df$SYMBOL != "", ]
  res_df <- res_df[!is.na(res_df$log2FoldChange) & !is.na(res_df$lfcSE), ]
  res_df <- dedup_first_by_padj(res_df)
  cat("Genes with a valid estimate after the refit:", nrow(res_df), "\n")

  # Pool again: the other three datasets stay unchanged, this dataset is replaced
  deg_iter <- deg_main[setdiff(all_datasets, ds)]
  deg_iter[[ds]] <- res_df[, c("SYMBOL", "log2FoldChange", "lfcSE")]
  pooled <- pool_tables(deg_iter)

  write.csv(pooled, file.path(out_dir, paste0("LOO_pooled_", ds, ".csv")), row.names = FALSE)

  n_sig <- sum(pooled$meta_padj < 0.05, na.rm = TRUE)
  get_padj <- function(g) pooled$meta_padj[pooled$SYMBOL == g]
  cat("Pooled genes:", nrow(pooled), "| meta-DEGs:", n_sig,
      "(main analysis:", main_baseline_degs, ")\n")
  cat("Cybb padj:", get_padj("Cybb"), "| Fcer1g padj:", get_padj("Fcer1g"),
      "| Btk padj:", get_padj("Btk"), "\n")

  summary_rows[[ds]] <- data.frame(
    dropped_dataset = ds, flagged_sample = flagged,
    global_meta_DEG_count = n_sig,
    Cybb_padj = get_padj("Cybb"), Fcer1g_padj = get_padj("Fcer1g"), Btk_padj = get_padj("Btk"),
    shrinkage_method = shrink_out$method, stringsAsFactors = FALSE)
  rm(dds_full, dds_loo, res_df, pooled)
}


# 3. Save the summary -----------------------------------------------------------
summary_df <- bind_rows(summary_rows)
write.csv(summary_df, file.path(out_dir, "LOO_global_summary.csv"), row.names = FALSE)
cat("\nSaved:", file.path(out_dir, "LOO_global_summary.csv"), "\n")
print(summary_df)

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "16_LOO_flagged_sample_sessionInfo.txt"))
