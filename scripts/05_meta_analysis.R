# =============================================================================
# 05_meta_analysis.R
# Random-effects meta-analysis of the four per-dataset DESeq2 results
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 5. Pools the outputs of scripts 01 to 04. Its result table is
#           the input of all later analyses (figures, network, enrichment,
#           WGCNA, GSEA, sensitivity analyses).
#
# Method
#   - Input per gene : apeglm-shrunk log2 fold change (log2FoldChange) and its
#                      standard error (lfcSE) from each dataset.
#   - Genes used     : genes with a symbol, a fold change and lfcSE > 0 in at
#                      least 2 of the 4 datasets (k >= 2). One row per gene and
#                      dataset (the first row is kept if a symbol repeats; the
#                      input tables are sorted by adjusted P value).
#   - Model          : random-effects model per gene, REML estimator
#                      (metafor::rma.uni). If REML returns an error, a
#                      fixed-effect model is used for that gene. If REML does
#                      not converge, metafor sets tau2 = 0 and the gene is
#                      kept. If both models fail, the gene is kept with NA
#                      estimates.
#   - Multiple tests : Benjamini-Hochberg correction across all pooled genes.
#   - Meta-DEG       : adjusted P value (meta_padj) < 0.05.
#   - meta_stat      : sign(meta_logFC) * -log10(meta_p), the ranking statistic
#                      used for GSEA.
#
# Input
#   results/deg_tables/GSE167459_7d.csv
#   results/deg_tables/GSE253476_7d.csv
#   results/deg_tables/GSE276647_7d.csv
#   results/deg_tables/GSE283560_7d.csv
#     Each needs the columns SYMBOL, log2FoldChange, lfcSE. ENSEMBL is used when
#     present (GSE276647 has no ENSEMBL column).
#
# Output
#   results/meta/Meta7dpi_all_genes_k2.csv       all pooled genes
#   results/meta/Meta7dpi_metaDEGs_FDR0.05.csv   meta-DEGs only
#   results/session_info/05_meta_analysis_sessionInfo.txt
#   Unshrunk arm (run through script 08): the same files with the suffix
#   _unshrunk, written to results/meta_unshrunk/ from results/deg_tables_unshrunk/.
#
# Usage   : run from the repository root, for example
#           Rscript scripts/05_meta_analysis.R
# Packages: dplyr, purrr, tibble, metafor
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(tibble)
  library(metafor)
})

# Two analysis arms use this same script:
#   "main"     : shrunk fold changes from scripts 01 to 04 (default)
#   "unshrunk" : unshrunk fold changes (sensitivity analysis, run by script 08,
#                which sets meta_arm before it calls this script)
arm <- if (exists("meta_arm")) meta_arm else "main"
stopifnot(arm %in% c("main", "unshrunk"))
cat("Analysis arm:", arm, "\n")

if (arm == "main") {
  input_dir  <- file.path("results", "deg_tables")
  result_dir <- file.path("results", "meta")
  suffix     <- ""
  out_prefix <- "05_meta_analysis"
} else {
  input_dir  <- file.path("results", "deg_tables_unshrunk")
  result_dir <- file.path("results", "meta_unshrunk")
  suffix     <- "_unshrunk"
  out_prefix <- "05_meta_analysis_unshrunk"
}
info_dir <- file.path("results", "session_info")

for (d in c(result_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

files <- paste0(c("GSE167459", "GSE253476", "GSE276647", "GSE283560"),
                "_7d", suffix, ".csv")

missing_files <- files[!file.exists(file.path(input_dir, files))]
if (length(missing_files) > 0)
  stop("Missing input file(s) in ", input_dir, ": ",
       paste(missing_files, collapse = ", "),
       ". Run scripts 01 to 04 first. Current folder: ", getwd())

# The unshrunk arm needs unshrunk fold changes, which have a 'stat' column.
if (arm == "unshrunk") {
  for (f in files)
    if (!"stat" %in% names(read.csv(file.path(input_dir, f), nrows = 1)))
      stop(f, " has no 'stat' column, so it is not an unshrunk result.")
  cat("All 4 input files are unshrunk (stat column present).\n")
}


# 1. Read the four result tables ---------------------------------------------
all_degs <- map_df(files, function(f) {
  df <- read.csv(file.path(input_dir, f))
  df$study <- tools::file_path_sans_ext(f)
  df
})

cat("Combined rows across the 4 datasets:", nrow(all_degs), "\n")
stopifnot(all(c("SYMBOL", "log2FoldChange", "lfcSE") %in% colnames(all_degs)))


# 2. Prepare the meta-analysis input -----------------------------------------
meta_input <- all_degs %>%
  filter(!is.na(SYMBOL),
         !is.na(log2FoldChange),
         !is.na(lfcSE),
         lfcSE > 0) %>%
  distinct(SYMBOL, study, .keep_all = TRUE)

# Keep genes observed in at least 2 datasets (needed for pooling)
meta_input_k2 <- meta_input %>%
  group_by(SYMBOL) %>%
  filter(n() >= 2) %>%
  ungroup()


# 3. Random-effects model per gene -------------------------------------------
pooled_row <- function(df_gene, fit, method_label, i2, tau2) {
  tibble(
    SYMBOL     = df_gene$SYMBOL[1],
    ENSEMBL    = {
      ens <- df_gene$ENSEMBL[!is.na(df_gene$ENSEMBL)]
      if (length(ens) == 0) NA_character_ else ens[1]
    },
    k          = fit$k,
    meta_logFC = as.numeric(fit$b),
    meta_SE    = fit$se,
    meta_z     = fit$zval,
    meta_p     = fit$pval,
    ci_lb      = fit$ci.lb,
    ci_ub      = fit$ci.ub,
    I2         = i2,
    tau2       = tau2,
    method     = method_label
  )
}

meta_per_gene <- function(df_gene) {
  # 1) random-effects, REML
  fit_re <- try(rma.uni(yi = df_gene$log2FoldChange, sei = df_gene$lfcSE,
                        method = "REML"), silent = TRUE)
  if (!inherits(fit_re, "try-error"))
    return(pooled_row(df_gene, fit_re, "REML", fit_re$I2, fit_re$tau2))

  # 2) fixed-effect fallback (I2 is not defined, tau2 is 0 by definition)
  fit_fe <- try(rma.uni(yi = df_gene$log2FoldChange, sei = df_gene$lfcSE,
                        method = "FE"), silent = TRUE)
  if (!inherits(fit_fe, "try-error"))
    return(pooled_row(df_gene, fit_fe, "FE", NA_real_, 0))

  # 3) neither model converged: keep the gene with NA estimates
  tibble(SYMBOL = df_gene$SYMBOL[1],
         ENSEMBL = {
           ens <- df_gene$ENSEMBL[!is.na(df_gene$ENSEMBL)]
           if (length(ens) == 0) NA_character_ else ens[1]
         },
         k = nrow(df_gene),
         meta_logFC = NA_real_, meta_SE = NA_real_, meta_z = NA_real_,
         meta_p = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
         I2 = NA_real_, tau2 = NA_real_, method = "FAILED")
}

meta_list <- meta_input_k2 %>%
  split(.$SYMBOL) %>%
  map_dfr(meta_per_gene)

meta_res <- meta_list %>%
  mutate(meta_padj = p.adjust(meta_p, method = "BH"),
         meta_stat = sign(meta_logFC) * -log10(meta_p)) %>%
  arrange(meta_padj)

cat("Model used per gene:\n")
print(table(meta_list$method))


# 4. Save ---------------------------------------------------------------------
meta_deg <- meta_res %>%
  filter(!is.na(meta_padj), meta_padj < 0.05)

write.csv(meta_res, file.path(result_dir, paste0("Meta7dpi_all_genes_k2", suffix, ".csv")),
          row.names = FALSE)
write.csv(meta_deg, file.path(result_dir, paste0("Meta7dpi_metaDEGs_FDR0.05", suffix, ".csv")),
          row.names = FALSE)


# 5. Summary and checks -------------------------------------------------------
cat("\nAnalyzable genes (k >= 2):", nrow(meta_res), "\n")
cat("Meta-DEGs (meta_padj < 0.05):", nrow(meta_deg), "\n")
cat("  Upregulated  :", sum(meta_deg$meta_logFC > 0), "\n")
cat("  Downregulated:", sum(meta_deg$meta_logFC < 0), "\n")
cat("Saved to:", result_dir, "\n")

stopifnot(nrow(meta_res) > 0, nrow(meta_deg) > 0)

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, paste0(out_prefix, "_sessionInfo.txt")))
