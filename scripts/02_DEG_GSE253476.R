# =============================================================================
# 02_DEG_GSE253476.R
# Per-dataset differential expression: GSE253476 (mouse microglia, TBI vs Sham,
# 7 days post-injury)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 2 of 4 per-dataset scripts. Their outputs are pooled by
#           05_meta_analysis.R
#
# Method
#   - Samples : 3 Sham (Sham_Mg_*) and 3 TBI at 7 dpi (WT_7d_Mg_*); no exclusions
#   - Counts  : salmon-merged gene counts (scaled, fractional); rounded to
#               integers for DESeq2
#   - Filter  : genes kept if rowSums(counts) >= 10, applied before DESeq()
#   - Model   : DESeq2, design ~ condition, Sham as the reference level
#   - LFC     : apeglm shrinkage (lfcShrink), contrast TBI_7d vs Sham
#   - Same filter, shrinkage and design are used for all four datasets.
#
# Input (download from GEO accession GSE253476, place in data/GSE253476/)
#   data/GSE253476/GSE253476_salmon.merged.gene_counts_scaled_GEO.txt
#   (a compressed copy ending in .txt.gz is also accepted; reading it needs
#   the R package R.utils)
#
# Output
#   results/deg_tables/GSE253476_7d.csv
#     Columns: baseMean, log2FoldChange, lfcSE, pvalue, padj, ENSEMBL, SYMBOL,
#     ENTREZID. No row names. No 'stat' column (apeglm does not return one).
#   results/deg_tables_unshrunk/GSE253476_7d_unshrunk.csv
#     Same genes and samples, fold change not shrunk (has a 'stat' column).
#     Used only by the shrinkage sensitivity analysis (script 08).
#   results/session_info/02_DEG_GSE253476_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/02_DEG_GSE253476.R
# Packages: data.table, DESeq2, apeglm, org.Mm.eg.db
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
})

# The count file may be plain text (.txt) or compressed (.txt.gz)
count_candidates <- file.path("data", "GSE253476",
                              c("GSE253476_salmon.merged.gene_counts_scaled_GEO.txt",
                                "GSE253476_salmon.merged.gene_counts_scaled_GEO.txt.gz"))
count_file <- count_candidates[file.exists(count_candidates)][1]
out_dir    <- file.path("results", "deg_tables")
info_dir   <- file.path("results", "session_info")

if (is.na(count_file))
  stop("Count file not found. Expected one of:\n  ",
       paste(count_candidates, collapse = "\n  "),
       "\nDownload it from GEO (GSE253476) and run this script from the ",
       "repository root. Current folder: ", getwd())

for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# 1. Load count matrix -------------------------------------------------------
dt <- fread(count_file, data.table = FALSE, check.names = FALSE)

gene_col <- which(names(dt) %in% c("Gene", "GeneID", "gene_id", "Geneid",
                                   "ENSEMBL", "Ensembl"))[1]
if (is.na(gene_col)) gene_col <- 1

genes <- as.character(dt[[gene_col]])
ok    <- !is.na(genes) & nzchar(genes)
dt    <- dt[ok, , drop = FALSE]
genes <- make.unique(genes[ok])

counts <- as.matrix(dt[, -gene_col, drop = FALSE])
mode(counts) <- "numeric"
stopifnot(nrow(counts) == length(genes))
rownames(counts) <- genes

cat("Count matrix: ", nrow(counts), " genes x ", ncol(counts),
    " samples\n", sep = "")


# 2. Select samples and assign groups ----------------------------------------
keep <- grep("^(Sham_Mg_|WT_7d_Mg_)", colnames(counts), value = TRUE)
cts  <- counts[, keep, drop = FALSE]

condition <- ifelse(grepl("^Sham_Mg_", keep), "Sham", "TBI_7d")
coldata   <- data.frame(row.names = keep,
                        condition = factor(condition,
                                           levels = c("Sham", "TBI_7d")))

cat("\nSamples used:\n")
print(table(coldata$condition))
stopifnot(sum(coldata$condition == "Sham")   == 3,
          sum(coldata$condition == "TBI_7d") == 3)

# Remove Ensembl version suffixes (for example ENSMUSG00000015340.12)
rownames(cts) <- sub("\\.\\d+$", "", rownames(cts))


# 3. DESeq2 and shrinkage ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(cts), colData = coldata,
                              design = ~ condition)

cat("\nGenes before filter:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after rowSums(counts) >= 10:", nrow(dds), "\n")

dds <- DESeq(dds)

# Save the fitted DESeq2 object. Later scripts (purity check, sample quality
# control, leave-one-sample-out) reuse it instead of rebuilding the dataset.
dds_dir <- file.path("results", "dds")
dir.create(dds_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds, file.path(dds_dir, "GSE253476_dds.rds"))

coef_name <- grep("condition_TBI_7d_vs_Sham", resultsNames(dds), value = TRUE)
stopifnot(length(coef_name) == 1)
cat("Shrinkage coefficient:", coef_name, "\n")

res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_unshrunk <- results(dds, name = coef_name)   # for script 08


# 4. Annotate and save -------------------------------------------------------
res_tbl <- as.data.frame(res)
res_tbl$ENSEMBL <- rownames(res_tbl)

res_tbl$SYMBOL   <- mapIds(org.Mm.eg.db, keys = res_tbl$ENSEMBL,
                           keytype = "ENSEMBL", column = "SYMBOL",
                           multiVals = "first")
res_tbl$ENTREZID <- mapIds(org.Mm.eg.db, keys = res_tbl$ENSEMBL,
                           keytype = "ENSEMBL", column = "ENTREZID",
                           multiVals = "first")
res_tbl <- res_tbl[order(res_tbl$padj), ]

# Pooling needs log2FoldChange and lfcSE; apeglm output has no 'stat' column.
stopifnot(all(c("log2FoldChange", "lfcSE", "SYMBOL", "ENSEMBL") %in% colnames(res_tbl)))
stopifnot(!"stat" %in% colnames(res_tbl))

out_path <- file.path(out_dir, "GSE253476_7d.csv")
write.csv(res_tbl, out_path, row.names = FALSE)

cat("\nRows written:", nrow(res_tbl),
    "| padj < 0.05:", sum(res_tbl$padj < 0.05, na.rm = TRUE), "\n")
cat("Saved:", out_path, "\n")

# Unshrunk (maximum likelihood) table for the sensitivity analysis (script 08).
# Same samples, filter, model and annotation; only the fold change is not shrunk.
unshrunk_dir <- file.path("results", "deg_tables_unshrunk")
dir.create(unshrunk_dir, recursive = TRUE, showWarnings = FALSE)
u_tab <- as.data.frame(res_unshrunk)[rownames(res_tbl), ]
res_u <- cbind(u_tab[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
               res_tbl[, c("ENSEMBL", "SYMBOL", "ENTREZID")])
stopifnot("stat" %in% colnames(res_u), nrow(res_u) == nrow(res_tbl))
u_path <- file.path(unshrunk_dir, "GSE253476_7d_unshrunk.csv")
write.csv(res_u, u_path, row.names = FALSE)
cat("Saved unshrunk table:", u_path, "| rows:", nrow(res_u), "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "02_DEG_GSE253476_sessionInfo.txt"))
