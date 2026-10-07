# =============================================================================
# 03_DEG_GSE276647.R
# Per-dataset differential expression: GSE276647 (mouse microglia, TBI vs Sham,
# 7 days post-injury)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 3 of 4 per-dataset scripts. Their outputs are pooled by
#           05_meta_analysis.R
#
# Method
#   - Samples : 4 Sham (Sham_1 to Sham_4) and 4 TBI (TBI_Iso_1 to TBI_Iso_4)
#   - Counts  : gene-level counts from the GEO supplementary table (7 dpi);
#               gene identifiers are mouse symbols or aliases, not Ensembl IDs.
#               Counts are rounded to integers for DESeq2.
#   - Filter  : genes kept if rowSums(counts) >= 10, applied before DESeq()
#   - Model   : DESeq2, design ~ condition, Sham as the reference level
#   - LFC     : apeglm shrinkage (lfcShrink), contrast TBI vs Sham
#   - Symbols : each identifier is matched to its official symbol through the
#               org.Mm.eg.db ALIAS table (first match); identifiers without a
#               match keep their original name.
#   - Same filter, shrinkage and design are used for all four datasets.
#     DESeq2 outlier replacement is not triggered here, because it applies
#     only to groups with 7 or more replicates and each group has 4.
#
# Input (download from GEO accession GSE276647, place in data/GSE276647/)
#   data/GSE276647/GSE276647_Supp_Table_2_7d.txt.gz
#   (an unzipped copy named ...Supp_Table_2_7d.txt is also accepted)
#
# Output
#   results/deg_tables/GSE276647_7d.csv
#     Columns: baseMean, log2FoldChange, lfcSE, pvalue, padj, SYMBOL, Gene.
#     No row names. No 'stat' column (apeglm does not return one).
#   results/deg_tables_unshrunk/GSE276647_7d_unshrunk.csv
#     Same genes and samples, fold change not shrunk (has a 'stat' column).
#     Used only by the shrinkage sensitivity analysis (script 08).
#   results/session_info/03_DEG_GSE276647_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/03_DEG_GSE276647.R
# Packages: data.table, DESeq2, apeglm, org.Mm.eg.db, AnnotationDbi
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
})

count_candidates <- file.path("data", "GSE276647",
                              c("GSE276647_Supp_Table_2_7d.txt.gz",
                                "GSE276647_Supp_Table_2_7d.txt"))
count_file <- count_candidates[file.exists(count_candidates)][1]
out_dir    <- file.path("results", "deg_tables")
info_dir   <- file.path("results", "session_info")

if (is.na(count_file))
  stop("Count file not found. Expected one of:\n  ",
       paste(count_candidates, collapse = "\n  "),
       "\nDownload it from GEO (GSE276647) and run this script from the ",
       "repository root. Current folder: ", getwd())

for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# 1. Load count matrix -------------------------------------------------------
deg_7d <- fread(count_file, skip = 1)
colnames(deg_7d)[1] <- "Gene"

counts <- as.data.frame(deg_7d)
rownames(counts) <- counts$Gene
counts <- counts[, -1]

cat("Count matrix: ", nrow(counts), " genes x ", ncol(counts),
    " samples\n", sep = "")


# 2. Sample table ------------------------------------------------------------
samples <- data.frame(
  row.names = c("Sham_1", "Sham_2", "Sham_3", "Sham_4",
                "TBI_Iso_1", "TBI_Iso_2", "TBI_Iso_3", "TBI_Iso_4"),
  condition = factor(c(rep("Sham", 4), rep("TBI", 4)),
                     levels = c("Sham", "TBI"))
)

stopifnot(all(rownames(samples) %in% colnames(counts)))
counts <- counts[, rownames(samples)]

cat("\nSamples used:\n")
print(table(samples$condition))
stopifnot(sum(samples$condition == "Sham") == 4,
          sum(samples$condition == "TBI")  == 4)


# 3. DESeq2 and shrinkage ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(counts),
                              colData   = samples,
                              design    = ~ condition)

cat("\nGenes before filter:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after rowSums(counts) >= 10:", nrow(dds), "\n")

dds <- DESeq(dds)

# Save the fitted DESeq2 object. Later scripts (purity check, sample quality
# control, leave-one-sample-out) reuse it instead of rebuilding the dataset.
dds_dir <- file.path("results", "dds")
dir.create(dds_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds, file.path(dds_dir, "GSE276647_dds.rds"))

res_coef <- grep("TBI_vs_Sham", resultsNames(dds), value = TRUE)
stopifnot(length(res_coef) == 1)
cat("Shrinkage coefficient:", res_coef, "\n")

res <- lfcShrink(dds, coef = res_coef, type = "apeglm")
res_unshrunk <- results(dds, name = res_coef)   # for script 08


# 4. Annotate and save -------------------------------------------------------
res_tbl <- as.data.frame(res)
res_tbl$Gene <- rownames(res_tbl)

# Match each identifier to its official symbol through the ALIAS table
mapped <- AnnotationDbi::select(org.Mm.eg.db,
                                keys    = res_tbl$Gene,
                                keytype = "ALIAS",
                                columns = "SYMBOL")
mapped <- mapped[!duplicated(mapped$ALIAS), ]

res_tbl$SYMBOL <- mapped$SYMBOL[match(res_tbl$Gene, mapped$ALIAS)]
res_tbl$SYMBOL[is.na(res_tbl$SYMBOL)] <- res_tbl$Gene[is.na(res_tbl$SYMBOL)]

res_tbl <- res_tbl[, c("baseMean", "log2FoldChange", "lfcSE", "pvalue",
                       "padj", "SYMBOL", "Gene")]
res_tbl <- res_tbl[order(res_tbl$padj), ]

# Pooling needs log2FoldChange and lfcSE; apeglm output has no 'stat' column.
stopifnot(all(c("log2FoldChange", "lfcSE", "SYMBOL", "Gene") %in% colnames(res_tbl)))
stopifnot(!"stat" %in% colnames(res_tbl))

out_path <- file.path(out_dir, "GSE276647_7d.csv")
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
               res_tbl[, c("SYMBOL", "Gene")])
stopifnot("stat" %in% colnames(res_u), nrow(res_u) == nrow(res_tbl))
u_path <- file.path(unshrunk_dir, "GSE276647_7d_unshrunk.csv")
write.csv(res_u, u_path, row.names = FALSE)
cat("Saved unshrunk table:", u_path, "| rows:", nrow(res_u), "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "03_DEG_GSE276647_sessionInfo.txt"))
