# =============================================================================
# 04_DEG_GSE283560.R
# Per-dataset differential expression: GSE283560 (mouse microglia, TBI vs
# Control, 7 days post-injury)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 4 of 4 per-dataset scripts. Their outputs are pooled by
#           05_meta_analysis.R
#
# Method
#   - Quantification : salmon quant.sf.gz files, 3 sequencing lanes per sample.
#                      Transcript counts are summarised to genes with tximport
#                      (transcript-to-gene table from the Ensembl GRCm39 GTF,
#                      transcript versions ignored).
#   - Lanes          : the 3 lanes of each sample are summed (counts) or
#                      averaged (abundance, length) into one column per sample.
#   - Samples        : microglia samples at 7 dpi, 3 Control and 8 TBI. The
#                      7 dpi subset is taken BEFORE the DESeq2 model is built,
#                      so size factors are estimated from these samples only.
#   - Filter         : genes kept if rowSums(counts) >= 10, applied before
#                      DESeq()
#   - Model          : DESeq2 on rounded gene counts, design ~ condition,
#                      Control as the reference level. Automatic outlier
#                      replacement is switched off (minReplicatesForReplace =
#                      Inf), because the TBI group has 8 samples and would
#                      otherwise trigger it. In the other three datasets every
#                      group has fewer than 7 samples, so it never applies.
#   - LFC            : apeglm shrinkage (lfcShrink), contrast TBI vs Control
#   - Symbols        : Ensembl IDs are mapped to symbols with org.Mm.eg.db;
#                      unmapped IDs keep the Ensembl ID. When several genes
#                      share a symbol, the row with the highest baseMean is kept.
#   - Same filter, shrinkage and design are used for all four datasets.
#
# Input
#   data/GSE283560/GSE283560_RAW.tar  GEO supplementary archive; it is extracted
#                                   automatically to data/GSE283560/GSE283560_RAW/
#                                   (182 quant.sf.gz files in total; the 60 files of
#                                   the 20 microglia samples, 3 lanes each, are used)
#   data/reference/Mus_musculus.GRCm39.110.gtf.gz
#                                   Ensembl release 110, mouse GRCm39 annotation
#   Sample metadata are retrieved from GEO at run time (internet required).
#
# Output
#   results/deg_tables/GSE283560_7d.csv
#     Columns: baseMean, log2FoldChange, lfcSE, pvalue, padj, ENSEMBL, SYMBOL.
#     No row names. No 'stat' column (apeglm does not return one).
#   results/deg_tables_unshrunk/GSE283560_7d_unshrunk.csv
#     Same genes and samples, fold change not shrunk (has a 'stat' column).
#     Used only by the shrinkage sensitivity analysis (script 08).
#   results/counts/GSE283560_7dpi_counts.csv
#     Gene counts of the 11 samples before the low-count filter (rows = Ensembl
#     IDs, 3 lanes summed, rounded). Input of the WGCNA script (21).
#   results/counts/GSE283560_7dpi_sample_table.csv   (sample, condition)
#   results/dds/GSE283560_dds.rds   fitted DESeq2 object (used by scripts 04b, 15-17)
#   results/session_info/04_DEG_GSE283560_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/04_DEG_GSE283560.R
# Packages: GEOquery, DESeq2, apeglm, tximport, GenomicFeatures,
#           AnnotationDbi, org.Mm.eg.db
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(apeglm)
  library(tximport)
  library(GenomicFeatures)
  library(AnnotationDbi)
  library(org.Mm.eg.db)
})

raw_dir   <- file.path("data", "GSE283560", "GSE283560_RAW")
gtf_dir   <- file.path("data", "reference")
cache_dir <- file.path("data", "geo_cache")
out_dir   <- file.path("results", "deg_tables")
info_dir  <- file.path("results", "session_info")

for (d in c(out_dir, info_dir, cache_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# 1. Sample metadata from GEO ------------------------------------------------
options(timeout = max(600, getOption("timeout")))
gse      <- getGEO("GSE283560", GSEMatrix = TRUE, getGPL = FALSE,
                   destdir = cache_dir)
metadata <- pData(gse[[1]])

metadata$celltype  <- factor(metadata[, "cell type:ch1"])
metadata$condition <- factor(ifelse(metadata[, "treatment:ch1"] == "TBI",
                                    "TBI", "Control"),
                             levels = c("Control", "TBI"))
metadata$time      <- factor(metadata[, "time:ch1"],
                             levels = c("0dpi", "7dpi", "30dpi"))


# 2. Transcript-to-gene table from the GTF -----------------------------------
gtf_candidates <- if (dir.exists(gtf_dir))
  list.files(gtf_dir, pattern = "GRCm39.*\\.gtf(\\.gz)?$",
             recursive = TRUE, full.names = TRUE, ignore.case = TRUE) else character(0)

if (length(gtf_candidates) == 0)
  stop("No GRCm39 GTF file found in ", gtf_dir, ". Download the Ensembl mouse ",
       "GRCm39 GTF (release 110, see README) and place it there, for example ",
       "data/reference/Mus_musculus.GRCm39.110.gtf.gz. Current folder: ", getwd())
if (length(gtf_candidates) > 1) {
  cat("Several GTF files found, using the first one:\n")
  print(gtf_candidates)
}

gtf_path <- gtf_candidates[1]
cat("GTF file:", gtf_path, "\n")

txdb    <- makeTxDbFromGFF(gtf_path)
k       <- AnnotationDbi::keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, keys = k, keytype = "TXNAME",
                                 columns = "GENEID")


# 3. Import all lane files of every microglia sample -------------------------
# If only the GEO archive (GSE283560_RAW.tar) is present, extract it first.
raw_tar <- file.path("data", "GSE283560", "GSE283560_RAW.tar")
if (!dir.exists(raw_dir) && file.exists(raw_tar)) {
  cat("Extracting", raw_tar, "...\n")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  untar(raw_tar, exdir = raw_dir)
}

files_raw <- list.files(raw_dir, pattern = "quant.sf.gz",
                        full.names = TRUE, recursive = TRUE)
if (length(files_raw) == 0)
  stop("No quant.sf.gz files found under ", raw_dir,
       ". Place GSE283560_RAW.tar (from GEO) in data/GSE283560/ or extract ",
       "it there. Current folder: ", getwd())

# Name each file by its GSM ID; every GSM appears once per lane.
names(files_raw) <- sub("_.*", "", basename(files_raw))

micro_metadata <- metadata[metadata$celltype == "Microglia", ]

# %in% keeps every lane file of a sample (match() would keep only the first).
micro_files <- files_raw[names(files_raw) %in% rownames(micro_metadata)]
cat("Files found:", length(files_raw), "| microglia files selected:",
    length(micro_files), "\n")
stopifnot(length(micro_files) > 0, length(micro_files) %% 3 == 0)

sample_table <- data.frame(files  = micro_files,
                           sample = names(micro_files))

txi_micro <- tximport(sample_table$files, type = "salmon",
                      tx2gene = tx2gene, ignoreTxVersion = TRUE)


# 4. Combine the lanes of each sample ----------------------------------------
group <- sample_table$sample
txi_micro$counts    <- as.matrix(sapply(unique(group), function(g)
  rowSums(txi_micro$counts[, group == g, drop = FALSE])))
txi_micro$abundance <- as.matrix(sapply(unique(group), function(g)
  rowMeans(txi_micro$abundance[, group == g, drop = FALSE])))
txi_micro$length    <- as.matrix(sapply(unique(group), function(g)
  rowMeans(txi_micro$length[, group == g, drop = FALSE])))

micro_metadata <- micro_metadata[unique(group), ]

cat("Cybb count in GSM8666315 (sum of 3 lanes):",
    round(txi_micro$counts["ENSMUSG00000015340", "GSM8666315"], 3), "\n")


# 5. Subset to 7 dpi, build the DESeq2 object, filter ------------------------
is_7dpi <- micro_metadata$time == "7dpi"

txi_7dpi           <- txi_micro
txi_7dpi$counts    <- txi_micro$counts[, is_7dpi, drop = FALSE]
txi_7dpi$abundance <- txi_micro$abundance[, is_7dpi, drop = FALSE]
txi_7dpi$length    <- txi_micro$length[, is_7dpi, drop = FALSE]

meta_7dpi           <- micro_metadata[is_7dpi, ]
meta_7dpi$condition <- droplevels(meta_7dpi$condition)
meta_7dpi$condition <- relevel(meta_7dpi$condition, "Control")

cat("\n7 dpi samples:", nrow(meta_7dpi), "\n")
print(table(meta_7dpi$condition))
stopifnot(nrow(meta_7dpi) == 11,
          sum(meta_7dpi$condition == "Control") == 3,
          sum(meta_7dpi$condition == "TBI")     == 8)

# Rounded gene counts (not DESeqDataSetFromTximport), as for the other datasets
counts_7dpi <- round(txi_7dpi$counts)

# Save the gene counts before the low-count filter (3 lanes summed, 11 samples)
# and the sample table. The WGCNA script (21) starts from these two files.
counts_dir <- file.path("results", "counts")
dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(counts_7dpi, file.path(counts_dir, "GSE283560_7dpi_counts.csv"))
write.csv(data.frame(sample    = rownames(meta_7dpi),
                     condition = as.character(meta_7dpi$condition)),
          file.path(counts_dir, "GSE283560_7dpi_sample_table.csv"),
          row.names = FALSE)

dds <- DESeqDataSetFromMatrix(countData = counts_7dpi,
                              colData   = meta_7dpi,
                              design    = ~ condition)

# Low-count filter before DESeq(), so size factors use the retained genes only
cat("\nGenes before filter:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after rowSums(counts) >= 10:", nrow(dds), "\n")


# 6. DESeq2 and shrinkage ----------------------------------------------------
dds <- DESeq(dds, minReplicatesForReplace = Inf)

# Save the fitted DESeq2 object. Later scripts (purity check, sample quality
# control, leave-one-sample-out) reuse it instead of rebuilding the dataset.
dds_dir <- file.path("results", "dds")
dir.create(dds_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds, file.path(dds_dir, "GSE283560_dds.rds"))

res_coef <- grep("TBI_vs_Control", resultsNames(dds), value = TRUE)
stopifnot(length(res_coef) == 1)
cat("Shrinkage coefficient:", res_coef, "\n")

res <- lfcShrink(dds, coef = res_coef, type = "apeglm")
res_unshrunk <- results(dds, name = res_coef)   # for script 08


# 7. Annotate, remove duplicate symbols, save --------------------------------
sym <- mapIds(org.Mm.eg.db, keys = rownames(res),
              keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")

res_df         <- as.data.frame(res)
res_df$ENSEMBL <- rownames(res_df)
res_df$SYMBOL  <- ifelse(is.na(sym), res_df$ENSEMBL, sym)
res_df         <- res_df[order(res_df$padj), ]

# Keep the row with the highest baseMean for each duplicated symbol
n_before <- nrow(res_df)
res_df   <- res_df[order(res_df$SYMBOL, -res_df$baseMean), ]
res_df   <- res_df[!duplicated(res_df$SYMBOL), ]
res_df   <- res_df[order(res_df$padj), ]
cat("Rows removed as duplicate symbols:", n_before - nrow(res_df), "\n")

# Pooling needs log2FoldChange and lfcSE; apeglm output has no 'stat' column.
stopifnot(all(c("log2FoldChange", "lfcSE", "SYMBOL", "ENSEMBL") %in% colnames(res_df)))
stopifnot(!"stat" %in% colnames(res_df))

out_path <- file.path(out_dir, "GSE283560_7d.csv")
write.csv(res_df, out_path, row.names = FALSE)

cat("\nRows written:", nrow(res_df),
    "| padj < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE),
    "| rows with Ensembl ID as symbol:", sum(res_df$SYMBOL == res_df$ENSEMBL), "\n")
cat("Saved:", out_path, "\n")

# Unshrunk (maximum likelihood) table for the sensitivity analysis (script 08).
# Same samples, filter, model and annotation; only the fold change is not shrunk.
unshrunk_dir <- file.path("results", "deg_tables_unshrunk")
dir.create(unshrunk_dir, recursive = TRUE, showWarnings = FALSE)
u_tab <- as.data.frame(res_unshrunk)[rownames(res_df), ]
res_u <- cbind(u_tab[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
               res_df[, c("ENSEMBL", "SYMBOL")])
stopifnot("stat" %in% colnames(res_u), nrow(res_u) == nrow(res_df))
u_path <- file.path(unshrunk_dir, "GSE283560_7d_unshrunk.csv")
write.csv(res_u, u_path, row.names = FALSE)
cat("Saved unshrunk table:", u_path, "| rows:", nrow(res_u), "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "04_DEG_GSE283560_sessionInfo.txt"))
