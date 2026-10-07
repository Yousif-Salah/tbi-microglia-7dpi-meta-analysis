# =============================================================================
# 01_DEG_GSE167459.R
# Per-dataset differential expression: GSE167459 (mouse microglia, FPI vs Sham,
# 7 days post-injury)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 1 of 4 per-dataset scripts. Their outputs are pooled by
#           05_meta_analysis.R
#
# Method
#   - Samples : 6 microglia samples (3 Sham, 3 FPI), selected by sample code:
#               the count-file columns whose code contains "M_S". Sample
#               1301M_S75 matches this pattern but is an astrocyte sample
#               according to its GEO record (GSM5105505), so it is excluded.
#               GEO lists 10 microglia samples for this series (5 Sham, 5 FPI).
#               The 4 samples that this pattern does not select (codes
#               308J-A, 310J-A, 309J-A, 315J-A) formed a separate batch on the
#               second principal component. The batch is inferred from the
#               sample codes, because GEO has no batch field. They are not part
#               of the primary analysis. The analysis with all 10 samples, with
#               and without batch in the design, is script 20 (Table S18).
#   - Filter  : genes kept if rowSums(counts) >= 10, applied before DESeq()
#   - Model   : DESeq2, design ~ group, Sham as the reference level
#   - LFC     : apeglm shrinkage (lfcShrink), contrast FPI vs Sham
#   - Same filter, shrinkage and design are used for all four datasets.
#
# Input (download from GEO accession GSE167459, place in data/GSE167459/)
#   data/GSE167459/GSE167459_annotated_combined.counts.txt.gz
#   (an unzipped copy named ...counts.txt is also accepted)
#   Sample metadata are retrieved from GEO at run time (internet required).
#
# Output
#   results/deg_tables/GSE167459_7d.csv
#     Row names = Ensembl ID. Columns: baseMean, log2FoldChange, lfcSE, pvalue,
#     padj, ENSEMBL, SYMBOL. No 'stat' column (apeglm does not return one).
#   results/deg_tables_unshrunk/GSE167459_7d_unshrunk.csv
#     Same genes and samples, fold change not shrunk (has a 'stat' column).
#     Used only by the shrinkage sensitivity analysis (script 08).
#   results/session_info/01_DEG_GSE167459_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/01_DEG_GSE167459.R
# Packages: GEOquery, DESeq2, apeglm, org.Mm.eg.db (Bioconductor)
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
})

# The count file may be compressed (.gz) or already unzipped (.txt)
count_candidates <- file.path("data", "GSE167459",
                              c("GSE167459_annotated_combined.counts.txt.gz",
                                "GSE167459_annotated_combined.counts.txt"))
count_file <- count_candidates[file.exists(count_candidates)][1]
out_dir    <- file.path("results", "deg_tables")
info_dir   <- file.path("results", "session_info")
cache_dir  <- file.path("data", "geo_cache")

if (is.na(count_file))
  stop("Count file not found. Expected one of:\n  ",
       paste(count_candidates, collapse = "\n  "),
       "\nDownload it from GEO (GSE167459) and run this script from the ",
       "repository root. Current folder: ", getwd())

for (d in c(out_dir, info_dir, cache_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# 1. Load count matrix -------------------------------------------------------
counts_raw <- read.table(count_file, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)

rownames(counts_raw) <- counts_raw$id
counts_all <- as.matrix(
  counts_raw[, !colnames(counts_raw) %in% c("id", "symbol")]
)
storage.mode(counts_all) <- "integer"

cat("Count matrix: ", nrow(counts_all), " genes x ", ncol(counts_all),
    " samples\n", sep = "")


# 2. Select microglia samples and assign groups ------------------------------
# Candidate microglia columns by name, minus the astrocyte sample (see header).
mg_cols   <- grep("M_S", colnames(counts_all), value = TRUE)
mg_cols   <- mg_cols[mg_cols != "1301M_S75"]
counts_mg <- counts_all[, mg_cols, drop = FALSE]

# Confirm each sample against its GEO record (title must contain "MG").
gset <- getGEO("GSE167459", GSEMatrix = TRUE, getGPL = FALSE,
               destdir = cache_dir)[[1]]
meta <- pData(gset)

find_hits <- function(s) {
  ix <- which(apply(meta, 1, function(r) any(grepl(s, r, fixed = TRUE))))
  ix[grepl("MG", meta$title[ix], ignore.case = TRUE)]
}

hits <- lapply(colnames(counts_mg), find_hits)

sample_map <- do.call(rbind, lapply(seq_along(hits), function(i) {
  h <- hits[[i]]
  if (length(h) == 1) {
    data.frame(sample = colnames(counts_mg)[i],
               gsm    = rownames(meta)[h],
               title  = meta$title[h],
               stringsAsFactors = FALSE)
  } else {
    data.frame(sample = colnames(counts_mg)[i],
               gsm = NA_character_, title = NA_character_,
               stringsAsFactors = FALSE)
  }
}))

sample_map$group <- ifelse(grepl("\\bFPI\\b", sample_map$title, ignore.case = TRUE), "FPI",
                    ifelse(grepl("\\bSham\\b", sample_map$title, ignore.case = TRUE), "Sham",
                           NA_character_))
sample_map <- sample_map[!is.na(sample_map$group), , drop = FALSE]

stopifnot(!any(grepl("Astro", sample_map$title, ignore.case = TRUE)))

counts_mg <- counts_mg[, sample_map$sample, drop = FALSE]
meta_mg   <- data.frame(row.names = sample_map$sample,
                        group = factor(sample_map$group,
                                       levels = c("Sham", "FPI")))
stopifnot(identical(colnames(counts_mg), rownames(meta_mg)))

cat("\nSamples used:\n")
print(sample_map[, c("sample", "gsm", "group")], row.names = FALSE)

# Expected design: 3 Sham and 3 FPI
stopifnot(all(table(meta_mg$group) == 3))


# 3. DESeq2 and shrinkage ----------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts_mg,
                              colData   = meta_mg,
                              design    = ~ group)

cat("\nGenes before filter:", nrow(dds), "\n")
dds <- dds[rowSums(counts(dds)) >= 10, ]
cat("Genes after rowSums(counts) >= 10:", nrow(dds), "\n")

dds <- DESeq(dds)

# Save the fitted DESeq2 object. Later scripts (purity check, sample quality
# control, leave-one-sample-out) reuse it instead of rebuilding the dataset.
dds_dir <- file.path("results", "dds")
dir.create(dds_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(dds, file.path(dds_dir, "GSE167459_dds.rds"))

coef_name <- resultsNames(dds)[grep("^group_.*_vs_", resultsNames(dds))][1]
cat("Shrinkage coefficient:", coef_name, "\n")
stopifnot(coef_name == "group_FPI_vs_Sham")

res <- lfcShrink(dds, coef = coef_name, type = "apeglm")
res_unshrunk <- results(dds, name = coef_name)   # for script 08


# 4. Annotate and save -------------------------------------------------------
ens_clean <- sub("\\..*$", "", rownames(res))
sym <- mapIds(org.Mm.eg.db, keys = ens_clean,
              keytype = "ENSEMBL", column = "SYMBOL",
              multiVals = "first")

res_df         <- as.data.frame(res)
res_df$ENSEMBL <- ens_clean
res_df$SYMBOL  <- unname(sym)
res_df         <- res_df[order(res_df$padj), ]

# Pooling needs log2FoldChange and lfcSE; apeglm output has no 'stat' column.
stopifnot(all(c("log2FoldChange", "lfcSE", "SYMBOL", "ENSEMBL") %in% colnames(res_df)))
stopifnot(!"stat" %in% colnames(res_df))

out_path <- file.path(out_dir, "GSE167459_7d.csv")
write.csv(res_df, out_path, row.names = TRUE)

cat("\nRows written:", nrow(res_df),
    "| padj < 0.05:", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("Saved:", out_path, "\n")

# Unshrunk (maximum likelihood) table for the sensitivity analysis (script 08).
# Same samples, filter, model and annotation; only the fold change is not shrunk.
unshrunk_dir <- file.path("results", "deg_tables_unshrunk")
dir.create(unshrunk_dir, recursive = TRUE, showWarnings = FALSE)
u_tab <- as.data.frame(res_unshrunk)[rownames(res_df), ]
res_u <- cbind(u_tab[, c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj")],
               res_df[, c("ENSEMBL", "SYMBOL")])
stopifnot("stat" %in% colnames(res_u), nrow(res_u) == nrow(res_df))
u_path <- file.path(unshrunk_dir, "GSE167459_7d_unshrunk.csv")
write.csv(res_u, u_path, row.names = FALSE)
cat("Saved unshrunk table:", u_path, "| rows:", nrow(res_u), "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "01_DEG_GSE167459_sessionInfo.txt"))
