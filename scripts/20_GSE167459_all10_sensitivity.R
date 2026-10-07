# =============================================================================
# 20_GSE167459_all10_sensitivity.R
# Sensitivity analysis: GSE167459 with all 10 microglia samples (5 Sham, 5 FPI)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 20. Tests whether the pooled results change when GSE167459 is
#           analysed with all 10 of its microglia samples instead of the 6
#           samples (3 Sham, 3 FPI) used in the primary analysis (script 01).
#           It does NOT change any primary file.
#
# Method
#   - Samples : all microglia samples of GSE167459 chosen from the GEO record
#               (cell type = Microglia; group from the sample title), not from
#               the column names. The astrocyte samples are not used.
#   - Batch   : the four samples with a "J-A" code (308, 309, 310, 315) were
#               prepared separately from the other six ("M" codes). GEO has no
#               batch field, so the batch label comes from the sample code.
#   - Filter  : genes kept if rowSums(counts) >= 10 over the 10 samples,
#               applied before DESeq() (same rule as scripts 01 to 04)
#   - Models  : DESeq2 with apeglm shrinkage (FPI vs Sham), two designs
#                 all10_group : ~ group
#                 all10_batch : ~ batch + group
#   - Pooling : exactly the rule of script 05 (one row per gene and dataset,
#               lfcSE > 0, k >= 2, random effects with REML, FE fallback,
#               Benjamini-Hochberg) with the three unchanged datasets from
#               results/deg_tables/.
#   - Check   : the same pooling code is first run on the primary tables and
#               must reproduce results/meta/Meta7dpi_all_genes_k2.csv.
#   - GSEA    : Hallmark and the 14-gene BTK/NOX2 panel, same settings as
#               scripts 13 and 14 (set.seed(1) before each run).
#
# Input
#   data/GSE167459/GSE167459_annotated_combined.counts.txt.gz (or unzipped .txt)
#   data/geo_cache/GSE167459_series_matrix.txt.gz (read offline if present,
#     otherwise downloaded from GEO)
#   results/deg_tables/GSE167459_7d.csv, GSE253476_7d.csv, GSE276647_7d.csv,
#   GSE283560_7d.csv, and results/meta/Meta7dpi_all_genes_k2.csv (scripts 01-05)
#
# Output (all in results/sensitivity_all10/)
#   GSE167459_7d_all10_group.csv, GSE167459_7d_all10_batch.csv
#   GSE167459_all10_sample_table.csv
#   Figure_S9_GSE167459_all10_PCA.tiff   (Figure S9, 85 mm wide, 300 dpi)
#   Meta7dpi_all_genes_k2_all10_group.csv, Meta7dpi_all_genes_k2_all10_batch.csv
#   Hallmark_all10_group.csv, Hallmark_all10_batch.csv
#   Summary_primary_vs_all10.csv         (one row per analysis)
#   Table_S18_GSE167459_all10_sensitivity.csv   (Table S18, same numbers, one
#                                        column per analysis)
#   results/session_info/20_GSE167459_all10_sensitivity_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/20_GSE167459_all10_sensitivity.R
# Packages: GEOquery, DESeq2, apeglm, org.Mm.eg.db, dplyr, purrr, tibble,
#           metafor, fgsea, msigdbr (25.1.1), clusterProfiler, ggplot2
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(GEOquery)
  library(DESeq2)
  library(apeglm)
  library(org.Mm.eg.db)
  library(dplyr)
  library(purrr)
  library(tibble)
  library(metafor)
  library(fgsea)
  library(msigdbr)
  library(clusterProfiler)
  library(ggplot2)
})

set.seed(1)

count_candidates <- file.path("data", "GSE167459",
                              c("GSE167459_annotated_combined.counts.txt.gz",
                                "GSE167459_annotated_combined.counts.txt"))
count_file  <- count_candidates[file.exists(count_candidates)][1]
cache_dir   <- file.path("data", "geo_cache")
matrix_file <- file.path(cache_dir, "GSE167459_series_matrix.txt.gz")
tab_dir     <- file.path("results", "deg_tables")
meta_file   <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
out_dir     <- file.path("results", "sensitivity_all10")
info_dir    <- file.path("results", "session_info")

if (is.na(count_file))
  stop("Count file not found. Expected one of:\n  ",
       paste(count_candidates, collapse = "\n  "),
       "\nRun this script from the repository root. Current folder: ", getwd())

other_sets <- c("GSE253476", "GSE276647", "GSE283560")
needed <- c(file.path(tab_dir, paste0(c("GSE167459", other_sets), "_7d.csv")), meta_file)
if (any(!file.exists(needed)))
  stop("Missing input file(s): ", paste(needed[!file.exists(needed)], collapse = ", "),
       ". Run scripts 01 to 05 first.")

for (d in c(out_dir, info_dir, cache_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)


# 1. Select microglia samples from the GEO record ----------------------------
counts_raw <- read.table(count_file, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)
rownames(counts_raw) <- counts_raw$id
counts_all <- as.matrix(counts_raw[, !colnames(counts_raw) %in% c("id", "symbol")])
storage.mode(counts_all) <- "integer"
cat("Count matrix:", nrow(counts_all), "genes x", ncol(counts_all), "samples\n")

gset <- if (file.exists(matrix_file)) {
  getGEO(filename = matrix_file, GSEMatrix = TRUE, getGPL = FALSE)
} else {
  getGEO("GSE167459", GSEMatrix = TRUE, getGPL = FALSE, destdir = cache_dir)[[1]]
}
geo <- pData(gset)

is_mg <- grepl("microglia", geo[["cell type:ch1"]], ignore.case = TRUE)
sel <- data.frame(sample = geo$description[is_mg],
                  gsm    = rownames(geo)[is_mg],
                  title  = geo$title[is_mg],
                  stringsAsFactors = FALSE)
sel$group <- ifelse(grepl("\\bFPI\\b",  sel$title, ignore.case = TRUE), "FPI",
             ifelse(grepl("\\bSham\\b", sel$title, ignore.case = TRUE), "Sham", NA))
sel$batch <- ifelse(grepl("J-A", sel$sample, fixed = TRUE), "J-A", "M")
sel <- sel[order(sel$group, sel$title), ]

# Safety checks: 10 samples, 5 per group, all present in the count file
stopifnot(nrow(sel) == 10, !any(is.na(sel$group)),
          all(table(sel$group) == 5),
          all(sel$sample %in% colnames(counts_all)),
          !any(grepl("Astro", sel$title, ignore.case = TRUE)),
          sum(sel$batch == "J-A") == 4,
          all(table(sel$group, sel$batch)["Sham", ] == table(sel$group, sel$batch)["FPI", ]))
cat("\nSamples used (all 10 microglia):\n")
print(sel, row.names = FALSE)
cat("\nGroup by batch:\n"); print(table(sel$group, sel$batch))
write.csv(sel, file.path(out_dir, "GSE167459_all10_sample_table.csv"), row.names = FALSE)

counts_mg <- counts_all[, sel$sample, drop = FALSE]
col_data  <- data.frame(row.names = sel$sample,
                        group = factor(sel$group, levels = c("Sham", "FPI")),
                        batch = factor(sel$batch, levels = c("M", "J-A")))
stopifnot(identical(colnames(counts_mg), rownames(col_data)))


# 2. DESeq2 with two designs --------------------------------------------------
run_deseq <- function(design) {
  dds <- DESeqDataSetFromMatrix(countData = counts_mg, colData = col_data, design = design)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  cat("Genes after rowSums(counts) >= 10:", nrow(dds), "\n")
  dds <- DESeq(dds)
  stopifnot("group_FPI_vs_Sham" %in% resultsNames(dds))
  res <- lfcShrink(dds, coef = "group_FPI_vs_Sham", type = "apeglm")
  ens <- sub("\\..*$", "", rownames(res))
  sym <- mapIds(org.Mm.eg.db, keys = ens, keytype = "ENSEMBL",
                column = "SYMBOL", multiVals = "first")
  df <- as.data.frame(res)
  df$ENSEMBL <- ens
  df$SYMBOL  <- unname(sym)
  # Same ordering as script 01 (by padj), so "first row per symbol" is the same rule
  df <- df[order(df$padj), ]
  stopifnot(all(c("log2FoldChange", "lfcSE", "SYMBOL", "ENSEMBL") %in% colnames(df)))
  list(dds = dds, table = df)
}

cat("\n--- Design ~ group ---\n")
fit_group <- run_deseq(~ group)
cat("\n--- Design ~ batch + group ---\n")
fit_batch <- run_deseq(~ batch + group)

write.csv(fit_group$table, file.path(out_dir, "GSE167459_7d_all10_group.csv"), row.names = TRUE)
write.csv(fit_batch$table, file.path(out_dir, "GSE167459_7d_all10_batch.csv"), row.names = TRUE)

# PCA to see the batch (report only)
vsd <- vst(fit_group$dds, blind = TRUE)
pca0 <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
pv   <- round(100 * attr(pca0, "percentVar"))
# Rebuild the table by hand: with two grouping columns plotPCA merges them into one
# "group" column ("FPI : M"), and the colours then do not match.
pca <- data.frame(PC1   = pca0$PC1, PC2 = pca0$PC2,
                  group = factor(as.character(colData(vsd)$group), levels = c("Sham", "FPI")),
                  batch = factor(as.character(colData(vsd)$batch), levels = c("M", "J-A")))
stopifnot(!anyNA(pca$group), !anyNA(pca$batch),
          identical(as.character(pca0$name), colnames(vsd)))
# Figure S9: 85 mm wide, Arial 9 pt, 300 dpi TIFF. Points are not labelled, so
# no text can overlap; sample codes and batch are in GSE167459_all10_sample_table.csv.
fig_font <- if (.Platform$OS.type == "windows") "Arial" else "sans"
p <- ggplot(pca, aes(PC1, PC2, colour = group, shape = batch)) +
  geom_point(size = 3, stroke = 1) +
  scale_colour_manual(values = c(Sham = "#2471A3", FPI = "#C0392B"), name = "Group") +
  scale_shape_manual(values = c(M = 16, `J-A` = 17),
                     labels = c(M = "Batch M (6 samples)", `J-A` = "Batch J-A (4 samples)"),
                     name = "Batch") +
  labs(x = paste0("PC1: ", pv[1], "% variance"), y = paste0("PC2: ", pv[2], "% variance")) +
  theme_bw(base_size = 9, base_family = fig_font) +
  theme(legend.position = "bottom", legend.box = "vertical",
        legend.spacing.y = unit(0, "mm"), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "Figure_S9_GSE167459_all10_PCA.tiff"), p,
       width = 85, height = 90, units = "mm", dpi = 300, compression = "lzw")

cat("\nGSE167459 alone, shrunk results (FPI vs Sham):\n")
show_gene <- function(df, g) {
  r <- df[which(df$SYMBOL == g)[1], c("log2FoldChange", "lfcSE", "padj")]
  cat(sprintf("  %-7s log2FC %6.3f  lfcSE %5.3f  padj %.3g\n", g,
              r$log2FoldChange, r$lfcSE, r$padj))
}
prim_gse <- read.csv(file.path(tab_dir, "GSE167459_7d.csv"), stringsAsFactors = FALSE)
for (nm in c("primary 3 vs 3", "all10 ~ group", "all10 ~ batch + group")) {
  cat(nm, "\n")
  d <- switch(nm, "primary 3 vs 3" = prim_gse,
                  "all10 ~ group" = fit_group$table,
                  "all10 ~ batch + group" = fit_batch$table)
  for (g in c("Cybb", "Fcer1g", "Btk")) show_gene(d, g)
}


# 3. Pooling (the rule of script 05) -----------------------------------------
pooled_row <- function(df_gene, fit, method_label, i2, tau2) {
  tibble(SYMBOL = df_gene$SYMBOL[1],
         ENSEMBL = { ens <- df_gene$ENSEMBL[!is.na(df_gene$ENSEMBL)]
                     if (length(ens) == 0) NA_character_ else ens[1] },
         k = fit$k, meta_logFC = as.numeric(fit$b), meta_SE = fit$se,
         meta_z = fit$zval, meta_p = fit$pval, ci_lb = fit$ci.lb, ci_ub = fit$ci.ub,
         I2 = i2, tau2 = tau2, method = method_label)
}

meta_per_gene <- function(df_gene) {
  fit_re <- try(rma.uni(yi = df_gene$log2FoldChange, sei = df_gene$lfcSE,
                        method = "REML"), silent = TRUE)
  if (!inherits(fit_re, "try-error"))
    return(pooled_row(df_gene, fit_re, "REML", fit_re$I2, fit_re$tau2))
  fit_fe <- try(rma.uni(yi = df_gene$log2FoldChange, sei = df_gene$lfcSE,
                        method = "FE"), silent = TRUE)
  if (!inherits(fit_fe, "try-error"))
    return(pooled_row(df_gene, fit_fe, "FE", NA_real_, 0))
  tibble(SYMBOL = df_gene$SYMBOL[1],
         ENSEMBL = { ens <- df_gene$ENSEMBL[!is.na(df_gene$ENSEMBL)]
                     if (length(ens) == 0) NA_character_ else ens[1] },
         k = nrow(df_gene), meta_logFC = NA_real_, meta_SE = NA_real_,
         meta_z = NA_real_, meta_p = NA_real_, ci_lb = NA_real_, ci_ub = NA_real_,
         I2 = NA_real_, tau2 = NA_real_, method = "FAILED")
}

# gse167459_tab: the GSE167459 table to use; the other three come from scripts 02-04
pool_with <- function(gse167459_tab) {
  others <- map_df(other_sets, function(s) {
    df <- read.csv(file.path(tab_dir, paste0(s, "_7d.csv")), stringsAsFactors = FALSE)
    df$study <- paste0(s, "_7d")
    df
  })
  first <- gse167459_tab
  first$study <- "GSE167459_7d"
  all_degs <- bind_rows(first, others)
  meta_input <- all_degs %>%
    filter(!is.na(SYMBOL), !is.na(log2FoldChange), !is.na(lfcSE), lfcSE > 0) %>%
    distinct(SYMBOL, study, .keep_all = TRUE) %>%
    group_by(SYMBOL) %>% filter(n() >= 2) %>% ungroup()
  res <- meta_input %>% split(.$SYMBOL) %>% map_dfr(meta_per_gene)
  res %>%
    mutate(meta_padj = p.adjust(meta_p, method = "BH"),
           meta_stat = sign(meta_logFC) * -log10(meta_p)) %>%
    arrange(meta_padj)
}

cat("\nPooling check: primary tables through this script's pooling code ...\n")
check <- pool_with(prim_gse)
ref   <- read.csv(meta_file, stringsAsFactors = FALSE)
mm    <- merge(check[, c("SYMBOL", "meta_logFC", "meta_p")],
               ref[, c("SYMBOL", "meta_logFC", "meta_p")], by = "SYMBOL")
cat("  genes:", nrow(check), "| reference genes:", nrow(ref),
    "| meta-DEGs:", sum(check$meta_padj < 0.05, na.rm = TRUE), "\n")
cat("  max |difference| meta_logFC:", max(abs(mm$meta_logFC.x - mm$meta_logFC.y)),
    "| meta_p:", max(abs(mm$meta_p.x - mm$meta_p.y)), "\n")
stopifnot(nrow(check) == nrow(ref), nrow(mm) == nrow(ref),
          max(abs(mm$meta_logFC.x - mm$meta_logFC.y)) < 1e-6,
          max(abs(mm$meta_p.x - mm$meta_p.y)) < 1e-6)
cat("  Pooling check passed: identical to script 05.\n")

cat("\nPooling with all10 ~ group ...\n")
meta_group <- pool_with(fit_group$table)
cat("Pooling with all10 ~ batch + group ...\n")
meta_batch <- pool_with(fit_batch$table)
write.csv(meta_group, file.path(out_dir, "Meta7dpi_all_genes_k2_all10_group.csv"), row.names = FALSE)
write.csv(meta_batch, file.path(out_dir, "Meta7dpi_all_genes_k2_all10_batch.csv"), row.names = FALSE)


# 4. GSEA: Hallmark and the 14-gene panel (settings of scripts 13 and 14) ------
gs_h <- msigdbr(species = "Mus musculus", collection = "H")
hallmark_sets <- split(as.character(gs_h$ncbi_gene), gs_h$gs_name)

btk_genes   <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1")
nox2_genes  <- c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
panel_genes <- c(btk_genes, nox2_genes)

make_ranks <- function(m) {
  m <- m[!is.na(m$meta_p), ]
  m$rank_stat <- sign(m$meta_logFC) * -log10(m$meta_p)
  if (any(is.infinite(m$rank_stat))) {
    small <- min(m$meta_p[m$meta_p > 0], na.rm = TRUE)
    ii <- is.infinite(m$rank_stat)
    m$rank_stat[ii] <- sign(m$meta_logFC[ii]) * -log10(small)
  }
  id_map <- bitr(m$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  mm <- merge(m, id_map, by = "SYMBOL") %>%
    group_by(ENTREZID) %>%
    slice_max(abs(rank_stat), n = 1, with_ties = FALSE) %>%
    ungroup()
  sort(setNames(mm$rank_stat, mm$ENTREZID), decreasing = TRUE)
}

gsea_arm <- function(m) {
  ranks <- make_ranks(m)
  set.seed(1)
  h <- fgsea(pathways = hallmark_sets, stats = ranks, minSize = 15, maxSize = 500,
             eps = 0, nPermSimple = 10000, BPPARAM = BiocParallel::SerialParam())
  panel_entrez <- unique(bitr(panel_genes, fromType = "SYMBOL", toType = "ENTREZID",
                              OrgDb = org.Mm.eg.db)$ENTREZID)
  stopifnot(length(panel_entrez) == 14, all(panel_entrez %in% names(ranks)))
  set.seed(1)
  p <- fgsea(pathways = list(BTK_NOX2 = panel_entrez), stats = ranks, minSize = 5,
             maxSize = 500, eps = 0, nPermSimple = 10000,
             BPPARAM = BiocParallel::SerialParam())
  list(hallmark = h[order(h$padj), ], panel = p)
}

cat("\nGSEA for primary, all10 ~ group and all10 ~ batch + group ...\n")
g_prim  <- gsea_arm(ref)
g_group <- gsea_arm(meta_group)
g_batch <- gsea_arm(meta_batch)

save_h <- function(g, f) {
  d <- as.data.frame(g$hallmark)
  d$leadingEdge <- sapply(d$leadingEdge, paste, collapse = "/")
  write.csv(d, file.path(out_dir, f), row.names = FALSE)
}
save_h(g_group, "Hallmark_all10_group.csv")
save_h(g_batch, "Hallmark_all10_batch.csv")


# 5. Compare with the primary analysis ---------------------------------------
primary_degs <- ref$SYMBOL[!is.na(ref$meta_padj) & ref$meta_padj < 0.05]

summ <- function(label, m, g) {
  deg <- m[!is.na(m$meta_padj) & m$meta_padj < 0.05, ]
  gene <- function(s, col) { v <- m[[col]][m$SYMBOL == s]; if (length(v)) v[1] else NA_real_ }
  hnes  <- function(nm, col) g$hallmark[[col]][g$hallmark$pathway == nm]
  data.frame(
    analysis = label,
    genes_pooled = nrow(m),
    meta_DEGs = nrow(deg),
    up = sum(deg$meta_logFC > 0),
    down = sum(deg$meta_logFC < 0),
    overlap_with_primary = sum(deg$SYMBOL %in% primary_degs),
    Cybb_log2FC = gene("Cybb", "meta_logFC"), Cybb_padj = gene("Cybb", "meta_padj"),
    Fcer1g_log2FC = gene("Fcer1g", "meta_logFC"), Fcer1g_padj = gene("Fcer1g", "meta_padj"),
    Btk_log2FC = gene("Btk", "meta_logFC"), Btk_padj = gene("Btk", "meta_padj"),
    panel_genes_padj_below_0.05 = sum(m$SYMBOL %in% panel_genes & m$meta_padj < 0.05, na.rm = TRUE),
    IFNG_NES = hnes("HALLMARK_INTERFERON_GAMMA_RESPONSE", "NES"),
    IFNG_padj = hnes("HALLMARK_INTERFERON_GAMMA_RESPONSE", "padj"),
    IFNA_NES = hnes("HALLMARK_INTERFERON_ALPHA_RESPONSE", "NES"),
    IFNA_padj = hnes("HALLMARK_INTERFERON_ALPHA_RESPONSE", "padj"),
    panel_GSEA_NES = g$panel$NES, panel_GSEA_padj = g$panel$padj,
    Cx3cr1_log2FC = gene("Cx3cr1", "meta_logFC"), Cx3cr1_padj = gene("Cx3cr1", "meta_padj"),
    Sall1_log2FC = gene("Sall1", "meta_logFC"),   Sall1_padj = gene("Sall1", "meta_padj"),
    Tjp1_log2FC = gene("Tjp1", "meta_logFC"),     Tjp1_padj = gene("Tjp1", "meta_padj"),
    stringsAsFactors = FALSE)
}

summary_tab <- rbind(
  summ("primary (GSE167459 3 vs 3)",          ref,        g_prim),
  summ("all10 ~ group (5 vs 5)",              meta_group, g_group),
  summ("all10 ~ batch + group (5 vs 5)",      meta_batch, g_batch))
write.csv(summary_tab, file.path(out_dir, "Summary_primary_vs_all10.csv"), row.names = FALSE)

# Table S18: the same numbers, one column per analysis (easier to read)
tab_s18 <- as.data.frame(t(summary_tab[, -1]))
colnames(tab_s18) <- summary_tab$analysis
# Round to 4 significant digits, but keep counts (genes, meta-DEGs, overlap) exact
is_count <- apply(tab_s18, 1, function(x) all(x == round(x)))
tab_s18[!is_count, ] <- signif(tab_s18[!is_count, ], 4)
tab_s18 <- cbind(measure = rownames(tab_s18), tab_s18)
write.csv(tab_s18, file.path(out_dir, "Table_S18_GSE167459_all10_sensitivity.csv"),
          row.names = FALSE)

cat("\n================ SUMMARY ================\n")
print(t(summary_tab), quote = FALSE)

cat("\nDirection and padj of selected genes in the pooled results:\n")
sel_genes <- c("Cybb", "Fcer1g", "Btk", "Cx3cr1", "Sall1", "Tjp1")
cmp <- do.call(rbind, lapply(sel_genes, function(s) {
  one <- function(m) { r <- m[m$SYMBOL == s, ]; if (nrow(r)) c(r$meta_logFC[1], r$meta_padj[1]) else c(NA, NA) }
  data.frame(gene = s,
             primary_log2FC = one(ref)[1], primary_padj = one(ref)[2],
             group_log2FC = one(meta_group)[1], group_padj = one(meta_group)[2],
             batch_log2FC = one(meta_batch)[1], batch_padj = one(meta_batch)[2])
}))
print(cmp, row.names = FALSE, digits = 3)

cat("\nPanel genes with padj < 0.05:\n")
for (nm in c("primary", "all10_group", "all10_batch")) {
  m <- switch(nm, primary = ref, all10_group = meta_group, all10_batch = meta_batch)
  cat(" ", nm, ":", paste(m$SYMBOL[m$SYMBOL %in% panel_genes & m$meta_padj < 0.05], collapse = ", "), "\n")
}

cat("\nSaved to:", out_dir, "\n")
writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "20_GSE167459_all10_sensitivity_sessionInfo.txt"))
