# =============================================================================
# 21_WGCNA_GSE283560.R
# Weighted gene co-expression network analysis (WGCNA) of GSE283560 at 7 dpi
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 21. Exploratory network analysis (Methods 2.6; Supplementary
#           Figures S7 and S8; Tables S12 to S15). It uses GSE283560 only,
#           the dataset with the most samples (11: 3 Control, 8 TBI).
#
# Method
#   1. Gene counts of the 11 samples (3 lanes summed, script 04) are filtered
#      to genes with CPM >= 1 in at least 3 samples.
#   2. Counts are transformed with vst (DESeq2, blind = TRUE,
#      minReplicatesForReplace = Inf), mapped to gene symbols (org.Mm.eg.db;
#      duplicated symbols keep the highest mean), and adjusted for log10
#      library size with limma::removeBatchEffect (Control/TBI term protected).
#   3. Signed network, blockwiseModules (minModuleSize 150, mergeCutHeight
#      0.30, deepSplit 1, maxBlockSize 5000). Power = lowest with signed
#      scale-free R2 >= 0.80, otherwise 18 (WGCNA FAQ default, n < 20).
#   4. Module eigengenes are correlated with TBI status (Pearson r), with
#      library size and with a vascular/fibroblast marker score. Exact
#      permutation p values (per module and family-wise) are computed.
#   5. Runs: primary (11 samples, depth adjusted), comparison (11 samples, no
#      depth adjustment), sensitivity (10 samples, without GSM8666331), and
#      leave-one-sample-out (11 networks). Modules are compared across runs by
#      gene overlap. GO Biological Process enrichment is run per module.
#
# Input
#   results/counts/GSE283560_7dpi_counts.csv        (script 04)
#   results/counts/GSE283560_7dpi_sample_table.csv  (script 04)
#
# Output (results/wgcna/)
#   depth_adjusted/, not_adjusted/, sensitivity_without_GSM8666331/
#       module tables, gene/kME tables, GO tables, figures of each run
#   overlap_*.csv, leave_one_out_results.csv, marker_expression_log2CPM.csv
#   console_log.txt
#   final/   the files used in the manuscript, named by their label:
#       Figure_S7a_sample_tree_11samples.png         Figure S7a
#       Figure_S7b_soft_threshold_11samples.png      Figure S7b
#       Figure_S7c_gene_dendrogram_11samples.tiff    Figure S7c
#       Figure_S7d_module_trait_heatmap_11samples.tiff   Figure S7d
#       Figure_S8a_module_trait_heatmap_10samples.tiff   Figure S8a
#       Figure_S8b_GO_dotplot_greenyellow_10samples.tiff Figure S8b
#       Table_S12_module_stats_11samples.csv         Table S12, sheet 1
#       Table_S12_module_stats_10samples.csv         Table S12, sheet 2
#       Table_S12_marker_expression_log2CPM.csv      Table S12, sheet 3
#       Table_S13_gene_module_kME_11samples.csv      Table S13, sheet 1
#       Table_S13_BTK_NOX2_panel_modules_11samples.csv   Table S13, sheet 2
#       Table_S14_GO_terms_11samples.csv             Table S14, sheet 1
#       Table_S14_GO_terms_10samples.csv             Table S14, sheet 2
#       Table_S15_leave_one_out_results.csv          Table S15
#   results/session_info/21_WGCNA_GSE283560_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/21_WGCNA_GSE283560.R
#           (takes several hours; the leave-one-out step takes most)
# Packages: WGCNA, DESeq2, limma, org.Mm.eg.db, AnnotationDbi, clusterProfiler,
#           dplyr, tibble, tidyr, ggplot2
# Note    : the figures use the font Arial (installed on most Windows and Mac
#           systems; on Linux ggplot2 falls back to another font).
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA); library(DESeq2); library(org.Mm.eg.db); library(AnnotationDbi)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2); library(clusterProfiler)
})
if (!requireNamespace("limma", quietly = TRUE))
  stop("Install limma first:  BiocManager::install('limma')")
select <- dplyr::select
cor    <- WGCNA::cor
options(stringsAsFactors = FALSE)
allowWGCNAThreads()

## ---------------------------------------------------------------------------
## SETTINGS (edit here only)
## ---------------------------------------------------------------------------
out_dir       <- file.path("results", "wgcna")
info_dir      <- file.path("results", "session_info")
count_file    <- file.path("results", "counts", "GSE283560_7dpi_counts.csv")
sample_file   <- file.path("results", "counts", "GSE283560_7dpi_sample_table.csv")
min_cpm       <- 1        # CPM = counts per million reads
min_samples   <- 3        # size of the smallest group (3 controls)
depth_flag_r  <- 0.60     # |r| >= 0.60 with n = 11 is nominally p < 0.05
run_loo       <- TRUE     # STEP 7: leave-one-sample-out, takes most of the run time
chemo_pat     <- "^Or[0-9]|^Olfr|^Vmn|^Taar|^Tas2r"   # smell, vomeronasal, trace amine, bitter taste receptors

# Order of the 11 samples (fixed, so that results do not depend on file order)
expected_cols <- c("GSM8666315", "GSM8666321", "GSM8666322", "GSM8666331", "GSM8666337",
                   "GSM8666338", "GSM8666346", "GSM8666352", "GSM8666353", "GSM8666365",
                   "GSM8666366")

while (sink.number() > 0) sink()          # close any log left open by an earlier failed run
if (!all(file.exists(c(count_file, sample_file))))
  stop("Input files not found. Run scripts 04 first (it writes results/counts/). Current folder: ", getwd())
for (d in c(out_dir, info_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
sink(file.path(out_dir, "console_log.txt"), split = TRUE)

## ---------------------------------------------------------------------------
## STEP 1. LOAD AND CHECK INPUT
## ---------------------------------------------------------------------------
raw_counts <- read.csv(count_file, row.names = 1, check.names = FALSE)
pheno <- read.csv(sample_file)
stopifnot(setequal(colnames(raw_counts), expected_cols))
raw_counts <- raw_counts[, expected_cols]

cybb_check <- raw_counts[rownames(raw_counts) == "ENSMUSG00000015340", "GSM8666315"]
cat("Fingerprint Cybb/GSM8666315 =", cybb_check, "(must be about 8158)\n")
stopifnot(length(cybb_check) == 1, abs(cybb_check - 8158.299) < 1, ncol(raw_counts) == 11)

id_col <- intersect(c("sample", "Sample", "sample_id", "SampleID", "GSM", "gsm", "X"), colnames(pheno))[1]
if (is.na(id_col)) id_col <- colnames(pheno)[1]
pheno <- pheno[match(colnames(raw_counts), pheno[[id_col]]), ]
stopifnot(!anyNA(pheno[[id_col]]), identical(colnames(raw_counts), as.character(pheno[[id_col]])))
rownames(pheno) <- as.character(pheno[[id_col]])
sample_ids <- rownames(pheno)
pheno$condition <- factor(pheno$condition, levels = c("Control", "TBI"))
stopifnot(sum(pheno$condition == "Control") == 3, sum(pheno$condition == "TBI") == 8)

libsize   <- colSums(raw_counts)
log_depth <- log10(libsize)
trait     <- as.numeric(pheno$condition) - 1
cat("\nLibrary size (million reads) by group:\n")
print(round(tapply(libsize / 1e6, pheno$condition, mean), 1))
print(round(sort(libsize / 1e6), 1))

## ---------------------------------------------------------------------------
## STEP 2. DEPTH-AWARE FILTER  (and evidence that it removes noise genes)
## ---------------------------------------------------------------------------
cpm_mat  <- sweep(as.matrix(raw_counts), 2, libsize / 1e6, "/")
keep_old <- rowSums(raw_counts >= 10) >= 2
keep_new <- rowSums(cpm_mat >= min_cpm) >= min_samples

symbol_all <- mapIds(org.Mm.eg.db, keys = rownames(raw_counts), keytype = "ENSEMBL",
                     column = "SYMBOL", multiVals = "first")
chemo_all <- grepl(chemo_pat, symbol_all)

cat("\nSTEP 2: filter comparison\n")
cat("  Old filter (>=10 counts in >=2 samples): ", sum(keep_old), "genes;",
    sum(keep_old & chemo_all), "chemoreceptor genes\n")
cat("  New filter (CPM >=", min_cpm, "in >=", min_samples, "samples):", sum(keep_new), "genes;",
    sum(keep_new & chemo_all), "chemoreceptor genes\n")

## ---------------------------------------------------------------------------
## STEP 3. NORMALIZE, VST, SYMBOLS, DEPTH ADJUSTMENT
## ---------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = round(raw_counts[keep_new, ]), colData = pheno, design = ~ condition)
dds <- DESeq(dds, minReplicatesForReplace = Inf)      # no automatic outlier replacement (as in Methods 2.2.4)
vsd <- vst(dds, blind = TRUE)
expr_mat <- assay(vsd)

sym <- symbol_all[rownames(expr_mat)]
ok  <- !is.na(sym) & sym != ""
cat("\nSTEP 3: genes without a symbol (dropped):", sum(!ok), "\n")
expr_df <- as.data.frame(expr_mat[ok, , drop = FALSE])
expr_df$SYMBOL    <- unname(sym[ok])
expr_df$mean_expr <- rowMeans(expr_mat[ok, , drop = FALSE])
expr_df <- expr_df[order(-expr_df$mean_expr), ]
n_dup   <- sum(duplicated(expr_df$SYMBOL))
expr_df <- expr_df[!duplicated(expr_df$SYMBOL), ]
cat("Duplicate symbols removed (kept highest mean expression):", n_dup, "\n")
expr_sym <- as.matrix(expr_df[, sample_ids]); rownames(expr_sym) <- expr_df$SYMBOL
cat("Genes entering the network:", nrow(expr_sym), "\n")

# Remove the depth effect; the Control/TBI term is kept in the design so it is protected.
design_mat <- model.matrix(~ condition, data = pheno)
expr_adj <- limma::removeBatchEffect(expr_sym,
                                     covariates = matrix(as.numeric(scale(log_depth)), ncol = 1),
                                     design = design_mat)

## ---------------------------------------------------------------------------
## STEP 4. MARKER CHECK (microglia vs vascular/fibroblast/immune), from log2 CPM
## ---------------------------------------------------------------------------
markers <- list(microglia = c("P2ry12", "Cx3cr1", "Tmem119", "Hexb"),
                vascular_fibroblast = c("Col1a1", "Col1a2", "Pdgfrb", "Acta2", "Dcn"),
                other_immune = c("Ptprc", "Cd3e"))
m_syms <- unlist(markers)
m_ens  <- mapIds(org.Mm.eg.db, keys = m_syms, keytype = "SYMBOL", column = "ENSEMBL", multiVals = "first")
m_ens  <- m_ens[!is.na(m_ens) & m_ens %in% rownames(cpm_mat)]
marker_log <- log2(cpm_mat[m_ens, , drop = FALSE] + 1)
rownames(marker_log) <- names(m_ens)
marker_set <- sapply(rownames(marker_log), function(g) names(markers)[sapply(markers, function(v) g %in% v)][1])
marker_tab <- data.frame(gene = rownames(marker_log), set = unname(marker_set),
                         round(marker_log, 2), check.names = FALSE)
write.csv(marker_tab, file.path(out_dir, "marker_expression_log2CPM.csv"), row.names = FALSE)
cat("\nSTEP 4: marker expression (log2 CPM + 1) per sample\n"); print(marker_tab, row.names = FALSE)

vasc_present <- intersect(markers$vascular_fibroblast, rownames(marker_log))
contam_score <- if (length(vasc_present) > 0)
  rowMeans(scale(t(marker_log[vasc_present, , drop = FALSE])), na.rm = TRUE) else rep(NA_real_, 11)
names(contam_score) <- sample_ids

## ---------------------------------------------------------------------------
## FUNCTION: one full WGCNA run on a genes x samples matrix
## ---------------------------------------------------------------------------
btk_genes  <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1")
nox2_genes <- c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
panel_genes <- c(btk_genes, nox2_genes)

run_wgcna <- function(expr_gs, tag, do_go = FALSE, keep = sample_ids) {
  # 'keep' lets a run use a subset of samples (sensitivity analysis); local copies shadow the global vectors
  ii <- match(keep, sample_ids)
  trait <- trait[ii]; log_depth <- log_depth[ii]; contam_score <- contam_score[ii]
  sample_ids <- keep
  expr_gs <- expr_gs[, keep, drop = FALSE]
  d <- file.path(out_dir, tag); dir.create(d, showWarnings = FALSE)
  cat("\n\n############ RUN:", tag, "(", length(keep), "samples ) ############\n")
  datExpr <- t(expr_gs)
  stopifnot(identical(rownames(datExpr), sample_ids))
  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  stopifnot(all(gsg$goodSamples))
  datExpr <- datExpr[, gsg$goodGenes]
  cat("Genes in network:", ncol(datExpr), "| chemoreceptor genes:",
      sum(grepl(chemo_pat, colnames(datExpr))), "\n")

  # sample tree with group and depth color bars
  tree <- hclust(dist(datExpr), method = "average")
  bars <- numbers2colors(cbind(TBI = trait, log10_depth = log_depth), signed = FALSE)
  png(file.path(d, "sample_tree.png"), width = 1000, height = 600)
  plotDendroAndColors(tree, bars, groupLabels = c("TBI (0/1)", "log10 depth"),
                      dendroLabels = sample_ids, main = paste("Sample clustering:", tag))
  dev.off()

  # soft threshold
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed", verbose = 0)
  fit <- sft$fitIndices; fit$signed_R2 <- -sign(fit$slope) * fit$SFT.R.sq
  print(round(fit[, c("Power", "signed_R2", "mean.k.")], 3))
  ok_p <- fit$Power[fit$signed_R2 >= 0.80]
  if (length(ok_p) > 0) {
    beta <- min(ok_p); beta_rule <- "lowest power with signed R2 >= 0.80"
  } else {
    # WGCNA FAQ: if R2 never reaches 0.80, use the default for signed networks by sample size (< 20 samples: 18)
    beta <- 18; beta_rule <- "FALLBACK: R2 never reached 0.80; WGCNA FAQ default for signed networks, n < 20"
    cat("WARNING: no power reached signed R2 >= 0.80 in run", tag, "-> using power 18 (WGCNA FAQ default)\n")
  }
  beta_R2 <- fit$signed_R2[fit$Power == beta]; mean_k <- fit$mean.k.[fit$Power == beta]
  cat("Chosen beta =", beta, "| signed R2 =", round(beta_R2, 3), "| mean connectivity =", round(mean_k), "\n")
  png(file.path(d, "soft_threshold.png"), width = 900, height = 500)
  par(mfrow = c(1, 2))
  plot(fit$Power, fit$signed_R2, type = "n", xlab = "Soft threshold (power)",
       ylab = "Scale-free topology fit, signed R^2", main = "Scale independence")
  text(fit$Power, fit$signed_R2, labels = fit$Power); abline(h = 0.80, col = "red")
  points(beta, beta_R2, pch = 1, cex = 3, col = "red")
  plot(fit$Power, fit$mean.k., type = "n", xlab = "Soft threshold (power)",
       ylab = "Mean connectivity", main = "Mean connectivity")
  text(fit$Power, fit$mean.k., labels = fit$Power)
  dev.off()

  # modules (settings unchanged from the manuscript)
  net <- blockwiseModules(datExpr, power = beta, networkType = "signed", TOMType = "signed",
                          minModuleSize = 150, mergeCutHeight = 0.30, deepSplit = 1,
                          maxBlockSize = 5000, numericLabels = TRUE, saveTOMs = FALSE, verbose = 2)
  mcol <- labels2colors(net$colors)
  cat("Blocks:", length(net$dendrograms), "| block sizes:", paste(lengths(net$blockGenes), collapse = ", "), "\n")
  cat("Real modules:", length(setdiff(unique(mcol), "grey")), "| grey genes:", sum(mcol == "grey"), "\n")
  print(sort(table(mcol), decreasing = TRUE))

  # module-trait, module-depth, module-marker correlations (p values are NOT adjusted)
  MEs <- orderMEs(moduleEigengenes(datExpr, colors = mcol)$eigengenes)
  nm  <- gsub("^ME", "", colnames(MEs))
  cc  <- function(v) {
    r <- cor(MEs, v, use = "p")[, 1]; p <- as.numeric(corPvalueStudent(r, nSamples = nrow(datExpr)))
    list(r = round(r, 3), p = round(p, 4))
  }
  a <- cc(trait); b <- cc(log_depth)

  # Exact permutation test. Modules were built without the TBI label, so the label can be shuffled.
  # All ways of choosing which samples are 'Control' are tried (n = 10: 120 ways; n = 11: 165 ways).
  # p_perm = chance of an |r| this large for THIS module; p_perm_familywise = chance that ANY module gets one this large.
  ME_use   <- as.matrix(MEs[, paste0("ME", setdiff(nm, "grey")), drop = FALSE])
  combs    <- combn(nrow(datExpr), sum(trait == 0))
  perm_abs <- sapply(seq_len(ncol(combs)), function(k) {
    tr <- rep(1, nrow(datExpr)); tr[combs[, k]] <- 0; abs(cor(ME_use, tr)[, 1]) })
  perm_max <- apply(perm_abs, 2, max)
  r_obs    <- cor(ME_use, trait)[, 1]
  p_perm   <- sapply(seq_along(r_obs), function(i) mean(perm_abs[i, ] >= abs(r_obs[i]) - 1e-12))
  p_fwer   <- sapply(seq_along(r_obs), function(i) mean(perm_max >= abs(r_obs[i]) - 1e-12))
  names(p_perm) <- names(p_fwer) <- setdiff(nm, "grey")
  cat("Permutation test:", ncol(combs), "label arrangements tried\n")
  cm <- if (all(is.finite(contam_score))) cc(contam_score) else list(r = NA, p = NA)
  mod_tab <- tibble(module = nm, n_genes = as.numeric(table(mcol)[nm]),
                    cor_TBI = a$r, p_TBI = a$p, p_perm = round(p_perm[nm], 4),
                    p_perm_familywise = round(p_fwer[nm], 4), cor_depth = b$r, p_depth = b$p,
                    cor_vascular_fibroblast = cm$r, p_vascular_fibroblast = cm$p) %>%
    filter(module != "grey") %>% arrange(desc(cor_TBI))
  write.csv(mod_tab, file.path(d, "module_trait_depth_marker_correlations.csv"), row.names = FALSE)
  cat("\nModule table (p values unadjusted):\n"); print(as.data.frame(mod_tab))
  cat("Modules with p_TBI < 0.05:", sum(mod_tab$p_TBI < 0.05), "\n")
  cat("Modules with family-wise permutation p < 0.05:", sum(mod_tab$p_perm_familywise < 0.05),
      if (any(mod_tab$p_perm_familywise < 0.05)) paste0("(", paste(mod_tab$module[mod_tab$p_perm_familywise < 0.05], collapse = ", "), ")") else "", "\n")
  dep_flag <- mod_tab$module[abs(mod_tab$cor_depth) >= depth_flag_r]
  cat("Modules with |r| >=", depth_flag_r, "with depth:", length(dep_flag),
      if (length(dep_flag)) paste0("(", paste(dep_flag, collapse = ", "), ")") else "", "\n")

  # genes: module + kME (supplementary table), panel summary
  kME <- signedKME(datExpr, MEs)
  idx <- match(paste0("kME", mcol), colnames(kME))
  gene_tab <- tibble(gene = colnames(datExpr), module = mcol,
                     kME_own_module = kME[cbind(seq_len(nrow(kME)), idx)])
  write.csv(gene_tab, file.path(d, "all_genes_module_kME.csv"), row.names = FALSE)
  panel_tab <- gene_tab %>% filter(gene %in% panel_genes) %>%
    mutate(pathway = ifelse(gene %in% btk_genes, "BTK", "NOX2")) %>% arrange(desc(kME_own_module))
  write.csv(panel_tab, file.path(d, "BTK_NOX2_panel_module_kME.csv"), row.names = FALSE)
  cat("\nBTK/NOX2 panel: ", nrow(panel_tab), "of 14 genes in the network, spread over",
      length(unique(panel_tab$module)), "modules\n"); print(as.data.frame(panel_tab))

  # figures
  long <- bind_rows(
    mod_tab %>% transmute(module, variable = "TBI status", r = cor_TBI, p = p_TBI),
    mod_tab %>% transmute(module, variable = "Sequencing depth", r = cor_depth, p = p_depth)) %>%
    mutate(module = factor(module, levels = rev(mod_tab$module)),
           label = paste0(round(r, 2), " (p=", round(p, 3), ")"))
  p_heat <- ggplot(long, aes(variable, module, fill = r)) + geom_tile(color = "white") +
    geom_text(aes(label = label), size = 3.2, family = "Arial") +
    scale_fill_gradient2(low = "#2471A3", mid = "white", high = "#C0392B", midpoint = 0,
                         limits = c(-1, 1), name = "r") +
    labs(x = NULL, y = NULL) + theme_minimal(base_size = 9, base_family = "Arial") +
    theme(axis.text = element_text(size = 9))
  ggsave(file.path(d, "module_trait_heatmap.tiff"), p_heat, width = 150, height = 150,
         units = "mm", dpi = 300, compression = "lzw")

  p_bar <- ggplot(mod_tab, aes(reorder(module, cor_TBI), cor_TBI, fill = cor_TBI > 0)) +
    geom_col() + coord_flip() +
    scale_fill_manual(values = c(`TRUE` = "#C0392B", `FALSE` = "#2471A3"), guide = "none") +
    labs(x = NULL, y = "Correlation with TBI status (r)") + theme_minimal(base_size = 9, base_family = "Arial")
  ggsave(file.path(d, "module_trait_barplot_supplement.tiff"), p_bar, width = 100, height = 140,
         units = "mm", dpi = 300, compression = "lzw")

  nb <- length(net$dendrograms)
  tiff(file.path(d, "gene_dendrogram_all_blocks.tiff"), width = 180, height = 60 * nb,
       units = "mm", res = 300, compression = "lzw", family = "sans")
  layout(matrix(seq_len(2 * nb), ncol = 1), heights = rep(c(3, 0.8), nb))   # own layout: one dendrogram + one color bar per block
  for (i in seq_len(nb))
    plotDendroAndColors(net$dendrograms[[i]], mcol[net$blockGenes[[i]]], "Module colors",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                        setLayout = FALSE, marAll = c(1, 5, 3, 1),
                        main = paste0("Block ", i, " (n=", length(net$blockGenes[[i]]), " genes)"))
  dev.off()
  for (i in seq_len(nb)) {
    png(file.path(d, paste0("gene_dendrogram_block", i, ".png")), width = 1000, height = 600)
    plotDendroAndColors(net$dendrograms[[i]], mcol[net$blockGenes[[i]]], "Module colors",
                        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                        main = paste0("Block ", i, " (n=", length(net$blockGenes[[i]]), " genes)"))
    dev.off()
  }

  p_kme <- ggplot(panel_tab, aes(reorder(gene, kME_own_module), kME_own_module, fill = module)) +
    geom_col() + scale_fill_identity() + coord_flip() +
    labs(x = NULL, y = "kME (own module)") + theme_minimal(base_size = 9, base_family = "Arial")
  ggsave(file.path(d, "BTK_NOX2_kME_supplement.tiff"), p_kme, width = 85, height = 100,
         units = "mm", dpi = 300, compression = "lzw")

  # GO for ALL modules (universe = all network genes; p cutoff 1 keeps null results visible)
  go_summary <- NULL
  if (do_go) {
    uni <- mapIds(org.Mm.eg.db, keys = colnames(datExpr), keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
    uni <- uni[!is.na(uni)]
    sig_list <- list(); ego_list <- list(); rows <- list()
    for (m in setdiff(unique(mcol), "grey")) {
      g   <- colnames(datExpr)[mcol == m]
      ent <- mapIds(org.Mm.eg.db, keys = g, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")
      ent <- ent[!is.na(ent)]
      ego <- enrichGO(gene = ent, universe = uni, OrgDb = org.Mm.eg.db, ont = "BP",
                      pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
      n_sig <- if (is.null(ego)) 0 else sum(ego@result$p.adjust < 0.05)
      rows[[m]] <- tibble(module = m, n_genes = length(g), n_mapped = length(ent), n_sig_terms = n_sig)
      if (!is.null(ego) && n_sig > 0) {
        ego_list[[m]] <- ego
        sig_list[[m]] <- ego@result %>% filter(p.adjust < 0.05) %>% mutate(module = m)
        cat("\n", m, "(", length(g), "genes,", n_sig, "significant terms) top terms:\n")
        print(head(as.data.frame(sig_list[[m]] %>% arrange(p.adjust) %>%
                                   dplyr::select(Description, p.adjust, Count)), 5))
      } else cat("\n", m, "(", length(g), "genes): 0 significant GO BP terms\n")
    }
    go_summary <- bind_rows(rows)
    write.csv(go_summary, file.path(d, "GO_summary_all_modules.csv"), row.names = FALSE)
    if (length(sig_list)) write.csv(bind_rows(sig_list), file.path(d, "GO_significant_terms_all_modules.csv"), row.names = FALSE)
    for (m in names(ego_list)) {
      p <- dotplot(ego_list[[m]], showCategory = 10) + theme_minimal(base_size = 9, base_family = "Arial")
      ggsave(file.path(d, paste0("GO_dotplot_", m, ".tiff")), p, width = 180, height = 110,
             units = "mm", dpi = 300, compression = "lzw")
    }
  }
  save(net, mcol, MEs, mod_tab, gene_tab, panel_tab, beta, datExpr, file = file.path(d, "WGCNA_main_objects.RData"))
  list(datExpr = datExpr, mcol = mcol, mod_tab = mod_tab, gene_tab = gene_tab, beta = beta,
       beta_R2 = beta_R2, beta_rule = beta_rule, mean_k = mean_k, nblocks = length(net$dendrograms), dep_flag = dep_flag)
}

## ---------------------------------------------------------------------------
## STEP 5. RUN BOTH VERSIONS
## ---------------------------------------------------------------------------
prim <- run_wgcna(expr_adj, "depth_adjusted", do_go = TRUE)   # primary
comp <- run_wgcna(expr_sym, "not_adjusted",   do_go = FALSE)  # comparison: same filter, no depth removal

# Sensitivity run: drop GSM8666331 (lowest depth, raised non-microglial markers, outlier in the sample tree).
# Depth adjustment is repeated on the remaining 10 samples. Filter and VST were computed on all 11 samples.
drop_sample <- "GSM8666331"
keep10 <- setdiff(sample_ids, drop_sample)
stopifnot(length(keep10) == 10, sum(pheno[keep10, "condition"] == "Control") == 3)
expr_adj10 <- limma::removeBatchEffect(expr_sym[, keep10],
                covariates = matrix(as.numeric(scale(log_depth[match(keep10, sample_ids)])), ncol = 1),
                design = model.matrix(~ condition, data = pheno[keep10, ]))
sens <- run_wgcna(expr_adj10, "sensitivity_without_GSM8666331", do_go = TRUE, keep = keep10)

## ---------------------------------------------------------------------------
## STEP 6. DO THE MODULES SURVIVE? (module overlap between versions)
## ---------------------------------------------------------------------------
cmp_labels <- function(g1, l1, g2, l2, file, name) {
  common <- intersect(g1, g2)
  ov <- overlapTable(l1[match(common, g1)], l2[match(common, g2)])
  write.csv(ov$countTable, file.path(out_dir, paste0("overlap_", file, "_counts.csv")))
  write.csv(signif(ov$pTable, 3), file.path(out_dir, paste0("overlap_", file, "_p.csv")))
  cat("\nOverlap table saved:", name, "(", length(common), "common genes )\n")
}
cmp_labels(colnames(prim$datExpr), prim$mcol, colnames(comp$datExpr), comp$mcol,
           "adjusted_vs_not_adjusted", "depth adjusted vs not adjusted")
cmp_labels(colnames(prim$datExpr), prim$mcol, colnames(sens$datExpr), sens$mcol,
           "adjusted_vs_without_GSM8666331", "depth adjusted (11 samples) vs without GSM8666331 (10 samples)")

## ---------------------------------------------------------------------------
## STEP 7. LEAVE-ONE-SAMPLE-OUT.
## Question: is the TBI-linked module of the 10-sample run real, or does it only appear when GSM8666331 is removed?
## Each of the 11 samples is dropped in turn; depth adjustment and WGCNA are repeated (filter and VST stay as in STEP 3).
## Read: best_p_familywise = permutation p that ANY module reaches the best |r|; frac_ref_in_best = share of the
## reference TBI-linked genes (sensitivity run) found in the best module of that run.
## ---------------------------------------------------------------------------
loo_tab <- NULL
if (run_loo) {
  cat("\n\n############ STEP 7: LEAVE-ONE-SAMPLE-OUT ############\n")
  ref_mods  <- sens$mod_tab$module[sens$mod_tab$p_TBI < 0.05]
  ref_genes <- colnames(sens$datExpr)[sens$mcol %in% ref_mods]
  nox2_chk  <- c("Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
  cat("Reference genes = modules", paste(ref_mods, collapse = ", "), "of the 10-sample run:", length(ref_genes), "genes\n")

  loo_one <- function(dr) {
    keep <- setdiff(sample_ids, dr); ii <- match(keep, sample_ids)
    tr <- trait[ii]; ld <- log_depth[ii]
    ea <- limma::removeBatchEffect(expr_sym[, keep, drop = FALSE],
                                   covariates = matrix(as.numeric(scale(ld)), ncol = 1),
                                   design = model.matrix(~ condition, data = pheno[keep, ]))
    dat <- t(ea); dat <- dat[, goodSamplesGenes(dat, verbose = 0)$goodGenes]
    sft <- pickSoftThreshold(dat, powerVector = 1:20, networkType = "signed", verbose = 0)
    fit <- sft$fitIndices; r2 <- -sign(fit$slope) * fit$SFT.R.sq
    okp <- fit$Power[r2 >= 0.80]; beta <- if (length(okp)) min(okp) else 18
    net  <- blockwiseModules(dat, power = beta, networkType = "signed", TOMType = "signed",
                             minModuleSize = 150, mergeCutHeight = 0.30, deepSplit = 1,
                             maxBlockSize = 5000, numericLabels = TRUE, saveTOMs = FALSE, verbose = 0)
    mc   <- labels2colors(net$colors)
    MEs_ <- moduleEigengenes(dat, colors = mc)$eigengenes
    nm_  <- gsub("^ME", "", colnames(MEs_)); use <- nm_ != "grey"
    MEm  <- as.matrix(MEs_[, use, drop = FALSE]); nm_ <- nm_[use]
    r    <- cor(MEm, tr)[, 1]
    cmb  <- combn(length(tr), sum(tr == 0))
    pmax_ <- apply(sapply(seq_len(ncol(cmb)), function(k) {
      t2 <- rep(1, length(tr)); t2[cmb[, k]] <- 0; abs(cor(MEm, t2)[, 1]) }), 2, max)
    b    <- which.max(abs(r)); gens <- colnames(dat)[mc == nm_[b]]
    data.frame(dropped = dr, group = as.character(pheno[dr, "condition"]), n_control = sum(tr == 0),
               beta = beta, R2_at_beta = round(r2[fit$Power == beta], 2), modules = length(nm_),
               grey = sum(mc == "grey"), best_module = nm_[b], best_size = length(gens),
               best_r = round(r[b], 3), best_p_familywise = round(mean(pmax_ >= abs(r[b]) - 1e-12), 3),
               frac_ref_in_best = round(mean(ref_genes %in% gens), 2),
               frac_best_in_ref = round(mean(gens %in% ref_genes), 2),
               NOX2_in_best = paste(intersect(nox2_chk, gens), collapse = ","), stringsAsFactors = FALSE)
  }
  loo_list <- list()
  for (dr in sample_ids) {
    cat("  dropping", dr, "...\n")
    loo_list[[dr]] <- tryCatch(loo_one(dr), error = function(e) {
      cat("   FAILED:", conditionMessage(e), "\n"); NULL })
  }
  loo_tab <- do.call(rbind, loo_list)
  write.csv(loo_tab, file.path(out_dir, "leave_one_out_results.csv"), row.names = FALSE)
  cat("\nLeave-one-out table saved.\n"); print(loo_tab, row.names = FALSE)
}

## ---------------------------------------------------------------------------
## SUMMARY
## ---------------------------------------------------------------------------
cat("\n================ SUMMARY ================\n")
cat("Filter: CPM >=", min_cpm, "in >=", min_samples, "samples | genes entering network:", nrow(expr_sym), "\n")
for (r in list(list("depth_adjusted (primary)", prim), list("not_adjusted (comparison)", comp),
                list("sensitivity: without GSM8666331 (10 samples)", sens))) {
  x <- r[[2]]
  cat("\n", r[[1]], "\n  genes:", ncol(x$datExpr), "| chemoreceptor genes:", sum(grepl(chemo_pat, colnames(x$datExpr))),
      "\n  beta:", x$beta, "| signed R2:", round(x$beta_R2, 3), "| mean connectivity:", round(x$mean_k), "\n  beta rule:", x$beta_rule,
      "\n  blocks:", x$nblocks, "| real modules:", nrow(x$mod_tab), "| grey genes:", sum(x$mcol == "grey"),
      "\n  modules with p_TBI < 0.05:", sum(x$mod_tab$p_TBI < 0.05),
      "\n  modules with family-wise permutation p < 0.05:", sum(x$mod_tab$p_perm_familywise < 0.05),
      "\n  modules with |r| >=", depth_flag_r, "with depth:", length(x$dep_flag), "\n")
}
if (!is.null(loo_tab)) { cat("\nLeave-one-sample-out (STEP 7):\n"); print(loo_tab, row.names = FALSE) }
cat("====================================================\n")

## ---------------------------------------------------------------------------
## STEP 8. COPY THE FILES USED IN THE MANUSCRIPT TO results/wgcna/final/
## ---------------------------------------------------------------------------
fin <- file.path(out_dir, "final"); dir.create(fin, showWarnings = FALSE)
p1 <- file.path(out_dir, "depth_adjusted")
p2 <- file.path(out_dir, "sensitivity_without_GSM8666331")
copy_map <- c(
  "Figure_S7a_sample_tree_11samples.png"               = file.path(p1, "sample_tree.png"),
  "Figure_S7b_soft_threshold_11samples.png"            = file.path(p1, "soft_threshold.png"),
  "Figure_S7c_gene_dendrogram_11samples.tiff"          = file.path(p1, "gene_dendrogram_all_blocks.tiff"),
  "Figure_S7d_module_trait_heatmap_11samples.tiff"     = file.path(p1, "module_trait_heatmap.tiff"),
  "Figure_S8a_module_trait_heatmap_10samples.tiff"     = file.path(p2, "module_trait_heatmap.tiff"),
  "Figure_S8b_GO_dotplot_greenyellow_10samples.tiff"   = file.path(p2, "GO_dotplot_greenyellow.tiff"),
  "Table_S13_gene_module_kME_11samples.csv"             = file.path(p1, "all_genes_module_kME.csv"),
  "Table_S13_BTK_NOX2_panel_modules_11samples.csv"      = file.path(p1, "BTK_NOX2_panel_module_kME.csv"),
  "Table_S12_module_stats_11samples.csv"               = file.path(p1, "module_trait_depth_marker_correlations.csv"),
  "Table_S12_module_stats_10samples.csv"               = file.path(p2, "module_trait_depth_marker_correlations.csv"),
  "Table_S12_marker_expression_log2CPM.csv"            = file.path(out_dir, "marker_expression_log2CPM.csv"),
  "Table_S15_leave_one_out_results.csv"                = file.path(out_dir, "leave_one_out_results.csv"),
  "Table_S14_GO_terms_11samples.csv"                   = file.path(p1, "GO_significant_terms_all_modules.csv"),
  "Table_S14_GO_terms_10samples.csv"                   = file.path(p2, "GO_significant_terms_all_modules.csv"))
for (nm_ in names(copy_map)) {
  if (file.exists(copy_map[[nm_]])) file.copy(copy_map[[nm_]], file.path(fin, nm_), overwrite = TRUE)
  else cat("NOT FOUND (check the module colour or run_loo):", copy_map[[nm_]], "\n")
}
cat("\nFiles in results/wgcna/final/:", length(list.files(fin)), "of", length(copy_map), "\n")

writeLines(capture.output(sessionInfo()), file.path(info_dir, "21_WGCNA_GSE283560_sessionInfo.txt"))
sink()
