# =============================================================================
# 12_Enrichment.R
# Functional enrichment of the meta-DEGs: overall (up and down separately) and
# for each MCODE network cluster
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 12a. Writes the enrichment tables. Script 12b draws the figures
#           from these tables.
#
# Method
#   - Test       : over-representation analysis (clusterProfiler) for GO
#                  Biological Process, KEGG and Reactome, with Benjamini-
#                  Hochberg (BH) adjustment.
#   - Gene lists : up-regulated and down-regulated meta-DEGs (sign of
#                  meta_logFC), and the genes of each MCODE cluster found with
#                  node score cutoff 0.2.
#   - Background : all analyzable genes of the meta-analysis (k >= 2 datasets),
#                  before any significance filter.
#   - All tested terms are saved (pvalueCutoff = 1). Significant terms are those
#     with p.adjust < 0.05; they are counted in the console and filtered in the
#     figure script.
#   - Cluster GO BP : also simplified with clusterProfiler::simplify()
#                  (similarity cutoff 0.7, keep the term with the lowest
#                  p.adjust). Both the full and the simplified table are saved.
#                  Simplification is not applied to the overall up/down tables.
#   - Cluster order : cluster 1, 2, 3 as listed in the MCODE table. The labels
#                  ISG, MHC1 and Ribo in the file names are the biological
#                  themes of these clusters.
#
# Input
#   results/meta/Meta7dpi_all_genes_k2.csv           (script 05)
#   results/meta/Meta7dpi_metaDEGs_FDR0.05.csv       (script 05)
#   results/network/ppi_nodes.csv                    (script 06)
#   data/mcode/MCODE_clusters_cutoff0.2.txt          (MCODE table exported from
#                                                     Cytoscape, see script 06b)
#   Internet access: KEGG is queried online, so KEGG results can change when
#   KEGG is updated.
#
# Output (results/enrichment/)
#   GO_BP_UP.csv, GO_BP_DOWN.csv, KEGG_UP.csv, KEGG_DOWN.csv,
#   Reactome_UP.csv, Reactome_DOWN.csv
#   Cluster<1-3>_<theme>_genes.csv
#   MCODE_Cluster<1-3>_<theme>_GO_BP_full.csv
#   MCODE_Cluster<1-3>_<theme>_GO_BP_simplified.csv
#   MCODE_Cluster<1-3>_<theme>_KEGG.csv, ..._Reactome.csv
#   results/session_info/12_Enrichment_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/12_Enrichment.R
# Packages: clusterProfiler, org.Mm.eg.db, ReactomePA
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ReactomePA)
})

universe_file <- file.path("results", "meta", "Meta7dpi_all_genes_k2.csv")
deg_file      <- file.path("results", "meta", "Meta7dpi_metaDEGs_FDR0.05.csv")
nodes_file    <- file.path("results", "network", "ppi_nodes.csv")
mcode_file    <- file.path("data", "mcode", "MCODE_clusters_cutoff0.2.txt")
out_dir       <- file.path("results", "enrichment")
info_dir      <- file.path("results", "session_info")

for (d in c(out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

for (f in c(universe_file, deg_file, nodes_file, mcode_file))
  if (!file.exists(f))
    stop("File not found: ", f, ". Run scripts 05 and 06 first, and keep the ",
         "MCODE table in data/mcode/. Current folder: ", getwd())

options(timeout = max(600, getOption("timeout")))


# 1. Enrichment function ------------------------------------------------------
run_enrichment_set <- function(entrez_ids, universe_entrez, label) {
  go <- enrichGO(gene = entrez_ids, universe = universe_entrez,
                 OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP",
                 pAdjustMethod = "BH", pvalueCutoff = 1, readable = TRUE)
  kegg <- enrichKEGG(gene = entrez_ids, universe = universe_entrez,
                     organism = "mmu", pAdjustMethod = "BH", pvalueCutoff = 1)
  reactome <- enrichPathway(gene = entrez_ids, universe = universe_entrez,
                            organism = "mouse", pAdjustMethod = "BH",
                            pvalueCutoff = 1, readable = TRUE)
  cat(label, "- terms tested  GO:", nrow(go@result), " KEGG:", nrow(kegg@result),
      " Reactome:", nrow(reactome@result), "\n")
  cat(label, "- significant (BH < 0.05)  GO:", sum(go@result$p.adjust < 0.05),
      " KEGG:", sum(kegg@result$p.adjust < 0.05),
      " Reactome:", sum(reactome@result$p.adjust < 0.05), "\n")
  list(go = go, kegg = kegg, reactome = reactome)
}


# 2. Background universe ------------------------------------------------------
universe_full <- read.csv(universe_file, stringsAsFactors = FALSE)
stopifnot("SYMBOL" %in% names(universe_full))
n_universe <- length(unique(universe_full$SYMBOL))
universe_ids <- bitr(unique(universe_full$SYMBOL), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Universe genes:", n_universe, "| mapped to Entrez:",
    length(unique(universe_ids$SYMBOL)), "\n")


# 3. Part A: overall enrichment, up versus down -------------------------------
deg <- read.csv(deg_file, stringsAsFactors = FALSE)
stopifnot(all(c("SYMBOL", "meta_logFC") %in% names(deg)))

up   <- deg[deg$meta_logFC > 0, ]
down <- deg[deg$meta_logFC < 0, ]
cat("Up:", nrow(up), "| down:", nrow(down), "\n")
stopifnot(nrow(up) > 0, nrow(down) > 0)

up_ids   <- bitr(up$SYMBOL,   fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
down_ids <- bitr(down$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
cat("Up mapped:", nrow(up_ids), "/", nrow(up),
    "| down mapped:", nrow(down_ids), "/", nrow(down), "\n")

# The meta-DEGs must be part of the background
stopifnot(all(up_ids$ENTREZID   %in% universe_ids$ENTREZID),
          all(down_ids$ENTREZID %in% universe_ids$ENTREZID))

res_up   <- run_enrichment_set(up_ids$ENTREZID,   universe_ids$ENTREZID, "UP")
res_down <- run_enrichment_set(down_ids$ENTREZID, universe_ids$ENTREZID, "DOWN")

write.csv(res_up$go@result,         file.path(out_dir, "GO_BP_UP.csv"),      row.names = FALSE)
write.csv(res_up$kegg@result,       file.path(out_dir, "KEGG_UP.csv"),       row.names = FALSE)
write.csv(res_up$reactome@result,   file.path(out_dir, "Reactome_UP.csv"),   row.names = FALSE)
write.csv(res_down$go@result,       file.path(out_dir, "GO_BP_DOWN.csv"),    row.names = FALSE)
write.csv(res_down$kegg@result,     file.path(out_dir, "KEGG_DOWN.csv"),     row.names = FALSE)
write.csv(res_down$reactome@result, file.path(out_dir, "Reactome_DOWN.csv"), row.names = FALSE)


# 4. Part B: enrichment per MCODE cluster -------------------------------------
# The MCODE table (exported from Cytoscape) has a header line that starts with
# "Cluster", then one line per cluster: cluster, score, nodes, edges, node ids.
mcode_lines <- readLines(mcode_file, warn = FALSE)
header_row  <- grep("^Cluster\t", mcode_lines)
stopifnot(length(header_row) == 1)
cluster_rows <- strsplit(mcode_lines[(header_row + 1):length(mcode_lines)], "\t")
cluster_rows <- cluster_rows[vapply(cluster_rows, length, integer(1)) >= 5]
cluster_ids  <- lapply(cluster_rows, function(x) trimws(strsplit(x[5], ",")[[1]]))
cat("MCODE clusters in the table:", length(cluster_ids),
    "| nodes per cluster:", paste(lengths(cluster_ids), collapse = ", "), "\n")
stopifnot(length(cluster_ids) >= 3)
cluster_ids <- cluster_ids[1:3]

cluster_names <- c("Cluster1_ISG", "Cluster2_MHC1", "Cluster3_Ribo")
names(cluster_ids) <- cluster_names

nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
cluster_results <- list()

for (i in 1:3) {
  cl_name <- cluster_names[i]
  genes <- nodes$gene[match(cluster_ids[[i]], nodes$id)]
  stopifnot(!any(is.na(genes)))

  write.csv(data.frame(gene = genes),
            file.path(out_dir, paste0(cl_name, "_genes.csv")), row.names = FALSE)

  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)
  cat("\n=== ", cl_name, " (genes: ", length(genes), ", mapped: ", nrow(ids), ") ===\n", sep = "")
  cat("Genes:", paste(genes, collapse = ", "), "\n")
  if (!all(ids$ENTREZID %in% universe_ids$ENTREZID))
    warning(cl_name, ": some genes are not in the universe.")

  cluster_results[[cl_name]] <- run_enrichment_set(ids$ENTREZID, universe_ids$ENTREZID, cl_name)

  go_simplified <- simplify(cluster_results[[cl_name]]$go, cutoff = 0.7,
                            by = "p.adjust", select_fun = min)
  cat(cl_name, "- GO BP terms before simplify():",
      nrow(cluster_results[[cl_name]]$go@result),
      "| after:", nrow(go_simplified@result), "\n")

  write.csv(cluster_results[[cl_name]]$go@result,
            file.path(out_dir, paste0("MCODE_", cl_name, "_GO_BP_full.csv")), row.names = FALSE)
  write.csv(go_simplified@result,
            file.path(out_dir, paste0("MCODE_", cl_name, "_GO_BP_simplified.csv")), row.names = FALSE)
  write.csv(cluster_results[[cl_name]]$kegg@result,
            file.path(out_dir, paste0("MCODE_", cl_name, "_KEGG.csv")), row.names = FALSE)
  write.csv(cluster_results[[cl_name]]$reactome@result,
            file.path(out_dir, paste0("MCODE_", cl_name, "_Reactome.csv")), row.names = FALSE)
}

cat("\nAll tables saved to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "12_Enrichment_sessionInfo.txt"))
