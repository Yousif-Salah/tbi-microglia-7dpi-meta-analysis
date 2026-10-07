# =============================================================================
# 06_STRING_network.R
# Protein-protein interaction network of the meta-DEGs (STRING)
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 6a. Builds the node and edge tables of the STRING network from
#           the meta-DEG list (script 05). Script 06b draws and clusters the
#           network in Cytoscape.
#
# Method
#   - Genes     : symbols of all meta-DEGs (meta_padj < 0.05).
#   - STRING    : version 12.0, mouse (taxonomy 10090), minimum combined
#                 confidence score 700 (0.70).
#   - Nodes     : genes mapped to a STRING protein id. One node per STRING id
#                 (two symbols can map to one id; the first symbol is kept).
#   - Edges     : STRING returns each pair twice (A-B and B-A). One row is kept
#                 per pair, and self loops are removed.
#   - Counts    : every count is printed and saved to a text file.
#
# Input
#   results/meta/Meta7dpi_metaDEGs_FDR0.05.csv     (script 05)
#   Internet access: STRING files (about 72 MB) are downloaded once to
#   data/STRING_cache/ and reused on later runs.
#
# Output (results/network/)
#   ppi_nodes.csv            columns id (STRING id), gene (symbol)
#   ppi_edges.csv            columns source, target, combined_score
#   ppi_counts_summary.txt   gene, node and edge counts, STRING version
#   results/session_info/06_STRING_network_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/06_STRING_network.R
# Packages: STRINGdb
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages(library(STRINGdb))

deg_file  <- file.path("results", "meta", "Meta7dpi_metaDEGs_FDR0.05.csv")
cache_dir <- file.path("data", "STRING_cache")
out_dir   <- file.path("results", "network")
info_dir  <- file.path("results", "session_info")

options(timeout = 1800)   # the mouse interaction file is large; default is 60 s
for (d in c(cache_dir, out_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(deg_file))
  stop("Meta-DEG file not found: ", deg_file,
       ". Run script 05 first. Current folder: ", getwd())


# 1. Input genes --------------------------------------------------------------
deg <- read.csv(deg_file, stringsAsFactors = FALSE)
stopifnot("SYMBOL" %in% names(deg))

genes <- unique(deg$SYMBOL)
genes <- genes[!is.na(genes) & genes != ""]
cat("Genes loaded:", length(genes), "\n")
stopifnot(length(genes) > 0)

# Information only: which genes of the BTK/NOX2 panel are meta-DEGs
panel <- c("Btk", "Lyn", "Syk", "Plcg2", "Blnk", "Hck", "Fcer1g", "Nfkb1",
           "Cybb", "Cyba", "Ncf1", "Ncf2", "Ncf4", "Rac2")
cat("Panel genes in the meta-DEG list:", paste(intersect(panel, genes), collapse = ", "), "\n")
cat("Panel genes not in the list     :", paste(setdiff(panel, genes), collapse = ", "), "\n")


# 2. STRING mapping -----------------------------------------------------------
string_db <- STRINGdb$new(version = "12.0", species = 10090,
                          score_threshold = 700, input_directory = cache_dir)
cat("STRING version used:", string_db$version, "\n")
stopifnot(string_db$version == "12.0")

mapped <- string_db$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
cat("Unique genes mapped:", length(unique(mapped$gene)), "of", length(genes), "\n")
unmapped <- setdiff(genes, mapped$gene)
cat("Unmapped genes:", paste(unmapped, collapse = ", "), "\n")

dup_ids <- mapped$STRING_id[duplicated(mapped$STRING_id)]
cat("STRING ids shared by more than one symbol:", length(unique(dup_ids)), "\n")

nodes <- mapped[!duplicated(mapped$STRING_id), c("STRING_id", "gene")]
names(nodes) <- c("id", "gene")
cat("Network nodes:", nrow(nodes), "\n")


# 3. Edges --------------------------------------------------------------------
interactions <- string_db$get_interactions(nodes$id)
cat("Raw edge rows returned:", nrow(interactions), "\n")

pair_key <- apply(interactions[, c("from", "to")], 1,
                  function(x) paste(sort(x), collapse = "_"))
edges <- interactions[!duplicated(pair_key), c("from", "to", "combined_score")]
names(edges) <- c("source", "target", "combined_score")
edges <- edges[edges$source != edges$target, ]
cat("Unique edges:", nrow(edges), "\n")


# 4. Save ---------------------------------------------------------------------
write.csv(nodes, file.path(out_dir, "ppi_nodes.csv"), row.names = FALSE)
write.csv(edges, file.path(out_dir, "ppi_edges.csv"), row.names = FALSE)
writeLines(c(
  paste("genes_input:", length(genes)),
  paste("unique_genes_mapped:", length(unique(mapped$gene))),
  paste("unmapped:", paste(unmapped, collapse = ";")),
  paste("nodes:", nrow(nodes)),
  paste("raw_edge_rows:", nrow(interactions)),
  paste("unique_edges:", nrow(edges)),
  paste("string_version:", string_db$version)
), file.path(out_dir, "ppi_counts_summary.txt"))
cat("Saved to:", out_dir, "\n")

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "06_STRING_network_sessionInfo.txt"))
