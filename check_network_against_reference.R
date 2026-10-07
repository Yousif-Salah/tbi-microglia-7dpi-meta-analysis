# check_network_against_reference.R
# Compares the STRING network of script 06 with the reference node and edge files.
# Run from the repository root after running script 06.
nodes_new <- read.csv("results/network/ppi_nodes.csv")
nodes_ref <- read.csv("reference/ppi_nodes_reference.csv")
edges_new <- read.csv("results/network/ppi_edges.csv")
edges_ref <- read.csv("reference/ppi_edges_reference.csv")

key <- function(e) paste(pmin(e$source, e$target), pmax(e$source, e$target), sep = "_")

cat("Nodes new / reference:", nrow(nodes_new), "/", nrow(nodes_ref), "\n")
nodes_ok <- setequal(paste(nodes_new$id, nodes_new$gene), paste(nodes_ref$id, nodes_ref$gene))
cat("Same nodes (id and gene):", nodes_ok, "\n")

cat("Edges new / reference:", nrow(edges_new), "/", nrow(edges_ref), "\n")
edges_same <- setequal(key(edges_new), key(edges_ref))
cat("Same edge pairs:", edges_same, "\n")
score_ok <- FALSE
if (edges_same) {
  m <- match(key(edges_new), key(edges_ref))
  score_ok <- all(edges_new$combined_score == edges_ref$combined_score[m])
}
cat("Same combined scores:", score_ok, "\n")

cat("\nCounts file (new):\n"); cat(readLines("results/network/ppi_counts_summary.txt"), sep = "\n")
cat("\nOverall:", if (nodes_ok && edges_same && score_ok) "PASS\n" else "FAIL, see the differences printed above\n")
