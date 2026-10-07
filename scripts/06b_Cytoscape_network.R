# =============================================================================
# 06b_Cytoscape_network.R
# Draw and cluster the STRING network in Cytoscape
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 6b. Needs the Cytoscape desktop program (version 3.9 or later)
#           to be open on the same computer. This script cannot run without it.
#
# What the script does (RCy3 sends commands to the open Cytoscape)
#   1. Creates the network "PPI_meta_7dpi" from results/network/ppi_nodes.csv
#      and ppi_edges.csv (script 06).
#   2. Colours the nodes by direction of the pooled log2 fold change:
#      up #C0392B, down #2471A3.
#   3. Removes isolated nodes (degree 0) and makes the sub-network
#      "PPI_meta_7dpi_main". This is the network shown in the figure.
#   4. Applies a force-directed layout weighted by the STRING combined score.
#
# Clustering with MCODE (done by hand in Cytoscape, not by this script)
#   Apps > MCODE > Analyze Current Network, on the full network
#   "PPI_meta_7dpi" (all 154 nodes).
#   Main setting : degree cutoff 2, node score cutoff 0.2, haircut on,
#                  fluff off, K-core 2, maximum depth 100.
#   Robustness   : repeated with node score cutoff 0.1 and 0.3, all other
#                  settings unchanged. Isolated nodes have no edges, so the
#                  clusters are the same on the full and on the main network.
#   Cluster tables were exported from the MCODE result panel.
#
# Input
#   results/meta/Meta7dpi_metaDEGs_FDR0.05.csv
#   results/network/ppi_nodes.csv, results/network/ppi_edges.csv
#
# Output
#   Networks inside the open Cytoscape session (save the session from Cytoscape).
#
# Usage   : open Cytoscape, then run from the repository root
#           source("scripts/06b_Cytoscape_network.R")
# Packages: RCy3
# =============================================================================


# 0. Packages and input -------------------------------------------------------
suppressPackageStartupMessages(library(RCy3))

deg_file   <- file.path("results", "meta", "Meta7dpi_metaDEGs_FDR0.05.csv")
nodes_file <- file.path("results", "network", "ppi_nodes.csv")
edges_file <- file.path("results", "network", "ppi_edges.csv")

for (f in c(deg_file, nodes_file, edges_file))
  if (!file.exists(f))
    stop("File not found: ", f, ". Run scripts 05 and 06 first. Current folder: ", getwd())

deg   <- read.csv(deg_file, stringsAsFactors = FALSE)
nodes <- read.csv(nodes_file, stringsAsFactors = FALSE)
edges <- read.csv(edges_file, stringsAsFactors = FALSE)
stopifnot(all(c("SYMBOL", "meta_logFC") %in% names(deg)))

net_all  <- "PPI_meta_7dpi"
net_main <- "PPI_meta_7dpi_main"


# 1. Create the network in Cytoscape -----------------------------------------
cytoscapePing()
createNetworkFromDataFrames(nodes = nodes, edges = edges,
                            title = net_all, collection = "TBI_meta_7dpi")
setNodeLabelMapping(table.column = "gene", network = net_all)


# 2. Colour nodes by direction ------------------------------------------------
deg$direction <- ifelse(deg$meta_logFC > 0, "Up", "Down")
cat("Meta-DEGs up:", sum(deg$direction == "Up"),
    "| down:", sum(deg$direction == "Down"), "\n")

dir_table <- data.frame(gene = deg$SYMBOL, direction = deg$direction,
                        stringsAsFactors = FALSE)
setCurrentNetwork(net_all)
loadTableData(dir_table, data.key.column = "gene", table.key.column = "gene")

setNodeColorMapping(table.column = "direction",
                    table.column.values = c("Up", "Down"),
                    colors = c("#C0392B", "#2471A3"),
                    mapping.type = "d", style.name = "default")


# 3. Keep connected nodes only ------------------------------------------------
analyzeNetwork(directed = FALSE)
node_table <- getTableColumns(table = "node", network = net_all)
stopifnot("Degree" %in% names(node_table))

n_isolated <- sum(node_table$Degree == 0)
cat("Total nodes:", nrow(node_table),
    "| isolated (degree 0):", n_isolated,
    "| connected:", nrow(node_table) - n_isolated, "\n")

connected_ids <- node_table$SUID[node_table$Degree >= 1]
selectNodes(as.character(connected_ids), by.col = "SUID", network = net_all)
createSubnetwork(nodes = "selected", subnetwork.name = net_main, network = net_all)


# 4. Layout --------------------------------------------------------------------
setCurrentNetwork(net_main)
layoutNetwork(paste("force-directed defaultSpringCoefficient=0.0000008",
                    "defaultSpringLength=50 defaultNodeMass=3",
                    "edgeAttribute=combined_score"))

cat("Network '", net_main, "' is ready. Run MCODE on it (see header).\n", sep = "")
cat("If nodes overlap, increase defaultSpringLength (80 to 100) and rerun the layout.\n")
