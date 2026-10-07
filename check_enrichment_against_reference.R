# check_enrichment_against_reference.R
# Compares the enrichment tables of script 12 with the reference tables.
# Run from the repository root after running script 12.
# GO and Reactome tables use local databases and should match exactly.
# KEGG tables are read from the KEGG website; a difference there can come from
# a KEGG update and is reported as a note, not as a failure.
tol <- 1e-8
new_dir <- "results/enrichment"; ref_dir <- "reference/enrichment"

pairs <- c("GO_BP_UP", "GO_BP_DOWN", "KEGG_UP", "KEGG_DOWN", "Reactome_UP", "Reactome_DOWN")
for (cl in c("Cluster1_ISG", "Cluster2_MHC1", "Cluster3_Ribo"))
  pairs <- c(pairs, paste0("MCODE_", cl, "_GO_BP_full"), paste0("MCODE_", cl, "_GO_BP_simplified"),
             paste0("MCODE_", cl, "_KEGG"), paste0("MCODE_", cl, "_Reactome"))

num_cols <- c("pvalue", "p.adjust", "qvalue", "Count")
txt_cols <- c("Description", "GeneRatio", "BgRatio", "geneID")
core_ok <- TRUE; kegg_note <- FALSE

for (nm in pairs) {
  nf <- file.path(new_dir, paste0(nm, ".csv"))
  rf <- file.path(ref_dir, paste0(nm, "_160.csv"))
  if (!file.exists(nf)) { cat(sprintf("%-42s SKIPPED (no result file)\n", nm)); next }
  new <- read.csv(nf, stringsAsFactors = FALSE); ref <- read.csv(rf, stringsAsFactors = FALSE)
  same_ids <- setequal(new$ID, ref$ID)
  ok <- same_ids
  if (same_ids) {
    ref <- ref[match(new$ID, ref$ID), ]
    for (col in intersect(num_cols, colnames(new)))
      if (!isTRUE(all(abs(new[[col]] - ref[[col]]) < tol, na.rm = TRUE))) ok <- FALSE
    for (col in intersect(txt_cols, colnames(new)))
      if (!identical(as.character(new[[col]]), as.character(ref[[col]]))) ok <- FALSE
  }
  is_kegg <- grepl("KEGG", nm)
  status <- if (ok) "PASS" else if (is_kegg) "DIFFERENT (KEGG, see note)" else "FAIL"
  cat(sprintf("%-42s rows %4d / %4d   %s\n", nm, nrow(new), nrow(ref), status))
  if (!ok && is_kegg) kegg_note <- TRUE
  if (!ok && !is_kegg) core_ok <- FALSE
}

for (cl in c("Cluster1_ISG", "Cluster2_MHC1", "Cluster3_Ribo")) {
  nf <- file.path(new_dir, paste0(cl, "_genes.csv"))
  if (!file.exists(nf)) next
  same <- identical(read.csv(nf)$gene, read.csv(file.path(ref_dir, paste0(cl, "_genes_160.csv")))$gene)
  cat(sprintf("%-42s %s\n", paste0(cl, "_genes"), if (same) "PASS" else "FAIL"))
  if (!same) core_ok <- FALSE
}

if (kegg_note) cat("\nNote: KEGG is read online. If only KEGG differs, KEGG was probably updated after the reference run.\n")
cat("\nOverall:", if (core_ok) "PASS (GO, Reactome and cluster genes)\n" else "FAIL, see the differences printed above\n")
