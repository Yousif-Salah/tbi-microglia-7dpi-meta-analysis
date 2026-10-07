# check_purity_loo_against_reference.R
# Compares the results of scripts 15, 16 and 17 with the reference files.
# Run from the repository root.
tol <- 1e-6
all_ok <- TRUE
rel <- function(a, b) max(abs(a - b) / pmax(abs(b), 1e-300), na.rm = TRUE)

cat("##### Script 15: purity check\n")
nf <- "results/purity/purity_check_markers.csv"
if (file.exists(nf)) {
  new <- read.csv(nf, stringsAsFactors = FALSE)
  ref <- read.csv("reference/purity/purity_check_markers_reference.csv", stringsAsFactors = FALSE)
  key_n <- paste(new$dataset, new$sample); key_r <- paste(ref$dataset, ref$sample)
  cat("Rows new / reference:", nrow(new), "/", nrow(ref), "| same samples:", setequal(key_n, key_r), "\n")
  ok <- setequal(key_n, key_r)
  if (ok) {
    ref <- ref[match(key_n, key_r), ]
    for (g in c("P2ry12", "Cx3cr1", "Ptprc", "Cd3e")) {
      na_same <- identical(is.na(new[[g]]), is.na(ref[[g]]))
      d <- rel(new[[g]], ref[[g]])
      cat(sprintf("  %-7s max relative difference = %.3g   NA pattern identical: %s\n", g, d, na_same))
      if (!na_same || d > tol) ok <- FALSE
    }
    same_cond <- identical(new$condition, ref$condition); cat("  condition identical:", same_cond, "\n")
    ok <- ok && same_cond
  }
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 15 first\n")

cat("\n##### Script 16: leave-one-sample-out\n")
sf <- "results/loo/LOO_global_summary.csv"
if (file.exists(sf)) {
  new <- read.csv(sf, stringsAsFactors = FALSE)
  ref <- read.csv("reference/loo/LOO_global_summary_reference.csv", stringsAsFactors = FALSE)
  m <- merge(new, ref, by = "dropped_dataset", suffixes = c("_new", "_reference"))
  print(m[, c("dropped_dataset", "global_meta_DEG_count_new", "global_meta_DEG_count_reference",
              "Cybb_padj_new", "Cybb_padj_reference", "Fcer1g_padj_new", "Fcer1g_padj_reference")],
        row.names = FALSE)
  ok <- nrow(m) == 4 &&
    all(m$global_meta_DEG_count_new == m$global_meta_DEG_count_reference) &&
    all(m$flagged_sample_new == m$flagged_sample_reference) &&
    all(m$shrinkage_method_new == m$shrinkage_method_reference) &&
    rel(m$Cybb_padj_new, m$Cybb_padj_reference) < tol &&
    rel(m$Fcer1g_padj_new, m$Fcer1g_padj_reference) < tol &&
    rel(m$Btk_padj_new, m$Btk_padj_reference) < tol
  cat("Summary identical to the reference:", ok, "\n")
  for (ds in c("GSE167459", "GSE253476", "GSE276647", "GSE283560")) {
    pf <- paste0("results/loo/LOO_pooled_", ds, ".csv")
    if (!file.exists(pf)) { cat(ds, ": pooled table missing\n"); ok <- FALSE; next }
    a <- read.csv(pf, stringsAsFactors = FALSE)
    b <- read.csv(paste0("reference/loo/LOO_pooled_", ds, "_reference.csv"), stringsAsFactors = FALSE)
    same_genes <- setequal(a$SYMBOL, b$SYMBOL)
    cat(sprintf("%s: genes new / reference = %d / %d | same genes: %s", ds, nrow(a), nrow(b), same_genes))
    if (same_genes) {
      b <- b[match(a$SYMBOL, b$SYMBOL), ]
      d1 <- rel(a$meta_logFC, b$meta_logFC); d2 <- rel(a$meta_p, b$meta_p); d3 <- rel(a$meta_padj, b$meta_padj)
      cat(sprintf(" | max rel. diff logFC %.2g, p %.2g, padj %.2g | k identical: %s\n",
                  d1, d2, d3, identical(a$k, b$k)))
      if (max(d1, d2, d3) > tol || !identical(a$k, b$k)) ok <- FALSE
    } else { cat("\n"); ok <- FALSE }
  }
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 16 first\n")

cat("\n##### Script 17: interferon sets and three genes\n")
gf <- "results/loo/LOO_interferon_GSEA.csv"
if (file.exists(gf)) {
  a <- read.csv(gf, stringsAsFactors = FALSE); b <- read.csv("reference/loo/LOO_interferon_GSEA_reference.csv", stringsAsFactors = FALSE)
  ka <- paste(a$run, a$pathway); kb <- paste(b$run, b$pathway)
  ok <- setequal(ka, kb)
  if (ok) {
    b <- b[match(ka, kb), ]
    for (col in c("size", "NES", "pval", "padj")) {
      d <- rel(a[[col]], b[[col]]); cat(sprintf("  %-5s max relative difference = %.3g\n", col, d)); if (d > tol) ok <- FALSE
    }
  }
  a2 <- read.csv("results/loo/LOO_Cx3cr1_Sall1_Tjp1.csv", stringsAsFactors = FALSE)
  b2 <- read.csv("reference/loo/LOO_Cx3cr1_Sall1_Tjp1_reference.csv", stringsAsFactors = FALSE)
  k2a <- paste(a2$run, a2$SYMBOL); k2b <- paste(b2$run, b2$SYMBOL)
  ok2 <- setequal(k2a, k2b)
  if (ok2) {
    b2 <- b2[match(k2a, k2b), ]
    for (col in c("meta_logFC", "meta_p", "meta_padj")) {
      d <- rel(a2[[col]], b2[[col]]); cat(sprintf("  %-10s max relative difference = %.3g\n", col, d)); if (d > tol) ok2 <- FALSE
    }
  }
  ok <- ok && ok2
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL, see the differences printed above\n"); if (!ok) all_ok <- FALSE
} else cat("SKIPPED: run script 17 first\n")

cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one difference, see the output above\n")
