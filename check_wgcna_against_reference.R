# check_wgcna_against_reference.R
# Compares the WGCNA output of script 21 (results/wgcna/final/) with the
# reference files in reference/wgcna/.
# Run from the repository root after running script 21.
#
# Notes
#   - Tables S12 and S15 were written with rounded numbers (3 digits), so numbers
#     are compared with tolerance 1e-3.
#   - kME values (Tables S13) are compared with tolerance 1e-4.
#   - Module colours are given by module size, so the same colour is expected on
#     every run. If a colour differs, the module assignment of the genes differs.
#   - Table S14 (GO terms): the number of terms per module must be identical and
#     the adjusted p values must agree (relative tolerance 1e-3).
tol_round <- 1e-3
tol_kme   <- 1e-4
all_ok <- TRUE
fin <- "results/wgcna/final"
ref <- "reference/wgcna"

num_same <- function(a, b, tol) {
  all(abs(as.numeric(a) - as.numeric(b)) <= tol, na.rm = TRUE) &&
    identical(is.na(a), is.na(b))
}

compare_table <- function(label, f_new, f_ref, key = NULL, tol = tol_round) {
  cat("\n#####", label, "\n")
  if (!file.exists(file.path(fin, f_new)) || !file.exists(file.path(ref, f_ref))) {
    cat("SKIPPED: run script 21 first, or reference file missing\n"); return(invisible())
  }
  new <- read.csv(file.path(fin, f_new), stringsAsFactors = FALSE)
  rf  <- read.csv(file.path(ref, f_ref), stringsAsFactors = FALSE)
  ok <- identical(names(new), names(rf)) && nrow(new) == nrow(rf)
  if (ok && !is.null(key)) {
    new <- new[order(new[[key]]), ]; rf <- rf[order(rf[[key]]), ]
    ok <- identical(new[[key]], rf[[key]])
  }
  if (ok) {
    for (col in names(new)) {
      if (is.numeric(new[[col]]) && is.numeric(rf[[col]])) {
        good <- num_same(new[[col]], rf[[col]], tol)
        if (!good) { ok <- FALSE; cat("  DIFFERENT numbers in column:", col, "\n") }
      } else {
        if (!identical(as.character(new[[col]]), as.character(rf[[col]]))) {
          ok <- FALSE; cat("  DIFFERENT text in column:", col, "\n")
        }
      }
    }
  } else cat("  different columns, rows or keys\n")
  cat("  rows:", nrow(new), "| columns:", ncol(new), "\n")
  cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n")
  if (!ok) all_ok <<- FALSE
}

compare_table("Table S12: module statistics, 11 samples (primary run)",
              "Table_S12_module_stats_11samples.csv",
              "Table_S12_module_stats_11samples_reference.csv", key = "module")
compare_table("Table S12: module statistics, 10 samples (without GSM8666331)",
              "Table_S12_module_stats_10samples.csv",
              "Table_S12_module_stats_10samples_reference.csv", key = "module")
compare_table("Table S13: BTK/NOX2 panel genes and their modules",
              "Table_S13_BTK_NOX2_panel_modules_11samples.csv",
              "Table_S13_BTK_NOX2_panel_modules_11samples_reference.csv",
              key = "gene", tol = tol_kme)
compare_table("Table S13: module and kME of the 13,388 network genes",
              "Table_S13_gene_module_kME_11samples.csv",
              "Table_S13_gene_module_kME_11samples_reference.csv",
              key = "gene", tol = tol_kme)
compare_table("Table S15: leave-one-sample-out",
              "Table_S15_leave_one_out_results.csv",
              "Table_S15_leave_one_out_results_reference.csv", key = "dropped")

for (g in c("11samples", "10samples")) {
  cat("\n##### Table S14: GO terms,", g, "\n")
  fn <- file.path(fin, paste0("Table_S14_GO_terms_", g, ".csv"))
  fr <- file.path(ref, paste0("Table_S14_GO_terms_", g, "_reference.csv"))
  if (file.exists(fn) && file.exists(fr)) {
    new <- read.csv(fn, stringsAsFactors = FALSE); rf <- read.csv(fr, stringsAsFactors = FALSE)
    same_n <- identical(table(new$module)[sort(names(table(new$module)))],
                        table(rf$module)[sort(names(table(rf$module)))])
    k1 <- paste(new$module, new$ID); k2 <- paste(rf$module, rf$ID)
    same_ids <- setequal(k1, k2)
    rel <- NA
    if (same_ids) {
      m <- match(k2, k1)
      rel <- max(abs(new$p.adjust[m] - rf$p.adjust) / pmax(abs(rf$p.adjust), 1e-300))
    }
    cat("  terms:", nrow(new), "(reference", nrow(rf), ") | same terms per module:", same_n,
        "| same module-term pairs:", same_ids, "| max relative difference in p.adjust:",
        signif(rel, 3), "\n")
    ok <- same_n && same_ids && is.finite(rel) && rel < 1e-3
    cat(if (ok) "RESULT: PASS\n" else "RESULT: FAIL\n"); if (!ok) all_ok <- FALSE
  } else cat("SKIPPED: run script 21 first, or reference file missing\n")
}

cat("\n##### Output files (14 expected)\n")
files <- list.files(fin)
cat("  files in results/wgcna/final/:", length(files), "\n")
need <- c(paste0("Figure_S7", c("a_sample_tree_11samples.png", "b_soft_threshold_11samples.png",
                               "c_gene_dendrogram_11samples.tiff", "d_module_trait_heatmap_11samples.tiff")),
          paste0("Figure_S8", c("a_module_trait_heatmap_10samples.tiff", "b_GO_dotplot_greenyellow_10samples.tiff")),
          "Table_S12_marker_expression_log2CPM.csv", "Table_S12_module_stats_10samples.csv",
          "Table_S12_module_stats_11samples.csv", "Table_S13_BTK_NOX2_panel_modules_11samples.csv",
          "Table_S13_gene_module_kME_11samples.csv", "Table_S14_GO_terms_10samples.csv",
          "Table_S14_GO_terms_11samples.csv", "Table_S15_leave_one_out_results.csv")
miss <- setdiff(need, files)
if (length(miss)) { cat("  MISSING:", paste(miss, collapse = ", "), "\n"); all_ok <- FALSE
} else cat("  all 14 expected files found\n")

cat("\nOverall:", if (all_ok) "all checked files pass\n" else "at least one difference, see the output above\n")
