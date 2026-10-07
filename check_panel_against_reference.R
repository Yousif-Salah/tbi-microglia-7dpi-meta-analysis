# Compare the Table S17 outputs of script 18 with the reference files.
# Run from the repository root after running script 18.
tol <- 1e-6
pairs <- list(
  c("results/panel/Table_S7_panel_per_gene.csv", "reference/panel/Table_S7_panel_per_gene.csv"),
  c("results/panel/Table_S7_formatted.csv",      "reference/panel/Table_S7_formatted.csv"))
for (p in pairs) {
  a <- read.csv(p[1], check.names = FALSE); b <- read.csv(p[2], check.names = FALSE)
  ok <- identical(names(a), names(b)) && identical(dim(a), dim(b))
  if (ok) {
    for (j in names(a)) {
      if (is.numeric(a[[j]])) ok <- ok && identical(is.na(a[[j]]), is.na(b[[j]])) &&
          all(abs(a[[j]] - b[[j]]) < tol, na.rm = TRUE)
      else ok <- ok && identical(a[[j]], b[[j]])
    }
  }
  cat(if (ok) "PASS" else "FAIL", "-", basename(p[1]), "\n")
}
