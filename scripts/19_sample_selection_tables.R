# =============================================================================
# 19_sample_selection_tables.R
# Sample-selection tables: every GEO sample of the four datasets, with the
# decision (used or excluded) and the reason
#
# Project : Transcriptomic meta-analysis of microglia at 7 days after
#           traumatic brain injury (four mouse RNA-seq datasets)
# Role    : Step 19. Documents which samples entered scripts 01 to 04.
#
# Rules (the same rules the scripts 01 to 04 apply)
#   GSE167459 : GEO lists 10 microglia samples (5 Sham, 5 FPI). The primary
#               analysis (script 01) uses the 6 samples whose code contains
#               "M_S" (3 Sham, 3 FPI). The 4 samples with a "J-A" code
#               (GSM5105489, GSM5105490, GSM5105494, GSM5105495) are microglia
#               too but were not selected by that code pattern. They are marked
#               "sensitivity_only": they are used in the all-10 sensitivity
#               analysis (script 20). The sample GSM5105505 (column
#               1301M_S75) is an astrocyte sample according to its GEO record
#               and is excluded.
#   GSE253476 : wild-type microglia: 3 Sham (Sham_Mg_*) and 3 TBI at 7 dpi
#               (WT_7d_Mg_*). Macrophages, the Fgg390-396A mutant mice and the
#               1 day samples are not used.
#   GSE276647 : isotype-treated Sham and TBI microglia at 7 dpi (4 + 4). The
#               anti-CD3 arm and all 1 month samples are not used.
#   GSE283560 : sorted resident microglia at 7 dpi (3 Control, 8 TBI).
#               Monocytes, macrophages and other time points are not used.
#
# Input   : sample records downloaded from GEO (internet required)
# Output  : data/sample_selection/<dataset>_samples.csv  (one file per dataset)
#           Columns: dataset, gsm, title, group, status, reason, plus the
#           characteristics recorded in GEO. Status is "used" (primary
#           analysis), "sensitivity_only" (GSE167459 only, script 20) or
#           "excluded".
#           results/session_info/19_sample_selection_tables_sessionInfo.txt
#
# Usage   : run from the repository root, for example
#           Rscript scripts/19_sample_selection_tables.R
# Packages: GEOquery
# =============================================================================


# 0. Packages and paths ------------------------------------------------------
suppressPackageStartupMessages(library(GEOquery))
options(timeout = max(600, getOption("timeout")))

out_dir   <- file.path("data", "sample_selection")
cache_dir <- file.path("data", "geo_cache")
info_dir  <- file.path("results", "session_info")
for (d in c(out_dir, cache_dir, info_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Read the sample records of one GEO series. A series matrix file saved in
# data/geo_cache/ is used when it exists (no internet needed); otherwise it is
# downloaded from GEO.
get_samples <- function(gse) {
  local_file <- file.path(cache_dir, paste0(gse, "_series_matrix.txt.gz"))
  g <- if (file.exists(local_file)) {
    list(getGEO(filename = local_file, getGPL = FALSE))
  } else {
    getGEO(gse, GSEMatrix = TRUE, getGPL = FALSE, destdir = cache_dir)
  }
  pd <- do.call(rbind, lapply(g, function(x) {
    p <- Biobase::pData(x)
    keep <- intersect(c("title", "geo_accession", "description"), names(p))
    keep <- c(keep, grep(":ch1$", names(p), value = TRUE))
    p[, keep, drop = FALSE]
  }))
  names(pd)[names(pd) == "geo_accession"] <- "gsm"
  if ("description" %in% names(pd)) names(pd)[names(pd) == "description"] <- "sample_code"
  rownames(pd) <- NULL
  pd
}

finish <- function(df, dataset, expected_used) {
  df$dataset <- dataset
  front <- c("dataset", "gsm", "title", "group", "status", "reason")
  df <- df[, c(front, setdiff(names(df), front))]
  df <- df[order(df$status != "used", df$group, df$title), ]
  write.csv(df, file.path(out_dir, paste0(dataset, "_samples.csv")), row.names = FALSE)
  n_used <- sum(df$status == "used")
  cat("\n", dataset, ": ", nrow(df), " samples in GEO | used: ", n_used,
      " (expected ", expected_used, ")\n", sep = "")
  print(table(df$group[df$status == "used"]))
  if (n_used != expected_used)
    warning(dataset, ": number of used samples differs from the expected ", expected_used)
  invisible(df)
}


# 1. GSE167459 ----------------------------------------------------------------
# The sample code (for example 1303M_S76) is the name of the column in the count
# file. Script 01 selects the columns whose code contains "M_S".
d <- get_samples("GSE167459")
is_mg <- grepl("MG", d$title, ignore.case = TRUE)
d$group <- ifelse(grepl("\\bFPI\\b", d$title, ignore.case = TRUE), "FPI",
           ifelse(grepl("\\bSham\\b", d$title, ignore.case = TRUE), "Sham", NA))
code_ok <- grepl("M_S", d$sample_code)
is_microglia <- is_mg & d$gsm != "GSM5105505" & !is.na(d$group)
d$status <- ifelse(is_microglia & code_ok, "used",
            ifelse(is_microglia, "sensitivity_only", "excluded"))
d$reason <- ifelse(d$status == "used", "Microglia sample used in the primary analysis",
            ifelse(d$status == "sensitivity_only",
                   paste0("Microglia sample (cell type Microglia, 7 days). Not used in the primary ",
                          "analysis, which selected the sample codes containing M_S; used in the ",
                          "all-10 sensitivity analysis (script 20). Sample code ", d$sample_code),
            ifelse(d$gsm == "GSM5105505",
                   "Astrocyte sample according to its GEO record (its sample code 1301M_S75 matches the microglia pattern)",
                   "Astrocyte sample, not microglia")))
stopifnot(sum(d$status == "sensitivity_only") == 4)
finish(d, "GSE167459", 6)


# 2. GSE253476 ----------------------------------------------------------------
d <- get_samples("GSE253476")
celltype <- sub("^celltype: ", "", d[["celltype:ch1"]])
genotype <- sub("^genotype: ", "", d[["genotype:ch1"]])
treat    <- sub("^treatment: ", "", d[["treatment:ch1"]])
d$group  <- ifelse(treat == "Sham control", "Sham",
            ifelse(treat == "7 days post TBI", "TBI_7d", "TBI_1d"))
ok <- grepl("^Microglia", celltype) & genotype == "wildtype" &
      treat %in% c("Sham control", "7 days post TBI")
d$status <- ifelse(ok, "used", "excluded")
d$reason <- ifelse(ok, "Wild-type microglia, Sham or 7 dpi TBI",
            ifelse(!grepl("^Microglia", celltype), "Macrophage sample, not microglia",
            ifelse(genotype != "wildtype", "Fgg390-396A mutant mouse, not wild type",
                   "1 day after TBI, not the 7 dpi time point")))
d$analysis_name <- d$title            # name of the column in the count file
finish(d, "GSE253476", 6)


# 3. GSE276647 ----------------------------------------------------------------
d <- get_samples("GSE276647")
is_7d  <- grepl("_7d$", d$title)
is_iso <- grepl("_Iso_", d$title)
d$group <- ifelse(grepl("^Sham", d$title), "Sham",
           ifelse(is_iso, "TBI", "TBI_antiCD3"))
ok <- is_7d & is_iso
d$status <- ifelse(ok, "used", "excluded")
d$reason <- ifelse(ok, "Isotype-treated microglia at 7 dpi (Sham or TBI)",
            ifelse(is_7d, "Nasal anti-CD3 treatment arm (treated group, not used)",
            ifelse(grepl("Female_Severe", d$title),
                   "1 month time point (female mice, severe TBI), not the 7 dpi time point",
            ifelse(grepl("delayed", d$title),
                   "1 month time point (delayed-treatment experiment), not the 7 dpi time point",
                   "1 month time point, not the 7 dpi time point"))))
# Name of the column in the count file: Sham_Iso_1_7d -> Sham_1, TBI_Iso_1_7d -> TBI_Iso_1
d$analysis_name <- ifelse(!ok, NA_character_,
                   ifelse(grepl("^Sham", d$title),
                          sub("^Sham_Iso_([0-9])_7d$", "Sham_\\1", d$title),
                          sub("_7d$", "", d$title)))
finish(d, "GSE276647", 8)


# 4. GSE283560 ----------------------------------------------------------------
d <- get_samples("GSE283560")
celltype <- d[["cell type:ch1"]]
treat    <- d[["treatment:ch1"]]
time     <- d[["time:ch1"]]
d$group  <- ifelse(treat == "TBI", "TBI", "Control")
ok <- celltype == "Microglia" & time == "7dpi"
d$status <- ifelse(ok, "used", "excluded")
d$reason <- ifelse(ok, "Sorted resident microglia at 7 dpi",
            ifelse(celltype != "Microglia",
                   paste0("Not resident microglia (", celltype, ")"),
                   paste0("Not the 7 dpi time point (", time, ")")))
finish(d, "GSE283560", 11)

writeLines(capture.output(sessionInfo()),
           file.path(info_dir, "19_sample_selection_tables_sessionInfo.txt"))
