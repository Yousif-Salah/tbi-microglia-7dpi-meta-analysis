# TBI Microglia 7-dpi Transcriptomic Meta-Analysis

**Analysis code for: "Reproducible Interferon-Centered Microglial Response at 7 Days after Traumatic Brain Injury: A Transcriptomic Meta-Analysis with Focus on the BTK/NOX2 Panel"**

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![R 4.5.1](https://img.shields.io/badge/R-4.5.1-276DC3.svg)

R code, sample-selection tables and parameters for a random-effects meta-analysis of four mouse microglia RNA-seq datasets (GEO: GSE167459, GSE253476, GSE276647, GSE283560) at 7 days after traumatic brain injury (TBI), with a secondary focus on a 14-gene BTK/NOX2 panel.

**Authors:** Yousif Salah, Ahmed Nadeem, Faleh Alqahtani\*  
Department of Pharmacology and Toxicology, College of Pharmacy, King Saud University, Riyadh, Saudi Arabia  
\*Corresponding author (contact details are in the manuscript)

**Manuscript:** submitted to *Biology* (MDPI). **Repository:** https://github.com/Yousif-Salah/tbi-microglia-7dpi-meta-analysis

## Contents

1. [Overview](#overview)
2. [Repository structure](#repository-structure)
3. [Quick start](#quick-start)
4. [Data](#data)
5. [Samples used](#samples-used)
6. [Run order](#run-order)
7. [Parameters](#parameters)
8. [Check that your run matches](#check-that-your-run-matches)
9. [Software](#software)
10. [How to cite](#how-to-cite)
11. [License](#license)

## Overview

The pipeline re-analyses raw counts of four independent datasets with one DESeq2 workflow (apeglm-shrunk log2 fold changes), pools the estimates gene by gene with a random-effects model (`metafor`, REML), and then characterises the result with a protein-protein interaction network (STRING, MCODE), over-representation analysis (GO, KEGG, Reactome), pre-ranked GSEA, and a focused test of the BTK/NOX2 panel. Sensitivity analyses cover unshrunk estimates, genes present in all four datasets, leave-one-dataset-out, leave-one-sample-out, and a ten-sample analysis of GSE167459. An exploratory WGCNA is run on GSE283560.

Main result of the primary analysis: 160 meta-DEGs (117 up, 43 down) among 16,636 analyzable genes (BH-adjusted p < 0.05).

## Repository structure

```text
.
├── scripts/                 Analysis scripts 01 to 21, numbered in run order
├── data/
│   ├── sample_selection/    One table per dataset: every GEO sample, used or excluded, with the reason
│   └── mcode/               MCODE cluster tables exported from Cytoscape (node score cutoff 0.1, 0.2, 0.3)
├── reference/               Expected output of key steps and REFERENCE_VALUES.md (used by the check scripts)
├── check_*.R                Check scripts: compare your output with reference/ and print PASS or FAIL
├── CITATION.cff             Citation metadata
└── LICENSE                  MIT license
```

Raw data are **not** stored in this repository. Download them from GEO (see [Data](#data)).

## Quick start

1. Clone the repository and open R (4.5.1 or later) in the repository root.
2. Install the packages listed under [Software](#software).
3. Download the input files (see [Data](#data)) into the folders shown.
4. Run the scripts in numbered order from the repository root, for example:

```bash
Rscript scripts/01_DEG_GSE167459.R
Rscript scripts/05_meta_analysis.R
```

Scripts 01 to 04 must run before all others. Script 21 takes several hours.

## Data

Download the files below and place them in the folders shown.

| Dataset | GEO | File(s) | Folder |
|---|---|---|---|
| GSE167459 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE167459 | `GSE167459_annotated_combined.counts.txt.gz` | `data/GSE167459/` |
| GSE253476 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE253476 | `GSE253476_salmon.merged.gene_counts_scaled_GEO.txt` (or `.gz`) | `data/GSE253476/` |
| GSE276647 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE276647 | `GSE276647_Supp_Table_2_7d.txt.gz` | `data/GSE276647/` |
| GSE283560 | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE283560 | `GSE283560_RAW.tar` | `data/GSE283560/` |
| Annotation | https://ftp.ensembl.org/pub/release-110/gtf/mus_musculus/ | `Mus_musculus.GRCm39.110.gtf.gz` (Ensembl release 110, GRCm39) | `data/reference/` |

Scripts 01, 04, 19 and 20 also read the GEO sample records at run time, so they need internet access. Scripts 06 (STRING) and 12 (KEGG) also need internet.

## Samples used

The full list is in `data/sample_selection/`. In short:

| Dataset | Primary analysis | Not used |
|---|---|---|
| GSE167459 | 6 microglia samples: 3 Sham, 3 FPI (selected by sample code, see script 01) | 4 microglia samples with codes 308J-A, 309J-A, 310J-A, 315J-A (`sensitivity_only`, analysed in script 20 and Table S18); 10 astrocyte samples (including 1301M_S75) |
| GSE253476 | 3 Sham, 3 TBI (wild-type microglia, 7 dpi) | Macrophages, mutant mice, 1 day samples |
| GSE276647 | 4 Sham, 4 TBI (isotype-treated, 7 dpi) | Anti-CD3 arm, 1 month samples |
| GSE283560 | 3 Control, 8 TBI (sorted microglia, 7 dpi) | Monocytes, macrophages, other time points |

For GSE167459, GEO has no batch field. The batch of the four extra samples is inferred from their sample codes ("J-A" codes versus "M" codes). On the PCA they form a separate group on PC2.

## Run order

Run each script from the repository root, for example `Rscript scripts/05_meta_analysis.R`. Each script writes a `sessionInfo` file to `results/session_info/`.

| Script | What it does | Main output | Manuscript item |
|---|---|---|---|
| `01` to `04` | DESeq2 per dataset, apeglm shrinkage | `results/deg_tables/<dataset>_7d.csv` (and `deg_tables_unshrunk/`) | Methods 2.2 |
| `04b_sample_QC.R` | Sample QC: library size, PCA, Spearman heatmap | `results/qc/` | Figures S1 to S4 |
| `05_meta_analysis.R` | Random-effects meta-analysis (REML), BH | `results/meta/Meta7dpi_all_genes_k2.csv`, `Meta7dpi_metaDEGs_FDR0.05.csv` | Table S9, Section 3.1 |
| `06_STRING_network.R`, `06b_Cytoscape_network.R` | STRING network; Cytoscape and MCODE (06b needs Cytoscape open) | `results/network/` | Figure 3 |
| `07_Figure2_volcano_forest.R` | Volcano and forest plots | `results/figures/` | Figure 2 |
| `08_Sensitivity_unshrunk.R` | Pooling with unshrunk fold changes | `results/meta_unshrunk/` | Table S4 |
| `09_k4_sensitivity.R` | Genes present in 2, 3 or 4 datasets | `results/sensitivity_k4/` | Table S5 |
| `10_LODO_Cybb_Fcer1g.R`, `11_LODO_genomewide.R` | Leave-one-dataset-out | `results/lodo/` | Table S6 |
| `12_Enrichment.R`, `12b_Figures_enrichment.R` | GO, KEGG, Reactome over-representation | `results/enrichment/`, `results/figures/` | Tables S10, S11; Figures 4, 5 |
| `13_GSEA.R` | Pre-ranked GSEA (GO BP, Hallmark, Reactome) | `results/gsea/` | Table S16, Figure 6 |
| `14_GSEA_BTK_NOX2_panel.R` | 14-gene panel: GSEA and gene-group t tests | `results/gsea_panel/` | Table 2 |
| `15_purity_check.R` | Microglia and immune-cell marker check | `results/purity/` | Table S3, Figure S5 |
| `16_LOO_flagged_sample.R`, `17_LOO_interferon_GSEA.R` | Leave-one-sample-out, interferon GSEA | `results/loo/` | Table S7 |
| `18_BTK_NOX2_forest_plots.R` | Panel table and forest plots | `results/panel/`, `figures/` | Table S17, Figures 7a, 7b, S6 |
| `19_sample_selection_tables.R` | Sample-selection tables | `data/sample_selection/` | Methods 2.2.1, Table 1 |
| `20_GSE167459_all10_sensitivity.R` | GSE167459 with all 10 microglia samples, with and without batch in the model | `results/sensitivity_all10/` | Table S18, Figure S9 |
| `21_WGCNA_GSE283560.R` | WGCNA of GSE283560 (11 samples and 10 samples, leave-one-sample-out) | `results/wgcna/final/` | Tables S12 to S15, Figures S7, S8 |

Scripts 01 to 04 must run before all others. Script 04 also writes the count table that script 21 reads, and it keeps the row with the highest baseMean when a gene symbol repeats. Script 21 takes several hours.

### Output file names and manuscript numbers

Some scripts were written before the final numbering of the manuscript. These output names differ from the manuscript labels:

| File written by the script | Manuscript label |
|---|---|
| `Table_S7_panel_per_gene.csv`, `Table_S7_formatted.csv` (script 18) | Table S17 |
| `Figure8a_panel_forest.tiff` (script 18) | Figure 7a |
| `Figure8b_Cybb_forest.tiff` (script 18) | Figure 7b |
| `FigureS7_Fcer1g_forest.tiff` (script 18) | Figure S6 |
| `FigureS6_purity_check.tiff` (script 15) | Figure S5 |
| `Figure7_GSEA.tiff` (script 13) | Figure 6 |
| `Enrichment_clusters.tiff` (script 12b) | Figure 4 |
| `Enrichment_overall.tiff` (script 12b) | Figure 5 |
| `QC_<dataset>_*.png` (script 04b) | Figures S1 (GSE253476), S2 (GSE283560), S3 (GSE276647), S4 (GSE167459) |

Script 21 already writes the final manuscript names into `results/wgcna/final/`. Figure 1 (workflow) was drawn in BioRender and is not produced by code. Tables S1 and S2 (screening records, study parameters) were compiled by hand.

## Parameters

| Step | Setting |
|---|---|
| Low-count filter | `rowSums(counts) >= 10`, applied before `DESeq()`, the same for all four datasets |
| DESeq2 model | `~ group` (or `~ condition`); Sham or Control is the reference level |
| Fold-change shrinkage | apeglm (`lfcShrink`, type `apeglm`) for all datasets. Unshrunk estimates are used only in script 08 |
| Outlier replacement | DESeq2 default (inactive for groups below 7 samples). GSE283560: `minReplicatesForReplace = Inf` |
| GSE283560 quantification | Salmon quant files (3 lanes per sample, summed), tximport, Ensembl release 110 (GRCm39) |
| Gene identifiers | Ensembl ID to symbol with `org.Mm.eg.db`; one row per symbol and dataset. In the pooling steps (scripts 05, 08, 09, 11, 16 and 20) the same rule is used: rows with a missing symbol, log2 fold change or standard error, or with a standard error of 0, are removed, and the first row per symbol is kept (tables are sorted by adjusted p value) |
| Meta-analysis | `metafor::rma.uni`, REML, genes in at least 2 datasets (k >= 2), `lfcSE > 0`; fixed-effect fallback if REML fails; Benjamini-Hochberg across genes. For a few genes metafor warns that the REML algorithm did not converge and sets the between-dataset variance (tau2) to 0; these genes are kept. Cybb and Fcer1g are not among them (tau2 = 0.51 and 0.035) |
| Meta-DEG | `meta_padj < 0.05` |
| STRING | v12.0, mouse (10090), combined score >= 700 |
| MCODE | degree cutoff 2, node score cutoff 0.2 (0.1 and 0.3 as robustness), K-core 2, maximum depth 100, haircut on, fluff off |
| Over-representation | `clusterProfiler` and `ReactomePA`, background = 16,636 analyzable genes, BH < 0.05; GO clusters simplified at similarity 0.7 |
| GSEA | `fgsea` 1.34.2, `msigdbr` 25.1.1 (human sets mapped to mouse), gene sets of 15 to 500 genes, 10,000 starting permutations, one core, `set.seed(1)`; ranking = sign(pooled log2FC) x -log10(pooled p); BTK/NOX2 panel as one set of 5 to 500 genes |
| Sample QC | `vst(blind = TRUE)`; PCA on the 500 most variable genes (scaled); Spearman correlation of samples; library size |
| Leave-one-sample-out | One sample per dataset removed, DESeq2 and pooling repeated: 1304M_S77 (GSE167459), Sham_Mg_R3 (GSE253476), TBI_Iso_1 (GSE276647), GSM8666331 (GSE283560) |
| All-ten-samples analysis | GSE167459 with 5 Sham and 5 FPI; designs `~ group` and `~ batch + group` |
| WGCNA (GSE283560 only) | Signed network, `blockwiseModules` (minModuleSize 150, mergeCutHeight 0.30, deepSplit 1, maxBlockSize 5000); soft power = lowest value with signed R2 >= 0.80, otherwise 18; library size removed with `limma::removeBatchEffect` before the network; exact permutation p values per module and family-wise; 11 leave-one-sample-out networks |

## Check that your run matches

After a script has run, use the check scripts (run from the repository root):

```r
source("check_against_reference.R")        # per-dataset tables
source("check_meta_against_reference.R")   # pooled results
source("check_qc_against_reference.R")     # sample QC (script 04b)
source("check_all10_against_reference.R")  # all-ten-samples analysis (script 20)
source("check_wgcna_against_reference.R")  # WGCNA (script 21)
```

The other check scripts match their step (`check_lodo_...`, `check_gsea_...`, `check_enrichment_...`, `check_k4_...`, `check_network_...`, `check_panel_...`, `check_purity_loo_...`, `check_unshrunk_...`). Expected values are in `reference/REFERENCE_VALUES.md`. Differences in the last digits can come from package versions.

## Software

R 4.5.1 (2025-06-13), Windows 10 x64. Bioconductor 3.21.

| Package | Version | Package | Version |
|---|---|---|---|
| DESeq2 | 1.48.1 | metafor | 4.8-0 |
| apeglm | 1.30.0 | fgsea | 1.34.2 |
| GEOquery | 2.76.0 | msigdbr | 25.1.1 |
| tximport | 1.36.1 | clusterProfiler | 4.16.0 |
| GenomicFeatures | 1.60.0 | ReactomePA | 1.52.0 |
| org.Mm.eg.db | 3.21.0 | STRINGdb | 2.20.0 |
| limma | 3.64.3 | RCy3 | 2.28.1 |
| data.table | 1.17.8 | ggplot2 | 4.0.0 |
| dplyr | 1.1.4 | ggrepel | 0.9.6 |
| purrr | 1.1.0 | pheatmap | 1.0.13 |
| tibble | 3.3.0 | patchwork | 1.3.2 |
| WGCNA | 1.73 | Cytoscape | 3.10.4 |
| AnnotationDbi | 1.70.0 | igraph | 2.1.4 |
| stringr | 1.5.2 | tidyr | 1.3.1 |

Figure text uses the Arial font. If Arial is not installed, R uses a default font and the figures look slightly different.

## How to cite

If you use this code, please cite the manuscript and this repository. GitHub shows a "Cite this repository" button based on `CITATION.cff`.

- Manuscript: Salah, Y.; Nadeem, A.; Alqahtani, F. *Biology* (MDPI), submitted.
- Archived release (v1.0.0): Zenodo, DOI [10.5281/zenodo.23212074](https://doi.org/10.5281/zenodo.23212074).

## License

Released under the [MIT License](LICENSE).
