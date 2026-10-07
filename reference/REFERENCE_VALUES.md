# Reference values

Values a user should see when they rerun the pipeline. Small differences in the
last digits can come from package versions. Source of each value is noted.

## Per-dataset scripts

| Script | Check | Expected value | Source |
|---|---|---|---|
| 01_DEG_GSE167459.R | Samples | 3 Sham, 3 FPI | script check |
| 01_DEG_GSE167459.R | Raw matrix | 39,017 genes x 20 samples | script output |
| 01_DEG_GSE167459.R | Genes after filter = rows written | 16,211 | script output, output matches the reference file |
| 01_DEG_GSE167459.R | Genes with padj < 0.05 | 469 | script output |
| 01_DEG_GSE167459.R | Cybb | log2FC 0.855243, lfcSE 0.283739, padj 0.00834 | script output |
| 01_DEG_GSE167459.R | Fcer1g | log2FC 0.243837, padj 0.417183 | script output |
| 02_DEG_GSE253476.R | Samples | 3 Sham, 3 TBI | script check |
| 02_DEG_GSE253476.R | Raw matrix | 45,706 genes x 21 samples | script output |
| 02_DEG_GSE253476.R | Genes after filter = rows written | 19,016 | script output |
| 02_DEG_GSE253476.R | Cybb | log2FC 1.651064, padj 0.0008529941 | script output |
| 02_DEG_GSE253476.R | Fcer1g | log2FC 0.5335984, padj 0.0001590768 | script output |
| 03_DEG_GSE276647.R | Samples | 4 Sham, 4 TBI | script check |
| 03_DEG_GSE276647.R | Raw matrix | 15,718 genes x 12 columns (8 used) | script output |
| 03_DEG_GSE276647.R | Genes after filter = rows written | 14,901 | script output, output matches the reference file |
| 03_DEG_GSE276647.R | Genes with padj < 0.05 | 255 | script output |
| 04_DEG_GSE283560.R | Samples (7 dpi, microglia) | 3 Control, 8 TBI | script check |
| 04_DEG_GSE283560.R | Files in archive / microglia files used | 182 / 60 (20 samples x 3 lanes) | script output |
| 04_DEG_GSE283560.R | Cybb count, GSM8666315 (3 lanes summed) | 8158.299 | script output |
| 04_DEG_GSE283560.R | Genes before / after filter | 56,848 / 46,201 | script output |
| 04_DEG_GSE283560.R | Rows removed as duplicate symbols | 40 | script output |
| 04_DEG_GSE283560.R | Rows written | 46,161 (17,153 keep the Ensembl ID as symbol) | script output, output matches the reference file |
| 04_DEG_GSE283560.R | Genes with padj < 0.05 | 462 | script output |
| 04_DEG_GSE283560.R | Cybb | log2FC 1.656493, padj 0.156030 | script output |
| 04_DEG_GSE283560.R | Fcer1g | log2FC -0.068946, padj 0.991924 | script output |

## Software and annotation

| Item | Version |
|---|---|
| Ensembl annotation (GTF) | release 110, mouse GRCm39 (`Mus_musculus.GRCm39.110.gtf.gz`) |

## Pooling script

| Script | Check | Expected value | Source |
|---|---|---|---|
| 05_meta_analysis.R | Combined rows across 4 datasets | 96,289 | script output |
| 05_meta_analysis.R | Model used per gene | 16,634 REML, 2 fixed-effect (FE: Egln3, Gm2237) | script output |
| 05_meta_analysis.R | Analyzable genes (k >= 2) | 16,636 | script output, output matches the reference file |
| 05_meta_analysis.R | Meta-DEGs (meta_padj < 0.05) | 160 (117 up, 43 down) | script output, output matches the reference file |
| 05_meta_analysis.R | Console note | exactly 47 metafor warnings ("Fisher scoring ... Setting tau^2 = 0"): REML did not converge for 47 genes, tau^2 set to 0, genes kept | script output |

## Network scripts

| Script | Check | Expected value | Source |
|---|---|---|---|
| 06_STRING_network.R | STRING version | 12.0, mouse (10090), confidence score 700 | script output |
| 06_STRING_network.R | Genes in / mapped / unmapped | 160 / 154 / 6 (Mx1, Mx2, Clec7a, Iigp1c, AU020206, Rpl14-ps1) | script output |
| 06_STRING_network.R | Nodes | 154 | script output, matches the reference file |
| 06_STRING_network.R | Raw edge rows / unique edges | 388 / 194 | script output, matches the reference file |
| 06_STRING_network.R | BTK/NOX2 panel genes among the meta-DEGs | Cybb only | script output |
| 06b_Cytoscape_network.R | MCODE clusters (score cutoff 0.2) | 3 clusters: 14 nodes / 82 edges, 5 / 10, 4 / 5 | Cytoscape MCODE table, run by hand |

## Enrichment script

| Script | Check | Expected value | Source |
|---|---|---|---|
| 12_Enrichment.R | Universe / mapped to Entrez | 16,636 / 16,636 | script output |
| 12_Enrichment.R | Up genes / down genes mapped | 117 / 117, 43 / 43 | script output |
| 12_Enrichment.R | Significant terms, up (GO / KEGG / Reactome) | 230 / 23 / 8 | script output |
| 12_Enrichment.R | Significant terms, down (GO / KEGG / Reactome) | 21 / 0 / 0 | script output |
| 12_Enrichment.R | Cluster 1 (14 genes): GO terms before / after simplify; significant KEGG / Reactome | 380 / 129; 16 / 29 | script output |
| 12_Enrichment.R | Cluster 2 (5 genes): GO terms before / after simplify; significant KEGG / Reactome | 65 / 23; 21 / 8 | script output |
| 12_Enrichment.R | Cluster 3 (4 genes): GO terms before / after simplify; significant KEGG / Reactome | 52 / 17; 2 / 17 | script output |
| 12_Enrichment.R | Note | GO and Reactome tables are matches the reference files. KEGG is read online, so adjusted P values can shift slightly (largest change 0.007) when KEGG is updated. The set of significant KEGG terms is the same. | script output |

## k = 4 sensitivity script

| Script | Check | Expected value | Source |
|---|---|---|---|
| 09_k4_sensitivity.R | Genes by k (k = 2 / 3 / 4) | 2,846 / 2,567 / 11,223 (17.1% / 15.4% / 67.5%) | script output (check_k4_against_reference.R passes) |
| 09_k4_sensitivity.R | Meta-DEGs in the k = 4 subset | 155 | script output (check_k4_against_reference.R passes) |
| 09_k4_sensitivity.R | Cybb | meta_padj 0.0188, meta_padj_k4 0.0143 | script output (check_k4_against_reference.R passes) |
| 09_k4_sensitivity.R | Fcer1g | meta_padj 0.9998, meta_padj_k4 0.742 | script output (check_k4_against_reference.R passes) |
| 09_k4_sensitivity.R | Btk | meta_padj 0.9998, meta_padj_k4 0.9998 | script output (check_k4_against_reference.R passes) |

## Sample QC script

| Script | Check | Expected value | Source |
|---|---|---|---|
| 04b_sample_QC.R | Samples x genes per dataset | GSE167459 6 x 16,211; GSE253476 6 x 19,016; GSE276647 8 x 14,901; GSE283560 11 x 46,201 | script output |
| 04b_sample_QC.R | PC1 / PC2 percent variance | GSE167459 58.9 / 17.6; GSE253476 60.9 / 12.5; GSE276647 20.9 / 18.1; GSE283560 46.0 / 16.9 | script output |
| 04b_sample_QC.R | Within-group Spearman, flagged sample | 1304M_S77 0.934; Sham_Mg_R3 0.927; TBI_Iso_1 0.756; GSM8666331 0.608 | script output |
| 04b_sample_QC.R | Flagged sample has the lowest within-group correlation | Yes in GSE167459, GSE276647 and GSE283560. No in GSE253476 (lowest is WT_7d_Mg_R2, 0.924; range 0.924 to 0.934) | script output |

## Unshrunk sensitivity script (script 08)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 08_Sensitivity_unshrunk.R | Genes pooled / meta-DEGs | 16,521 / 1,411 (915 up, 496 down) | script output |
| 08_Sensitivity_unshrunk.R | Cybb | pooled log2FC 1.885, meta_padj 0.00033, I2 74.6% | script output |
| 08_Sensitivity_unshrunk.R | Fcer1g | pooled log2FC 0.513, 95% CI 0.307 to 0.719, meta_padj 5.9e-05, I2 0% | script output |

## Leave-one-dataset-out (scripts 10 and 11)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 11_LODO_genomewide.R | Meta-DEGs when GSE167459 / GSE253476 / GSE276647 / GSE283560 is dropped | 122 / 77 / 311 / 174 (range 77 to 311; baseline 160) | script output |
| 11_LODO_genomewide.R | Genes tested, same order | 15,693 / 15,368 / 15,946 / 13,845 | script output |
| 11_LODO_genomewide.R | Cybb meta_padj, same order | 7.8e-06 / 0.700 / 0.025 / 0.086 | script output |
| 10_LODO_Cybb_Fcer1g.R | Cybb pooled log2FC, same order | 2.11 / 1.68 / 1.19 / 1.67 (nominal p 3.0e-08 to 0.0059) | script output |
| 10_LODO_Cybb_Fcer1g.R | Fcer1g pooled log2FC, same order | 0.29 / 0.13 / 0.41 / 0.32 | script output |

## Purity script (script 15)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 15_purity_check.R | Cd3e | not detected after filtering in GSE253476 and GSE167459 | script output |
| 15_purity_check.R | GSM8666331 (GSE283560) | lowest P2ry12 (17.4, others 19.8 to 21.3), highest Cd3e (7.8) | script output |

## GSEA scripts (scripts 13 and 14)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 13_GSEA.R | Sets tested / significant (padj < 0.05): GO BP | 3,748 / 54 (37 positive, 17 negative) | script output |
| 13_GSEA.R | Hallmark | 50 / 4 (all positive) | script output |
| 13_GSEA.R | Reactome | 1,032 / 26 (24 positive, 2 negative) | script output |
| 13_GSEA.R | Interferon Gamma Response | NES 2.22, padj 4.2e-10 | script output |
| 13_GSEA.R | Interferon Alpha Response | NES 2.19, padj 5.0e-07 | script output |
| 13_GSEA.R | Cholesterol Homeostasis / Angiogenesis | NES 2.00, padj 0.0050 / NES 1.95, padj 0.0053 | script output |
| 14_GSEA_BTK_NOX2_panel.R | 14-gene panel, fgsea | NES 1.04, p 0.417 | script output |
| 14_GSEA_BTK_NOX2_panel.R | One-sample t tests of pooled log2FC | BTK (8 genes): mean -0.025, t(7) -0.50, p 0.631; NOX2 (6 genes): mean 0.334, t(5) 1.25, p 0.266 | script output |

## Leave-one-sample-out scripts (scripts 16 and 17)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 16_LOO_flagged_sample.R | Meta-DEGs after dropping the flagged sample of GSE167459 / GSE253476 / GSE276647 / GSE283560 | 136 / 158 / 228 / 105 (baseline 160) | script output |
| 16_LOO_flagged_sample.R | Cybb meta_padj, same order | 0.0747 / 0.00866 / 0.0170 / 0.0125 | script output |
| 16_LOO_flagged_sample.R | Cx3cr1, Sall1, Tjp1 | log2FC -0.46 to -0.65, meta_padj at most 0.0059 in every run | script output |
| 17_LOO_interferon_GSEA.R | Interferon Gamma Response across the four runs | NES 2.17 to 2.30, padj at most 1.1e-07 | script output |
| 17_LOO_interferon_GSEA.R | Interferon Alpha Response across the four runs | NES 2.16 to 2.26, padj at most 3.0e-06 (GSE283560 run) | script output |
| 16_LOO_flagged_sample.R | Safety check | stops if the unmodified pooled table does not give 160 meta-DEGs | script |

## Panel script (script 18)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 18_BTK_NOX2_forest_plots.R | Cybb | pooled log2FC 1.662, 95% CI 0.802 to 2.523, p 1.54e-04, meta_padj 0.0188, I2 72.8%, k 4 | script output |
| 18_BTK_NOX2_forest_plots.R | Other 13 genes | meta_padj 0.9998 | script output |
| 18_BTK_NOX2_forest_plots.R | Fcer1g | pooled log2FC 0.296, 95% CI 0.012 to 0.579, p 0.041, I2 43.6% | script output |
| 18_BTK_NOX2_forest_plots.R | Btk | pooled log2FC -0.102 | script output |

## All-ten-samples sensitivity script (GSE167459)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 20_GSE167459_all10_sensitivity.R | Samples | 10 (5 Sham, 5 FPI); 4 with batch J-A, 2 per group | script output |
| 20_GSE167459_all10_sensitivity.R | Meta-DEGs: primary / all10 ~ group / all10 ~ batch + group | 160 / 148 / 156 | script output |
| 20_GSE167459_all10_sensitivity.R | Overlap with the 160 primary meta-DEGs | 106 / 108 | script output |
| 20_GSE167459_all10_sensitivity.R | Cybb pooled padj: primary / ~ group / ~ batch + group | 0.0188 / 0.763 / 7.5e-06 | script output |
| 20_GSE167459_all10_sensitivity.R | Interferon Gamma NES (primary / ~ group / ~ batch + group) | 2.22 / 2.19 / 2.19 | script output |
| 20_GSE167459_all10_sensitivity.R | Interferon Alpha NES (primary / ~ group / ~ batch + group) | 2.19 / 2.14 / 2.14 | script output |
| 20_GSE167459_all10_sensitivity.R | Genes pooled: primary / ~ group / ~ batch + group | 16,636 / 16,870 / 16,870 | script output |
| 20_GSE167459_all10_sensitivity.R | Cx3cr1 pooled padj: primary / ~ group / ~ batch + group | 1.1e-05 / 0.118 / 0.071 | script output |
| 20_GSE167459_all10_sensitivity.R | Sall1 pooled padj: primary / ~ group / ~ batch + group | 9.7e-04 / 0.225 / 0.113 | script output |
| 20_GSE167459_all10_sensitivity.R | Tjp1 pooled padj: primary / ~ group / ~ batch + group | 1.8e-05 / 0.017 / 0.019 | script output |
| 20_GSE167459_all10_sensitivity.R | GSE167459 alone, Cybb (shrunk): primary / ~ group / ~ batch + group | log2FC 0.855, padj 0.0083 / log2FC 0.057, padj 0.47 / log2FC 1.41, padj 8.2e-05 | script output |
| 20_GSE167459_all10_sensitivity.R | Output files | Figure S9 (PCA, 85 x 90 mm, 300 dpi TIFF), Table S18 | script output |
| 20_GSE167459_all10_sensitivity.R | PCA of the ten samples | PC1 48%, PC2 24% (injury along PC1, batch along PC2) | script output |

## WGCNA script (GSE283560)

| Script | Check | Expected value | Source |
|---|---|---|---|
| 21_WGCNA_GSE283560.R | Fingerprint, Cybb count of GSM8666315 | 8158 | script output |
| 21_WGCNA_GSE283560.R | Genes after filter (CPM >= 1 in >= 3 samples) / genes in the network | 15,635 / 13,388 (2,234 without symbol; 13 duplicate symbols) | script output |
| 21_WGCNA_GSE283560.R | Primary run (11 samples, depth adjusted) | power 16, signed R2 0.825, mean connectivity 322; 4 blocks; 12 modules; 82 grey genes | script output |
| 21_WGCNA_GSE283560.R | Strongest module with TBI, primary run | purple (338 genes), r 0.598, family-wise permutation p 0.267; no module with family-wise p < 0.05 | script output |
| 21_WGCNA_GSE283560.R | BTK/NOX2 panel, primary run | 14 of 14 genes in the network, spread over 5 modules; Cybb in purple (kME 0.668) | script output |
| 21_WGCNA_GSE283560.R | Run without GSM8666331 (10 samples) | no power reached R2 0.80, so power 18 (R2 0.53); 12 modules; greenyellow (210 genes) r 0.919, family-wise permutation p 0.0083; Cybb in greenyellow | script output |
| 21_WGCNA_GSE283560.R | Leave-one-sample-out (11 networks) | family-wise p < 0.05 only after dropping GSM8666331 (greenyellow, 210 genes, p 0.0083) or GSM8666365 (magenta, 252 genes, r 0.961, p 0.0083); the other 9 runs: p 0.18 to 0.98 | script output |
| 21_WGCNA_GSE283560.R | GO terms (significant, BH < 0.05): 11 samples / 10 samples | 346 / 323 rows in the tables | script output |
| 21_WGCNA_GSE283560.R | Output | 14 files in `results/wgcna/final/` (Figures S7a-d and S8a-b, Tables S12 to S15) | script output |
| 21_WGCNA_GSE283560.R | Run time | several hours on a laptop (the leave-one-sample-out step builds 11 networks) | script output |
