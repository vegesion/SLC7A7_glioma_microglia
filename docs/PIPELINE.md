# Pipeline

Scripts run in the order below. Each one sources `R/setup.R`, which in turn
loads `config/config.R`; nothing else needs to be configured inside a script.

Steps marked **slow** read the full single-cell atlas or the Visium objects and
take tens of minutes to hours depending on the machine.

---

## 1. Bulk cohorts and survival

| Script | Does | Writes |
|---|---|---|
| `01a_tcga_download_prepare.R` | Downloads TCGA-LGG and TCGA-GBM through the GDC and builds normalised count matrices | `data/processed/*_TCGA_norm.csv` |
| `01b_cgga_data_prepare.R` | Prepares CGGA-693 and CGGA-325 from the files downloaded by hand | `data/processed/CGGA_*` |
| `02a_curate_transporter_geneset.R` | Curates the amino acid transport gene set used for the module score | transport gene set |
| `02b_slc_cox_screen.R` | Multivariable Cox screen of the SLC family in both TCGA cohorts | `*_all_cox_zscore.csv` |
| `03_lasso_cox_signature.R` | LASSO-penalised Cox over the curated transport genes | signature objects |
| `04_multicohort_validation.R` | The survival models reported in the paper, across four cohorts, with IDH strata, time-dependent AUC and decision curves | `Fig3_*` |
| `05a_bulk_clinical_gsea_decon.R` | Deconvolution, clinical associations, and bulk correlation GSEA | `Fig4_*`, `Fig5_*` |
| `05b_infiltration_specificity.R` | Head-to-head comparison of SLC7A7 with pan-myeloid markers | head-to-head table |
| `13_collect_cox_event_counts.R` | Recovers n and event count for every Cox model without redrawing figures | event-count table |

## 2. Single-cell atlas -- slow

| Script | Does | Writes |
|---|---|---|
| `06_sc_preprocess_gbmap.R` | Reads the GBmap `.h5ad` directly and builds the on-disk matrix | BPCells matrix |
| `07a_donor_pseudobulk_icc.R` | The batch-structure gate: intraclass correlation of donor-level SLC7A7 attributable to dataset | `ICC_gate_summary.csv` |
| `07b_sc_atlas_donor.R` | Donor-level expression across all annotated cell types | atlas tables |
| `07c_sc_atlas_donor_microglia.R` | The same restricted to microglia | microglia tables |
| `07d_ucell_module_axes.R` | UCell module scores and the module screen against the earlier matched null | `UC_*.csv` |
| `08a_sc_myeloid_donor.R` | Donor quality control, variance components, myeloid subtype proportions | `Supp_predictor_QC.csv`, `Supp_microglia_ICC_depth.csv`, `Fig6_donor_subtype_proportion_regression.csv` |
| `08b_sc_within_between_models.R` | The Mundlak within-between mixed models over the 14 modules | `Fig6_withinbetween_modules.csv` |
| `08c_sc_within_between_gsea.R` | Gene-level within and between statistics, and Hallmark enrichment on each | `Fig6_wb_gene_stats.csv`, `Fig6_wb_gsea_*.csv` |
| `08d_sc_donor_gsea_heatmap.R` | Donor-wise enrichment and the depth-only negative control | donor-wise tables |
| `09_trajectory_exploratory.R` | Trajectory screen. Descriptive only; no test is applied and no claim rests on it | Fig 6 panels |

## 3. Spatial -- slow

| Script | Does | Writes |
|---|---|---|
| `10_spatial.R` | Section- and tumor-level analysis of 97 Visium sections from 12 tumors: composition after centered log-ratio transformation, partial correlation, histopathological regions, pseudobulk differential expression and enrichment, and the pan-myeloid compositional control | `TAMcomposition_*`, `partial_cor_*`, `SLC7A7_histopath.csv`, `spatial_*` |

## 4. Specificity and robustness controls

| Script | Does | Writes |
|---|---|---|
| `11a_matched_null_donor_level.R` | Matched-control gene null for the donor-level contrast, with two caliper widths and with and without excluding co-regulated genes | `Supp_null_*` |
| `11b_null_diagnostics.R` | Diagnostics on the matching: pool sizes, covariate distances | null diagnostics |
| `11c_camera_correlation_adjusted.R` | Composition-adjusted enrichment | `Supp_GSEA_composition_adj_*` |
| `12_matched_null_within_donor.R` | The same null built for the within-donor contrast | `Supp_wdnull_*` |

## 5. Figures

Everything in `figures/` reads results and draws. No figure script computes a
reported statistic.

| Script | Figure |
|---|---|
| `fig3_stratified_meta_forest.R` | Fig 3 strata and the random-effects meta-analysis |
| `fig4CD_specificity_panels.R` | Fig 4C, 4D |
| `fig6A_hallmark.R` | Fig 6A |
| `fig6_ucell_module_panels.R` | Fig 6 module panels |
| `fig7_gsea_all.R` | Fig 7A cross-layer panel |
| `supp_ucell_graph.R` | Supplementary UCell panel |
| `figures_main.R` | Remaining main-figure assembly |

---

## Seeds and determinism

A single seed, `SEED <- 42`, is set in `config/config.R` and used for LASSO
cross-validation, UCell scoring, gene set enrichment and every permutation.
`fgsea` is called with `eps = 0` and `nPermSimple = 100000` throughout. Both
matter: with the default `eps`, the most significant pathways return `NA`
rather than a p value, which truncates the tail of any null built on top of
them.
