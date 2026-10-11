# Output map

Which script produces which additional file and which figure. This is the
fastest way to check a number in the manuscript against the code.

## Additional files

| # | Content | Produced by | Source table |
|---|---|---|---|
| 1 | Solute carrier survival screen | `analysis/02b_slc_cox_screen.R` | `LGG_all_cox_zscore.csv`, `GBM_all_cox_zscore.csv` |
| 2 | Atlas quality control and variance components | `analysis/08a_sc_myeloid_donor.R`, `analysis/07a_donor_pseudobulk_icc.R` | `Supp_predictor_QC.csv`, `ICC_gate_summary.csv`, `Supp_microglia_ICC_depth.csv` |
| 3 | BV2 experiment | not computational | -- |
| 4 | Cox models and survival sensitivity | `analysis/04_multicohort_validation.R`, `analysis/13_collect_cox_event_counts.R`, `figures/fig3_stratified_meta_forest.R` | `Fig3_multivariable_cox.csv`, `Fig3_multivariable_full_*.csv`, `Fig3_stratified_cox_idh.csv`, `Fig3_meta_summary.csv`, `Fig3_timeROC_AUC.csv`, `Fig3_DCA.csv` |
| 5 | SLC3A2 donor-level co-expression | `analysis/11a_matched_null_donor_level.R` | `Supp_null_SLC3A2_partialcor_microglia.csv` |
| 6 | Supplementary figures | `analysis/08d_sc_donor_gsea_heatmap.R` (Fig. S1A), `analysis/08b_sc_within_between_models.R` (Fig. S1B) | -- |
| 7 | Within- and between-donor module coefficients | `analysis/08b_sc_within_between_models.R` | `Fig6_withinbetween_modules.csv` |
| 8 | Within-donor Hallmark enrichment | `analysis/08c_sc_within_between_gsea.R` | `Fig6_wb_gsea_within.csv` |
| 9 | Donor-level gene-level results | `analysis/08c_sc_within_between_gsea.R` | `Fig6_wb_gene_stats.csv` |
| 10 | Specificity and robustness controls | `analysis/11a_matched_null_donor_level.R`, `analysis/12_matched_null_within_donor.R`, `analysis/07d_ucell_module_axes.R`, `analysis/11c_camera_correlation_adjusted.R` | `Supp_null_*.csv`, `Supp_wdnull_*.csv`, `UC_empirical_p.csv`, `Supp_GSEA_composition_adj_microglia.csv` |
| 11 | Myeloid subtype proportion regression | `analysis/08a_sc_myeloid_donor.R` | `Fig6_donor_subtype_proportion_regression.csv` |
| 12 | Tumor-level spatial composition | `analysis/10_spatial.R` | `TAMcomposition_matched_tumorlevel.csv`, `partial_cor_SLC7A7_cellstate.csv`, `SLC7A7_histopath.csv` |

## Main figures

| Figure | Produced by |
|---|---|
| Fig 1 | schematic, not computational |
| Fig 2 | `analysis/02b_slc_cox_screen.R` + `figures/figures_main.R` |
| Fig 3 | `analysis/04_multicohort_validation.R` + `figures/fig3_stratified_meta_forest.R` |
| Fig 4 | `analysis/05a_bulk_clinical_gsea_decon.R` + `figures/fig4CD_specificity_panels.R` |
| Fig 5 | `analysis/07b_sc_atlas_donor.R`, `analysis/07c_sc_atlas_donor_microglia.R` |
| Fig 6 | `analysis/09_trajectory_exploratory.R`, `figures/fig6A_hallmark.R`, `figures/fig6_ucell_module_panels.R` |
| Fig 7 | `analysis/08c_sc_within_between_gsea.R`, `analysis/08d_sc_donor_gsea_heatmap.R`, `figures/fig7_gsea_all.R` |
| Fig 8 | `analysis/10_spatial.R` |
| Fig 9 | experimental, not computational |

## Key reported quantities

| Quantity | Script | Note |
|---|---|---|
| Intraclass correlation attributable to dataset (0.5620, 77 donors, 15 datasets) | `07a_donor_pseudobulk_icc.R` | `ICC_gate_summary.csv`. The gate that forced the donor-level design |
| Within-donor module coefficients | `08b_sc_within_between_models.R` | 14 modules; between-donor terms are adjustment terms |
| Depth-only negative control | `08d_sc_donor_gsea_heatmap.R` | The sign reversal that the specificity argument rests on |
| Matched-control gene null, both contrasts | `11a`, `12` | Both negative; reported as such |
| Pan-myeloid compositional control (spatial) | `10_spatial.R` | The second leg of the specificity argument |
