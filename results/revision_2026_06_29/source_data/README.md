# Revision source data manifest

Date: 2026-06-29

This directory stages the source data needed for the active revised manuscript:

`results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.md`

`source_data_file_manifest.csv` records file size, SHA256 checksum, row count, column count and leading column names for all staged source-data CSV files.

## Main figures

- `Source_Data_Fig_1.csv`: Figure 1 quantitative discovery and robustness panels.
- `Source_Data_Fig_2.csv`: Figure 2 GSE272362 site-level and random-core panels.
- `Source_Data_Fig_3A.csv`, `Source_Data_Fig_3B.csv`, `Source_Data_Fig_3C.csv` and `Source_Data_Fig_3_candidate_NC_style_panel_index.csv`: Figure 3 ecotype, candidate-axis and panel-index source data.
- `Source_Data_Fig_4A.csv`, `Source_Data_Fig_4A_multiresolution_scale.csv`, `Source_Data_Fig_4B.csv`, `Source_Data_Fig_4B_GSE274557.csv`, `Source_Data_Fig_4C_GSE274673.csv`, `Source_Data_Fig_4D_GSE274673.csv` and `Source_Data_Fig_4_NC_style_selected_xenium_cell_scores.csv`: Figure 4 multi-resolution validation and selected Xenium cell-map source data.

## Regenerated supplementary modules

- `Source_Data_supplementary_module1_spatial_specificity_robustness.csv`: Extended Data Figure 1 module.
- `Source_Data_supplementary_module2_metastatic_immune_decoupling.csv`: Extended Data Figure 2 module.
- `Source_Data_supplementary_module3_cell_state_multiresolution_validation.csv`: Extended Data Figure 3 module.
- `Source_Data_supplementary_module4_mechanism_interface_priority.csv`: Extended Data Figure 4 module.
- `Source_Data_supplementary_module5_pathology_tcga_tls_boundaries.csv`: Extended Data Figure 5 module.
- `Source_Data_supplementary_module6_spatial_architecture_mechanism_deepening.csv`: Extended Data Figure 6 module.

## Revision-specific Extended Data figures

- `Source_Data_Extended_Data_Fig_7_Specificity_Sensitivity.csv`: Extended Data Figure 7 contiguous-null, anchor-component, gene-overlap, LN leave-one-out and NMF-rank source data.
- `Supplementary_Table_8_Mechanism_Triangulation_Scoring.csv`: mechanism-triangulation scoring thresholds and assigned evidence scores.
- `Source_Data_Extended_Data_Fig_8_LN_Individual_Spatial_Maps.csv`: Extended Data Figure 8 all-sample LN metastasis spatial-map source data.

## Additional legacy main-figure source tables

The directory also retains older numbered source tables (`Source_Data_Fig_5A.csv`, `Source_Data_Fig_5B.csv`, `Source_Data_Fig_5C.csv`, `Source_Data_Fig_5D.csv`, `Source_Data_Fig_6A.csv`, `Source_Data_Fig_7A_B.csv` and `Source_Data_Fig_7C.csv`) because several current main-figure and regenerated-module panels were assembled from the same earlier analysis products.

## Notes

- Source data were copied from `results/source_data/` and the revision-specific analysis outputs without changing values.
- Raw spatial images, count matrices and RDS/H5AD objects remain in their original public-data locations or local raw-data mirrors; these large primary data files are not duplicated here.
