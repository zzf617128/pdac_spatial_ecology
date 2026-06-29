# References, Software Versions, and Reproducibility Commands

Last updated: 2026-06-28

## Core Dataset References

1. Lyubetskaya A, Rabe B, Kavran A, Bai Y, Fisher A, Font-Tello A, et al. In situ multi-modal characterization of pancreatic cancer reveals tumor cell identity as a defining factor of the surrounding microenvironment. *Cell Reports*. 2026;45(1):116827. doi: `10.1016/j.celrep.2025.116827`. PMID: `41533516`.
   - GEO: `GSE282302`, Spatial Transcriptomics subseries.
   - GEO SuperSeries: `GSE310353`.
   - Relevant GEO sample metadata state post-neoadjuvant resected PDAC tumor and neoadjuvant FOLFIRINOX treatment.

2. Liu Y, Sinjab A, Min J, Han G, Paradiso F, Zhang Y, et al. Conserved spatial subtypes and cellular neighborhoods of cancer-associated fibroblasts revealed by single-cell spatial multi-omics. *Cancer Cell*. 2025;43(5):905-924.e6. doi: `10.1016/j.ccell.2025.03.004`. PMID: `40154487`. PMCID: `PMC12074878`.
   - GEO: `GSE274103`.
   - GEO title: "Spatial transcriptomics on treatment-naive pancreatic ductal adenocarcinoma (PDAC) patients".
   - GEO design: Visium spatial transcriptomics on FFPE tissues from five treatment-naive PDAC patient tissues.

3. Khaliq AM, Rajamohan M, Saeed O, Mansouri K, Adil A, Zhang C, et al. Spatial transcriptomic analysis of primary and metastatic pancreatic cancers highlights tumor microenvironmental heterogeneity. *Nature Genetics*. 2024;56(11):2455-2465. doi: `10.1038/s41588-024-01914-4`. PMID: `39294496`.
   - GEO: `GSE272362`.
   - Zenodo: `10.5281/zenodo.10712047`.
   - Local file used: `PDAC_Updated.rds`.

4. Pei S, et al. Metastatic pancreatic ductal adenocarcinoma lineage plasticity atlas. *Nature*. 2025. PMID: `40269162`.
   - GEO SuperSeries / spatial resources inspected: `GSE274557`, `GSE277782`.
   - Local GSE274557 files used: GEO series matrix, RAW file list, per-sample filtered feature-barcode HDF5 matrices and Space Ranger spatial archives.
   - Local GSE277782 files downloaded for feasibility and possible follow-up: CosMx SCT matrix, cell metadata, slide expression matrices, slide metadata and FOV positions.

5. GSE274673 Xenium PDAC cell-resolution dataset.
   - Local files used: GEO series family metadata, RAW file list, six per-sample Xenium archives, extracted `cells.csv.gz`, `cell_feature_matrix.h5`, feature/barcode files, metrics, clustering and UMAP files.
   - Samples analyzed: three treatment-naive and three chemoradiotherapy-treated PDAC sections.

## External Resource Links

- GSE282302 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE282302`
- GSE310353 GEO SuperSeries: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE310353`
- GSE274103 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274103`
- GSE272362 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272362`
- Zenodo 10712047: `https://zenodo.org/records/10712047`
- GSE274557 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274557`
- GSE277782 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277782`
- GSE274673 GEO: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274673`
- TCGA PAAD legacy Xena HiSeqV2: `https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/HiSeqV2.gz`
- TCGA PAAD legacy Xena clinical matrix: `https://tcga.xenahubs.net/download/TCGA.PAAD.sampleMap/PAAD_clinicalMatrix`

## Software Versions

R:

- R version 4.6.0 (2026-04-24 ucrt)
- Seurat 5.5.0
- SeuratObject 5.4.0
- data.table 1.18.4
- FNN 1.1.4.1
- RANN 2.6.2
- ggplot2 4.0.3
- png 0.1.9
- gridExtra 2.3
- patchwork 1.3.2

Python:

- Python 3.11.0
- numpy 2.4.3
- pandas 3.0.1
- matplotlib 3.10.8
- Pillow/PIL 12.1.1
- scipy 1.17.1
- scikit-learn 1.8.0
- h5py 3.16.0
- PyMuPDF 1.26.6

## Command-Level Reproducibility

Set the project root:

```powershell
$PROJECT = "E:\PDAC_TLS\pdac_spatial_ecology"
$PYTHON = "C:\Users\zzf61\AppData\Local\Programs\Python\Python311\python.exe"
$RSCRIPT = "C:\Program Files\R\R-4.6.0\bin\Rscript.exe"
```

Main processing and scoring:

```powershell
& $PYTHON "$PROJECT\scripts\03_mvp_score_visium.py"
& $PYTHON "$PROJECT\scripts\07_edge_background_qc.py"
& $PYTHON "$PROJECT\scripts\09_compute_caf_myeloid_niche_gradients.py"
& $PYTHON "$PROJECT\scripts\10_test_caf_myeloid_niche_gradients.py"
& $RSCRIPT "$PROJECT\scripts\11_score_gse272362_pdac_updated_rds.R" $PROJECT
```

Random-core controls and robustness:

```powershell
& $RSCRIPT "$PROJECT\scripts\15_mvp_random_core_permutation.R" $PROJECT
& $RSCRIPT "$PROJECT\scripts\14_gse272362_random_core_permutation.R" $PROJECT
& $RSCRIPT "$PROJECT\scripts\18_caf_core_threshold_sensitivity.R" $PROJECT
```

H&E overlays and manuscript figures:

```powershell
& $PYTHON "$PROJECT\scripts\17_make_gse272362_rds_he_overlays.py"
& $RSCRIPT "$PROJECT\scripts\19_make_manuscript_figure_drafts.R" $PROJECT
& $RSCRIPT "$PROJECT\scripts\20_gse272362_patient_matched_decomposition.R" $PROJECT
& $PYTHON "$PROJECT\scripts\21_mvp_he_morphology_screen.py"
& $PYTHON "$PROJECT\scripts\22_spatial_ecotype_deep_dive.py"
& $PYTHON "$PROJECT\scripts\23_score_gse235315_anchor.py"
& $PYTHON "$PROJECT\scripts\24_gse235315_random_core_anchor.py"
& $PYTHON "$PROJECT\scripts\25_mechanism_candidate_axes.py"
& $PYTHON "$PROJECT\scripts\26_make_submission_figure_suite.py"
& $RSCRIPT "$PROJECT\scripts\27_export_gse272362_target_gene_expression.R" $PROJECT
& $PYTHON "$PROJECT\scripts\27_targeted_gene_axis_validation.py"
& $PYTHON "$PROJECT\scripts\28_make_extended_data_robustness_figure.py"
& $PYTHON "$PROJECT\scripts\29_make_extended_data_he_morphology_figure.py"
& $PYTHON "$PROJECT\scripts\30_make_submission_cohort_summary.py"
& $PYTHON "$PROJECT\scripts\31_export_submission_docx.py"
& $PYTHON "$PROJECT\scripts\32_cxcl9_spp1_polarity_gap1.py"
& $PYTHON "$PROJECT\scripts\33_focused_ligand_receptor_gap3.py"
& $PYTHON "$PROJECT\scripts\34_cell_state_marker_attribution_gap2.py"
& $PYTHON "$PROJECT\scripts\35_prepare_gse111672_reference_summary.py"
& $PYTHON "$PROJECT\scripts\36_reference_projection_deconvolution_gap2.py"
& $PYTHON "$PROJECT\scripts\36_reference_projection_deconvolution_gap2.py" --reference full
& $PYTHON "$PROJECT\scripts\37_compare_gse202051_reference_projection.py"
& $PYTHON "$PROJECT\scripts\38_tcga_paad_bulk_context.py"
& $PYTHON "$PROJECT\scripts\39_make_extended_data_figure6_cell_state_reference_support.py"
& $PYTHON "$PROJECT\scripts\40_make_extended_data_figures8_9_mechanism_support.py"
& $PYTHON "$PROJECT\scripts\41_external_metastatic_pdac_feasibility.py"
& $PYTHON "$PROJECT\scripts\42_gse274557_pilot_caf_core_validation.py" --input-dir "$PROJECT\data\external\GSE274557\all_h5_spatial" --output-prefix gse274557_full --n-random 1000
& $PYTHON "$PROJECT\scripts\43_make_extended_data_figure10_gse274557_external_validation.py"
& $PYTHON "$PROJECT\scripts\44_make_extended_data_figure10_source_data.py"
& $PYTHON "$PROJECT\scripts\45_gse277782_cosmx_cell_neighborhood_validation.py"
& $PYTHON "$PROJECT\scripts\46_gse274673_xenium_pilot_expression_domain.py"
& $PYTHON "$PROJECT\scripts\47_gse274673_xenium_two_sample_anchor_sensitivity.py"
& $PYTHON "$PROJECT\scripts\48_gse274673_xenium_fixed_anchor_validation.py"
& $PYTHON "$PROJECT\scripts\49_make_candidate_figure4_multiresolution_validation.py"
& $PYTHON "$PROJECT\scripts\50_make_reference_style_supplement_figures.py"
& $PYTHON "$PROJECT\scripts\51_make_deep_reference_style_figures.py"
& $PYTHON "$PROJECT\scripts\52_make_ecotype_flow_figure.py"
& $PYTHON "$PROJECT\scripts\53_make_candidate_main_figure3_v2.py"
& $PYTHON "$PROJECT\scripts\54_make_candidate_main_figure3_nc_style.py"
& $PYTHON "$PROJECT\scripts\55_make_nc_style_main_figures.py"
```

## Main Reproducibility Outputs

- `results/tables/mvp_spot_level_scores_with_edge_qc.csv`
- `results/tables/gse272362_rds_spot_level_scores.csv`
- `results/tables/mvp_random_core_permutation_summary.csv`
- `results/tables/gse272362_rds_random_core_permutation_summary.csv`
- `results/tables/caf_core_threshold_sensitivity_summary.csv`
- `results/tables/gse272362_rds_overlay_manifest.csv`
- `results/figures/main/figure1_main.pdf`
- `results/figures/main/figure2_main.pdf`
- `results/figures/main/figure3_gse272362_matched_decomposition.pdf`
- `results/figures/main/figure4_he_morphology_screen.pdf`
- `results/figures/main/figure5_spatial_ecotype_deep_dive.pdf`
- `results/figures/main/figure6_gse235315_anchor.pdf`
- `results/figures/main/figure7_mechanism_candidate_axes.pdf`
- `results/figures/submission/figure1_submission_spatial_specificity.pdf`
- `results/figures/submission/figure1_submission_spatial_specificity_nc_style.pdf`
- `results/figures/submission/figure1_supplement_submission_post_nact_spatial_example.pdf`
- `results/figures/submission/figure2_submission_metastatic_decoupling.pdf`
- `results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.pdf`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes.pdf`
- `results/figures/submission/figure3_supplement_targeted_gene_axis_validation.pdf`
- `results/figures/submission/extended_data_figure4_he_morphology_bridge.pdf`
- `results/figures/submission/extended_data_figure5_external_anchor_robustness.pdf`
- `results/figures/submission/extended_data_gap1_cxcl9_spp1_polarity.pdf`
- `results/figures/submission/extended_data_gap2_cell_state_attribution.pdf`
- `results/figures/submission/extended_data_gap3_focused_lr_interface.pdf`
- `results/figures/submission/extended_data_figure10_gse274557_external_validation.pdf`
- `results/figures/submission/figure2_supplement_submission_spatial_examples.pdf`
- `results/tables/targeted_gene_axis_validation_summary.csv`
- `results/tables/submission_cohort_summary.csv`
- `results/tables/gap1_cxcl9_spp1_polarity_per_sample.csv`
- `results/tables/gap1_cxcl9_spp1_polarity_context_summary.csv`
- `results/tables/gap1_cxcl9_spp1_polarity_decoupling_correlations.csv`
- `results/tables/gap2_cell_state_marker_attribution_per_sample.csv`
- `results/tables/gap2_cell_state_marker_attribution_context_summary.csv`
- `results/tables/gap2_cell_state_marker_attribution_correlations.csv`
- `results/tables/gse111672_reference_cell_label_summary.csv`
- `results/tables/gse111672_reference_marker_means.csv`
- `results/tables/gse202051_reference_projection_signature_matrix.csv`
- `results/tables/gap2_reference_projection_deconvolution_per_sample.csv`
- `results/tables/gap2_reference_projection_deconvolution_context_summary.csv`
- `results/tables/gap2_reference_projection_deconvolution_correlations.csv`
- `results/tables/gse202051_full_reference_projection_signature_matrix.csv`
- `results/tables/gap2_full_reference_projection_deconvolution_per_sample.csv`
- `results/tables/gap2_full_reference_projection_deconvolution_context_summary.csv`
- `results/tables/gap2_full_reference_projection_deconvolution_correlations.csv`
- `results/tables/gap2_reference_projection_small_vs_full_comparison.csv`
- `results/tables/gap3_focused_lr_interface_per_sample.csv`
- `results/tables/gap3_focused_lr_interface_context_summary.csv`
- `results/tables/gap3_focused_lr_interface_correlations.csv`
- `results/tables/gse274557_series_metadata.csv`
- `results/tables/gse274557_raw_file_counts_by_sample.csv`
- `results/tables/gse277782_cosmx_file_feasibility.csv`
- `results/tables/gse277782_metadata_column_audit.csv`
- `results/tables/gse274557_full_caf_core_gradients.csv`
- `results/tables/gse274557_full_caf_core_enrichment.csv`
- `results/tables/gse274557_full_caf_core_context_summary.csv`
- `results/tables/gse274557_full_signature_coverage.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_10A.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_10B.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_10_Sample_Metadata.csv`
- `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.pdf`
- `results/figures/submission/extended_data_figure12_distance_to_caf_core_dynamics.pdf`
- `results/figures/submission/extended_data_figure13_xenium_cell_domain_maps.pdf`
- `results/figures/submission/extended_data_figure14_spatial_atlas_overview.pdf`
- `results/figures/submission/extended_data_figure15_local_spatial_program_maps.pdf`
- `results/figures/submission/extended_data_figure16_interface_compartment_maps.pdf`
- `results/figures/submission/extended_data_figure17_he_patch_examples.pdf`
- `results/figures/submission/extended_data_figure18_xenium_program_neighborhoods.pdf`
- `results/figures/submission/extended_data_figure19_random_core_null_diagnostics.pdf`
- `results/figures/submission/extended_data_figure20_ecotype_context_flow.pdf`
- `results/figures/submission/figure3_candidate_v2_ecotype_interface_story.pdf`
- `results/figures/submission/figure3_candidate_nc_style_ecotype_interface_story.pdf`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.pdf`
- `results/figures/submission/figure4_submission_multiresolution_validation_nc_style.pdf`
- `results/tables/gse274673_xenium_fixed_anchor_gradients.csv`
- `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv`
- `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv`
- `results/tables/gse274673_xenium_fixed_anchor_signature_coverage.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_11A_B.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_11C.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_11D.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_12D_core_to_far.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_13_selected_xenium_cell_scores.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_14A_atlas_counts.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_15_local_program_maps.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_16_interface_compartments.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_17_he_patch_examples.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_18_xenium_program_neighborhoods.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_19_random_core_null_diagnostics.csv`
- `results/source_data/Source_Data_Extended_Data_Fig_20_ecotype_context_flow.csv`
- `results/tables/gse277782_cosmx_annotation_composition.csv`
- `results/tables/gse277782_cosmx_caf_neighborhood_proximity.csv`
- `results/tables/gse277782_cosmx_caf_neighborhood_context_summary.csv`
- `results/reports/gse277782_cosmx_cell_neighborhood_validation.md`
- `results/tables/tcga_paad_bulk_context_scores.csv`
- `results/tables/tcga_paad_bulk_context_gene_coverage.csv`
- `results/tables/tcga_paad_bulk_context_axis_correlations.csv`
- `results/tables/tcga_paad_bulk_context_clinical_exploratory.csv`
- `results/reports/submission_cohort_summary.md`
- `results/reports/gap1_cxcl9_spp1_polarity_report.md`
- `results/reports/gap2_cell_state_marker_attribution_report.md`
- `results/reports/gse111672_reference_download_summary.md`
- `results/reports/gap2_reference_projection_deconvolution_report.md`
- `results/reports/gap2_full_reference_projection_deconvolution_report.md`
- `results/reports/gap2_reference_projection_small_vs_full_comparison.md`
- `results/reports/gap3_focused_ligand_receptor_interface_report.md`
- `results/reports/tcga_paad_bulk_context_report.md`
- `results/manuscript/pdac_caf_myeloid_spatial_niche_submission_v2_full.docx`
- `results/source_data/README.md`

## Reproducibility Caveats

- GSE272362 RDS-derived analyses do not yet use the same edge/background-risk filtering as GSE282302/GSE274103.
- Patient-level or outcome-level claims require additional metadata curation and are not part of the current reproducible claim.
- The H&E morphology screen uses simple patch-level color and texture features with grouped cross-validation; it is exploratory and not a clinical-grade H&E-only prediction model.
- GSE235315 was downloaded and scored as an external paired-ST anchor. Its sample metadata must be audited before any patient-level, treatment-level, or clinical-context claim.
- Candidate mechanism-axis analyses nominate pathway-level biology from module scores; they do not establish direct ligand-receptor signaling or causality.
- Marker-level cell-state attribution supports CAF/TAM and immune-state interpretation, but it is not formal reference-based deconvolution, segmentation, or orthogonal immunostaining.
- `GSE111672` is currently downloaded as a small reference prototype only. It is useful for method development but should not be the sole final reference for top-journal-grade deconvolution without comparison to a larger PDAC sc/snRNA atlas.
- `GSE202051` marker-constrained reference projection is a prototype abundance estimate. It strengthens cell-state interpretation but should not be described as final validated deconvolution or single-cell-resolved cell abundance.
- The full `GSE202051_totaldata-final-toshare.h5ad` validation preserves the main direction of the small-reference projection, but both analyses remain marker-constrained projection prototypes rather than orthogonal cell abundance measurements.
- TCGA PAAD bulk context supports broad axis relationships but is non-spatial and should not be used to validate CAF-core localization or clinical outcomes.
- GSE274557 validates broad CAF-core-centered organization across primary, liver, lung and peritoneal PDAC contexts, but it does not validate lymph-node-specific immune decoupling because lymph-node metastases are not present in the analyzed non-PDX Visium set.
- The first GSE277782 CosMx annotation-level nearest-neighbor stress test does not provide positive CAF-neighborhood validation and should not be promoted without a refined expression-domain or region-level model.
- GSE274673 Xenium supports cell-resolution antigen-presenting and SPP1/TAM-linked CAF-domain organization of immune/myeloid programs, but tumor epithelial and SPP1-tumor-like programs are not CAF-domain centered in this fixed-anchor analysis.
- Program-defined local-program, CAF/tumor/interface and Xenium neighborhood labels are expression-derived visualizations, not pathologist-annotated tissue compartments or formal cell-type abundance calls.
- The Extended Data Figure 20 ecotype flow summarizes sample-level dominant ecotype composition and should not be described as temporal transition or lineage progression.
- Representative H&E patch panels are visual examples from the exploratory morphology bridge and should not be described as clinical-grade pathology prediction.
- Random-core controls test spatial specificity against arbitrary same-size regions; they do not establish causality.
- The exported Word manuscript draft can be opened and structurally inspected with `python-docx`, but LibreOffice headless rendering hung during local PNG/PDF visual QA on 2026-06-25; visually inspect the `.docx` in Word before target-journal upload.
