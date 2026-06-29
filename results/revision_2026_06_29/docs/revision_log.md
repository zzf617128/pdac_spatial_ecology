# Revision log

Initialized: 2026-06-29

## Completed

- Created revision workspace.
- Copied current manuscript, supplementary information and staged PDF display items.
- Generated Supplementary Table 1: datasets.
- Generated Supplementary Table 2: gene modules.
- Generated Supplementary Table 3: random-core results.
- Completed spatially contiguous random-core null across 205 spatial samples with 1,000 contiguous random cores per sample.
- Generated Supplementary Table 4: stronger null sensitivity.
- Added `analysis_outputs/stronger_null_contiguous_random_core_per_sample.csv`.
- Added `analysis_outputs/stronger_null_contiguous_random_core_summary.csv`.
- Added `docs/stronger_null_contiguous_random_core_report.md`.
- Completed CAF-only / myeloid-only / CAF-myeloid combined anchor comparison across 205 spatial samples with 1,000 same-size random anchors per sample.
- Generated Supplementary Table 4B: anchor component comparison.
- Added `analysis_outputs/caf_myeloid_component_anchor_comparison_per_sample.csv`.
- Added `analysis_outputs/caf_myeloid_component_anchor_comparison_summary.csv`.
- Added `docs/caf_myeloid_component_anchor_comparison_report.md`.
- Completed gene-module overlap sensitivity across 205 spatial samples.
- Generated Supplementary Table 5: module overlap sensitivity.
- Added `analysis_outputs/module_overlap_matrix_all_pairs.csv`.
- Added `analysis_outputs/caf_myeloid_target_module_overlap_matrix.csv`.
- Added `analysis_outputs/gene_module_overlap_sensitivity_per_sample.csv`.
- Added `analysis_outputs/gene_module_overlap_sensitivity_summary.csv`.
- Added `docs/gene_module_overlap_sensitivity_report.md`.
- Completed LN metastasis leave-one-out analysis across the 5 lymph-node metastasis samples.
- Generated Supplementary Table 6: LN leave-one-out sensitivity.
- Added `analysis_outputs/ln_metastasis_individual_sample_summary.csv`.
- Added `analysis_outputs/ln_metastasis_leave_one_out_summary.csv`.
- Added `docs/ln_metastasis_leave_one_out_report.md`.
- Completed NMF rank stability and rank-4 justification analysis for the 143-sample Stage 22 ecotype matrix.
- Generated Supplementary Table 7: NMF rank stability.
- Added `analysis_outputs/nmf_rank_stability_summary.csv`.
- Added `analysis_outputs/nmf_rank_stability_run_level.csv`.
- Added `analysis_outputs/nmf_rank_stability_component_labels.csv`.
- Added `analysis_outputs/nmf_rank_nndsvda_reference.csv`.
- Added NMF rank-stability outputs and later archived the standalone NMF rank figure under `figures/archive_supporting/` to avoid Extended Data numbering conflicts.
- Added `docs/nmf_rank_stability_report.md`.
- Added `docs/research_modification_tracker_p1_completed.md` mapping P1 manuscript-improvement gaps to completed analyses, evidence files and manuscript actions.
- Added `manuscript/Manuscript_NatureSubjournal_revised.md` with P1 analyses integrated into title, abstract, results, methods, discussion and figure legends.
- Added `manuscript/Manuscript_NatureSubjournal_revised.docx` exported from the revised Markdown source.
- Added `scripts/85_export_revised_manuscript_docx.py`.
- Added `docs/docx_export_QA_note.md`; DOCX structural QA passed, but LibreOffice visual render QA stalled and should be repeated before final submission.
- Added `figures/Extended_Data_Figure_7_Specificity_Sensitivity.*`, a multi-panel specificity and sensitivity suite.
- Added `analysis_outputs/extended_data_figure7_specificity_sensitivity_source_data.csv`.
- Updated revised manuscript and figure-panel map to cite Extended Data Figure 7.
- Added `figures/Extended_Data_Figure_8_LN_Individual_Spatial_Maps.*`, showing all five lymph-node metastasis samples as H&E-anchored maps for CAF core, IFN/MHC, immune core and tumor-aggressive programs.
- Added `analysis_outputs/extended_data_figure8_ln_individual_spatial_maps_source_data.csv`.
- Updated revised manuscript and figure-panel map to cite Extended Data Figure 8.
- Reframed this workspace as manuscript modification rather than response-letter preparation; added `docs/research_modification_tracker_p1_completed.md`.
- Updated manuscript declaration sections with formal pre-submission placeholders for acknowledgements, author contributions, funding, competing interests, ethics statement and supplementary information.
- Updated Data availability and Code availability language to include public dataset accessions and repository/publication-ready code availability expectations.
- Completed a terminology pass to remove obvious abbreviation inconsistencies such as `SPP1-TAM`, `immune core` and bracket-style author placeholders from the revised manuscript.
- Added `docs/p0_manuscript_completeness_audit.md` documenting completed P0 manuscript-completeness items and remaining author-provided information.
- Updated `docs/claim_boundary_checklist.md` and added `docs/claim_language_audit_2026_06_29.md` after causal, clinical, TLS, mechanism and cell-resolution language scans.
- Reorganized the Results section from six analysis-forward subsections into four story modules aligned with the manuscript-modification task list.
- Added `docs/results_reorganization_audit_2026_06_29.md`.
- Rewrote the main Figure 1-4 legends to align with the four Results modules and added direct Extended Data Figure 7/8 support references where appropriate.
- Rebuilt `docs/figure_panel_map.md` as a module-level evidence map and added `docs/figure_legend_alignment_audit_2026_06_29.md`.
- Added revision-local `source_data/` staging with 28 CSV files, `source_data/README.md` and `source_data/source_data_file_manifest.csv` for file sizes, SHA256 checksums and row/column counts.
- Updated `results/reports/submission_package_index.md` and refreshed `docs/revision_workspace_inventory.csv` to include the revision-local source-data staging.
- Added `docs/revision_task_completion_audit_2026_06_29.md` to map the original manuscript-modification task list to completed evidence, remaining author-input items and journal-formatting items.
- Added `Supplementary_Table_8_Mechanism_Triangulation_Scoring.csv` and cleaned Extended Data Figure 7/NMF rank-stability file naming so each official Extended Data figure number maps to one display item.
- Revised the active manuscript according to the latest remaining-issues brief: compressed the abstract, added a spatial-core definition, softened lymph-node wording, expanded contiguous-null implementation details, added random-seed reporting, shortened Code availability language, added a conceptual-contribution Discussion paragraph and cleaned internal figure-legend language.
- Regenerated Figure 1 with a main-panel spatially contiguous random-core specificity summary and synchronized `figures/Figure_1_current.pdf`.
- Added `figure_plot_code_archive/` with per-figure plotting-script copies, `panel_code_index.csv`, `plot_code_file_manifest.csv` and update rules for synchronizing future figure edits with source data and manuscript legends.
- Updated the manuscript title block, author affiliations, corresponding-author details, equal-contribution note, acknowledgements, author contributions, competing-interest statement and Code availability GitHub-reviewer-link placeholder using author-supplied metadata.

## Pending analyses

- None in the initial P1 manuscript-strengthening batch.
