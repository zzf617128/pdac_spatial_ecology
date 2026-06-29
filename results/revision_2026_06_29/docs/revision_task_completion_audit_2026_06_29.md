# Revision task completion audit

Date: 2026-06-29

Source task file:

`E:\PDAC_TLS\PDAC_CAF_myeloid_manuscript_revision_tasks_for_Codex.md`

## Executive conclusion

The scientific and manuscript-revision tasks are substantially complete. The core P0/P1 requirements have been implemented in the revised manuscript, supplementary tables, Extended Data figures, source data and provenance documents.

The remaining items are not scientific-analysis gaps. They are author- or journal-dependent finalization items: author metadata, funding, acknowledgements, author contributions, competing-interest confirmation, public repository/DOI, target-journal figure-format conversion and final visual render QA.

## Task-by-task status

| Task area | Status | Evidence | Residual item |
|---|---|---|---|
| Four-module Results structure | Complete | `Manuscript_NatureSubjournal_revised.md` Results now contains four modules: CAF-myeloid architecture, metastatic-site remodeling, cell-state/interface axes, independent validation/pathology context. | None. |
| Reduce analysis stacking | Complete | Results reorganized from analysis-forward sections into four story modules; audit in `results_reorganization_audit_2026_06_29.md`. | None. |
| Causal/mechanism language downgrade | Complete | `claim_language_audit_2026_06_29.md`; active manuscript explicitly states that matrix-integrin/SPP1-CD44 axes are candidates, not causal proof. | Continue guarding wording during later edits. |
| P0 declarations | Mostly complete | Acknowledgements, author contributions, funding, competing interests, ethics, data availability, code availability and supplementary information sections are present. | Requires author-provided information and repository URL before final submission. |
| Statistical analysis section | Complete | Methods lines under `Statistical analysis` define sample-level inference, random-core delta, empirical p-value formula, multiple testing, NMF, NNLS/reference and TCGA boundaries. | None. |
| Gene-module supplementary table | Complete | `Supplementary_Table_2_Gene_Modules.csv`: 315 rows across 34 modules with source/reference fields. | Source references can be polished later for journal style. |
| Abbreviation/terminology pass | Complete | No active-manuscript hits for obvious mixed forms such as `SPP1-TAM`, `IFN-MHC`, `TGFb`, `post NACT` or local-only code availability language. | Continue guarding during future edits. |
| P1-1 stronger null | Complete | `Supplementary_Table_4_Stronger_Null_Sensitivity.csv`; `stronger_null_contiguous_random_core_report.md`; ED7 panel. | Edge/density or autocorrelation-preserving nulls were optional and not required because contiguous null was completed. |
| P1-2 CAF-only/myeloid-only/combined anchor comparison | Complete | `Supplementary_Table_4B_Anchor_Component_Comparison.csv`; `caf_myeloid_component_anchor_comparison_report.md`; ED7 panel. | None. |
| P1-3 gene-overlap sensitivity | Complete | `Supplementary_Table_5_Module_Overlap_Sensitivity.csv`; overlap matrices and overlap-sensitivity report; ED7 panel. | None. |
| P1-4 LN leave-one-out and individual maps | Complete | `Supplementary_Table_6_LN_Leave_One_Out.csv`; `ln_metastasis_*` outputs; `Extended_Data_Figure_8_LN_Individual_Spatial_Maps.*`. | The LN result remains deliberately hypothesis-generating. |
| P1-5 NMF rank stability | Complete | `Supplementary_Table_7_NMF_Rank_Stability.csv`; `nmf_rank_stability_report.md`; ED7 panel. | None. |
| Figure 1 updates | Substantially complete | Figure 1 legend defines negative delta interpretation and points to ED7 for contiguous null, anchor comparison and overlap sensitivity. | Current main figure PDFs were staged; if journal requires direct in-panel relabeling rather than legend support, figure artwork can be regenerated. |
| Figure 2 updates | Complete in evidence, mostly complete in placement | Figure 2 legend frames LN as five-sample/hypothesis-generating; ED8 shows all LN maps; Supplementary Table 6 provides leave-one-out. | LN leave-one-out is in ED8/ED7/Table 6 rather than literally embedded into original ED2. This is scientifically adequate but can be rearranged if a target journal prefers all LN evidence in ED2. |
| Figure 3 updates | Complete | Figure 3 legend states ecotype annotations are loading-based and candidate axes are not causal ligand-receptor inference; NMF rank support in ED7/Table 7. | If desired, later figure-artwork refinement can convert dense heatmap panels to dot plots. |
| Figure 4 updates | Complete | Figure 4 legend states GSE274557 does not validate LN-specific decoupling and Xenium is targeted-panel/cell-resolution support, not causal or full-transcriptome validation. | None. |
| Extended Data modules | Complete | ED1-ED8 staged; ED7/ED8 added for specificity suite and all-LN maps. | Placement can be optimized for final journal figure limits. |
| Supplementary tables 1-7 | Complete | Supplementary Tables 1, 2, 3, 4, 4B, 5, 6 and 7 exist and have appropriate fields. | Could convert CSV to XLSX if requested by journal. |
| Source data staging | Complete | `results/revision_2026_06_29/source_data/` plus `source_data_file_manifest.csv`. | None. |
| Submission package | Complete as working package | `results/REVISION_SUBMISSION_PACKAGE_2026_06_29/` and `.zip` created with manifest. | Journal-specific upload conversion still pending. |

## Items not complete because they require author or target-journal input

- Final author list, affiliations, corresponding author and ORCID details.
- Funding bodies and grant numbers.
- Acknowledgements and resource acknowledgements.
- Author contributions.
- Competing-interest confirmation.
- Public repository, DOI or reviewer-access link for code/source tables.
- Target-journal figure-format conversion, if TIFF/EPS/separate panels are required.
- Final visual render QA for DOCX/PDF, because a previous LibreOffice render attempt stalled.

## Bottom line

All core scientific review/modification tasks have been addressed. The manuscript is not yet a final click-to-submit package only because administrative metadata, repository deposition, final visual QA and target-journal formatting remain outside the available author-independent analysis workflow.

