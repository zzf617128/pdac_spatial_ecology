# Revision progress summary

Date: 2026-06-29

## Completed in this revision round

Created `results/revision_2026_06_29/` with manuscript copies, staged current figures, submission-ready supplementary tables, analysis outputs and traceable revision documents.

## Supplementary tables

- `Supplementary_Table_1_Datasets.csv`: 8 dataset/context rows.
- `Supplementary_Table_2_Gene_Modules.csv`: 315 module-gene rows.
- `Supplementary_Table_3_Random_Core_Results.csv`: 1,002 random-core result rows.
- `Supplementary_Table_4_Stronger_Null_Sensitivity.csv`: 1,435 contiguous-null result rows.
- `Supplementary_Table_4B_Anchor_Component_Comparison.csv`: 8,610 anchor-comparison rows.
- `Supplementary_Table_5_Module_Overlap_Sensitivity.csv`: module-overlap and overlap-removal sensitivity summary.
- `Supplementary_Table_6_LN_Leave_One_Out.csv`: LN metastasis individual and leave-one-out summary.
- `Supplementary_Table_7_NMF_Rank_Stability.csv`: NMF rank 2-8 stability summary.

## P1-1 Stronger Spatial Null

Completed a spatially contiguous random-core null across 205 spatial samples, with 1,000 contiguous random cores per sample.

Key result: CAF-core centering remains stronger than arbitrary contiguous tissue cores in the major primary and liver-metastasis settings, while LN metastases retain tumor-aggressive, SPP1/TAM and TGF-beta/EMT centering more consistently than IFN/MHC or immune-core centering.

Reviewer-safe use: strengthens tissue-architecture specificity, but remains observational.

## P1-2 Anchor Component Comparison

Completed CAF-only, myeloid-only, CAF-myeloid combined and tumor-high anchor comparisons across 205 spatial samples.

Key result: the combined CAF-myeloid anchor is a broad balanced anchor, but myeloid-only better explains IFN/MHC, immune-core and SPP1/TAM gradients, whereas CAF-only/combined better explains TGF-beta/EMT.

Reviewer-safe use: frame the core as a CAF-dominant, myeloid-enriched stromal-myeloid architecture with target-specific component contributions, not as a universally superior combined marker.

## P1-3 Gene Module Overlap

Completed overlap matrix and overlap-removal sensitivity.

Key result: CAF-myeloid core markers do not overlap with IFN/MHC, immune-core, tumor-aggressive or TGF-beta/EMT target modules, supporting that these gradients are not a direct marker-overlap artifact. SPP1/TAM and myCAF/matrix are expected core-component interpretation modules and should not be used as independent non-overlap target modules.

Reviewer-safe use: use non-overlapping target modules for independent downstream claims; use SPP1/TAM and myCAF/matrix as cell-state/core-component interpretation.

## P1-4 LN Leave-One-Out

Completed individual-sample and leave-one-out analysis across the 5 LN metastasis samples.

Key result: LN metastases consistently retain tumor-aggressive CAF-core coupling, while IFN/MHC and immune-core gradients are weak or sample-sensitive. Leave-one-out median tumor-aggressive delta stays negative, whereas immune metrics vary more widely.

Reviewer-safe use: present LN metastasis as a hypothesis-generating extension and a candidate immune-decoupled state, not as a definitive subtype.

## P1-5 NMF Rank Stability

Completed rank 2-8 NMF stability on the 143-sample Stage 22 ecotype matrix, using 50 randomized NNDSVDar starts per rank plus one NNDSVDa reference fit matching the original implementation.

Key result: rank 4 explains 0.980 of the non-negative CAF-core feature matrix, with ARI 0.993, PAC 0.012 and component cosine reproducibility 1.000 across randomized starts. Rank 5-8 add smaller reconstruction gains and begin producing very small or empty dominant components.

Reviewer-safe use: state that rank 4 resolves four recurrent, reproducible CAF-core ecological axes used for interpretation. Do not claim four ecotypes are the unique mathematical optimum.

## Next Step

Continue converting these completed analyses into manuscript edits: Methods additions, Results claim refinements, Extended Data/Supplementary figure placement and source-data/provenance mapping.

## Manuscript Integration Completed

Generated `manuscript/Manuscript_NatureSubjournal_revised.md` and exported `manuscript/Manuscript_NatureSubjournal_revised.docx`.

Integrated the P1 evidence into:

- title and abstract;
- Results sections for stronger spatial nulls, LN leave-one-out, NMF rank stability, anchor-component comparison and gene-overlap sensitivity;
- Methods sections for contiguous nulls, alternative anchors, overlap sensitivity, LN leave-one-out, NMF rank stability and unified statistical analysis;
- Discussion claim boundaries;
- Figure and Extended Data figure legends.

DOCX structural QA passed, but LibreOffice visual render QA stalled; see `docs/docx_export_QA_note.md`.

## Figure Integration Completed

Generated the official Extended Data Figure 7 specificity and sensitivity suite (`figures/Extended_Data_Figure_7_Specificity_Sensitivity.pdf/.svg/.png`) with five panels:

- contiguous random-core null;
- CAF-only / myeloid-only / combined anchor dissection;
- marker-overlap sensitivity;
- LN individual and leave-one-out sensitivity;
- NMF rank stability.

The source data are in `analysis_outputs/extended_data_figure7_specificity_sensitivity_source_data.csv`.

Generated `figures/Extended_Data_Figure_8_LN_Individual_Spatial_Maps.pdf/.svg/.png` with all five GSE272362 lymph-node metastasis samples:

- H&E plus CAF-myeloid core rings;
- IFN/MHC, immune-core and tumor-aggressive program maps on the same H&E backgrounds;
- per-sample deltas, immune-decoupling index, dominant CAF-core ecotype, spot count and image-match count.

The source data are in `analysis_outputs/extended_data_figure8_ln_individual_spatial_maps_source_data.csv`.
