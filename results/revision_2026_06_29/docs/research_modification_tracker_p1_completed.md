# P1 research-modification tracker after added analyses

Date: 2026-06-29

Scope: maps the completed P1 manuscript-improvement analyses from `PDAC_CAF_myeloid_manuscript_revision_tasks_for_Codex.md` to evidence files, claim changes and manuscript insertion points. This is not a response-letter document.

## Executive status

The initial specificity and sensitivity suite is complete. The core story is now supported by stronger spatial nulls, anchor-component dissection, gene-overlap sensitivity, LN leave-one-out, all-LN H&E spatial maps and NMF rank stability. The next manuscript task is to keep integrating these results into Methods, Results, figure legends and source-data/provenance files.

## Tracker

| Task ID | Manuscript gap | Added analysis | Key evidence | Manuscript action | Claim boundary |
|---|---|---|---|---|---|
| P1-1 | Same-size random cores may be too weak because CAF cores are spatially contiguous tissue regions. | Spatially contiguous random-core null across 205 spatial samples, 1,000 null cores per sample. | `Supplementary_Table_4_Stronger_Null_Sensitivity.csv`; `stronger_null_contiguous_random_core_report.md`; `stronger_null_contiguous_random_core_summary.csv`. | Add Methods subsection "Stronger spatial null sensitivity analysis"; add Results sentence that CAF-core centering remains stronger than contiguous random tissue cores in dominant primary/liver contexts. | Supports tissue-architecture specificity beyond arbitrary contiguous regions; does not prove causality. |
| P1-2 | CAF-myeloid signal may be driven by CAF-only, myeloid-only or tumor-high anchors. | CAF-only, myeloid-only, CAF-myeloid combined and control-anchor comparison across 205 spatial samples. | `Supplementary_Table_4B_Anchor_Component_Comparison.csv`; `caf_myeloid_component_anchor_comparison_report.md`. | Add Results text that myeloid-only better centers IFN/MHC, immune-core and SPP1/TAM, whereas CAF-only/combined better centers TGF-beta/EMT. | Do not claim combined CAF-myeloid anchor is universally strongest. |
| P1-3 | Target gradients may be marker-overlap artifacts. | Module overlap matrix and overlap-removal sensitivity. | `Supplementary_Table_5_Module_Overlap_Sensitivity.csv`; `caf_myeloid_target_module_overlap_matrix.csv`; `gene_module_overlap_sensitivity_report.md`. | Add Methods subsection "Gene-module overlap sensitivity"; Results should emphasize no overlap for IFN/MHC, immune-core, tumor-aggressive and TGF-beta/EMT target modules. | Use SPP1/TAM and myCAF/matrix as core-component interpretation modules, not independent non-overlap targets. |
| P1-4 | LN metastasis finding may be driven by one of five samples. | Individual LN sample summary, leave-one-out analysis and all five LN H&E-anchored spatial maps. | `Supplementary_Table_6_LN_Leave_One_Out.csv`; `ln_metastasis_individual_sample_summary.csv`; `ln_metastasis_leave_one_out_report.md`; `Extended_Data_Figure_8_LN_Individual_Spatial_Maps.*`. | Keep Figure 2 language at "five-sample lymph-node metastasis subset suggests selective immune decoupling"; cite ED8 for all-sample visual evidence. | LN finding remains hypothesis-generating; do not define a definitive LN subtype. |
| P1-5 | NMF rank 4 may be arbitrary. | Rank 2-8 NMF stability on the 143-sample Stage 22 ecotype matrix, using 50 randomized NNDSVDar starts per rank and one NNDSVDa reference sweep. | `Supplementary_Table_7_NMF_Rank_Stability.csv`; `nmf_rank_stability_report.md`; ED7 panel E in `Extended_Data_Figure_7_Specificity_Sensitivity.*`; standalone NMF plots archived under `figures/archive_supporting/`. | Add Methods description of NMF input matrix, rank sweep and label annotation; cite ED7 for rank-sensitivity integration. | Rank 4 is a reproducible and interpretable working resolution, not a unique mathematical optimum. |

## Results language already integrated

### Stronger null

`To test whether CAF-core associations simply reflected spatial contiguity, we repeated the random-core analysis using same-size spatially contiguous null regions generated within each section. CAF-myeloid cores retained stronger target-program centering than these contiguous null regions across the dominant primary and liver-metastasis contexts, supporting a spatial-architecture signal beyond arbitrary contiguous tissue regions.`

### Anchor components

`Alternative-anchor analysis showed target-specific contributions of the CAF and myeloid components. Myeloid-only anchors most strongly centered IFN/MHC, immune-core and SPP1/TAM programs, whereas CAF-only or combined CAF-myeloid anchors better captured TGF-beta/EMT-associated gradients. These results support a CAF-dominant, myeloid-enriched stromal-myeloid architecture rather than a generic high-expression anchor.`

### Gene overlap

`The IFN/MHC, immune-core, tumor-aggressive and TGF-beta/EMT target modules shared no genes with the CAF-myeloid core definition, and their CAF-core gradients were directionally preserved in overlap-sensitivity analyses. In contrast, SPP1/TAM and myCAF/matrix modules are best interpreted as core-state annotation modules because they share marker content with the core definition.`

### LN leave-one-out and all-sample maps

`In the five-sample lymph-node metastasis subset, leave-one-out analysis preserved tumor-aggressive CAF-core coupling but showed weaker and more sample-sensitive IFN/MHC and immune-core gradients. H&E-anchored maps for all five LN samples show retained tumor-aggressive CAF-core coupling with heterogeneous immune/IFN organization. The lymph-node result is therefore treated as a hypothesis-generating immune-decoupling lead rather than a definitive metastatic subtype.`

### NMF rank stability

`NMF rank sensitivity supported rank 4 as a stable and interpretable working resolution of the CAF-core ecotype matrix. Across 50 randomized NNDSVDar starts, rank 4 explained 0.980 of the non-negative feature matrix with ARI 0.993, PAC 0.012 and component cosine reproducibility 1.000; higher ranks provided smaller reconstruction gains and began producing very small or empty dominant components.`

## Remaining manuscript-improvement work

- P0 manuscript declarations still require author-provided administrative details or formal submission placeholders.
- P0 consolidated Statistical analysis section should be checked once more against every displayed figure.
- Final language pass should search for strong causal, clinical, TLS and mechanism terms.
