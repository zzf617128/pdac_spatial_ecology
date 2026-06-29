# Submission Figure Captions and Source Map

Last updated: 2026-06-29

## 2026-06-29 Revision Addendum

The active manuscript-revision draft is staged under `results/revision_2026_06_29/`. In that revision-specific figure numbering, two additional Extended Data figures have been added:

| display item | content | source |
|---|---|---|
| Extended Data Figure 7 | specificity suite integrating contiguous nulls, anchor-component comparison, gene-overlap sensitivity, LN leave-one-out and NMF rank stability | `results/revision_2026_06_29/analysis_outputs/extended_data_figure7_p1_specificity_source_data.csv` |
| Extended Data Figure 8 | all five GSE272362 lymph-node metastasis H&E-anchored spatial maps for CAF core, IFN/MHC, immune core and tumor-aggressive programs | `results/revision_2026_06_29/analysis_outputs/extended_data_figure8_ln_individual_spatial_maps_source_data.csv`; `results/tables/gse272362_rds_spot_level_scores.csv`; matched `data/raw/GSE272362` H&E images |

**Revision figure files:** `results/revision_2026_06_29/figures/Extended_Data_Figure_7_P1_Specificity_Suite.*` and `results/revision_2026_06_29/figures/Extended_Data_Figure_8_LN_Individual_Spatial_Maps.*`.

This addendum preserves the older full-archive figure numbering below while documenting the current revision-specific numbering used by `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.md`.

## Figure 1

**CAF-myeloid cores define reproducible spatial organizing regions in PDAC.**

**Take-home message:** CAF-myeloid-rich regions behave as spatial organizing cores for inflammatory, immune-core and tumor-aggressive programs, and this claim is robust to random-core controls, threshold perturbation and an external paired-ST anchor.

| panel | content | source |
|---|---|---|
| A | Cohort scale and evidence roles | manual cohort counts; `results/source_data/Source_Data_Main_Figures_NC_style_panel_index.csv` |
| B | Observed distance-to-core correlations versus 1,000 same-size random-core nulls across GSE282302, GSE274103 and GSE235315 | `results/source_data/Source_Data_Fig_1.csv`; `results/tables/gse235315_random_core_anchor_summary.csv` |
| C | CAF-core threshold sensitivity across top 15%, 10% and 5% definitions | `results/source_data/Source_Data_Fig_1.csv` |
| D | Representative random-core null intervals | `results/source_data/Source_Data_Extended_Data_Fig_19_random_core_null_diagnostics.csv` |
| E | Program decay from CAF-core to far regions | `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv` |
| F-J | Representative GSE282302 post-neoadjuvant H&E, CAF-core and spatial program maps | `results/tables/mvp_spot_level_scores_with_edge_qc.csv`; `results/tables/mvp_overlay_manifest.csv`; raw H&E image and scalefactor files |

**Generated figure:** `results/figures/submission/figure1_submission_spatial_specificity_nc_style.pdf`.

**Claim boundary:** This figure supports spatial organization around CAF-myeloid cores, not mature TLS formation, causality or clinical outcome prediction.

## Extended Figure 1

**Representative post-NACT PDAC section.**

**Take-home message:** A representative GSE282302 post-neoadjuvant section shows that CAF-myeloid cores occupy spatially coherent tissue regions and co-occur with IFN/MHC, tumor-aggressive and immune-core gradients.

Source: `results/tables/mvp_spot_level_scores_with_edge_qc.csv`, `results/tables/mvp_overlay_manifest.csv`, `data/raw/GSE282302/GSM8641105_C3_D8_ROI3_tissue_hires_image.png`, and `data/raw/GSE282302/GSM8641105_C3_D8_ROI3_scalefactors_json.json`.

## Figure 2

**Lymph-node metastases uncouple immune programs from CAF-myeloid cores.**

**Take-home message:** Primary tumors and liver metastases validate CAF-core-centered inflammatory and tumor-aggressive spatial organization, whereas lymph-node metastases retain stromal-tumor CAF coupling but decouple immune/IFN organization.

| panel | content | source |
|---|---|---|
| A | GSE272362 validation atlas composition by specimen site | `results/source_data/Source_Data_Fig_2.csv` |
| B | Random-core validation across primary tumors, liver metastases and lymph-node metastases | `results/source_data/Source_Data_Fig_2.csv` |
| C | CAF-core subprogram decomposition across specimen sites | `results/source_data/Source_Data_Fig_3A.csv` |
| D | Immune-decoupling index across contexts | `results/source_data/Source_Data_Fig_5C.csv` |
| E-H | Representative primary-tumor H&E, IFN/MHC, immune-core and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |
| I-L | Representative liver-metastasis H&E, IFN/MHC, immune-core and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |
| M-P | Representative lymph-node-metastasis H&E, IFN/MHC, immune-core and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |

**Generated figure:** `results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.pdf`.

**Claim boundary:** The lymph-node result is a spatial decoupling contrast from five samples; it should be framed as a strong biological lead, not a definitive clinical subtype.

## Extended Figure 2

**Representative GSE272362 spatial programs around CAF-myeloid cores.**

**Take-home message:** Large tissue-level overlays show the primary/liver positive pattern and the lymph-node immune-decoupling contrast at inspectable scale.

Source: `results/tables/gse272362_rds_spot_level_scores.csv`, `results/tables/gse272362_rds_overlay_manifest.csv`, matched raw H&E images and scalefactor files.

## Figure 3

**CAF-core ecotypes expose invasive-interface and immune-coupling axes.**

**Take-home message:** CAF cores are not one uniform compartment; they separate into recurrent spatial ecotypes and expose SPP1-TAM/matrix and TGF-beta/EMT invasive-interface axes, with immune/APC and lymphoid axes inversely related to immune decoupling. The NC-style main figure combines cohort-level ecotype evidence, candidate-axis statistics and representative interface maps in one A-R panel.

| panel | content | source |
|---|---|---|
| A | CAF-core NMF ecotype loadings | `results/source_data/Source_Data_Fig_5A.csv` |
| B | Context-to-CAF-core ecotype architecture flow | `results/source_data/Source_Data_Extended_Data_Fig_20_ecotype_context_flow.csv`; `results/tables/spatial_ecotype_context_counts.csv` |
| C | Dominant CAF-core ecotypes by context | `results/source_data/Source_Data_Fig_5B.csv`; `results/tables/spatial_ecotype_context_counts.csv` |
| D | Immune-decoupling index by dominant CAF-core ecotype | `results/tables/spatial_ecotype_sample_summary.csv` |
| E | Candidate-axis CAF-core and interface enrichment | `results/source_data/Source_Data_Fig_7A_B.csv` |
| F | Candidate-axis correlations with immune-decoupling index | `results/source_data/Source_Data_Fig_7C.csv` |
| G-J | Representative primary tumor compartment, SPP1/TAM, TGF-beta and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |
| K-N | Representative liver-metastasis compartment, SPP1/TAM, TGF-beta and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |
| O-R | Representative lymph-node-metastasis compartment, SPP1/TAM, TGF-beta and tumor-aggressive spatial maps | `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; matched raw H&E images and scalefactor files |

**Generated figure:** `results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.pdf`.

**Panel index:** `results/source_data/Source_Data_Fig_3_candidate_NC_style_panel_index.csv`.

**Claim boundary:** These are pathway/module nominations and representative spatial examples for follow-up validation, not direct ligand-receptor or perturbational mechanism proof. Spatial maps are representative visual anchors; quantitative support is provided by the cohort-level panels and companion source tables.

## Figure 4

**Independent Visium and Xenium data validate CAF-domain organization at complementary resolutions.**

**Take-home message:** Independent external Visium and Xenium datasets support the CAF-domain model at complementary resolutions: GSE274557 validates broad CAF-core-centered organization across primary and metastatic organ contexts, while GSE274673 supports cell-resolution CAF-APC and CAF-SPP1/TAM immune/myeloid domains.

| panel | content | source |
|---|---|---|
| A | External validation scale across GSE274557 Visium tissue contexts and GSE274673 Xenium treatment contexts | `results/source_data/Source_Data_Fig_4A_multiresolution_scale.csv`; `results/tables/gse274557_full_caf_core_gradients.csv`; `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv` |
| B | GSE274557 Visium CAF-core validation across primary, liver, lung and peritoneal contexts | `results/source_data/Source_Data_Fig_4B_GSE274557.csv`; `results/tables/gse274557_full_caf_core_context_summary.csv` |
| C | GSE274673 Xenium fixed-anchor cohort summary for CAF-APC and CAF-SPP1/TAM anchors | `results/source_data/Source_Data_Fig_4C_GSE274673.csv`; `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv` |
| D | GSE274673 CAF-SPP1/TAM anchor deltas across all six Xenium sections | `results/source_data/Source_Data_Fig_4D_GSE274673.csv`; `results/tables/gse274673_xenium_fixed_anchor_gradients.csv` |
| E-L | Representative treatment-naive and CRT Xenium cell maps for CAF-SPP1/TAM anchor, SPP1/TAM, IFN/APC and tumor epithelial programs | `results/source_data/Source_Data_Fig_4_NC_style_selected_xenium_cell_scores.csv`; `results/tables/gse274673_xenium_fixed_anchor_cell_scores.csv` |

**Generated figure:** `results/figures/submission/figure4_submission_multiresolution_validation_nc_style.pdf`.

**Claim boundary:** This main figure supports independent multi-resolution spatial validation. It does not test lymph-node-specific immune decoupling in GSE274557, does not support direct CAF-to-tumor epithelial proximity in GSE274673 and does not establish causality.

## Extended Figure 3

**Targeted gene-level support for candidate CAF-core axes.**

**Take-home message:** SPP1-TAM/matrix and TGF-beta/EMT invasive genes are recurrently enriched in CAF cores across tumor contexts, while IFN/APC antigen and T-cell/checkpoint genes are attenuated in lymph-node metastasis CAF cores.

Source: `results/tables/gse272362_rds_target_gene_expression.csv`, `results/tables/targeted_gene_axis_validation_per_sample.csv`, and `results/tables/targeted_gene_axis_validation_summary.csv`.

**Claim boundary:** This supports candidate-axis nomination at the targeted-gene level; it is not ligand-receptor, perturbational or causal validation.

## Extended Figure 4

**H&E morphology provides an exploratory pathology bridge.**

**Take-home message:** Simple H&E patch color/texture features recover the strongest signal for the CAF-myeloid program, with weaker but above-shuffle signal for tumor-aggressive, IFN/MHC and immune-core programs.

| panel | content | source |
|---|---|---|
| A | Median within-sample H&E feature correlations with spatial programs | `results/source_data/Source_Data_Fig_4A.csv`; `results/tables/mvp_he_patch_feature_correlation_summary.csv` |
| B | Grouped held-out ridge-model performance and target-shuffle controls | `results/source_data/Source_Data_Fig_4B.csv`; `results/tables/mvp_he_patch_grouped_cv_metrics.csv` |
| C | Patch-feature workflow schematic for exploratory H&E bridge analysis | `results/tables/mvp_he_patch_morphology_features.csv`; `results/tables/mvp_he_patch_grouped_cv_metrics.csv` |

**Claim boundary:** This supports a pathology-visible stromal correlate, not a clinical-grade H&E-only predictor or independent histopathology model.

## Extended Figure 5

**External paired-ST anchor and robustness analyses.**

**Take-home message:** Cross-cohort random-core support, CAF-core threshold sensitivity and GSE235315 per-sample deltas reinforce the core spatial claim while showing that immune-maturity-like signals are weaker.

| panel | content | source |
|---|---|---|
| A | Random-core support across discovery, validation, metastasis and external-anchor contexts | `results/tables/mvp_random_core_permutation_summary.csv`; `results/tables/gse272362_rds_random_core_permutation_summary.csv`; `results/tables/gse235315_random_core_anchor_summary.csv` |
| B | CAF-core threshold sensitivity across 15%, 10% and 5% core definitions | `results/tables/caf_core_threshold_sensitivity_summary.csv` |
| C | GSE235315 per-sample external-anchor deltas against random same-size cores | `results/tables/gse235315_random_core_anchor_summary.csv` |
| D | Cross-context support fraction for each spatial program | `results/source_data/Source_Data_ED_Fig_5D_cross_context_support_fraction.csv` |

**Claim boundary:** This figure supports robustness and external anchoring; it does not support mature TLS, clinical outcome or causal-mechanism claims.

## Extended Data Figure 8

**SPP1/TAM enrichment and CXCL9/IFN attenuation around CAF-myeloid cores.**

**Take-home message:** CAF cores are consistently SPP1/TAM-high across tumor contexts, while immune-decoupled samples are better distinguished by attenuated CXCL9/IFN CAF-core coupling than by a stronger combined SPP1-high/CXCL9-low ratio alone.

| panel | content | source |
|---|---|---|
| A | CAF-core enrichment of SPP1, CXCL9-family genes and TAM markers/programs across contexts | `results/tables/gap1_cxcl9_spp1_polarity_context_summary.csv` |
| B | Context-level SPP1-high/CXCL9-low polarity enrichment with sample points | `results/tables/gap1_cxcl9_spp1_polarity_per_sample.csv` |
| C | CAF-core polarity enrichment versus immune-decoupling index | `results/tables/gap1_cxcl9_spp1_polarity_decoupling_correlations.csv`; `results/tables/mechanism_candidate_axis_sample_summary.csv` |
| D | Sample-level SPP1-high/CXCL9-low polarity enrichment versus immune-decoupling index | `results/tables/gap1_cxcl9_spp1_polarity_per_sample.csv`; `results/tables/gap1_cxcl9_spp1_polarity_decoupling_correlations.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure8_cxcl9_spp1_polarity.pdf`.

**Claim boundary:** This is targeted gene-level spatial evidence. It does not establish macrophage-state causality, ligand-receptor signaling or outcome prediction.

## Extended Data Figure 6

**Cell-state support for CAF-myeloid core biology.**

**Take-home message:** Marker-level attribution and GSE202051 reference projection support the cellular interpretation of CAF-myeloid cores. myCAF/matrix and SPP1/TAM signals remain CAF-core enriched across contexts, whereas immune-decoupled contexts show weaker DC/APC, T/NK and B/plasma coupling. Small-reference and full-reference projections show consistent directional support for the key immune-decoupling result.

| panel | content | source |
|---|---|---|
| A | Marker-level CAF-core enrichment across tumor contexts | `results/tables/gap2_cell_state_marker_attribution_context_summary.csv` |
| B | Full-reference GSE202051 CAF-core projection enrichment across tumor contexts | `results/tables/gap2_full_reference_projection_deconvolution_context_summary.csv` |
| C | Full-reference projection association with immune-decoupling index | `results/tables/gap2_full_reference_projection_deconvolution_correlations.csv` |
| D | Small-reference versus full-reference projection stability | `results/tables/gap2_reference_projection_small_vs_full_comparison.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure6_cell_state_reference_support.pdf`.

**Claim boundary:** Marker-level attribution and reference projection support cell-state interpretation, but they are not image segmentation, immunostaining, final validated deconvolution or single-cell-resolved spatial abundance.

## Gap 2 Candidate Supplement

**Marker-level cell-state attribution of CAF-myeloid cores.**

**Take-home message:** CAF cores are consistently enriched for myCAF/matrix and SPP1/TAM marker states, while immune-decoupled contexts show weaker DC/APC and T/NK marker coupling around CAF cores, especially in lymph-node metastases.

| panel | content | source |
|---|---|---|
| A | CAF-core marker enrichment heatmap for CAF, TAM, immune, epithelial, endothelial, neural and normal/acinar marker states | `results/tables/gap2_cell_state_marker_attribution_context_summary.csv` |
| B | Spot-level correlation between CAF-myeloid score and marker-state scores | `results/tables/gap2_cell_state_marker_attribution_per_sample.csv` |
| C | Tumor-context contrast showing retained CAF/TAM core identity and attenuated immune marker coupling in lymph-node metastasis | `results/tables/gap2_cell_state_marker_attribution_context_summary.csv`; `results/tables/gap2_cell_state_marker_attribution_correlations.csv` |

**Claim boundary:** This is marker-level attribution from spatial expression data. It is not formal reference-based deconvolution, image segmentation, immunostaining or causal cell-state validation.

## Gap 2 Reference-Projection Supplement

**GSE202051 reference-projection deconvolution prototype.**

**Take-home message:** A marker-constrained projection against downloaded GSE202051 PDAC h5ad references supports recurrent SPP1/TAM projection enrichment in CAF cores and shows lymph-node attenuation of DC/APC and T/NK projection coupling.

| panel | content | source |
|---|---|---|
| A | CAF-core enrichment of GSE202051-derived marker-constrained cell-state projections across contexts | `results/tables/gap2_reference_projection_deconvolution_context_summary.csv` |
| B | Spot-level correlation between reference-projected cell-state fractions and CAF-myeloid score | `results/tables/gap2_reference_projection_deconvolution_context_summary.csv` |
| C | Primary/liver/LN metastatic-site contrast for selected reference-projected states | `results/tables/gap2_reference_projection_deconvolution_per_sample.csv` |

**Claim boundary:** This is a reference-projection prototype, not final validated deconvolution. It uses selected marker genes and clipped least-squares projection against GSE202051 signatures; it should support cell-state interpretation but not be presented as single-cell-resolved abundance.

## Gap 2 Full-Reference Validation Supplement

**GSE202051 total-reference validation of CAF-core cell-state projection.**

**Take-home message:** The full GSE202051 total h5ad reference preserves and strengthens the key Gap 2 result: immune decoupling tracks weaker DC/APC, T/NK and B/plasma CAF-core projection and stronger myCAF/matrix projection.

Source: `results/tables/gap2_full_reference_projection_deconvolution_context_summary.csv`, `results/tables/gap2_full_reference_projection_deconvolution_correlations.csv`.

**Claim boundary:** This validates the direction of the marker-constrained reference projection against a larger reference, but it remains a projection-based abundance prototype rather than final single-cell-resolved deconvolution.

## Gap 2 Small-vs-Full Reference Comparison

**Stability of GSE202051 projection results across reference size.**

**Take-home message:** Small-reference and full-reference projections agree in the direction of immune-decoupling correlations for all key states, and per-sample core-enrichment agreement is high for myCAF/matrix, B/plasma, epithelial/tumor, DC/APC and T/NK projections.

Source: `results/tables/gap2_reference_projection_small_vs_full_comparison.csv`.

**Claim boundary:** This is a robustness comparison between two projection references, not a benchmark against orthogonal histology, flow cytometry or perturbation.

## Extended Data Figure 9

**Focused CAF-core/interface candidate communication axes.**

**Take-home message:** CAF-core/interface regions nominate matrix-integrin, SPP1-CD44/integrin and TGF-beta response axes as spatially organized candidate communication programs. Matrix-integrin shows the clearest association with immune decoupling, whereas SPP1-CD44/integrin and TGF-beta/TGFBR more strongly track stromal-tumor coupling.

| panel | content | source |
|---|---|---|
| A | Ligand enrichment in CAF cores and receptor/response enrichment at tumor-stroma interfaces | `results/tables/gap3_focused_lr_interface_context_summary.csv` |
| B | Combined directional nomination scores across tumor contexts | `results/tables/gap3_focused_lr_interface_per_sample.csv` |
| C | Matrix-integrin directional nomination score versus immune-decoupling index | `results/tables/gap3_focused_lr_interface_correlations.csv`; `results/tables/mechanism_candidate_axis_sample_summary.csv` |
| D | Schematic of ligand-core, receptor-interface and response-interface scoring | same source tables plus Methods definition |

**Generated figure:** `results/figures/submission/extended_data_figure9_focused_interface_axes.pdf`.

**Claim boundary:** This nominates candidate communication biology. It is not causal ligand-receptor proof and should not be written as perturbational evidence.

## Extended Data Figure 7

**TCGA PAAD bulk RNA-seq context for nominated CAF/TAM, matrix and immune axes.**

**Take-home message:** TCGA PAAD bulk RNA-seq supports a broader stromal-myeloid/matrix expression continuum and a decoupling-like axis in which higher stromal-myeloid/matrix signal is accompanied by lower immune/APC/T-cell bulk signals.

| panel | content | source |
|---|---|---|
| A | Spearman correlation structure among nominated bulk gene-set axes | `results/tables/tcga_paad_bulk_context_axis_correlations.csv` |
| B | Bulk stromal-myeloid index versus immune projection index, colored by bulk decoupling-like index | `results/tables/tcga_paad_bulk_context_scores.csv` |
| C | Rolling gene-set score trends across samples ranked by bulk decoupling-like index | `results/tables/tcga_paad_bulk_context_scores.csv` |

**Generated figure:** `results/figures/submission/extended_data_tcga_paad_bulk_context.pdf`.

**Claim boundary:** TCGA bulk RNA-seq is non-spatial and should be used only as external biological context. It does not validate CAF-core localization, metastatic-site decoupling, clinical response or survival prediction.

## Extended Data Figure 10

**Independent GSE274557 metastatic PDAC validation.**

**Take-home message:** A Nature 2025 metastatic PDAC Visium atlas independently validates broad CAF-core-centered organization across primary, liver, lung and peritoneal PDAC contexts, while defining the boundary that this dataset cannot test lymph-node-specific immune decoupling.

| panel | content | source |
|---|---|---|
| A | GSE274557 cohort composition across primary PDAC and metastatic organ sites | `results/source_data/Source_Data_Extended_Data_Fig_10A.csv`; `results/source_data/Source_Data_Extended_Data_Fig_10_Sample_Metadata.csv` |
| B | Median observed-minus-random distance-to-CAF-core Spearman rho after 1,000 same-size random cores per sample | `results/source_data/Source_Data_Extended_Data_Fig_10B.csv`; `results/tables/gse274557_full_caf_core_gradients.csv` |
| C | Per-sample tumor-aggressive observed-minus-random deltas by tissue context | `results/tables/gse274557_full_caf_core_gradients.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure10_gse274557_external_validation.pdf`.

**Claim boundary:** This figure validates broad CAF-core organization in an independent primary/metastatic Visium resource. It does not validate lymph-node-specific immune decoupling and does not establish causal signaling.

## Extended Data Figure 11

**GSE274673 Xenium cell-resolution validation of CAF-domain immune/myeloid organization.**

**Take-home message:** Cell-resolution Xenium data support antigen-presenting and SPP1/TAM-linked CAF domains as organizers of SPP1/TAM, IFN/APC, T/NK and TGF-beta/EMT programs across treatment-naive and chemoradiotherapy-treated PDAC sections, while tumor epithelial and SPP1-tumor-like programs are not CAF-domain centered.

| panel | content | source |
|---|---|---|
| A | Observed-minus-random target-program gradients around fixed CAF-APC anchors | `results/source_data/Source_Data_Extended_Data_Fig_11A_B.csv`; `results/tables/gse274673_xenium_fixed_anchor_gradients.csv` |
| B | Observed-minus-random target-program gradients around fixed CAF-SPP1/TAM anchors | `results/source_data/Source_Data_Extended_Data_Fig_11A_B.csv`; `results/tables/gse274673_xenium_fixed_anchor_gradients.csv` |
| C | Cohort-level support for immune/myeloid target programs | `results/source_data/Source_Data_Extended_Data_Fig_11C.csv`; `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv` |
| D | Xenium sample cell counts and treatment context | `results/source_data/Source_Data_Extended_Data_Fig_11D.csv`; `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.pdf`.

**Claim boundary:** This figure addresses the mixed-spot criticism by using cell-level Xenium coordinates and expression-domain anchors. It does not establish ligand-receptor causality, perturbational mechanism, lymph-node-specific biology, clinical outcome prediction or direct CAF-to-tumor epithelial proximity.

## Extended Data Figure 12

**Distance-to-CAF-core dynamics of spatial programs.**

**Take-home message:** Explicit distance-bin dynamics show that IFN/MHC, immune-core and tumor-aggressive programs decline from CAF-core-proximal to distal regions in discovery/support and primary/liver validation contexts, whereas lymph-node metastases retain tumor-aggressive CAF-core proximity but lose immune-core/IFN proximity.

| panel | content | source |
|---|---|---|
| A | Distance-bin program dynamics in discovery/support cohorts | `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv`; `results/tables/caf_myeloid_niche_distance_bins.csv` |
| B | Distance-bin program dynamics in GSE272362 primary and liver metastasis contexts | `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv`; `results/tables/gse272362_rds_caf_myeloid_distance_bins.csv` |
| C | Distance-bin program dynamics in GSE272362 lymph-node metastasis contexts | `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv`; `results/tables/gse272362_rds_caf_myeloid_distance_bins.csv` |
| D | Core-minus-far enrichment by cohort/context and target program | `results/source_data/Source_Data_Extended_Data_Fig_12D_core_to_far.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure12_distance_to_caf_core_dynamics.pdf`.

**Claim boundary:** These panels show spatial gradients relative to CAF cores, not temporal dynamics or causality.

## Extended Data Figure 13

**Representative GSE274673 Xenium cell-domain maps.**

**Take-home message:** Cell-level maps make the Figure 4/Extended Data Figure 11 Xenium result visually inspectable, showing CAF-SPP1/TAM anchor domains alongside SPP1/TAM, IFN/APC and tumor epithelial programs in representative treatment-naive and chemoradiotherapy-treated sections.

| panel | content | source |
|---|---|---|
| A-D | Cell-level maps of CAF-SPP1/TAM anchor, SPP1/TAM, IFN/APC and tumor epithelial scores in selected treatment-naive and chemoradiotherapy-treated Xenium sections | `results/source_data/Source_Data_Extended_Data_Fig_13_selected_xenium_cell_scores.csv`; `results/tables/gse274673_xenium_fixed_anchor_cell_scores.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure13_xenium_cell_domain_maps.pdf`.

**Claim boundary:** This is representative visualization; the six-section quantitative validation remains Figure 4 and Extended Data Figure 11.

## Extended Data Figure 14

**PDAC spatial ecology atlas overview.**

**Take-home message:** An atlas-style overview summarizes cohort scale, evidence roles, CAF-core ecotype composition and representative H&E sections, matching the visual grammar of spatial atlas papers while preserving claim boundaries.

| panel | content | source |
|---|---|---|
| A | Cohort scale by analyzed sections/samples | `results/source_data/Source_Data_Extended_Data_Fig_14A_atlas_counts.csv`; `results/tables/submission_cohort_summary.csv`; `results/tables/gse274557_full_spot_scores.csv`; `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv` |
| B | Evidence role by cohort | manual evidence matrix in `scripts/50_make_reference_style_supplement_figures.py`; `results/reports/submission_claim_evidence_matrix.md` |
| C | CAF-core ecotype composition by context | `results/tables/spatial_ecotype_context_counts.csv` |
| D | Representative H&E sections across post-NACT, primary, liver-metastasis and lymph-node-metastasis contexts | `results/tables/mvp_overlay_manifest.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; raw H&E images |

**Generated figure:** `results/figures/submission/extended_data_figure14_spatial_atlas_overview.pdf`.

**Claim boundary:** The evidence-role matrix is a visual guide to the study design, not a statistical test.

## Extended Data Figure 15

**Local spatial program and CAF/tumor compartment maps.**

**Take-home message:** Representative discovery, primary, liver-metastasis and lymph-node-metastasis sections show that the CAF-core model can be inspected as tissue maps of local dominant programs and program-defined CAF/tumor/interface compartments, not only as aggregate statistics.

| panel | content | source |
|---|---|---|
| A-D | H&E, CAF-core mask, dominant local program and program-defined CAF/tumor/interface compartment maps for representative contexts | `results/source_data/Source_Data_Extended_Data_Fig_15_local_program_maps.csv`; `results/tables/mvp_spot_level_scores_with_edge_qc.csv`; `results/tables/gse272362_rds_spot_level_scores.csv`; overlay manifests and raw H&E images |

**Generated figure:** `results/figures/submission/extended_data_figure15_local_spatial_program_maps.pdf`.

**Claim boundary:** Dominant programs and compartments are computed from expression-program scores and spatial distances. They are not pathologist-annotated tissue regions or cell-type calls.

## Extended Data Figure 16

**Tumor-stroma interface compartments and candidate programs.**

**Take-home message:** GSE272362 representative maps make the candidate interface model visible by overlaying program-defined CAF-core, tumor-high and interface compartments with SPP1/TAM, TGF-beta and tumor-aggressive program gradients.

| panel | content | source |
|---|---|---|
| A-C | Program-defined compartment, SPP1/TAM, TGF-beta and tumor-aggressive maps in representative primary, liver-metastasis and lymph-node-metastasis sections | `results/source_data/Source_Data_Extended_Data_Fig_16_interface_compartments.csv`; `results/tables/gse272362_rds_spot_level_scores.csv`; `results/tables/gse272362_rds_overlay_manifest.csv`; raw H&E images |

**Generated figure:** `results/figures/submission/extended_data_figure16_interface_compartment_maps.pdf`.

**Claim boundary:** This figure visualizes expression-defined interface candidates. It does not establish histologic interface annotation, direct ligand-receptor signaling or perturbational causality.

## Extended Data Figure 17

**Representative H&E patches from program-high and program-low regions.**

**Take-home message:** Patch-level examples provide a pathology-style visual bridge for the H&E morphology analysis by showing regions with high CAF-myeloid, tumor-aggressive and IFN/MHC programs and CAF-myeloid-low controls.

| panel | content | source |
|---|---|---|
| A-D | Representative H&E patches selected from analysis-eligible spots by program score strata | `results/source_data/Source_Data_Extended_Data_Fig_17_he_patch_examples.csv`; `results/tables/mvp_he_patch_morphology_features.csv`; `metadata/dataset_manifest_curated.csv`; raw H&E images |

**Generated figure:** `results/figures/submission/extended_data_figure17_he_patch_examples.pdf`.

**Claim boundary:** Patch examples are illustrative and exploratory. They do not constitute a validated clinical-grade histology predictor.

## Extended Data Figure 18

**Program-defined neighborhoods around Xenium CAF-SPP1/TAM domains.**

**Take-home message:** GSE274673 Xenium neighborhoods around CAF-SPP1/TAM expression-domain anchors show modest but interpretable enrichment of CAF/matrix and SPP1/TAM-like program labels relative to matched random cells, complementing the distance-gradient validation in Figure 4 and Extended Data Figure 11.

| panel | content | source |
|---|---|---|
| A | Stacked program-defined neighbor composition around CAF-SPP1/TAM anchors and random cells | `results/source_data/Source_Data_Extended_Data_Fig_18_xenium_program_neighborhoods.csv`; `results/tables/gse274673_xenium_fixed_anchor_cell_scores.csv` |
| B | CAF-anchor minus random-neighbor fraction by program label and sample | same as panel A |
| C | Treatment-context summary of median neighbor fractions | same as panel A |

**Generated figure:** `results/figures/submission/extended_data_figure18_xenium_program_neighborhoods.pdf`.

**Claim boundary:** Labels are program-defined from the Xenium panel and should not be written as formal cell-type abundance, ligand-receptor interaction or treatment-response evidence.

## Extended Data Figure 19

**Representative same-size random-core null distributions.**

**Take-home message:** Representative null histograms show that observed CAF-core gradients for IFN/MHC, immune-core and tumor-aggressive programs sit outside same-size random-core distributions in discovery and validation contexts.

| panel | content | source |
|---|---|---|
| A | Observed distance-to-CAF-core Spearman rho compared with 1,000 same-size random-core rho values in representative GSE282302 and GSE272362 samples | `results/source_data/Source_Data_Extended_Data_Fig_19_random_core_null_diagnostics.csv`; `results/tables/mvp_random_core_permutation_null_rhos.csv`; `results/tables/gse272362_rds_random_core_permutation_null_rhos.csv`; corresponding sample-stat tables |

**Generated figure:** `results/figures/submission/extended_data_figure19_random_core_null_diagnostics.pdf`.

**Claim boundary:** This is a diagnostic visualization of the random-core control, not an additional biological validation cohort.

## Extended Data Figure 20

**CAF-core ecotype architecture across PDAC contexts.**

**Take-home message:** Context-to-ecotype flow and composition panels show that CAF-core ecotypes vary across post-neoadjuvant, primary and metastatic contexts, while EMT/myCAF and basal/tumor-dominant ecotypes align with higher immune-decoupling values.

| panel | content | source |
|---|---|---|
| A | Flow from cohort/context to dominant CAF-core NMF ecotype | `results/source_data/Source_Data_Extended_Data_Fig_20_ecotype_context_flow.csv`; `results/tables/spatial_ecotype_context_counts.csv` |
| B | Sample-level immune-decoupling index by dominant CAF-core ecotype | `results/tables/spatial_ecotype_sample_summary.csv` |
| C | CAF-core ecotype composition by cohort/context | `results/source_data/Source_Data_Extended_Data_Fig_20_ecotype_context_flow.csv`; `results/tables/spatial_ecotype_context_counts.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure20_ecotype_context_flow.pdf`.

**Claim boundary:** The flow summarizes sample-level dominant ecotype composition. It does not imply temporal transition, lineage progression or histologic cell-state validation.

## Extended Data Figure 21

**Mechanism triangulation prioritizes perturbation-ready interface axes.**

**Take-home message:** Multi-layer triangulation moves the mechanism discussion from a generic limitation to an actionable priority list: matrix-integrin and SPP1-CD44/integrin have the strongest combined support for future blockade or perturbation experiments.

| panel | content | source |
|---|---|---|
| A | Seven-layer mechanism evidence matrix across focused CAF-core/interface candidates | `results/source_data/Source_Data_Extended_Data_Fig_21_mechanism_triangulation.csv`; `results/tables/mechanism_triangulation_priority_matrix.csv` |
| B | Median directional core-to-interface score and sample support fraction | `results/tables/gap3_focused_lr_interface_context_summary.csv`; `results/tables/mechanism_triangulation_priority_matrix.csv` |
| C | Spearman association between directional candidate scores and immune-decoupling index | `results/tables/gap3_focused_lr_interface_correlations.csv`; `results/tables/mechanism_triangulation_priority_matrix.csv` |
| D | Triangulated perturbation-priority ranking | `results/tables/mechanism_triangulation_priority_matrix.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure21_mechanism_triangulation_priority.pdf`.

**Claim boundary:** This figure ranks mechanism candidates for experimental follow-up. It does not establish causal signaling without perturbation, lineage tracing or orthogonal protein-level validation.

## Extended Data Figure 22

**TCGA PAAD survival context for nominated stromal-myeloid axes.**

**Take-home message:** TCGA bulk RNA-seq provides exploratory clinical-context support for the matrix/stromal-myeloid side of the model, with matrix-integrin showing the strongest adverse univariable Cox association, but this is not spatial or clinical-grade validation.

| panel | content | source |
|---|---|---|
| A | Univariable Cox hazard ratios per one standard-deviation increase in nominated axes | `results/tables/tcga_paad_survival_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_22_TCGA_survival_context.csv` |
| B | Kaplan-Meier visual stress test for matrix-integrin high versus low median split | `results/tables/tcga_paad_bulk_context_scores.csv`; `data/external/TCGA_PAAD/TCGA-PAAD.survival.tsv.gz` |
| C | Kaplan-Meier visual stress test for bulk decoupling-like high versus low median split | `results/tables/tcga_paad_bulk_context_scores.csv`; `data/external/TCGA_PAAD/TCGA-PAAD.survival.tsv.gz` |

**Generated figure:** `results/figures/submission/extended_data_figure22_tcga_survival_context.pdf`.

**Claim boundary:** This is a non-spatial TCGA bulk clinical-context analysis. It does not establish spatial localization, treatment response, independent clinical validation or a prognostic model.

## Extended Data Figure 23

**TLS-maturity stress test of CAF-myeloid cores.**

**Take-home message:** CAF-core regions can contain immune-hub features, but only a minority of samples pass a stringent multi-compartment TLS-maturity gate; the central framing should remain CAF-myeloid spatial niche rather than mature TLS.

| panel | content | source |
|---|---|---|
| A | CAF-core enrichment of TLS chemokine, lymphoid/plasma, FDC/GC-like, immune-maturity and CAF-myeloid modules across contexts | `results/tables/tls_maturity_stress_test_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_23_TLS_maturity_stress_test.csv` |
| B | Fraction and count of samples passing the stringent TLS-maturity gate | `results/tables/tls_maturity_stress_test_context_summary.csv` |
| C | Sample-level TLS three-compartment score | `results/tables/tls_maturity_stress_test_per_sample.csv` |
| D | CAF-myeloid minus immune-maturity enrichment balance | `results/tables/tls_maturity_stress_test_context_summary.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure23_tls_maturity_stress_test.pdf`.

**Claim boundary:** This is an expression-derived TLS-maturity stress test. It does not replace histologic TLS annotation, FDC network validation or germinal-center maturation assays.

## Extended Data Figure 24

**Reviewer-risk resolution module for the CAF-myeloid spatial niche model.**

**Take-home message:** The manuscript's vulnerable claims are now either supported, ranked as perturbation candidates or explicitly bounded: spatial organization and metastatic immune decoupling are supported, mechanism is prioritized rather than proven causal, TCGA is clinical context only and mature TLS is not the central framing.

| panel | content | source |
|---|---|---|
| A | Evidence spine from spatial specificity to claim-boundary testing | `results/reports/extended_data_figure24_review_risk_resolution_notes.md`; `results/source_data/Source_Data_Extended_Data_Fig_24_review_risk_resolution.csv` |
| B | Mechanism evidence matrix across core, interface, directionality, decoupling, targeted genes, Xenium and TCGA support | `results/tables/mechanism_triangulation_priority_matrix.csv`; `results/source_data/Source_Data_Extended_Data_Fig_21_mechanism_triangulation.csv`; `results/source_data/Source_Data_Extended_Data_Fig_24_review_risk_resolution.csv` |
| C | Prioritized interface axes versus immune decoupling | `results/tables/mechanism_triangulation_priority_matrix.csv` |
| D | TCGA PAAD bulk Cox survival-context forest plot | `results/tables/tcga_paad_survival_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_22_TCGA_survival_context.csv` |
| E | Stringent TLS-maturity gate fractions by context | `results/tables/tls_maturity_stress_test_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_23_TLS_maturity_stress_test.csv` |
| F | Claim-scope map after stress testing | `results/reports/top_journal_integrated_story_audit_2026_06_28.md`; `results/reports/final_story_package.md` |

**Generated figure:** `results/figures/submission/extended_data_figure24_review_risk_resolution_nc_style.pdf`.

**Claim boundary:** This is a synthesis and audit-readiness figure. It consolidates existing analyses and does not add causal, prognostic or mature-TLS claims.

## Extended Data Figure 25

**Spatial robustness module for the CAF-myeloid niche model.**

**Take-home message:** The organizing-core signal is robust to random-core controls, CAF-core threshold choices, distance-gradient summaries and independent Visium/Xenium validation layers.

| panel | content | source |
|---|---|---|
| A | Cross-cohort random-core specificity heatmap with median deltas and supported sample counts | `results/source_data/Source_Data_Extended_Random_Core_MVP.csv`; `results/source_data/Source_Data_Extended_Random_Core_GSE272362.csv`; `results/source_data/Source_Data_Fig_6A.csv`; `results/source_data/Source_Data_Extended_Data_Fig_25_spatial_robustness_module.csv` |
| B | CAF-core top 15%, top 10% and top 5% threshold sensitivity | `results/source_data/Source_Data_Extended_Threshold_Sensitivity.csv`; `results/tables/caf_core_threshold_sensitivity_summary.csv` |
| C | Discovery/support distance-gradient dynamics from core to far regions | `results/source_data/Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv`; `results/tables/caf_myeloid_niche_distance_bins.csv` |
| D | Context-level core-to-far effect sizes | `results/source_data/Source_Data_Extended_Data_Fig_12D_core_to_far.csv` |
| E | Independent GSE274557 Visium summary across tissue contexts | `results/source_data/Source_Data_Fig_4B_GSE274557.csv`; `results/tables/gse274557_full_caf_core_context_summary.csv` |
| F | GSE274673 Xenium cell-resolution CAF-domain support | `results/source_data/Source_Data_Fig_4C_GSE274673.csv`; `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv` |
| G | Cross-context support fraction for main programs | `results/source_data/Source_Data_ED_Fig_5D_cross_context_support_fraction.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure25_spatial_robustness_module_nc_style.pdf`.

**Claim boundary:** This module supports spatial robustness and external consistency. It does not establish perturbational causality, clinical prediction or histology-annotated compartment identity.

## Extended Data Figure 26

**Cell-state, reference-projection and Xenium support module.**

**Take-home message:** Marker-state enrichment, GSE202051 reference projection and GSE274673 Xenium cell-resolution analyses converge on a cellular interpretation of the CAF-myeloid spatial niche: myCAF/matrix and SPP1/TAM states are CAF-core aligned, immune-state coupling is attenuated in immune-decoupled contexts and Xenium CAF domains center immune/myeloid programs rather than tumor epithelial programs.

| panel | content | source |
|---|---|---|
| A | Marker-state CAF-core enrichment across tumor contexts | `results/tables/gap2_cell_state_marker_attribution_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_26_cell_state_reference_xenium_module.csv` |
| B | Marker-state spot-level correlation with CAF-myeloid score | `results/tables/gap2_cell_state_marker_attribution_context_summary.csv` |
| C | Marker-level and full-reference associations with immune-decoupling index | `results/tables/gap2_cell_state_marker_attribution_correlations.csv`; `results/tables/gap2_full_reference_projection_deconvolution_correlations.csv` |
| D | Full GSE202051 reference-projection CAF-core enrichment across contexts | `results/tables/gap2_full_reference_projection_deconvolution_context_summary.csv` |
| E | Small-reference versus full-reference per-sample projection stability | `results/tables/gap2_reference_projection_small_vs_full_comparison.csv` |
| F | GSE274673 Xenium CAF-domain target-program centering around CAF-APC and CAF-SPP1/TAM anchors | `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv`; `results/source_data/Source_Data_Fig_4C_GSE274673.csv` |
| G | Xenium sample scale and signature gene coverage | `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv`; `results/tables/gse274673_xenium_fixed_anchor_signature_coverage.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure26_cell_state_reference_xenium_module_nc_style.pdf`.

**Claim boundary:** This module supports cell-state interpretation but is not formal image segmentation, immunostaining, final validated deconvolution or causal cell-cell interaction evidence.

## Extended Data Figure 27

**Metastatic-site immune-decoupling module.**

**Take-home message:** GSE272362 supports a metastatic-site contrast in which primary tumors and liver metastases preserve CAF-core-centered immune and tumor-aggressive organization, whereas lymph-node metastases retain stromal-tumor coupling but selectively weaken immune/IFN coupling to the CAF core.

| panel | content | source |
|---|---|---|
| A | GSE272362 primary, liver-metastasis and lymph-node-metastasis sample and spot scale | `results/tables/gse272362_rds_sample_specimen_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_27_metastatic_immune_decoupling_module.csv` |
| B | Random-core support for IFN/MHC, immune-core and tumor-aggressive programs by site | `results/source_data/Source_Data_Fig_2.csv`; `results/tables/gse272362_rds_random_core_permutation_summary.csv` |
| C | CAF-core subprogram distance gradients across sites | `results/tables/gse272362_caf_core_subprogram_gradient_summary.csv`; `results/source_data/Source_Data_Fig_3A.csv` |
| D | Stromal-tumor coupling, immune-core coupling and immune-decoupling index | `results/tables/immune_decoupling_context_summary.csv`; `results/source_data/Source_Data_Fig_5C.csv` |
| E | Patient-matched primary-to-metastasis median deltas | `results/tables/gse272362_patient_matched_site_delta_summary.csv`; `results/source_data/Source_Data_Fig_3B.csv` |
| F | Matched-patient delta distributions for immune-core, IFN/MHC and tumor-aggressive programs | `results/tables/gse272362_patient_matched_site_deltas.csv`; `results/source_data/Source_Data_Fig_3C.csv` |
| G | Interpretation boundary for the lymph-node result | `results/reports/extended_data_figure27_metastatic_immune_decoupling_module_notes.md` |

**Generated figure:** `results/figures/submission/extended_data_figure27_metastatic_immune_decoupling_module_nc_style.pdf`.

**Claim boundary:** The lymph-node result is a spatial decoupling contrast from five lymph-node metastasis samples. It supports a metastatic immune-remodeling lead, not a clinical subtype or causal mechanism.

## Extended Data Figure 28

**Strict NNLS reference-deconvolution sensitivity analysis.**

**Take-home message:** Per-spot non-negative least-squares deconvolution using full GSE202051-derived reference signatures preserves the key cell-state direction: immune decoupling is associated with weaker T/NK, DC/APC and B/plasma CAF-core enrichment and stronger myCAF/matrix enrichment.

| panel | content | source |
|---|---|---|
| A | NNLS-derived CAF-core enrichment of key cell-state fractions across contexts | `results/tables/strict_nnls_reference_deconvolution_context_summary.csv`; `results/source_data/Source_Data_Extended_Data_Fig_28_strict_nnls_deconvolution.csv` |
| B | Median spot-level Spearman correlation between NNLS state fraction and CAF-myeloid score | `results/tables/strict_nnls_reference_deconvolution_context_summary.csv` |
| C | Association between NNLS CAF-core cell-state enrichment and immune-decoupling index | `results/tables/strict_nnls_reference_deconvolution_correlations.csv` |
| D | Sample-level agreement between strict NNLS enrichment and prior full-reference projection | `results/tables/strict_nnls_vs_projection_comparison.csv`; `results/tables/gap2_full_reference_projection_deconvolution_per_sample.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure28_strict_nnls_deconvolution_sensitivity.pdf`.

**Claim boundary:** This is a stricter reference-dependent computational deconvolution sensitivity analysis. It is not immunostaining, image segmentation or single-cell-resolved ground truth.

## Extended Data Figure 29

**Mechanism gene/interface module for perturbation-ready axes.**

**Take-home message:** Matrix-integrin and SPP1-CD44/integrin have the broadest triangulated support across CAF-core ligand signal, interface response, directional structure, targeted genes, Xenium support and TCGA context; TGF-beta/TGFBR remains a strong invasive-interface follow-up axis.

| panel | content | source |
|---|---|---|
| A | Triangulated evidence layers for candidate axes | `results/tables/mechanism_triangulation_priority_matrix.csv`; `results/source_data/Source_Data_Extended_Data_Fig_29_mechanism_gene_interface_module.csv` |
| B | Targeted-gene CAF-core enrichment across contexts | `results/tables/targeted_gene_axis_validation_summary.csv` |
| C | Targeted-gene tumor-stroma-interface enrichment across contexts | `results/tables/targeted_gene_axis_validation_summary.csv` |
| D | Focused ligand-core, receptor-interface, response-interface and directional metrics | `results/tables/gap3_focused_lr_interface_context_summary.csv` |
| E | Directional-score associations with immune decoupling and stromal-tumor coupling | `results/tables/gap3_focused_lr_interface_correlations.csv` |
| F | Perturbation-priority ranking | `results/tables/mechanism_triangulation_priority_matrix.csv`; `results/source_data/Source_Data_Extended_Data_Fig_21_mechanism_triangulation.csv` |
| G | Compact candidate mechanism model for follow-up | `results/reports/extended_data_figure29_mechanism_gene_interface_module_notes.md` |

**Generated figure:** `results/figures/submission/extended_data_figure29_mechanism_gene_interface_module_nc_style.pdf`.

**Claim boundary:** This module converts observational spatial evidence into a ranked perturbation agenda. It does not establish causal ligand-receptor signaling.

## Extended Data Figure 30

**Alternative biological-anchor specificity of CAF-myeloid cores.**

**Take-home message:** CAF-myeloid cores organize multiple cross-program spatial gradients beyond same-size random anchors and beyond several alternative biological anchors. The analysis also shows useful specificity: immune-core, tumor-aggressive, SPP1/TAM and tumor-epithelial anchors most strongly center their own target programs, whereas CAF-myeloid cores preferentially coordinate TGFb/EMT, SPP1/TAM, IFN/MHC, immune-core and tumor-aggressive programs but not tumor epithelial signal.

| panel | content | source |
|---|---|---|
| A | Anchor-target observed-minus-random specificity matrix across CAF-myeloid, tumor-aggressive, tumor-epithelial, immune-core, panCAF and SPP1/TAM anchors | `results/source_data/Source_Data_Extended_Data_Fig_30_alternative_anchor_specificity.csv`; `results/tables/alternative_biological_anchor_specificity_summary.csv` |
| B | CAF-myeloid core versus strongest alternative anchor for IFN/MHC, immune-core and tumor-aggressive programs | `results/tables/alternative_biological_anchor_specificity_summary.csv` |
| C | CAF-core specificity by tissue/specimen context | `results/tables/alternative_biological_anchor_specificity_per_sample.csv` |
| D,E | Distance-bin profiles from selected biological anchors | `results/source_data/Source_Data_Extended_Data_Fig_30_distance_curves.csv`; `results/tables/alternative_biological_anchor_distance_curves.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure30_alternative_biological_anchor_specificity.pdf`.

**Claim boundary:** This is a biological-anchor specificity test from observational spatial transcriptomics. It strengthens the organizing-core interpretation but does not prove CAF-to-immune or CAF-to-tumor causality.

## Extended Data Figure 31

**Xenium cell-neighborhood networks resolve CAF-domain spatial organization.**

**Take-home message:** GSE274673 cell-resolution Xenium data show that CAF-domain neighborhoods are enriched for CAF/matrix, SPP1/TAM, IFN/APC and T/NK states and depleted for tumor epithelial states relative to matched random anchor neighborhoods. A k-nearest-neighbor cell-state graph highlights IFN/APC-T/NK, SPP1/TAM-immune and CAF/matrix-TGFb/EMT spatial adjacency layers.

| panel | content | source |
|---|---|---|
| A | Representative cell-level neighborhood map with dominant program states | `results/tables/gse274673_xenium_fixed_anchor_cell_scores.csv` |
| B | Enriched cell-state adjacency network | `results/source_data/Source_Data_Extended_Data_Fig_31_xenium_neighborhood_network.csv`; `results/tables/gse274673_xenium_cell_state_adjacency_summary.csv` |
| C | Observed/expected cell-state adjacency matrix | `results/tables/gse274673_xenium_cell_state_adjacency_summary.csv` |
| D | Cell-state composition by treatment context | `results/source_data/Source_Data_Extended_Data_Fig_31_cell_state_composition.csv`; `results/tables/gse274673_xenium_cell_state_composition.csv` |
| E,F | CAF-SPP1/TAM and CAF-APC neighborhood composition versus matched random anchors | `results/source_data/Source_Data_Extended_Data_Fig_31_anchor_neighborhood_summary.csv`; `results/tables/gse274673_xenium_anchor_neighborhood_summary.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure31_xenium_cell_neighborhood_network.pdf`.

**Claim boundary:** This is cell-state adjacency enrichment in observational Xenium data. It supports a spatial neighborhood grammar but not direct ligand-receptor signaling or perturbational causality.

## Extended Data Figure 32

**Core-to-interface transition trajectories deepen the CAF-domain model.**

**Take-home message:** A pseudo-spatial coordinate from CAF-myeloid cores to tumor-high regions separates CAF-core-proximal stromal/myeloid and immune programs from tumor-high epithelial programs. Primary tumors and liver metastases show strong CAF-core-proximal SPP1/TAM, TGFb/EMT and IFN/APC trajectories, whereas lymph-node metastases show a shifted tumor-aggressive peak and weaker CAF/myCAF core-to-tumor gradient.

| panel | content | source |
|---|---|---|
| A | Core-to-interface coordinate schematic | `results/reports/core_to_interface_transition_model_report.md` |
| B | Representative primary, liver-metastasis and lymph-node-metastasis spatial maps colored by core-to-tumor coordinate | `results/source_data/Source_Data_Extended_Data_Fig_32_selected_spot_maps.csv`; `results/tables/core_to_interface_transition_selected_spot_maps.csv` |
| C | Program trajectories for SPP1/TAM, TGFb/EMT, IFN/APC, T/NK, tumor-aggressive and tumor-epithelial programs | `results/source_data/Source_Data_Extended_Data_Fig_32_core_to_interface_transition.csv`; `results/tables/core_to_interface_transition_context_summary.csv` |
| D | Median peak coordinate by program and metastatic context | `results/source_data/Source_Data_Extended_Data_Fig_32_sample_summary.csv`; `results/tables/core_to_interface_transition_sample_summary.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure32_core_to_interface_transition_model.pdf`.

**Claim boundary:** This is a spatial trajectory/gradient model, not temporal ordering and not causal progression. It supports spatial process structure along CAF-core to tumor-high axes.

## Extended Data Figure 33

**CAF-core geometry reveals tissue architecture behind spatial niche organization.**

**Take-home message:** CAF-myeloid cores can be quantified as tissue structures rather than only high-score spots. Less fragmented CAF cores and larger dominant connected components associate with stronger stromal-tumor coupling, while tumor-contact geometry tracks SPP1/TAM and TGFb/EMT transition behavior. Interface-core fraction does not directly explain immune decoupling, which helps bound the geometry claim.

| panel | content | source |
|---|---|---|
| A | CAF-core architecture metric schematic | `results/reports/caf_core_geometry_tissue_architecture_report.md` |
| B | Representative primary, liver-metastasis and lymph-node-metastasis maps showing CAF-core, tumor-high and core-tumor interface spots | `results/source_data/Source_Data_Extended_Data_Fig_33_selected_spot_maps.csv`; `results/tables/caf_core_geometry_selected_spot_maps.csv` |
| C-F | Largest connected component fraction, fragmentation, core-tumor interface fraction and core-to-tumor distance across contexts | `results/source_data/Source_Data_Extended_Data_Fig_33_caf_core_geometry.csv`; `results/tables/caf_core_geometry_metrics_per_sample.csv` |
| G | Interface geometry versus immune-decoupling index | `results/tables/caf_core_geometry_metrics_per_sample.csv`; `results/tables/mechanism_candidate_axis_sample_summary.csv` |
| H | Geometry-biological readout association matrix | `results/source_data/Source_Data_Extended_Data_Fig_33_biological_correlations.csv`; `results/tables/caf_core_geometry_biological_correlations.csv` |

**Generated figure:** `results/figures/submission/extended_data_figure33_caf_core_geometry_tissue_architecture.pdf`.

**Claim boundary:** This is expression-defined tissue architecture, not pathologist-annotated histologic compartment segmentation. It supports a spatial-geometry layer for the CAF-myeloid niche model but does not establish causal remodeling.

## Regenerated Supplementary Module 1

**Spatial specificity and robustness of CAF-myeloid cores.**

**Take-home message:** The CAF-core organizing signal is robust to same-size random-core controls, CAF-core threshold choices, distance-gradient summaries and independent Visium/Xenium validation layers.

| panel | content | source |
|---|---|---|
| A-G | Cross-cohort random-core specificity, threshold sensitivity, distance-gradient summaries, independent Visium/Xenium support and cross-context support fraction | `results/source_data/Source_Data_supplementary_module1_spatial_specificity_robustness.csv`; `results/reports/supplementary_module1_spatial_specificity_robustness_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module1_spatial_specificity_robustness.pdf`.

**Claim boundary:** This regenerated module supports spatial specificity and robustness. It remains observational and program-defined.

## Regenerated Supplementary Module 2

**Metastatic-site remodeling and lymph-node immune decoupling.**

**Take-home message:** GSE272362 lymph-node metastases retain stromal-tumor coupling while weakening immune/IFN coupling to CAF cores.

| panel | content | source |
|---|---|---|
| A-G | Site scale, random-core support, subprogram distance gradients, immune-decoupling metrics, matched-site deltas and lymph-node interpretation boundary | `results/source_data/Source_Data_supplementary_module2_metastatic_immune_decoupling.csv`; `results/reports/supplementary_module2_metastatic_immune_decoupling_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module2_metastatic_immune_decoupling.pdf`.

**Claim boundary:** This regenerated module supports a metastatic-site immune-remodeling lead, not a definitive clinical subtype or causal mechanism.

## Regenerated Supplementary Module 3

**Cell-state interpretation and multi-resolution validation.**

**Take-home message:** Marker, reference-projection, strict NNLS and Xenium layers support a cellular interpretation of CAF-myeloid cores while preserving boundaries against overclaiming deconvolution or ground truth.

| panel | content | source |
|---|---|---|
| A-H | Marker-state enrichment, marker/full-reference/strict-NNLS associations, reference stability, Xenium support, sample/signature coverage and cell-state claim ladder | `results/source_data/Source_Data_supplementary_module3_cell_state_multiresolution_validation.csv`; `results/reports/supplementary_module3_cell_state_multiresolution_validation_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module3_cell_state_multiresolution_validation.pdf`.

**Claim boundary:** This regenerated module supports cell-state interpretation but is not immunostaining, image segmentation, single-cell-resolved abundance or causal interaction evidence.

## Regenerated Supplementary Module 4

**Mechanism, interface biology and perturbation-priority axes.**

**Take-home message:** Matrix-integrin and SPP1-CD44/integrin carry the broadest triangulated support and form the leading perturbation-ready candidates.

| panel | content | source |
|---|---|---|
| A-H | Mechanism evidence matrix, perturbation-priority ranking, targeted-gene core/interface enrichment, focused ligand-response metrics, directional associations, external-context anchors and claim ladder | `results/source_data/Source_Data_supplementary_module4_mechanism_interface_priority.csv`; `results/reports/supplementary_module4_mechanism_interface_priority_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module4_mechanism_interface_priority.pdf`.

**Claim boundary:** This regenerated module ranks perturbation candidates and does not establish causal ligand-receptor signaling.

## Regenerated Supplementary Module 5

**Pathology bridge, TCGA context and TLS/claim boundaries.**

**Take-home message:** H&E, TCGA and TLS analyses provide bounded translational context and explicitly delimit mature TLS, clinical outcome and causal-signaling claims.

| panel | content | source |
|---|---|---|
| A-G | H&E model performance, H&E feature directions, TCGA bulk correlation, exploratory TCGA survival context, TLS compatibility, immune-compartment stress testing and claim-boundary map | `results/source_data/Source_Data_supplementary_module5_pathology_tcga_tls_boundaries.csv`; `results/reports/supplementary_module5_pathology_tcga_tls_boundaries_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module5_pathology_tcga_tls_boundaries.pdf`.

**Claim boundary:** This regenerated module strengthens translational plausibility and claim control; it does not establish mature TLS biology, clinical-grade pathology prediction, prognosis, therapy response, spatial TCGA localization or causal signaling.

## Regenerated Supplementary Module 6

**Deep spatial architecture and mechanism-deepening evidence.**

**Take-home message:** ED30-ED33 converge on a deeper spatial-architecture model: CAF-myeloid cores are biologically specific, supported at Xenium cell-neighborhood resolution, organized along a core-to-interface axis and modulated by CAF-core geometry.

| panel | content | source |
|---|---|---|
| A-H | evidence chain, alternative-anchor specificity, Xenium CAF-domain neighborhoods, cell-state adjacency network, core-to-interface transition curves, CAF-core geometry, geometry-biological associations and claim boundary | `results/source_data/Source_Data_supplementary_module6_spatial_architecture_mechanism_deepening.csv`; `results/reports/supplementary_module6_spatial_architecture_mechanism_deepening_notes.md` |

**Generated figure:** `results/figures/submission/supplementary_module6_spatial_architecture_mechanism_deepening.pdf`.

**Claim boundary:** This regenerated module strengthens spatial architecture and perturbation-priority logic; it should not be written as direct causal proof without perturbation experiments.

## Generated Files

- `results/figures/submission/figure1_submission_spatial_specificity.pdf`
- `results/figures/submission/figure1_submission_spatial_specificity.svg`
- `results/figures/submission/figure1_submission_spatial_specificity.png`
- `results/figures/submission/figure1_submission_spatial_specificity_nc_style.pdf`
- `results/figures/submission/figure1_submission_spatial_specificity_nc_style.svg`
- `results/figures/submission/figure1_submission_spatial_specificity_nc_style.png`
- `results/figures/submission/figure1_supplement_submission_post_nact_spatial_example.pdf`
- `results/figures/submission/figure1_supplement_submission_post_nact_spatial_example.svg`
- `results/figures/submission/figure1_supplement_submission_post_nact_spatial_example.png`
- `results/figures/submission/figure2_submission_metastatic_decoupling.pdf`
- `results/figures/submission/figure2_submission_metastatic_decoupling.svg`
- `results/figures/submission/figure2_submission_metastatic_decoupling.png`
- `results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.pdf`
- `results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.svg`
- `results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.png`
- `results/figures/submission/figure2_supplement_submission_spatial_examples.pdf`
- `results/figures/submission/figure2_supplement_submission_spatial_examples.svg`
- `results/figures/submission/figure2_supplement_submission_spatial_examples.png`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes.pdf`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes.svg`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes.png`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.pdf`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.svg`
- `results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.png`
- `results/figures/submission/figure3_supplement_targeted_gene_axis_validation.pdf`
- `results/figures/submission/figure3_supplement_targeted_gene_axis_validation.svg`
- `results/figures/submission/figure3_supplement_targeted_gene_axis_validation.png`
- `results/figures/submission/figure4_submission_multiresolution_validation_nc_style.pdf`
- `results/figures/submission/figure4_submission_multiresolution_validation_nc_style.svg`
- `results/figures/submission/figure4_submission_multiresolution_validation_nc_style.png`
- `results/figures/submission/supplementary_module6_spatial_architecture_mechanism_deepening.pdf`
- `results/figures/submission/supplementary_module6_spatial_architecture_mechanism_deepening.svg`
- `results/figures/submission/supplementary_module6_spatial_architecture_mechanism_deepening.png`
- `results/figures/submission/extended_data_figure4_he_morphology_bridge.pdf`
- `results/figures/submission/extended_data_figure4_he_morphology_bridge.svg`
- `results/figures/submission/extended_data_figure4_he_morphology_bridge.png`
- `results/figures/submission/extended_data_figure5_external_anchor_robustness.pdf`
- `results/figures/submission/extended_data_figure5_external_anchor_robustness.svg`
- `results/figures/submission/extended_data_figure5_external_anchor_robustness.png`
- `results/figures/submission/extended_data_figure6_cell_state_reference_support.pdf`
- `results/figures/submission/extended_data_figure6_cell_state_reference_support.svg`
- `results/figures/submission/extended_data_figure6_cell_state_reference_support.png`
- `results/figures/submission/extended_data_figure8_cxcl9_spp1_polarity.pdf`
- `results/figures/submission/extended_data_figure8_cxcl9_spp1_polarity.svg`
- `results/figures/submission/extended_data_figure8_cxcl9_spp1_polarity.png`
- `results/figures/submission/extended_data_figure9_focused_interface_axes.pdf`
- `results/figures/submission/extended_data_figure9_focused_interface_axes.svg`
- `results/figures/submission/extended_data_figure9_focused_interface_axes.png`
- `results/figures/submission/extended_data_figure10_gse274557_external_validation.pdf`
- `results/figures/submission/extended_data_figure10_gse274557_external_validation.svg`
- `results/figures/submission/extended_data_figure10_gse274557_external_validation.png`
- `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.pdf`
- `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.svg`
- `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.png`
- `results/figures/submission/extended_data_figure12_distance_to_caf_core_dynamics.pdf`
- `results/figures/submission/extended_data_figure12_distance_to_caf_core_dynamics.svg`
- `results/figures/submission/extended_data_figure12_distance_to_caf_core_dynamics.png`
- `results/figures/submission/extended_data_figure13_xenium_cell_domain_maps.pdf`
- `results/figures/submission/extended_data_figure13_xenium_cell_domain_maps.svg`
- `results/figures/submission/extended_data_figure13_xenium_cell_domain_maps.png`
- `results/figures/submission/extended_data_figure14_spatial_atlas_overview.pdf`
- `results/figures/submission/extended_data_figure14_spatial_atlas_overview.svg`
- `results/figures/submission/extended_data_figure14_spatial_atlas_overview.png`
- `results/figures/submission/extended_data_figure15_local_spatial_program_maps.pdf`
- `results/figures/submission/extended_data_figure15_local_spatial_program_maps.svg`
- `results/figures/submission/extended_data_figure15_local_spatial_program_maps.png`
- `results/figures/submission/extended_data_figure16_interface_compartment_maps.pdf`
- `results/figures/submission/extended_data_figure16_interface_compartment_maps.svg`
- `results/figures/submission/extended_data_figure16_interface_compartment_maps.png`
- `results/figures/submission/extended_data_figure17_he_patch_examples.pdf`
- `results/figures/submission/extended_data_figure17_he_patch_examples.svg`
- `results/figures/submission/extended_data_figure17_he_patch_examples.png`
- `results/figures/submission/extended_data_figure18_xenium_program_neighborhoods.pdf`
- `results/figures/submission/extended_data_figure18_xenium_program_neighborhoods.svg`
- `results/figures/submission/extended_data_figure18_xenium_program_neighborhoods.png`
- `results/figures/submission/extended_data_figure19_random_core_null_diagnostics.pdf`
- `results/figures/submission/extended_data_figure19_random_core_null_diagnostics.svg`
- `results/figures/submission/extended_data_figure19_random_core_null_diagnostics.png`
- `results/figures/submission/extended_data_figure20_ecotype_context_flow.pdf`
- `results/figures/submission/extended_data_figure20_ecotype_context_flow.svg`
- `results/figures/submission/extended_data_figure20_ecotype_context_flow.png`
- `results/figures/submission/extended_data_figure21_mechanism_triangulation_priority.pdf`
- `results/figures/submission/extended_data_figure21_mechanism_triangulation_priority.svg`
- `results/figures/submission/extended_data_figure21_mechanism_triangulation_priority.png`
- `results/figures/submission/extended_data_figure22_tcga_survival_context.pdf`
- `results/figures/submission/extended_data_figure22_tcga_survival_context.svg`
- `results/figures/submission/extended_data_figure22_tcga_survival_context.png`
- `results/figures/submission/extended_data_figure23_tls_maturity_stress_test.pdf`
- `results/figures/submission/extended_data_figure23_tls_maturity_stress_test.svg`
- `results/figures/submission/extended_data_figure23_tls_maturity_stress_test.png`
- `results/figures/submission/extended_data_figure24_review_risk_resolution_nc_style.pdf`
- `results/figures/submission/extended_data_figure24_review_risk_resolution_nc_style.svg`
- `results/figures/submission/extended_data_figure24_review_risk_resolution_nc_style.png`
- `results/figures/submission/extended_data_figure25_spatial_robustness_module_nc_style.pdf`
- `results/figures/submission/extended_data_figure25_spatial_robustness_module_nc_style.svg`
- `results/figures/submission/extended_data_figure25_spatial_robustness_module_nc_style.png`
- `results/figures/submission/extended_data_figure26_cell_state_reference_xenium_module_nc_style.pdf`
- `results/figures/submission/extended_data_figure26_cell_state_reference_xenium_module_nc_style.svg`
- `results/figures/submission/extended_data_figure26_cell_state_reference_xenium_module_nc_style.png`
- `results/figures/submission/extended_data_figure27_metastatic_immune_decoupling_module_nc_style.pdf`
- `results/figures/submission/extended_data_figure27_metastatic_immune_decoupling_module_nc_style.svg`
- `results/figures/submission/extended_data_figure27_metastatic_immune_decoupling_module_nc_style.png`
- `results/figures/submission/extended_data_figure28_strict_nnls_deconvolution_sensitivity.pdf`
- `results/figures/submission/extended_data_figure28_strict_nnls_deconvolution_sensitivity.svg`
- `results/figures/submission/extended_data_figure28_strict_nnls_deconvolution_sensitivity.png`
- `results/figures/submission/extended_data_figure29_mechanism_gene_interface_module_nc_style.pdf`
- `results/figures/submission/extended_data_figure29_mechanism_gene_interface_module_nc_style.svg`
- `results/figures/submission/extended_data_figure29_mechanism_gene_interface_module_nc_style.png`
- `results/figures/submission/extended_data_figure30_alternative_biological_anchor_specificity.pdf`
- `results/figures/submission/extended_data_figure30_alternative_biological_anchor_specificity.svg`
- `results/figures/submission/extended_data_figure30_alternative_biological_anchor_specificity.png`
- `results/figures/submission/extended_data_figure31_xenium_cell_neighborhood_network.pdf`
- `results/figures/submission/extended_data_figure31_xenium_cell_neighborhood_network.svg`
- `results/figures/submission/extended_data_figure31_xenium_cell_neighborhood_network.png`
- `results/figures/submission/extended_data_figure32_core_to_interface_transition_model.pdf`
- `results/figures/submission/extended_data_figure32_core_to_interface_transition_model.svg`
- `results/figures/submission/extended_data_figure32_core_to_interface_transition_model.png`
- `results/figures/submission/extended_data_figure33_caf_core_geometry_tissue_architecture.pdf`
- `results/figures/submission/extended_data_figure33_caf_core_geometry_tissue_architecture.svg`
- `results/figures/submission/extended_data_figure33_caf_core_geometry_tissue_architecture.png`
- `results/figures/submission/extended_data_gap1_cxcl9_spp1_polarity.pdf`
- `results/figures/submission/extended_data_gap1_cxcl9_spp1_polarity.svg`
- `results/figures/submission/extended_data_gap1_cxcl9_spp1_polarity.png`
- `results/figures/submission/extended_data_gap2_cell_state_attribution.pdf`
- `results/figures/submission/extended_data_gap2_cell_state_attribution.svg`
- `results/figures/submission/extended_data_gap2_cell_state_attribution.png`
- `results/figures/submission/extended_data_gap2_reference_projection_deconvolution.pdf`
- `results/figures/submission/extended_data_gap2_reference_projection_deconvolution.svg`
- `results/figures/submission/extended_data_gap2_reference_projection_deconvolution.png`
- `results/figures/submission/extended_data_gap2_full_reference_projection_deconvolution.pdf`
- `results/figures/submission/extended_data_gap2_full_reference_projection_deconvolution.svg`
- `results/figures/submission/extended_data_gap2_full_reference_projection_deconvolution.png`
- `results/figures/submission/extended_data_gap2_reference_projection_small_vs_full_comparison.pdf`
- `results/figures/submission/extended_data_gap2_reference_projection_small_vs_full_comparison.svg`
- `results/figures/submission/extended_data_gap2_reference_projection_small_vs_full_comparison.png`
- `results/figures/submission/extended_data_gap3_focused_lr_interface.pdf`
- `results/figures/submission/extended_data_gap3_focused_lr_interface.svg`
- `results/figures/submission/extended_data_gap3_focused_lr_interface.png`
- `results/figures/submission/extended_data_tcga_paad_bulk_context.pdf`
- `results/figures/submission/extended_data_tcga_paad_bulk_context.svg`
- `results/figures/submission/extended_data_tcga_paad_bulk_context.png`
