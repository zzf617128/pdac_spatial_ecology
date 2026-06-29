# MVP Decision Report

Last updated UTC: 2026-06-24T07:30:00+00:00

## Current Data Status

- Directly scored Visium samples: 113
- GSE282302 scored samples: 108
- GSE274103 scored samples: 5
- GSE272362 Zenodo `PDAC_Updated.rds` scored samples: 30; spots: 91,496.
- GSE272362 tissue-site groups from RDS `Origin`: 10 primary tumors, 12 liver metastases, 5 lymph node metastases, 3 normal pancreas.
- Random-core permutation controls completed with 1,000 random same-size cores per sample for both main MVP cohorts and GSE272362.
- Six GSE272362 representative H&E overlays completed with exact RDS-to-GEO coordinate matching.
- CAF-core threshold sensitivity completed across top 15%, top 10%, and top 5% definitions.
- Manuscript-style Figure 1/2 drafts and source tables generated.
- GSE272362 specimen-group provenance audited against GEO, the linked Nature Genetics article, Zenodo, and local RDS metadata.

## Preliminary Recommendation

Current top MVP direction: **CAF-myeloid inflammatory stromal niche with limited immune organization, validated across primary and metastatic PDAC contexts**.

This direction is now supported by manual H&E/spatial overlay review, specimen-group metadata provenance, random-core controls, and CAF-core threshold sensitivity. It should still be framed as spatial association rather than causality.

## Key Metrics

- n_samples: 113
- n_gse282302: 108
- n_gse274103: 5
- median_barrier_fraction_z_gt1: 0.0535179640718562
- median_immune_core_fraction_z_gt1: 0.0777327935222672
- median_tumor_aggressive_fraction_z_gt1: 0.0480167014613778
- median_rho_barrier_vs_immune_core: 0.3710413099562499
- median_rho_barrier_vs_ifn_mhc: 0.41676439674012594
- median_rho_barrier_vs_tumor_aggressive: 0.3700604692432501
- median_rho_immune_core_vs_maturity: 0.5694026102650378
- median_safe_rho_barrier_vs_immune_core: 0.34164116736265465
- median_safe_rho_barrier_vs_ifn_mhc: 0.4342439453116909
- median_safe_rho_barrier_vs_tumor_aggressive: 0.38648023208456933
- n_detectable_immune_hub: 0
- n_mature_like_hub: 0

## Candidate Story Ranking

| candidate_story | mvp_score_0_to_3 | evidence | caveat | next_action |
|---|---:|---|---|---|
| CAF-myeloid inflammatory stromal niche with limited immune organization | 3 | median safe rho barrier~immune=0.342; median safe rho barrier~tumor_aggressive=0.386; random-core controls show true CAF-myeloid cores beat random same-size cores in GSE282302 and GSE272362 primary/liver metastasis samples; GSE272362 overlays show aligned internal tissue compartments. | Current MVP shows CAF-myeloid co-localization with immune/IFN rather than simple immune exclusion. Mature TLS is not supported. | Move to manuscript-style Figure 1/2 draft; keep TLS as a negative/secondary claim. |
| Post-neoadjuvant residual spatial ecology | 2 | GSE282302 contributes 108 directly scored ST-H&E samples. | Treatment context and response metadata must be resolved before any therapy-response claim. | Curate GSE282302 sample/patient metadata and compare residual ecotype axes. |
| Immune hub maturation arrest | 1 | detectable immune-hub-like signal in 0/113 samples; mature-like signal in 0/113 samples; median rho immune_core~maturity=0.569 | MVP detection uses signatures only; do not call mature TLS without B/T/FDC/GC morphology review. | Keep as secondary unless manual overlays show true lymphoid aggregates. |
| Primary-to-metastatic spatial ecology remodeling | 2 | GSE272362 RDS adds 30 samples across primary tumor, liver metastasis, lymph node metastasis, and normal pancreas; random-core controls support strong CAF-core-centered IFN/MHC and immune-core gradients in liver metastases and tumor-aggressive gradients in primary/liver/lymph node metastases. | Tissue-site labels are derived from RDS `Origin` and specimen-group provenance is audited; lymph node metastasis immune gradients diverge from primary/liver patterns. | Use as validation and biological divergence layer. Keep patient-level claims separate unless metadata are further audited. |
| H&E-readable PDAC spatial ecotype | 1 | All directly scored samples have paired tissue images; model not trained yet. | No claim until patient/sample split H&E feature model is validated. | Train a small exploratory ResNet/timm feature model only after story axis is chosen. |
| Neural niche-associated immune exclusion | 0 | Neural niche has not been targeted in MVP; GSE202740 not yet added. | Only add if neural_schwann signal is strong in current scored samples or if project pivots. | Defer. |

## Manuscript-Level Decision

Proceed with the CAF-myeloid inflammatory stromal niche as the main manuscript story. The strongest claim is no longer simply "co-localization"; it is now a random-core-controlled and threshold-stable spatial-gradient result, with independent primary/liver-metastasis validation and lymph-node-specific divergence.

Current Figure 1/2 artifacts:

- `results/figures/main/figure1_draft.pdf`
- `results/figures/main/figure1_draft.png`
- `results/tables/figure1_source.csv`
- `results/figures/main/figure2_draft.pdf`
- `results/figures/main/figure2_draft.png`
- `results/tables/figure2_source.csv`
- `results/reports/manuscript_story_outline.md`
- `metadata/gse272362_provenance_audit.md`

Do not promote mature TLS, survival, therapy response, causality, or H&E-only prediction claims unless new evidence is added.

## Immediate Next Step

1. Upgrade Figure 1/2 drafts into final panel-ready figures and write captions.
2. Write a full reproducible Methods section from the current notes.
3. Decide whether patient-level metadata are sufficient for matched primary-metastasis analysis; otherwise keep claims at specimen/site level.
4. Decide whether to stop at a complete spatial transcriptomics story or add H&E feature modeling as a secondary extension.
