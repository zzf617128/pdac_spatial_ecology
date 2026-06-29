# Manual Review Log

## 2026-06-24 MVP H&E Overlay Review

Reviewed representative overlays from:

- `GSM8641033_C2_D11_ROI2_s2`
- `GSM8641031_C2_D11_ROI1_s2`

Preliminary observations:

1. Spot coordinates and H&E image alignment appear broadly correct after applying `tissue_hires_scalef`.
2. High CAF-myeloid, immune-core, tumor-aggressive, and IFN/MHC signals often appear at tissue edges, interfaces, or broad tissue compartments.
3. Current overlays do not yet show obvious mature TLS-like aggregates. Immune/IFN-high regions look more like diffuse or compartment-level inflammatory programs.
4. The current MVP supports a cautious direction: CAF-myeloid / inflammatory stromal niche with limited immune organization.
5. Main risk: edge effects, tissue folds, necrosis, low-tissue spots, or broad tissue composition may drive some high-score regions.

Required before main claims:

- Add tissue/edge/background QC to overlays.
- Review more top and low-control samples.
- Confirm high CAF-myeloid regions correspond to plausible stromal/fibrotic histology.
- Confirm immune-core high regions are not simply tissue edge or low-quality regions.
- Do not call mature TLS without convincing B/T/FDC/GC/plasma and morphology support.

## 2026-06-24 Edge-Aware Overlay And Gradient Review

Reviewed updated overlays with red rings marking edge/background-risk spots.

Additional observations:

1. Edge/background-risk spots are visible and nontrivial in GSE282302, but high CAF-myeloid signal is not globally dominated by risk spots.
2. A cleaner representative sample is `GSM8641105_C3_D8_ROI3`, where CAF-myeloid, immune core, IFN/MHC, and tumor-aggressive programs occupy a broad tissue region rather than only the tissue edge.
3. Edge-QC-safe gradient analysis shows immune core, IFN/MHC, tumor-aggressive, and immune-maturity scores are generally higher near high CAF-myeloid cores and decline with distance.
4. The current strongest story should be framed as a **CAF-myeloid inflammatory stromal niche with limited immune organization in post-neoadjuvant PDAC**, not as mature TLS abundance or simple immune exclusion.

Recommended main result direction:

- Define high CAF-myeloid cores in each section.
- Show spatial gradients of IFN/MHC, immune core, tumor-aggressive, and maturity-like signals around those cores.
- Use H&E overlays to show these cores correspond to real tissue regions.
- Treat immune hub/TLS as secondary and mostly negative unless manual morphology review proves otherwise.

## 2026-06-24 GSE272362 RDS Integration Review

Processed Zenodo `PDAC_Updated.rds` for `GSE272362`.

Summary:

1. The RDS contains 91,496 spots and 30 samples.
2. Tissue-site labels are available in the RDS metadata column `Origin`.
3. Interpreted groups are: 10 primary tumor, 12 liver metastasis, 5 lymph node metastasis, and 3 normal pancreas samples.
4. CAF-myeloid distance gradients replicate in primary tumor and liver metastasis samples:
   - primary tumors: IFN/MHC and tumor-aggressive gradients are negative in 10/10 samples.
   - liver metastases: immune core and IFN/MHC gradients are negative in 12/12 samples.
   - lymph node metastases show more heterogeneous immune gradients, but tumor-aggressive gradients are negative in 5/5 samples.
5. This strengthens the story as a conserved CAF-myeloid inflammatory niche rather than a dataset-specific artifact.

Status before manuscript-level claims:

- Completed: documented specimen-group provenance for the RDS `Origin` field in `metadata/gse272362_provenance_audit.md`.
- Completed: generated representative GSE272362 H&E overlays for primary tumor, liver metastasis, and lymph node metastasis.
- Completed: added random-core controls for distance-gradient specificity.
- Avoid patient-level survival or treatment claims until patient/sample metadata are audited.

## 2026-06-24 Random-Core Permutation Review

Completed 1,000-iteration random-core controls for:

- `GSE282302` and `GSE274103` using edge-QC-safe spots.
- `GSE272362` using RDS-derived spot scores and tissue-site labels from `Origin`.

Interpretation:

1. The CAF-myeloid niche signal is not explained by random spatial cores.
2. In `GSE282302`, true CAF-myeloid cores beat random same-size cores in 93/108 samples for IFN/MHC, 94/108 for tumor-aggressive, and 87/108 for immune core.
3. In `GSE272362` liver metastases, true CAF-myeloid cores beat random cores in 12/12 samples for IFN/MHC and immune core.
4. In `GSE272362` primary tumors, true CAF-myeloid cores beat random cores in 9/10 samples for IFN/MHC and immune core, and 10/10 for tumor-aggressive.
5. Lymph node metastases appear biologically divergent: tumor-aggressive signal remains CAF-core associated, while immune/IFN programs are not consistently CAF-core centered.

Updated manual priority:

- Inspect representative overlays from `GSE272362` primary tumor, liver metastasis, and lymph node metastasis.
- Confirm whether liver metastasis CAF-core/IFN regions correspond to plausible fibrotic-inflammatory tissue compartments.
- For lymph node metastases, specifically check whether immune-rich regions are organized by lymphoid architecture rather than CAF-myeloid cores.

## 2026-06-24 GSE272362 H&E Overlay Review

Generated six representative overlays with exact RDS-to-GEO coordinate matching:

- `IU_PDA_T3` primary tumor, selected for tumor-aggressive CAF-core association.
- `IU_PDA_T1` primary tumor, selected for IFN/MHC CAF-core association.
- `IU_PDA_HM10` liver metastasis, selected for IFN/MHC CAF-core association.
- `IU_PDA_HM5` liver metastasis, selected for tumor-aggressive CAF-core association.
- `IU_PDA_LNM12` lymph node metastasis, selected for tumor-aggressive CAF-core association.
- `IU_PDA_LNM7` lymph node metastasis, selected as immune-divergent.

Technical check:

- All six overlays had 400/400 exact coordinate matches between RDS spot coordinates and GEO `tissue_positions.csv`.
- This supports reliable H&E/spot alignment for the selected GSE272362 examples.

Visual observations:

1. `IU_PDA_HM10` is currently the strongest visual example: CAF-myeloid core, IFN/MHC, tumor-aggressive, and immune-core signals co-localize within a broad internal tissue compartment rather than only at slide/background edges.
2. `IU_PDA_HM5` independently supports a liver-metastasis CAF-core region with overlapping IFN/MHC and tumor-aggressive programs.
3. `IU_PDA_T1` and `IU_PDA_T3` support the same general primary-tumor pattern, with CAF-core regions spanning internal tissue regions and overlapping IFN/MHC or tumor-aggressive signals.
4. `IU_PDA_LNM12` supports CAF-core-associated tumor-aggressive signal in lymph node metastasis.
5. `IU_PDA_LNM7` illustrates the lymph-node divergence: immune/IFN-rich regions are not simply centered on CAF-myeloid cores and may reflect lymphoid tissue architecture.

Recommended figure use:

- Use `IU_PDA_HM10` as the main GSE272362 visual validation panel.
- Use `IU_PDA_T1` or `IU_PDA_T3` as primary-tumor companion examples.
- Use `IU_PDA_LNM7` as a divergence/contrast example if the manuscript includes lymph-node-specific interpretation.

## 2026-06-24 Manuscript Figure and Robustness Review

Added a top-journal-oriented figure and robustness pass:

- CAF-core threshold sensitivity was run for top 15%, top 10%, and top 5% CAF-myeloid cores.
- GSE272362 specimen labels were corrected so normal pancreas remains separate from primary tumor in sensitivity summaries.
- Figure 1 draft now tests the main-cohort claim with cohort design, random-core controls, threshold sensitivity, and a representative GSE282302 overlay.
- Figure 2 draft now tests the independent validation/divergence claim with site counts, site-specific random-core controls, and full-width liver, primary, and lymph-node overlays.
- Source tables were written for both figure drafts.

Artifacts reviewed:

- `results/tables/caf_core_threshold_sensitivity_summary.csv`
- `results/figures/mvp/sensitivity/caf_core_threshold_sensitivity.png`
- `results/figures/main/figure1_draft.png`
- `results/figures/main/figure2_draft.png`
- `results/tables/figure1_source.csv`
- `results/tables/figure2_source.csv`
- `results/reports/manuscript_story_outline.md`

Top-journal claim boundary:

1. Supported: CAF-myeloid spatial cores organize IFN/MHC, immune-core, and tumor-aggressive programs beyond random same-size cores.
2. Supported: the core-gradient result is not dependent on a single CAF-core threshold.
3. Supported: primary tumors and liver metastases validate the CAF-core-centered inflammatory/tumor-aggressive pattern.
4. Supported: lymph node metastases retain tumor-aggressive CAF association but diverge for immune/IFN organization.
5. Not supported: mature TLS, clinical response, survival, causal mechanism, or H&E-only prediction.

Next manual checks:

- Keep GSE272362 claims at specimen/site level unless patient-level clinical metadata are separately audited.
- Write final captions with exact `n`, null model, threshold definitions, and coordinate-matching details.

## 2026-06-24 Figure 1d Representative Overlay Review

Reviewed the available GSE282302 H&E overlay contact sheet and inspected candidate examples individually.

Decision:

- Keep `GSM8641105_C3_D8_ROI3` as the Figure 1d representative GSE282302 overlay.

Rationale:

1. The H&E tissue area is large and intact enough for a main-text visual example.
2. CAF-myeloid, immune-core, tumor-aggressive and IFN/MHC panels show a coherent shared spatial gradient.
3. Compared with `GSM8641070_C3_D12_ROI2_s2`, the selected example has fewer large central holes and less boundary-dominated interpretation risk.
4. Compared with `GSM8641067_C3_D12_ROI1_s1`, the selected example communicates the spatial-gradient claim more directly.

Remaining figure work:

- Final journal typography and source-data naming should be adjusted only after a target journal format is chosen.
