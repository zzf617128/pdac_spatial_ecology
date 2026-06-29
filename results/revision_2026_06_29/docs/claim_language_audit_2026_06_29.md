# Claim language audit

Date: 2026-06-29

## Scope

This audit checks whether the revised manuscript uses causal, clinical, TLS, mechanism and cell-resolution language within the evidence boundaries defined by the manuscript-modification task list.

## Search Terms Audited

- `causal`, `causality`
- `clinical validation`, `prognostic`, `prediction`, `predictor`, `response`, `therapy`
- `TLS`, `tertiary lymphoid`
- `ligand-receptor`, `mechanism`
- `single-cell-resolved`, `ground truth`
- `subtype`
- `drive`, `mediate`
- `validate`, `validation`

## Findings

- Causal language is bounded by explicit statements such as "does not establish causality", "not causal ligand-receptor inference" and "perturbation is required for causal inference".
- Clinical language is bounded by "non-spatial context", "exploratory univariable survival context", "not clinical validation" and "not a clinical prediction model".
- TLS language is bounded by the TLS-maturity stress test and the statement that immune-hub features should not be reframed as mature TLS without histologic, FDC or germinal-center validation.
- LN language is bounded by "five-sample lymph-node metastasis subset", "hypothesis-generating lead" and "not a definitive metastatic subtype".
- Xenium language is bounded by targeted-panel and cell-resolution-support phrasing, not whole-transcriptome cell-state abundance.
- Reference projection and NNLS are described as marker-constrained computational projections, not single-cell-resolved spatial abundance or ground truth.
- `validate` is used for external dataset support and is paired with explicit limitations where needed, especially for GSE274557 and GSE274673.

## Residual Risk

The title and abstract still use compact phrasing such as "define a reproducible CAF-myeloid spatial architecture". This is acceptable because the Methods and Results define the architecture statistically, but the wording should remain paired with "observational", "associated with" and "requiring perturbational validation" language.

## Outcome

No unbounded causal, mature-TLS, clinical-grade prediction, prognostic-model, treatment-response, spatial-TCGA, direct ligand-receptor or single-cell-ground-truth claims were detected in the active revised manuscript after the current pass.
