# P0 manuscript completeness audit

Date: 2026-06-29

## Status

P0 manuscript-completeness items from `PDAC_CAF_myeloid_manuscript_revision_tasks_for_Codex.md` have been implemented where they can be completed without author-specific administrative information.

## Completed

- Added formal manuscript sections for Acknowledgements, Author contributions, Funding, Competing interests, Ethics statement and Supplementary information.
- Replaced bracket-style author placeholders with submission-ready pre-submission text.
- Added author byline, affiliations, corresponding-author details and equal-contribution note from author-supplied metadata.
- Replaced the Acknowledgements, Author contributions and Competing interests placeholders with manuscript-ready text.
- Updated Data availability to list public accessions: GSE282302, GSE274103, GSE235315, GSE272362, GSE274557, GSE274673, GSE202051 and TCGA-PAAD.
- Updated Code availability to state that analysis code and processed source tables will be deposited in a persistent public repository before publication, with controlled-access/private-review availability during journal evaluation.
- Confirmed `Statistical analysis` is present and includes section/specimen-level inference, spot/cell pseudo-replication boundary, random-core delta definition, empirical p-value formula, multiple-testing language, NMF rank-sensitivity interpretation, reference projection/NNLS boundary and TCGA non-spatial context boundary.
- Confirmed Supplementary Table 1 has dataset/accession/platform/sample/role fields for all eight core data resources.
- Confirmed Supplementary Table 2 has `module_name`, `gene_symbol`, `source_reference`, `source_type`, `used_in_figures` and `notes`, with 315 rows across 34 modules.
- Confirmed Supplementary Table 5 provides module-overlap sensitivity with shared-gene count, Jaccard index, overlap-sensitive delta and support summaries.
- Added Supplementary Table 8 with mechanism-triangulation scoring metrics, thresholds, observed values and assigned `+`/`++` evidence scores.
- Completed terminology scan for obvious mixed forms including `SPP1-TAM`, `SPP1 TAM`, `immune core`, `IFN-MHC`, `IFN MHC`, `TGFb`, `TGF-beta EMT`, `treat-naive`, `post NACT`, bracket-style author placeholders and local-only code availability language.

## Still Requires Author Input

- Persistent code repository or DOI: the private GitHub repository `https://github.com/zzf617128/pdac_spatial_ecology` has been created and pushed. Reviewer access settings may still need to be enabled according to journal instructions, and a public DOI/release should be created before publication.
- Optional ORCID details if required by the target journal.

## Verification

- Revised manuscript source: `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.md`.
- Revised Word export: `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.docx`.
- Current DOCX structural check confirms the Funding, Ethics statement, Supplementary information, ED7 and ED8 legends are present.
