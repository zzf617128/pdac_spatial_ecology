# Reproducibility Lock 2026-06-30 With ED10 v1

Enhanced submission candidate archive. This directory is additive and does not overwrite the base lock.

## Core Candidate
- Manuscript: `manuscript/Manuscript_with_ED10_v1_submission_safe_slidelevel_patch.docx`
- ED10 PDF: `ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf`
- ED10 SVG: `ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg`
- ED10 source data: `ed10/source_data/Source_Data_Extended_Data_Figure_10_v1.csv`

## Interpretation Boundaries
- GSE310352 is slide/FOV-level support only.
- GSE310352 cell states are rule-based because author annotations were unavailable.
- TGF/EMT is interpreted as stromal-interface, not tumor-intrinsic EMT.
- No causal signaling or direct SPP1-CD44 validation is claimed.
- Cho IMC is source-only and not part of ED10 v1.

## Manifests
- `enhanced_submission_file_manifest.csv`
- `enhanced_ed10_panel_map.csv`
- `enhanced_source_data_manifest.csv`
- `enhanced_script_manifest.csv`
- `enhanced_dataset_manifest.csv`
- `enhanced_parameter_manifest.csv`
- `checksums_sha256_enhanced_lock.txt`

## Optional backfill status
- GSE240078 metadata/module-score tables were found in `outputs/orthogonal_validation_2026_06_30/tables/` and copied.
- Final and cross-platform evidence summary docs/tables were created inside this lock from archived outputs.
