# PDAC CAF-Myeloid Spatial Ecology

This repository contains the analysis code, figure-generation scripts, processed source tables and revision package for a public-data spatial transcriptomics study of CAF-myeloid stromal cores in pancreatic ductal adenocarcinoma.

## Active Revision Package

The current manuscript-revision workspace is:

- `results/revision_2026_06_29/`

The current submission-ready working package is:

- `results/REVISION_SUBMISSION_PACKAGE_2026_06_29/`

Key files:

- Revised manuscript: `results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.docx`
- Figure panel map: `results/revision_2026_06_29/docs/figure_panel_map.md`
- Source-data manifest: `results/revision_2026_06_29/source_data/source_data_file_manifest.csv`
- Plotting-code archive: `results/revision_2026_06_29/figure_plot_code_archive/`
- Panel-level plotting-code index: `results/revision_2026_06_29/figure_plot_code_archive/panel_code_index.csv`

## Scope

The repository is intended for peer-review access to:

- spatial preprocessing and module-scoring scripts;
- same-size random-core and spatially contiguous random-core analyses;
- CAF-only, myeloid-only and CAF-myeloid anchor sensitivity analyses;
- marker-overlap sensitivity analyses;
- lymph-node metastasis leave-one-out analyses;
- NMF rank-stability analyses;
- Xenium, TCGA and H&E context analyses;
- manuscript display-item source data and plotting scripts.

## Data Policy

Raw public datasets are not committed to this repository because they are large and remain available from GEO, TCGA, UCSC Xena, Zenodo or linked public resources. The active revision package includes processed source data for plotted quantitative panels and manifests that identify file sizes, checksums and provenance.

## Reproducibility Notes

The revision workspace includes:

- `results/revision_2026_06_29/docs/revision_log.md`
- `results/revision_2026_06_29/docs/claim_language_audit_2026_06_29.md`
- `results/revision_2026_06_29/docs/revision_task_completion_audit_2026_06_29.md`
- `results/revision_2026_06_29/figure_plot_code_archive/README.md`

These documents describe the active manuscript state, claim boundaries, figure mapping and figure-code synchronization rules.
