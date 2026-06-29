# Figure plotting-code archive

Date: 2026-06-29

This directory preserves the plotting code used for each main and Extended Data figure in the active revised manuscript.

Active manuscript:

`results/revision_2026_06_29/manuscript/Manuscript_NatureSubjournal_revised.md`

## Contents

- `Figure_1/` to `Figure_4/`: copied plotting script used to assemble the current main figures.
- `Extended_Data_Figure_1/` to `Extended_Data_Figure_8/`: copied plotting scripts used to assemble the current Extended Data figures.
- `shared_scripts/`: supporting historical figure-generation scripts that contributed to the current figure suite or source-data generation.
- `panel_code_index.csv`: panel-level map from figure panel to output file, plotting script, source data and notes.
- `plot_code_file_manifest.csv`: SHA256 checksum manifest for archived plotting scripts and index files.

## Update rule

Every time a figure panel is changed, update all of the following in the same revision:

1. The plotting script in the relevant figure folder.
2. The source data file used by the panel, if values changed.
3. `panel_code_index.csv`, if the script, panel label, output file or source data changed.
4. `plot_code_file_manifest.csv`, after any plotting-code or index update.
5. The manuscript legend and `docs/figure_panel_map.md`, if the panel meaning or numbering changed.
6. The submission package manifest and zip, if a new package is being prepared.

## Notes

- Several current figure scripts generate a full multi-panel figure rather than one script per panel. The panel-level index therefore maps each panel to the full plotting script and the relevant source-data file.
- Raw H&E images, Visium objects, Xenium matrices and other large public-data mirrors are not duplicated in this archive. They remain in the project raw-data locations or public repositories and are referenced through the source-data/provenance documents.
