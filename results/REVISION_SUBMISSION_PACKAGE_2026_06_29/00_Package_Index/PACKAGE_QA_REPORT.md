# Package QA report

Date: 2026-06-29

## Package contents checked

- Manuscript files: 4 files in `01_Manuscript/`.
- Main figures: 4 PDF files in `02_Main_Figures/`.
- Extended Data figures: 12 official files in `03_Extended_Data_Figures/`.
- Supplementary tables: 9 CSV files in `04_Supplementary_Tables/`.
- Source data: 31 files in `05_Source_Data/`, including `source_data_file_manifest.csv`.
- Provenance and QA documents: 21 files in `06_Methods_Provenance_QA/`.
- Revision analysis outputs: 16 CSV files in `07_Analysis_Outputs/`.
- Script index: 4 scripts in `08_Scripts_Index/`.

## Automated checks completed

- The revised manuscript DOCX was regenerated from the revised Markdown source.
- DOCX structural checks confirmed the updated Figure 1, Figure 2 and Extended Data Figure 7 headings.
- Figure 1 was regenerated with a contiguous-null panel and rendered successfully for QA preview.
- ED7/NMF file naming conflict was resolved by keeping ED7 as the integrated specificity/sensitivity suite and removing the standalone NMF rank plot from the official Extended Data folder.
- Source data and package manifests record SHA256 checksums.

## Known limitations before final upload

- Author metadata, funding, acknowledgements, author contributions and competing interests still require author confirmation.
- Code availability still needs a real reviewer-accessible repository link.
- Main figure and Extended Data formats may need conversion after the target journal is chosen.
- Full visual DOCX/PDF render QA should be repeated before upload.

