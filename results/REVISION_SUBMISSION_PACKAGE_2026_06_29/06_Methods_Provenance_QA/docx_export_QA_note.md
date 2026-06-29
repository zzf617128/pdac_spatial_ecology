# DOCX export QA note

Date: 2026-06-29

## Exported file

- `manuscript/Manuscript_NatureSubjournal_revised.md`
- `manuscript/Manuscript_NatureSubjournal_revised.docx`

## Structural QA

The revised DOCX was generated from the revised Markdown source using `scripts/85_export_revised_manuscript_docx.py`.

Structural open check with `python-docx` passed:

- Non-empty paragraphs: 150
- Heading paragraphs: 44
- File size: 58,568 bytes
- First paragraph/title: `CAF-myeloid stromal cores mark inflammatory and tumor-aggressive spatial programs in pancreatic cancer`

## Render QA status

The required LibreOffice render step was attempted with the documents skill `render_docx.py`, but LibreOffice stalled during DOCX-to-PDF conversion and produced no PNG/PDF output in `manuscript/docx_render_QA/`. The stalled render processes were stopped after checking the output directory.

Therefore, this DOCX has passed structural QA but has not passed visual render QA. Before submission, open the DOCX in Word or LibreOffice and inspect page layout, spacing and figure legend wrapping.
