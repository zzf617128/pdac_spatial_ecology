# Negative Or Weak Results Log

## Zenodo 10246315 PDAC vascular-niche IMC

- Status: loaded and analyzed in MVP.
- Positive: processed Seurat object is usable; cell type/neighborhood annotations, coordinates, CD44, Collagen, CD68, CD11b, CD14, HLA-DR, pan-Keratin are available.
- Weak point: corrected ROI-level kNN analysis showed CAF/stromal to myeloid/macrophage adjacency near random (median log2 observed/expected around 0.005).
- Decision: do not include in candidate strong ED figure. Keep as internal feasibility and marker availability support.

## GSE240078 GeoMx caveats

- Strong: stroma/TME AOIs enrich CAF/myCAF, panCAF/matrix, myeloid/macrophage, immune-core, TGF/EMT and IFN/APC modules.
- Boundary: SPP1/TAM module is carcinoma-enriched in MVP and must not be used as stromal SPP1/TAM support.
- Boundary: B/plasma module has insufficient targeted-panel coverage.

## Lower-Priority/Excluded Sources

- 10x Xenium pancreatic cancer demo: useful for pipeline testing, not strong manuscript validation due weak clinical metadata.
- Existing GSE274673 Xenium: already part of current project; do not repeat unless used as benchmark.
- Non-spatial scRNA references: can support module origin, not ED10 spatial validation.

## GSE310352 CosMx automated deep-analysis caveats

- SPP1 and CD44 are both available, but SPP1+ macrophage to CD44+ tumor/stromal adjacency was negative in the automated FOV-level analysis (median slide log2 O/E = -0.509 among analyzable slides). Do not claim direct SPP1-CD44 validation from this dataset.
- CAF/matrix to tumor-aggressive epithelial adjacency was consistently negative (median slide log2 O/E = -1.305). Do not display this as supportive evidence.
- CAF/matrix to myeloid/macrophage adjacency was positive where myeloid-like populations were called, but only 4/8 slides were analyzable for this rule-based myeloid state. Treat as secondary/source-level support unless patient/annotation mapping resolves the missing myeloid calls.
- The cleanest automated positive signal is CAF/matrix to TGF/EMT-interface adjacency across 8/8 slides.
