# GSE310352 Cell-State Definition Transparency Report

## Scope
This QC analysis evaluates whether the existing rule-based GSE310352 CosMx cell-state labels are marker- and context-consistent. It does not redefine labels, change adjacency results, modify ED10, or make patient-level claims.

## Data Boundary
- Public metadata did not allow reliable recovery of patient, specimen or tissue-block identifiers for the processed CosMx slides.
- Processed CSV files did not include author cell-type annotations.
- All support remains slide-level/FOV-level and rule-based.

## CAF/matrix-like Label
- CAF/matrix-like cells show positive overall CAF/matrix module enrichment versus the rest of cells (delta=1.37).
- The fraction above each slide's 80th percentile is low for PanCK (0), CD45 (0) and CD3 (0.029), consistent with a non-epithelial, non-immune stromal-like assignment.

## TGF/EMT Stromal-Interface Label
- TGF/EMT stromal-interface cells show positive TGF/EMT module enrichment versus the rest of cells (delta=0.795).
- Among TGF/EMT interface cells, the fraction overlapping CAF/matrix-like cells is 0.415 and the fraction overlapping tumor epithelial-like cells is 0.585.
- The fraction overlapping immune-like cells is 0.012, supporting the interpretation that this label is not primarily an immune state.
- Because the TGF/EMT interface state contains both CAF/matrix-like and epithelial-proximal components, the safest label remains stromal-interface rather than tumor-intrinsic EMT.

## Tumor And Immune IF Checks
- Tumor epithelial-like cells show high PanCK enrichment: fraction above slide q80=0.585.
- Immune-like cells show elevated CD45/CD3 context: CD45 fraction above slide q80=0.545; CD3 fraction above slide q80=0.355.

## Spatial Identity Context
Nearest-distance QC was computed within FOVs using sampled TGF/EMT interface query cells where needed for tractability. Overall summary:
```
                        target_type median_actual_distance
                             <char>                  <num>
1:                  CAF/matrix-like               77.94983
2: PanCK-high tumor epithelial-like               90.69788
   median_random_query_distance median_random_target_distance
                          <num>                         <num>
1:                     147.4955                     107.59583
2:                     140.8976                      94.69563
   ratio_actual_to_random_query ratio_actual_to_random_target
                          <num>                         <num>
1:                    0.5256403                     0.7121360
2:                    0.6859623                     0.9870283
```
These distances are QC context only and are not a new biological claim or a replacement for the original adjacency analysis.

## CAF/TGF Gene-Overlap Boundary
The direct CAF versus TGF/EMT gene overlap was previously limited to MMP2. Prior robustness tables retained positive CAF/matrix to TGF/EMT adjacency after removing MMP2 and after removing MMP2 plus ITGA5:
```
                      tgfb_module n_slides positive_slides median_slide_log2_oe
                           <char>    <int>           <int>                <num>
1:                       tgfb_emt        8               8            0.6781272
2:             tgfb_emt_no_shared        8               8            0.5025719
3: tgfb_emt_no_shared_no_integrin        8               8            0.4491418
```
This supports the view that the displayed adjacency is not solely a single shared-gene artifact, while still remaining observational and rule-based.

## Safest Current Claim
GSE310352 CosMx provides slide-level/FOV-level transparency support that rule-based CAF/matrix-like and TGF/EMT stromal-interface labels are marker-consistent enough for use as orthogonal stromal-interface support.

## Claims Not Supported
- Do not claim patient-level or specimen-level validation.
- Do not claim author-annotated cell types.
- Do not claim tumor-intrinsic EMT.
- Do not claim causal EMT induction or CAF-to-tumor signaling.
- Do not claim direct Visium distance-gradient reconstruction.

## Outputs
- `tables/gse310352_cell_state_definition_marker_enrichment.csv`
- `tables/gse310352_cell_state_if_marker_qc.csv`
- `tables/gse310352_cell_state_overlap_audit.csv`
- `tables/gse310352_tgfemt_spatial_identity_context.csv`
- `figures/gse310352_cell_state_marker_enrichment_heatmap.pdf`
- `figures/gse310352_cell_state_if_marker_qc.pdf`
- `figures/gse310352_cell_state_overlap_upset_or_barplot.pdf`
- `figures/gse310352_tgfemt_spatial_identity_context.pdf`
- `figures/gse310352_cell_state_definition_transparency_qc.pdf`
- `figures/gse310352_cell_state_definition_transparency_qc.svg`
