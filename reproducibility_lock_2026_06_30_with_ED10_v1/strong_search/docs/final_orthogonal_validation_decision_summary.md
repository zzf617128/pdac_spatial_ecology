# Final Orthogonal Validation Decision Summary

This summary records the final display decision for the enhanced ED10 v1 submission candidate. It is a provenance document only and does not add new biological analysis.

## Included In ED10 v1

### GSE240078 GeoMx DSP
- Decision: included in ED10 v1.
- Role: compartment-level support for CAF/matrix and immune/TME programs in TME/stromal AOIs versus carcinoma AOIs.
- Boundary: compartment-level support only; SPP1/TAM is not treated as stromal support; no causal signaling or Visium distance-gradient reconstruction.

### GSE199102 GeoMx WTA
- Decision: included in ED10 v1.
- Role: independent GeoMx replication/concordance using CAF-plus-immune versus epithelial segment comparisons and paired ROI support.
- Boundary: compartment-level support only; no direct SPP1-CD44 validation and no causal signaling.

### GSE310352 CosMx
- Decision: included in ED10 v1.
- Role: slide/FOV-level CAF/matrix-associated TGF/EMT stromal-interface support.
- Boundary: public metadata did not allow reliable recovery of patient, specimen or tissue-block identifiers; processed CSV files lacked author cell-type annotations; cell states are rule-based; TGF/EMT is interpreted as stromal-interface rather than tumor-intrinsic EMT.

## Not Included In ED10 v1

### Cho IMC Zenodo 15596960
- Decision: source-only archive.
- Role: modest protein-level CAF/stromal-to-myeloid adjacency context.
- Reason not included: main adjacency was positive and null-supported but modest; stromal-neighborhood macrophage enrichment and tumor-stroma-myeloid interface metrics were weak/borderline; HLA-DR interface context was negative.

### Vascular-niche IMC Zenodo 10246315
- Decision: internal feasibility only.
- Role: marker/processed-object feasibility.
- Reason not included: corrected ROI-level CAF/stromal-to-myeloid adjacency was near random.

## Final Figure Decision

ED10 v1 remains the final manuscript-facing enhanced orthogonal validation figure. No ED10 v2 was generated.

## Key Claim Boundaries

- GeoMx supports compartment-level CAF/matrix and immune/TME programs.
- CosMx supports CAF/matrix-associated TGF/EMT stromal-interface organization at slide/FOV level.
- No patient-level or specimen-level validation is claimed for GSE310352.
- No tumor-intrinsic EMT claim is made.
- No causal signaling is claimed.
- No direct SPP1-CD44 validation is claimed.
- No Visium distance-gradient reconstruction is claimed.
- No lymph-node immune-uncoupling validation is claimed from ED10, GSE310352 or Cho IMC.
