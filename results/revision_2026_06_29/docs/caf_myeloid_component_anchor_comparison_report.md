# CAF-only, myeloid-only and CAF-myeloid anchor comparison

n_perm per sample and anchor: 1000

This analysis asks whether the combined CAF-myeloid core adds value beyond CAF-only or myeloid-only anchors.
Negative delta indicates stronger target-program centering around the biological anchor than around same-size random anchors.

## Key anchor summary

- CAF-myeloid combined -> IFN/MHC: median delta -0.285, support 184/205.
- CAF-myeloid combined -> SPP1/TAM: median delta -0.377, support 203/205.
- CAF-myeloid combined -> TGF-beta/EMT: median delta -0.478, support 203/205.
- CAF-myeloid combined -> immune-core: median delta -0.261, support 177/205.
- CAF-myeloid combined -> tumor epithelial: median delta 0.018, support 97/205.
- CAF-myeloid combined -> tumor-aggressive: median delta -0.277, support 175/205.
- CAF-only -> IFN/MHC: median delta -0.104, support 145/205.
- CAF-only -> SPP1/TAM: median delta -0.113, support 161/205.
- CAF-only -> TGF-beta/EMT: median delta -0.427, support 198/205.
- CAF-only -> immune-core: median delta -0.103, support 140/205.
- CAF-only -> tumor epithelial: median delta 0.096, support 83/205.
- CAF-only -> tumor-aggressive: median delta -0.202, support 158/205.
- immune-high -> IFN/MHC: median delta -0.383, support 205/205.
- immune-high -> SPP1/TAM: median delta -0.091, support 146/205.
- immune-high -> TGF-beta/EMT: median delta -0.021, support 109/205.
- immune-high -> immune-core: median delta -0.539, support 205/205.
- immune-high -> tumor epithelial: median delta 0.128, support 55/205.
- immune-high -> tumor-aggressive: median delta 0.092, support 74/205.
- myeloid-only -> IFN/MHC: median delta -0.328, support 197/205.
- myeloid-only -> SPP1/TAM: median delta -0.475, support 205/205.
- myeloid-only -> TGF-beta/EMT: median delta -0.275, support 181/205.
- myeloid-only -> immune-core: median delta -0.284, support 193/205.
- myeloid-only -> tumor epithelial: median delta -0.044, support 120/205.
- myeloid-only -> tumor-aggressive: median delta -0.210, support 167/205.
- tumor-high -> IFN/MHC: median delta -0.071, support 126/205.
- tumor-high -> SPP1/TAM: median delta -0.124, support 149/205.
- tumor-high -> TGF-beta/EMT: median delta -0.039, support 115/205.
- tumor-high -> immune-core: median delta 0.064, support 85/205.
- tumor-high -> tumor epithelial: median delta -0.632, support 205/205.
- tumor-high -> tumor-aggressive: median delta -0.304, support 183/205.

## Interpretation rule

If CAF-only and CAF-myeloid combined anchors perform similarly, the manuscript should frame the architecture as a CAF-dominant stromal core with myeloid enrichment rather than claiming that the combined score is uniquely superior.
