# Cho IMC Final Decision Report

Decision: B. Use Cho IMC as secondary/source-only support

## Decision summary
                                  metric n_roi n_patients  median_roi
                                  <char> <int>      <int>       <num>
1:         CAF_stromal_to_myeloid_log2OE    64         12  0.18523527
2: stromal_neighborhood_macrophage_delta    64         12  0.00370224
3: tumor_stroma_myeloid_interface_log2OE    64         12  0.06087448
4:        HLA_DR_interface_myeloid_delta    64         12 -0.49857557
   roi_support_fraction_positive patient_median
                           <num>          <num>
1:                    0.79687500    0.203076568
2:                    0.51562500    0.005223962
3:                    0.53125000    0.005064981
4:                    0.09677419   -0.425445520
   patient_support_fraction_positive median_random_label_p median_abundance_p
                               <num>                 <num>              <num>
1:                         0.9166667            0.00990099         0.00990099
2:                         0.5833333                    NA                 NA
3:                         0.5000000                    NA                 NA
4:                         0.1818182                    NA                 NA
   median_density_p
              <num>
1:       0.00990099
2:               NA
3:               NA
4:               NA

## Comparison with existing orthogonal evidence
- GSE240078 GeoMx DSP: strong compartment-level TME-vs-carcinoma support for CAF/matrix and immune/TME programs.
- GSE199102 GeoMx WTA: strong independent compartment-level concordance with GSE240078 using CAF+Immune versus epithelial segments.
- GSE310352 CosMx: robust cell-level CAF/matrix-associated TGF/EMT stromal-interface organization; not direct SPP1-CD44 support.
- Cho rapid-autopsy IMC: protein-level spatial-neighborhood/context support as summarized above; interpretation depends on ROI/patient support and robustness.
- Prior vascular-niche IMC Zenodo 10246315: previously logged as weak/near-random adjacency and not suitable for ED10 main support.

## Claim boundary
Cho IMC can support protein-level neighborhood or spatial-protein context only. It does not establish causality, direct SPP1-CD44 validation, Visium gradient reconstruction, lymph-node immune uncoupling, clinical prediction or mature TLS validation.
