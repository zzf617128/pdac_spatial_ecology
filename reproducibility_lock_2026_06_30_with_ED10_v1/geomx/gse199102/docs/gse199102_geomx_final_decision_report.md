# GSE199102 GeoMx Final Decision Report

Decision: A. Strong ED support

## Summary
Key: <module>
                module n_patients median_stromal_minus_tumor
                <char>      <int>                      <num>
 1:           b_plasma         20                 0.12224771
 2:          caf_mycaf         20                 1.87429470
 3:            ifn_apc         20                -0.01538592
 4:        immune_core         20                 0.17750256
 5: myeloid_macrophage         20                 0.17971336
 6:      pancaf_matrix         20                 1.72813537
 7:           spp1_tam         20                -0.36218014
 8:             t_cell         20                 0.37100246
 9:           tgfb_emt         20                 0.70679046
10:   tumor_aggressive         20                -0.97821966
11:   tumor_epithelial         20                -1.05164005
    patient_support_fraction_stromal_gt_tumor median_immune_minus_tumor
                                        <num>                     <num>
 1:                                      0.70                 0.1458852
 2:                                      1.00                 0.3223619
 3:                                      0.40                 0.6759012
 4:                                      0.95                 0.6036377
 5:                                      0.90                 1.1589378
 6:                                      1.00                 0.4977632
 7:                                      0.00                 0.1020123
 8:                                      0.95                 1.0379396
 9:                                      0.90                 0.2856156
10:                                      0.00                -0.9969362
11:                                      0.00                -1.1128108
    patient_support_fraction_immune_gt_tumor
                                       <num>
 1:                                     0.65
 2:                                     1.00
 3:                                     0.95
 4:                                     1.00
 5:                                     1.00
 6:                                     1.00
 7:                                     0.50
 8:                                     1.00
 9:                                     0.70
10:                                     0.00
11:                                     0.00
    median_fibroblast_minus_nonfibroblast
                                    <num>
 1:                            0.06348929
 2:                            1.58255255
 3:                           -0.16804846
 4:                            0.05639317
 5:                           -0.12006101
 6:                            1.35906492
 7:                           -0.36262298
 8:                            0.13615271
 9:                            0.63822666
10:                           -0.35815952
11:                           -0.45344678
    patient_support_fraction_fibroblast_gt_nonfibroblast median_tme_minus_tumor
                                                   <num>                  <num>
 1:                                                 0.60              0.1184957
 2:                                                 1.00              0.8075536
 3:                                                 0.30              0.2262469
 4:                                                 0.60              0.3279197
 5:                                                 0.25              0.4779628
 6:                                                 1.00              0.8235727
 7:                                                 0.00             -0.2775013
 8:                                                 0.70              0.5776957
 9:                                                 0.90              0.4443280
10:                                                 0.00             -0.9752634
11:                                                 0.00             -1.0934998
    patient_support_fraction_tme_gt_tumor expected_stromal_direction
                                    <num>                     <char>
 1:                                  0.70                immune_high
 2:                                  1.00               stromal_high
 3:                                  0.80     immune_or_stromal_high
 4:                                  1.00     immune_or_stromal_high
 5:                                  1.00     stromal_or_immune_high
 6:                                  1.00               stromal_high
 7:                                  0.10          report_separately
 8:                                  1.00                immune_high
 9:                                  0.85               stromal_high
10:                                  0.00                 tumor_high
11:                                  0.00                 tumor_high
    control_pass
          <lgcl>
 1:         TRUE
 2:         TRUE
 3:         TRUE
 4:         TRUE
 5:         TRUE
 6:         TRUE
 7:           NA
 8:         TRUE
 9:         TRUE
10:         TRUE
11:         TRUE

## Interpretation
GSE199102 provides independent GeoMx WTA compartment-level testing using CAF, Immune and Epithelial segment labels. Interpret only at compartment level, not as cell-cell spatial proximity.

## Claim Boundaries
- No causality.
- No direct SPP1-CD44 validation.
- No Visium gradient reconstruction.
- No LN immune uncoupling.
