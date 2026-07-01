# GSE199102 and GSE240078 Direction Concordance Interpretation

GSE199102 was analyzed as CAF/Immune/Epithelial GeoMx WTA segment data and compared with the previous GSE240078 stroma-vs-carcinoma DSP result.

Key: <module>
                module gse199102_median_tme_minus_tumor
                <char>                            <num>
 1:           b_plasma                        0.1184957
 2:          caf_mycaf                        0.8075536
 3:            ifn_apc                        0.2262469
 4:        immune_core                        0.3279197
 5: myeloid_macrophage                        0.4779628
 6:      pancaf_matrix                        0.8235727
 7:           spp1_tam                       -0.2775013
 8:             t_cell                        0.5776957
 9:           tgfb_emt                        0.4443280
10:   tumor_aggressive                       -0.9752634
11:   tumor_epithelial                       -1.0934998
    gse199102_tme_patient_support gse199102_median_stromal_minus_tumor
                            <num>                                <num>
 1:                          0.70                           0.12224771
 2:                          1.00                           1.87429470
 3:                          0.80                          -0.01538592
 4:                          1.00                           0.17750256
 5:                          1.00                           0.17971336
 6:                          1.00                           1.72813537
 7:                          0.10                          -0.36218014
 8:                          1.00                           0.37100246
 9:                          0.85                           0.70679046
10:                          0.00                          -0.97821966
11:                          0.00                          -1.05164005
    gse199102_median_immune_minus_tumor gse240078_delta_stroma_minus_tumor
                                  <num>                              <num>
 1:                           0.1458852                                 NA
 2:                           0.3223619                          1.0763115
 3:                           0.6759012                          0.1660568
 4:                           0.6036377                          0.7091361
 5:                           1.1589378                          0.9262532
 6:                           0.4977632                          1.0656854
 7:                           0.1020123                         -0.3123146
 8:                           1.0379396                          0.9668998
 9:                           0.2856156                          0.3786529
10:                          -0.9969362                         -0.7448272
11:                          -1.1128108                         -1.1990072
    gse240078_fdr          gse199102_direction gse240078_direction
            <num>                       <char>              <char>
 1:            NA                     tme_high                <NA>
 2:  0.000000e+00                     tme_high        stromal_high
 3:  1.440586e-02                     tme_high        stromal_high
 4:  0.000000e+00                     tme_high        stromal_high
 5:  0.000000e+00                     tme_high        stromal_high
 6:  0.000000e+00                     tme_high        stromal_high
 7:  5.490414e-11 tumor_high_report_separately          tumor_high
 8:  0.000000e+00                     tme_high        stromal_high
 9:  8.912315e-13                     tme_high        stromal_high
10:  1.130179e-26                   tumor_high          tumor_high
11:  5.429516e-35                   tumor_high          tumor_high
    direction_concordant
                  <lgcl>
 1:                   NA
 2:                 TRUE
 3:                 TRUE
 4:                 TRUE
 5:                 TRUE
 6:                 TRUE
 7:                   NA
 8:                 TRUE
 9:                 TRUE
10:                 TRUE
11:                 TRUE

SPP1/TAM is reported separately and should not be forced into the stromal-high group.
