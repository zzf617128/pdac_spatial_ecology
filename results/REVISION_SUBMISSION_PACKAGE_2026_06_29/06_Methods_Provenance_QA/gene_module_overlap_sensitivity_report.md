# Gene-module overlap sensitivity

n_perm per sample: 1000

This analysis quantifies marker overlap between CAF-myeloid core genes and target modules.
For target modules without shared genes, the original target score is unchanged and overlap sensitivity is equivalent to the core-gradient analysis.
For component modules that become empty after removing shared genes, the target-removed score is not estimable; instead, an anchor-side sensitivity was run by excluding the overlapping component from the CAF-myeloid core score where possible.

## CAF-myeloid overlap summary

- CAF-myeloid vs IFN/MHC antigen-presentation: shared 0/9 target genes (Jaccard 0.000); shared genes: none.
- CAF-myeloid vs SPP1/TAM: shared 6/6 target genes (Jaccard 0.240); shared genes: APOE;CCL18;CTSB;LGALS3;SPP1;TREM2.
- CAF-myeloid vs TGF-beta/EMT: shared 0/13 target genes (Jaccard 0.000); shared genes: none.
- CAF-myeloid vs immune-core: shared 0/32 target genes (Jaccard 0.000); shared genes: none.
- CAF-myeloid vs myCAF/matrix: shared 11/11 target genes (Jaccard 0.440); shared genes: ACTA2;COL1A1;COL1A2;COL3A1;DCN;FAP;LUM;MYL9;PDPN;TAGLN;THY1.
- CAF-myeloid vs tumor-aggressive: shared 0/20 target genes (Jaccard 0.000); shared genes: none.

## Sensitivity summary

- IFN/MHC: mode=no_shared_genes_original_core; status=target_score_unchanged_no_overlap; median delta -0.285; support 184/205.
- SPP1/TAM: mode=core_excluding_spp1_tam_component; status=target_empty_after_overlap_removal; median delta -0.237; support 187/205.
- TGF-beta/EMT: mode=no_shared_genes_original_core; status=target_score_unchanged_no_overlap; median delta -0.478; support 203/205.
- immune-core: mode=no_shared_genes_original_core; status=target_score_unchanged_no_overlap; median delta -0.261; support 177/205.
- myCAF/matrix: mode=core_excluding_caf_matrix_component; status=target_empty_after_overlap_removal; median delta -0.253; support 174/205.
- tumor-aggressive: mode=no_shared_genes_original_core; status=target_score_unchanged_no_overlap; median delta -0.277; support 175/205.

## Interpretation

The main IFN/MHC, immune-core and tumor-aggressive target modules do not share marker genes with the current CAF-myeloid core definition, so their CAF-core gradients are not explained by direct marker overlap.
SPP1/TAM and myCAF/matrix are core-component interpretation modules. They should be described as components of the CAF-myeloid architecture rather than independent non-overlapping targets.
