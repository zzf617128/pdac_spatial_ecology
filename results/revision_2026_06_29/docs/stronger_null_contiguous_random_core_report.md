# Spatially contiguous random-core null

n_perm per sample: 1000
datasets: mvp, gse272362, gse235315, gse274557

This analysis tests whether CAF-core gradients exceed random contiguous tissue regions of the same size.
Negative observed-minus-null delta indicates stronger enrichment near the observed CAF-myeloid core than near a matched contiguous random core.

## Summary

- GSE235315 primary_tumor IFN/MHC: median delta -0.188, support 6/7.
- GSE235315 primary_tumor SPP1/TAM: median delta -0.336, support 6/7.
- GSE235315 primary_tumor TGF-beta/EMT: median delta -0.355, support 6/7.
- GSE235315 primary_tumor immune-core: median delta -0.195, support 5/7.
- GSE235315 primary_tumor tumor-aggressive: median delta -0.162, support 5/7.
- GSE272362 liver_metastasis IFN/MHC: median delta -0.325, support 12/12.
- GSE272362 liver_metastasis SPP1/TAM: median delta -0.276, support 11/12.
- GSE272362 liver_metastasis TGF-beta/EMT: median delta -0.393, support 12/12.
- GSE272362 liver_metastasis immune-core: median delta -0.296, support 12/12.
- GSE272362 liver_metastasis tumor-aggressive: median delta -0.297, support 11/12.
- GSE272362 lymph_node_metastasis IFN/MHC: median delta 0.003, support 2/5.
- GSE272362 lymph_node_metastasis SPP1/TAM: median delta -0.187, support 5/5.
- GSE272362 lymph_node_metastasis TGF-beta/EMT: median delta -0.426, support 5/5.
- GSE272362 lymph_node_metastasis immune-core: median delta 0.110, support 2/5.
- GSE272362 lymph_node_metastasis tumor-aggressive: median delta -0.218, support 4/5.
- GSE272362 normal_pancreas IFN/MHC: median delta -0.236, support 3/3.
- GSE272362 normal_pancreas SPP1/TAM: median delta -0.175, support 3/3.
- GSE272362 normal_pancreas TGF-beta/EMT: median delta -0.185, support 3/3.
- GSE272362 normal_pancreas immune-core: median delta -0.241, support 3/3.
- GSE272362 normal_pancreas tumor-aggressive: median delta -0.096, support 3/3.
- GSE272362 primary_tumor IFN/MHC: median delta -0.165, support 9/10.
- GSE272362 primary_tumor SPP1/TAM: median delta -0.272, support 10/10.
- GSE272362 primary_tumor TGF-beta/EMT: median delta -0.409, support 10/10.
- GSE272362 primary_tumor immune-core: median delta -0.147, support 9/10.
- GSE272362 primary_tumor tumor-aggressive: median delta -0.179, support 8/10.
- GSE274103 metadata_required IFN/MHC: median delta -0.237, support 5/5.
- GSE274103 metadata_required SPP1/TAM: median delta -0.281, support 5/5.
- GSE274103 metadata_required TGF-beta/EMT: median delta -0.317, support 5/5.
- GSE274103 metadata_required immune-core: median delta -0.236, support 5/5.
- GSE274103 metadata_required tumor-aggressive: median delta -0.167, support 4/5.
- GSE274557 liver_metastasis IFN/MHC: median delta -0.269, support 15/16.
- GSE274557 liver_metastasis SPP1/TAM: median delta -0.316, support 16/16.
- GSE274557 liver_metastasis TGF-beta/EMT: median delta -0.404, support 16/16.
- GSE274557 liver_metastasis immune-core: median delta -0.283, support 15/16.
- GSE274557 liver_metastasis tumor-aggressive: median delta -0.167, support 10/16.
- GSE274557 lung_metastasis IFN/MHC: median delta -0.399, support 5/6.
- GSE274557 lung_metastasis SPP1/TAM: median delta -0.376, support 6/6.
- GSE274557 lung_metastasis TGF-beta/EMT: median delta -0.301, support 6/6.
- GSE274557 lung_metastasis immune-core: median delta -0.391, support 6/6.
- GSE274557 lung_metastasis tumor-aggressive: median delta -0.017, support 4/6.
- GSE274557 peritoneal_metastasis IFN/MHC: median delta -0.221, support 14/14.
- GSE274557 peritoneal_metastasis SPP1/TAM: median delta -0.349, support 14/14.
- GSE274557 peritoneal_metastasis TGF-beta/EMT: median delta -0.346, support 13/14.
- GSE274557 peritoneal_metastasis immune-core: median delta -0.277, support 14/14.
- GSE274557 peritoneal_metastasis tumor-aggressive: median delta -0.202, support 11/14.
- GSE274557 primary_tumor IFN/MHC: median delta -0.246, support 19/19.
- GSE274557 primary_tumor SPP1/TAM: median delta -0.280, support 19/19.
- GSE274557 primary_tumor TGF-beta/EMT: median delta -0.338, support 19/19.
- GSE274557 primary_tumor immune-core: median delta -0.292, support 19/19.
- GSE274557 primary_tumor tumor-aggressive: median delta -0.058, support 13/19.
- GSE282302 metadata_required IFN/MHC: median delta -0.261, support 92/108.
- GSE282302 metadata_required SPP1/TAM: median delta -0.368, support 102/108.
- GSE282302 metadata_required TGF-beta/EMT: median delta -0.526, support 104/108.
- GSE282302 metadata_required immune-core: median delta -0.175, support 89/108.
- GSE282302 metadata_required tumor-aggressive: median delta -0.275, support 91/108.

## Claim boundary

This is still an observational spatial null model. It controls for arbitrary contiguous tissue regions but does not establish causal CAF-to-immune or CAF-to-tumor signaling.
