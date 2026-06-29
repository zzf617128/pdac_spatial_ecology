# GSE272362 Provenance Audit

Last updated: 2026-06-24

## External Sources Checked

1. GEO accession page: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE272362`
   - Status: public.
   - Title: "Spatial Transcriptomic Analysis Of Primary And Metastatic Pancreatic Cancers Highlights Tumor Microenvironmental Heterogeneity".
   - GEO overall design states an atlas of 30 specimens: 10 primary PDAC tumors, 3 normal pancreata, 12 matched hepatic metastases, and 5 matched lymph node metastases, totaling 91,496 spot transcriptomes.
   - GEO links the series to Nature Genetics article `s41588-024-01914-4`.

2. Nature Genetics article: `https://www.nature.com/articles/s41588-024-01914-4`
   - Article title: "Spatial transcriptomic analysis of primary and metastatic pancreatic cancers highlights tumor microenvironmental heterogeneity".
   - Published: 18 September 2024.
   - Data availability states that ST and scRNA-seq analysis data are available through Zenodo record `10712047`, and generated ST data are deposited in GEO as `GSE272362`.

3. Zenodo record: `https://zenodo.org/records/10712047`
   - Record title: "Spatially Resolved Transcriptomics Atlas of Matched Primary and Metastatic Pancreatic Cancer Reveal Principles of Ecological Adaptation".
   - DOI: `10.5281/zenodo.10712047`.
   - File used in this project: `PDAC_Updated.rds`.
   - Zenodo page lists `PDAC_Updated.rds` with md5 `d0f0b12e0fb013f3def1a62d0f925cbf`.

## Local Data Checked

Local RDS-derived metadata table:

- `metadata/gse272362_rds_sample_mapping.csv`

Local scored spot table:

- `results/tables/gse272362_rds_spot_level_scores.csv`

Observed local sample groups:

| group | samples | spots |
|---|---:|---:|
| primary_tumor | 10 | 35,458 |
| liver_metastasis | 12 | 28,520 |
| lymph_node_metastasis | 5 | 17,698 |
| normal_pancreas | 3 | 9,820 |
| total | 30 | 91,496 |

These local totals match the GEO overall design at the level of specimen counts and total spot transcriptomes.

## Origin Field Interpretation

The local RDS metadata column `Origin` contains the tissue-origin labels:

- `Pancreas` -> `primary_tumor`
- `Liver` -> `liver_metastasis`
- `Lymph node` -> `lymph_node_metastasis`
- `Normal Pancreas` -> `normal_pancreas`

This mapping is consistent with the GEO overall design count structure:

- 10 primary PDAC tumors correspond to local `Origin == "Pancreas"`.
- 12 hepatic metastases correspond to local `Origin == "Liver"`.
- 5 lymph node metastases correspond to local `Origin == "Lymph node"`.
- 3 normal pancreata correspond to local `Origin == "Normal Pancreas"`.

## Coordinate-Matching Evidence

Representative overlay manifest:

- `results/tables/gse272362_rds_overlay_manifest.csv`

Selected overlays:

- `IU_PDA_T3`: 400/400 exact RDS-to-GEO coordinate matches.
- `IU_PDA_T1`: 400/400 exact RDS-to-GEO coordinate matches.
- `IU_PDA_HM10`: 400/400 exact RDS-to-GEO coordinate matches.
- `IU_PDA_HM5`: 400/400 exact RDS-to-GEO coordinate matches.
- `IU_PDA_LNM12`: 400/400 exact RDS-to-GEO coordinate matches.
- `IU_PDA_LNM7`: 400/400 exact RDS-to-GEO coordinate matches.

## Manuscript-Ready Statement

GSE272362 was used as an independent spatial transcriptomics validation cohort. The GEO overall design describes 30 specimens comprising 10 primary PDAC tumors, three normal pancreata, 12 matched hepatic metastases, and five matched lymph node metastases, totaling 91,496 spot transcriptomes. We used the Zenodo `PDAC_Updated.rds` object linked from the associated Nature Genetics article and mapped its `Origin` metadata field to tissue-site labels. Local sample and spot counts exactly matched the GEO overall design, and selected H&E overlays were generated after exact coordinate matching between RDS spot coordinates and GEO tissue-position files.

## Residual Caveat

This audit confirms tissue-site provenance at the specimen-group level. It does not by itself establish clinical outcome labels, treatment response, survival, or any patient-level claim.
