# PDAC CAF-Myeloid Spatial Ecology

This repository contains analysis code, figure-generation scripts, processed source tables and reproducibility manifests for a public-data spatial transcriptomics study of CAF-myeloid stromal architecture in pancreatic ductal adenocarcinoma.

## Current Submission Candidate

The current manuscript-facing package is:

- `results/submission_package_v3_2026_07_01/`

The submission-ready subset is:

- `results/submission_package_v3_2026_07_01/submission_ready/`

The Zenodo-ready source-data archive is:

- `results/PDAC_spatial_ecology_source_data.zip`

The enhanced reproducibility lock is:

- `reproducibility_lock_2026_06_30_with_ED10_v1/`

This version adds Extended Data Fig. 10 as an orthogonal validation figure while preserving the original frozen base lock:

- `reproducibility_lock_2026_06_30/`

## Planned Release

Recommended release tag:

- `v2026.07.01-ed10-submission`

Recommended GitHub release title:

- `PDAC spatial ecology ED10 submission package (2026-07-01)`

The release should point to the exact commit used for submission and should be paired with a Zenodo archive for the submission source-data package.

## ED10 Scope

Extended Data Fig. 10 provides orthogonal support for CAF/matrix compartment and TGF/EMT stromal-interface organization:

- GSE240078 GeoMx: compartment-level support.
- GSE199102 GeoMx: independent compartment-level replication.
- GSE310352 CosMx: slide/FOV-level CAF/matrix-associated TGF/EMT stromal-interface support.
- Cho IMC: source-only archive, not included in ED10 v1.

Claim boundaries:

- No causal signaling is claimed.
- No direct SPP1-CD44 validation is claimed.
- No tumor-intrinsic EMT claim is made from GSE310352.
- No Visium distance-gradient reconstruction is claimed from ED10.
- GSE310352 is interpreted at slide/FOV level because public metadata did not allow reliable patient/specimen recovery for processed slides.

## Code Availability

Analysis scripts, figure-generation code and reproducibility manifests are maintained in this public repository:

- `https://github.com/zzf617128/pdac_spatial_ecology`

For submission, cite the fixed release:

- `v2026.07.01-ed10-submission`

If a Zenodo DOI is minted for the GitHub release or data package, add the DOI to the manuscript availability statement before final submission.

Source-data DOI:

- `10.5281/zenodo.21092084`
- `https://zenodo.org/records/21092084`

## Data Policy

Raw public datasets are not committed to this repository because they are large and remain available from GEO, TCGA, UCSC Xena, Zenodo or linked public resources.

The source-data package contains the plotted quantitative source data and is intended to be archived separately on Zenodo:

- `results/PDAC_spatial_ecology_source_data.zip`
- `https://doi.org/10.5281/zenodo.21092084`

Large derived files should remain in the Zenodo archive or journal source-data package rather than being committed directly to the Git repository.

## Reproducibility Notes

Key v3 files:

- `results/submission_package_v3_2026_07_01/README_submission_package_v3.md`
- `results/submission_package_v3_2026_07_01/docs/submission_v3_QA_report.md`
- `results/submission_package_v3_2026_07_01/docs/local_path_audit_submission_v3.md`
- `results/submission_package_v3_2026_07_01/manifest/checksums_sha256_submission_v3.txt`
- `reproducibility_lock_2026_06_30_with_ED10_v1/checksums_sha256_enhanced_lock.txt`
