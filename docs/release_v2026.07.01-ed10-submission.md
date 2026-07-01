# Release notes: v2026.07.01-ed10-submission

## Title

PDAC spatial ecology ED10 submission package (2026-07-01)

## Summary

This release fixes the code, figure-generation scripts, reproducibility manifests and submission package references for the ED10-enhanced manuscript candidate.

## Manuscript-facing additions

- Adds Extended Data Fig. 10 as an orthogonal validation figure.
- Integrates public GeoMx compartment-level validation from GSE240078 and GSE199102.
- Integrates public CosMx GSE310352 slide/FOV-level support for CAF/matrix-associated TGF/EMT stromal-interface organization.
- Preserves the original base reproducibility lock at `reproducibility_lock_2026_06_30/`.
- Adds the enhanced ED10 reproducibility lock at `reproducibility_lock_2026_06_30_with_ED10_v1/`.

## Submission and source-data package

Current package:

- `results/submission_package_v3_2026_07_01/`

Zenodo source-data archive:

- `results/PDAC_spatial_ecology_source_data.zip`
- `https://doi.org/10.5281/zenodo.21092084`

SHA256:

```text
7e2b08e368c7008610dbae51de743fa379af3b4d742a92dfc4b2a70d896ca215  results/PDAC_spatial_ecology_source_data.zip
```

Full internal package:

- `results/submission_package_v3_2026_07_01.zip`

SHA256:

```text
b2fc0518dd74bbe7348da09a65a44fc4e404e2cd53b3704bacf4e4fbc9ea0253  results/submission_package_v3_2026_07_01.zip
```

## Claim boundaries

- GeoMx datasets support compartment-level CAF/matrix and immune/TME programs.
- GSE310352 CosMx supports slide/FOV-level CAF/matrix-associated TGF/EMT stromal-interface organization.
- GSE310352 is not interpreted as patient-level or specimen-level validation.
- GSE310352 cell states are rule-based because public processed files lacked author cell-type annotations.
- TGF/EMT is interpreted as a stromal-interface state, not tumor-intrinsic EMT.
- No causal signaling is claimed.
- No direct SPP1-CD44 validation is claimed.
- No Visium distance-gradient reconstruction is claimed from ED10.
- Cho IMC remains source-only and is not included in ED10 v1.

## Recommended GitHub release command

Run after committing and tagging the release state:

```powershell
gh release create v2026.07.01-ed10-submission `
  --repo zzf617128/pdac_spatial_ecology `
  --title "PDAC spatial ecology ED10 submission package (2026-07-01)" `
  --notes-file docs/release_v2026.07.01-ed10-submission.md
```

The source-data archive should be deposited on Zenodo for DOI-backed data availability rather than committed directly to the Git repository.

Published Zenodo source-data DOI:

- `10.5281/zenodo.21092084`
- `https://zenodo.org/records/21092084`
