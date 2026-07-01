# How To Rebuild ED10 v1

This enhanced lock archives the ED10 v1 figure/source outputs and the scripts needed to regenerate or audit the components available in this workspace.

## Inputs
- GSE240078 GeoMx DSP quantities are preserved in `ed10/source_data/Source_Data_Extended_Data_Figure_10_v1.csv` and concordance outputs.
- GSE199102 GeoMx WTA source tables and script are archived under `geomx/gse199102/` and `scripts/analyze_gse199102_geomx.R`.
- GSE310352 CosMx source/robustness/transparency outputs are archived under `gse310352/` with scripts under `scripts/`.
- Cho IMC is archived under `cho_imc_source_only/` but is not part of ED10 v1.

## Rebuild Steps
1. Re-run upstream public-dataset scripts only if the public input data are available locally.
2. Re-run GSE199102: `Rscript scripts/analyze_gse199102_geomx.R` from the project root.
3. Re-run GSE310352 initial/robustness/transparency scripts from the project root if needed.
4. Re-run ED10 figure generation with `scripts/make_extended_data_fig10_v1.R` if present and inputs are available.
5. Compare outputs to `checksums_sha256_enhanced_lock.txt`.

## Claim Boundaries
Keep GSE310352 slide/FOV-level, rule-based and stromal-interface. Do not convert it into patient-level validation, tumor-intrinsic EMT or causal signaling.
