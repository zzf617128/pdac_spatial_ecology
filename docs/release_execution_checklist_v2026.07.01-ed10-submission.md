# Release execution checklist: v2026.07.01-ed10-submission

## 1. Secret safety

- Revoke the previously exposed Zenodo token.
- Generate a fresh Zenodo token with `deposit:write` and, only if publishing by script, `deposit:actions`.
- Store the fresh token only in the local environment:

```powershell
$env:ZENODO_TOKEN="new_token_here"
```

Do not commit tokens or paste them into notes.

## 2. GitHub code release

The GitHub release should contain code, scripts, README/release notes and lightweight ED10 provenance. Do not commit the large v3 submission zip or large base lock.

Suggested commands:

```powershell
git add README.md .gitignore `
  docs/release_v2026.07.01-ed10-submission.md `
  docs/availability_wording_v2026.07.01-ed10-submission.md `
  docs/zenodo_metadata_v2026.07.01-ed10-submission.json `
  docs/release_execution_checklist_v2026.07.01-ed10-submission.md `
  scripts/create_submission_package_v3.py `
  scripts/upload_zenodo_submission_package.py `
  reproducibility_lock_2026_06_30_with_ED10_v1

git commit -m "Prepare ED10 submission release"
git tag -a v2026.07.01-ed10-submission -m "PDAC spatial ecology ED10 submission package"
git push origin main
git push origin v2026.07.01-ed10-submission
gh release create v2026.07.01-ed10-submission `
  --repo zzf617128/pdac_spatial_ecology `
  --title "PDAC spatial ecology ED10 submission package (2026-07-01)" `
  --notes-file docs/release_v2026.07.01-ed10-submission.md
```

## 3. Zenodo data archive

Create a draft deposition and upload the submission-ready package:

```powershell
py -3 scripts/upload_zenodo_submission_package.py
```

Default uploaded file:

- `results/submission_package_v3_2026_07_01_submission_ready_only.zip`

Default behavior:

- creates a Zenodo draft;
- uploads the package;
- applies metadata from `docs/zenodo_metadata_v2026.07.01-ed10-submission.json`;
- leaves the deposition unpublished for manual review.

Publish only after checking the Zenodo draft:

```powershell
py -3 scripts/upload_zenodo_submission_package.py --publish
```

If using the sandbox first:

```powershell
py -3 scripts/upload_zenodo_submission_package.py --sandbox
```

## 4. Manuscript availability update

After Zenodo publication, replace `[DOI to be added]` in the manuscript and release notes with the final DOI.

Recommended wording is in:

- `docs/availability_wording_v2026.07.01-ed10-submission.md`

## 5. Final checks

- Confirm GitHub release page is public.
- Confirm Zenodo DOI resolves.
- Confirm uploaded archive SHA256 matches:

```text
51adc30a844f780c8a94d7bb49389786f2b7408150be2e3d6b1cda4e1e4e3bc0
```

- Confirm no patient/specimen-level wording is used for GSE310352.
- Confirm Cho IMC remains source-only and is not described as part of ED10 v1.

