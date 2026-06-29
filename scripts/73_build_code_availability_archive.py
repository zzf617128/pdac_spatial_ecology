from __future__ import annotations

from datetime import date
from pathlib import Path
import zipfile


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "code_availability"
ZIP_PATH = OUT / "pdac_spatial_ecology_code_archive_preview.zip"
README = OUT / "README_code_availability.md"


def iter_archive_files() -> list[Path]:
    files: list[Path] = []
    for pattern in ["*.py", "*.R"]:
        files.extend(sorted((ROOT / "scripts").glob(pattern)))

    report_names = [
        "references_versions_reproducibility.md",
        "submission_package_index.md",
        "submission_figure_captions_and_source_map.md",
        "submission_claim_evidence_matrix.md",
        "final_story_package.md",
        "final_claim_language_audit_2026_06_28.md",
        "complete_workflow_audit_and_optimization_plan_2026_06_28.md",
        "nature_subjournal_editorial_highlights_2026_06_28.md",
        "nature_subjournal_submission_readiness_2026_06_28.md",
        "nature_subjournal_reporting_summary_draft_2026_06_28.md",
        "strict_single_cell_deconvolution_decision_2026_06_28.md",
        "strict_nnls_reference_deconvolution_report.md",
        "supplementary_module_manuscript_integration_audit_2026_06_28.md",
        "supplementary_module6_spatial_architecture_mechanism_deepening_notes.md",
    ]
    for name in report_names:
        p = ROOT / "results" / "reports" / name
        if p.exists():
            files.append(p)

    table_names = [
        "submission_final_figure_manifest.csv",
        "supplement_module_reorganization_manifest.csv",
        "supplementary_modules_1_5_export_manifest.csv",
        "nature_subjournal_display_item_renumbering.csv",
        "submission_cohort_summary.csv",
    ]
    for name in table_names:
        p = ROOT / "results" / "tables" / name
        if p.exists():
            files.append(p)

    metadata_dir = ROOT / "metadata"
    if metadata_dir.exists():
        files.extend(sorted(metadata_dir.glob("*.md")))

    return sorted(set(files), key=lambda p: p.as_posix())


def write_readme(archive_files: list[Path]) -> None:
    script_count = sum(1 for p in archive_files if p.parent.name == "scripts")
    readme = f"""# PDAC Spatial Ecology Code Availability Preview

Last updated: {date.today().isoformat()}

This folder contains a reviewer-ready preview archive of custom analysis code and reproducibility metadata for the manuscript:

**CAF-myeloid stromal cores organize inflammatory and tumor-aggressive spatial programs in pancreatic cancer**

## Contents

- Code archive: `pdac_spatial_ecology_code_archive_preview.zip`
- Scripts included: {script_count}
- Reports, manifest tables and metadata files included: {len(archive_files) - script_count}

## Scope

The archive contains custom Python and R scripts under `scripts/`, selected reproducibility reports, figure/source-data maps, package indexes, cohort summaries and provenance notes.

The archive does not contain raw public expression matrices, large image files, Visium/Xenium raw outputs or generated figure/source-data files. Those are handled through public accessions and the separate staged Source Data package.

## Submission Use

For initial submission, this preview archive can support reviewer access if a journal allows code upload through the manuscript system. For final publication, deposit this archive or an equivalent repository snapshot in a persistent repository that mints a DOI, such as Zenodo or Code Ocean, and update the manuscript Code availability section with the DOI.

## Reproducibility Pointers

- Command and software provenance: `results/reports/references_versions_reproducibility.md`
- Package index: `results/reports/submission_package_index.md`
- Figure/source-data mapping: `results/reports/submission_figure_captions_and_source_map.md`
- Nature-subjournal readiness: `results/reports/nature_subjournal_submission_readiness_2026_06_28.md`
"""
    OUT.mkdir(parents=True, exist_ok=True)
    README.write_text(readme, encoding="utf-8", newline="\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    archive_files = iter_archive_files()
    write_readme(archive_files)
    archive_files.append(README)

    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(archive_files, key=lambda p: p.as_posix()):
            zf.write(path, arcname=path.relative_to(ROOT).as_posix())

    print(f"Wrote {README}")
    print(f"Wrote {ZIP_PATH}")
    print(f"Archived files: {len(archive_files)}")


if __name__ == "__main__":
    main()
