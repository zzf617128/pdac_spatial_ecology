from __future__ import annotations

from pathlib import Path
import csv
import shutil
import zipfile


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "nature_subjournal_submission"

DIRS = {
    "manuscript": OUT / "manuscript",
    "figures": OUT / "figures",
    "source_data": OUT / "source_data",
    "supplementary_information": OUT / "supplementary_information",
    "reports": OUT / "reports",
    "code": OUT / "code",
}


def copy(src: Path, dst: Path, role: str, rows: list[dict[str, str]]) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    rows.append({"role": role, "source": str(src.relative_to(ROOT)), "staged_file": str(dst.relative_to(ROOT))})


def main() -> None:
    for d in DIRS.values():
        d.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, str]] = []

    copy(
        ROOT / "results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal.docx",
        DIRS["manuscript"] / "Manuscript_NatureSubjournal.docx",
        "main manuscript Word",
        rows,
    )
    copy(
        ROOT / "results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal.md",
        DIRS["manuscript"] / "Manuscript_NatureSubjournal.md",
        "main manuscript Markdown provenance",
        rows,
    )
    copy(
        ROOT / "results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal_cover_letter_draft.md",
        DIRS["manuscript"] / "Cover_Letter_Draft.md",
        "cover letter draft",
        rows,
    )
    copy(
        ROOT / "results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.docx",
        DIRS["supplementary_information"] / "Supplementary_Information.docx",
        "supplementary information Word",
        rows,
    )
    copy(
        ROOT / "results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.md",
        DIRS["supplementary_information"] / "Supplementary_Information.md",
        "supplementary information Markdown provenance",
        rows,
    )

    figure_map = [
        ("Figure_1", "figure1_submission_spatial_specificity_nc_style"),
        ("Figure_2", "figure2_submission_metastatic_decoupling_nc_style"),
        ("Figure_3", "figure3_submission_ecotypes_mechanism_axes_nc_style"),
        ("Figure_4", "figure4_submission_multiresolution_validation_nc_style"),
        ("Extended_Data_Figure_1", "supplementary_module1_spatial_specificity_robustness"),
        ("Extended_Data_Figure_2", "supplementary_module2_metastatic_immune_decoupling"),
        ("Extended_Data_Figure_3", "supplementary_module3_cell_state_multiresolution_validation"),
        ("Extended_Data_Figure_4", "supplementary_module4_mechanism_interface_priority"),
        ("Extended_Data_Figure_5", "supplementary_module5_pathology_tcga_tls_boundaries"),
        ("Extended_Data_Figure_6", "supplementary_module6_spatial_architecture_mechanism_deepening"),
    ]
    for staged_base, source_base in figure_map:
        for ext in [".pdf", ".svg", ".png"]:
            copy(
                ROOT / "results/figures/submission" / f"{source_base}{ext}",
                DIRS["figures"] / f"{staged_base}{ext}",
                "display figure",
                rows,
            )

    source_files = [
        "Source_Data_Fig_1.csv",
        "Source_Data_Fig_2.csv",
        "Source_Data_Fig_3A.csv",
        "Source_Data_Fig_3_candidate_NC_style_panel_index.csv",
        "Source_Data_Fig_4A_multiresolution_scale.csv",
        "Source_Data_Fig_4B_GSE274557.csv",
        "Source_Data_Fig_4C_GSE274673.csv",
        "Source_Data_Fig_4D_GSE274673.csv",
        "Source_Data_supplementary_module1_spatial_specificity_robustness.csv",
        "Source_Data_supplementary_module2_metastatic_immune_decoupling.csv",
        "Source_Data_supplementary_module3_cell_state_multiresolution_validation.csv",
        "Source_Data_supplementary_module4_mechanism_interface_priority.csv",
        "Source_Data_supplementary_module5_pathology_tcga_tls_boundaries.csv",
        "Source_Data_supplementary_module6_spatial_architecture_mechanism_deepening.csv",
    ]
    for name in source_files:
        copy(ROOT / "results/source_data" / name, DIRS["source_data"] / name, "source data", rows)

    copy(
        ROOT / "results/tables/nature_subjournal_display_item_renumbering.csv",
        DIRS["reports"] / "nature_subjournal_display_item_renumbering.csv",
        "display item renumbering",
        rows,
    )
    copy(
        ROOT / "results/reports/nature_subjournal_submission_readiness_2026_06_28.md",
        DIRS["reports"] / "nature_subjournal_submission_readiness_2026_06_28.md",
        "submission readiness report",
        rows,
    )
    copy(
        ROOT / "results/reports/nature_subjournal_reporting_summary_draft_2026_06_28.md",
        DIRS["reports"] / "nature_subjournal_reporting_summary_draft_2026_06_28.md",
        "reporting summary draft",
        rows,
    )
    copy(
        ROOT / "results/reports/nature_subjournal_editorial_highlights_2026_06_28.md",
        DIRS["reports"] / "nature_subjournal_editorial_highlights_2026_06_28.md",
        "editorial highlights",
        rows,
    )
    copy(
        ROOT / "results/reports/final_claim_language_audit_2026_06_28.md",
        DIRS["reports"] / "final_claim_language_audit_2026_06_28.md",
        "claim-language audit",
        rows,
    )
    copy(
        ROOT / "results/reports/complete_workflow_audit_and_optimization_plan_2026_06_28.md",
        DIRS["reports"] / "complete_workflow_audit_and_optimization_plan_2026_06_28.md",
        "complete workflow audit",
        rows,
    )

    code_readme = ROOT / "results/code_availability/README_code_availability.md"
    code_zip = ROOT / "results/code_availability/pdac_spatial_ecology_code_archive_preview.zip"
    if code_readme.exists():
        copy(code_readme, DIRS["code"] / code_readme.name, "code availability README", rows)
    if code_zip.exists():
        copy(code_zip, DIRS["code"] / code_zip.name, "code availability archive preview", rows)

    zip_path = OUT / "source_data_nature_subjournal.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in sorted(DIRS["source_data"].glob("*")):
            zf.write(path, arcname=path.name)
    rows.append({"role": "source data zip", "source": str(DIRS["source_data"].relative_to(ROOT)), "staged_file": str(zip_path.relative_to(ROOT))})

    manifest = OUT / "upload_manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["role", "source", "staged_file"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {manifest}")
    print(f"Wrote {zip_path}")


if __name__ == "__main__":
    main()
