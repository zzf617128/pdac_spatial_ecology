from __future__ import annotations

import csv
import hashlib
import os
import re
import shutil
import zipfile
from datetime import date
from pathlib import Path
from xml.etree import ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
BASE_LOCK = ROOT / "reproducibility_lock_2026_06_30"
ENHANCED_LOCK = ROOT / "reproducibility_lock_2026_06_30_with_ED10_v1"
PACKAGE = ROOT / "results" / f"submission_package_v3_{date.today().isoformat().replace('-', '_')}"
ZIP_PATH = PACKAGE.with_suffix(".zip")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def copy_file(src: Path, dst: Path, rows: list[dict[str, str]], component: str, note: str = "") -> None:
    if not src.exists():
        raise FileNotFoundError(src)
    ensure_dir(dst.parent)
    shutil.copy2(src, dst)
    rows.append(
        {
            "component": component,
            "package_path": str(dst.relative_to(PACKAGE)).replace("\\", "/"),
            "source_lock": "base" if BASE_LOCK in src.parents else "enhanced",
            "source_relative_path": str(
                src.relative_to(BASE_LOCK if BASE_LOCK in src.parents else ENHANCED_LOCK)
            ).replace("\\", "/"),
            "status": "copied",
            "note": note,
        }
    )


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    ensure_dir(path.parent)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_text(path: Path, text: str) -> None:
    ensure_dir(path.parent)
    path.write_text(text, encoding="utf-8", newline="\n")


def extract_docx_text(path: Path) -> str:
    pieces: list[str] = []
    with zipfile.ZipFile(path) as zf:
        for name in zf.namelist():
            if name.startswith("word/") and name.endswith(".xml"):
                try:
                    root = ET.fromstring(zf.read(name))
                except ET.ParseError:
                    continue
                for node in root.iter():
                    if node.tag.endswith("}t") and node.text:
                        pieces.append(node.text)
    return "\n".join(pieces)


def audit_local_paths(paths: list[Path]) -> list[dict[str, str]]:
    patterns = [
        r"E:/",
        r"E:\\",
        r"results/",
        r"results\\",
        r"reproducibility_lock_2026_06_30",
        r"orthogonal_validation_strong_search",
        r"revision_2026_06_30_with_ED10_v1",
        r"pdac_spatial_ecology",
    ]
    rows: list[dict[str, str]] = []
    for path in paths:
        if path.suffix.lower() == ".docx":
            text = extract_docx_text(path)
        elif path.suffix.lower() in {".md", ".txt", ".csv"}:
            text = path.read_text(encoding="utf-8", errors="ignore")
        else:
            continue
        for pat in patterns:
            matches = list(re.finditer(pat, text, flags=re.IGNORECASE))
            rows.append(
                {
                    "file": str(path.relative_to(PACKAGE)).replace("\\", "/"),
                    "pattern": pat,
                    "match_count": str(len(matches)),
                    "status": "pass" if not matches else "review",
                }
            )
    return rows


def main() -> None:
    if PACKAGE.exists() or ZIP_PATH.exists():
        raise SystemExit(f"Refusing to overwrite existing package or zip: {PACKAGE}")

    rows: list[dict[str, str]] = []

    submission = PACKAGE / "submission_ready"
    provenance = PACKAGE / "provenance_internal"

    # Manuscript
    copy_file(
        ENHANCED_LOCK / "manuscript" / "Manuscript_with_ED10_v1_submission_safe_slidelevel_patch.docx",
        submission / "manuscript" / "Manuscript_submission_v3_with_ED10.docx",
        rows,
        "manuscript",
        "Submission-safe ED10 v1 manuscript candidate with GSE310352 slide/FOV-level patch.",
    )

    # Main and Extended Data figures
    base_fig_dir = BASE_LOCK / "frozen_submission_candidate" / "figures"
    figure_map = {
        "Figure_1_current.pdf": "main_figures/Figure_1.pdf",
        "Figure_2_current.pdf": "main_figures/Figure_2.pdf",
        "Figure_3_current.pdf": "main_figures/Figure_3.pdf",
        "Figure_4_current.pdf": "main_figures/Figure_4.pdf",
        "Extended_Data_Figure_1_current.pdf": "extended_data_figures/Extended_Data_Figure_1.pdf",
        "Extended_Data_Figure_2_current.pdf": "extended_data_figures/Extended_Data_Figure_2.pdf",
        "Extended_Data_Figure_3_current.pdf": "extended_data_figures/Extended_Data_Figure_3.pdf",
        "Extended_Data_Figure_4_current.pdf": "extended_data_figures/Extended_Data_Figure_4.pdf",
        "Extended_Data_Figure_5_current.pdf": "extended_data_figures/Extended_Data_Figure_5.pdf",
        "Extended_Data_Figure_6_current.pdf": "extended_data_figures/Extended_Data_Figure_6.pdf",
        "Extended_Data_Figure_7_Specificity_Sensitivity.pdf": "extended_data_figures/Extended_Data_Figure_7.pdf",
        "Extended_Data_Figure_8_LN_Individual_Spatial_Maps.pdf": "extended_data_figures/Extended_Data_Figure_8.pdf",
        "Extended_Data_Figure_9_Matched_Contiguous_Null.pdf": "extended_data_figures/Extended_Data_Figure_9.pdf",
    }
    for src_name, dst_rel in figure_map.items():
        copy_file(base_fig_dir / src_name, submission / "figures" / dst_rel, rows, "figure")

    copy_file(
        ENHANCED_LOCK / "ed10" / "figures" / "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf",
        submission
        / "figures"
        / "extended_data_figures"
        / "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf",
        rows,
        "figure",
        "ED10 v1 manuscript-facing figure.",
    )
    copy_file(
        ENHANCED_LOCK / "ed10" / "figures" / "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg",
        provenance
        / "editable_figures"
        / "Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg",
        rows,
        "editable_figure",
        "Editable ED10 SVG, stored in internal provenance rather than submission-ready figures.",
    )

    # Supplementary tables
    supp_dir = BASE_LOCK / "frozen_submission_candidate" / "supplementary_tables"
    for src in sorted(supp_dir.glob("*")):
        if src.is_file():
            copy_file(src, submission / "supplementary_tables" / src.name, rows, "supplementary_table")

    # Source data package: base source data plus ED10 source data.
    source_dir = BASE_LOCK / "frozen_submission_candidate" / "source_data_package"
    for src in sorted(source_dir.glob("*")):
        if src.is_file():
            copy_file(src, submission / "source_data" / src.name, rows, "source_data")
    copy_file(
        ENHANCED_LOCK / "ed10" / "source_data" / "Source_Data_Extended_Data_Figure_10_v1.csv",
        submission / "source_data" / "Source_Data_Extended_Data_Figure_10_v1.csv",
        rows,
        "source_data",
        "Source data for Extended Data Fig. 10.",
    )

    # Internal provenance and reproducibility support.
    provenance_files = [
        (ENHANCED_LOCK / "README.md", provenance / "enhanced_lock_README.md", "provenance"),
        (
            ENHANCED_LOCK / "enhanced_claim_boundary_checklist.md",
            provenance / "claim_boundary" / "enhanced_claim_boundary_checklist.md",
            "claim_boundary",
        ),
        (
            ENHANCED_LOCK / "enhanced_unresolved_items.md",
            provenance / "review_notes" / "enhanced_unresolved_items.md",
            "review_note",
        ),
        (
            ENHANCED_LOCK / "docs" / "enhanced_submission_summary.md",
            provenance / "review_notes" / "enhanced_submission_summary.md",
            "review_note",
        ),
        (
            ENHANCED_LOCK / "docs" / "what_changed_from_base_lock.md",
            provenance / "review_notes" / "what_changed_from_base_lock.md",
            "review_note",
        ),
        (
            ENHANCED_LOCK / "docs" / "how_to_rebuild_ED10.md",
            provenance / "rebuild" / "how_to_rebuild_ED10.md",
            "rebuild_note",
        ),
        (
            ENHANCED_LOCK / "ed10" / "docs" / "ed10_v1_figure_legend_draft.md",
            provenance / "ed10_docs" / "ed10_v1_figure_legend_draft.md",
            "ed10_doc",
        ),
        (
            ENHANCED_LOCK / "ed10" / "docs" / "ed10_v1_methods_draft.md",
            provenance / "ed10_docs" / "ed10_v1_methods_draft.md",
            "ed10_doc",
        ),
        (
            ENHANCED_LOCK / "ed10" / "docs" / "ed10_v1_results_paragraph_draft.md",
            provenance / "ed10_docs" / "ed10_v1_results_paragraph_draft.md",
            "ed10_doc",
        ),
        (
            ENHANCED_LOCK / "ed10" / "docs" / "ed10_v1_claim_boundary_notes.md",
            provenance / "ed10_docs" / "ed10_v1_claim_boundary_notes.md",
            "ed10_doc",
        ),
        (
            ENHANCED_LOCK / "gse310352" / "docs" / "gse310352_patient_mapping_recovery_report.md",
            provenance / "gse310352" / "gse310352_patient_mapping_recovery_report.md",
            "gse310352_boundary",
        ),
        (
            ENHANCED_LOCK / "gse310352" / "docs" / "gse310352_cell_state_definition_transparency_report.md",
            provenance / "gse310352" / "gse310352_cell_state_definition_transparency_report.md",
            "gse310352_qc",
        ),
        (
            ENHANCED_LOCK / "strong_search" / "docs" / "final_orthogonal_validation_decision_summary.md",
            provenance / "strong_search" / "final_orthogonal_validation_decision_summary.md",
            "decision_summary",
        ),
        (
            ENHANCED_LOCK / "strong_search" / "docs" / "cross_platform_evidence_summary.md",
            provenance / "strong_search" / "cross_platform_evidence_summary.md",
            "evidence_summary",
        ),
        (
            ENHANCED_LOCK / "strong_search" / "tables" / "cross_platform_evidence_summary.csv",
            provenance / "strong_search" / "cross_platform_evidence_summary.csv",
            "evidence_summary",
        ),
        (
            ENHANCED_LOCK / "cho_imc_source_only" / "docs" / "cho_imc_final_decision_report.md",
            provenance / "cho_imc_source_only" / "cho_imc_final_decision_report.md",
            "source_only_archive",
        ),
        (
            ENHANCED_LOCK / "cho_imc_source_only" / "source_data" / "Source_Data_Cho_IMC.csv",
            provenance / "cho_imc_source_only" / "Source_Data_Cho_IMC.csv",
            "source_only_archive",
        ),
    ]
    for src, dst, component in provenance_files:
        copy_file(src, dst, rows, component)

    for src in sorted((ENHANCED_LOCK / "manifest").glob("*")):
        if src.is_file():
            copy_file(src, provenance / "enhanced_lock_manifest" / src.name, rows, "manifest")

    root_manifest_files = [
        "enhanced_submission_file_manifest.csv",
        "enhanced_source_data_manifest.csv",
        "enhanced_script_manifest.csv",
        "enhanced_dataset_manifest.csv",
        "enhanced_parameter_manifest.csv",
        "checksums_sha256_enhanced_lock.txt",
    ]
    for name in root_manifest_files:
        copy_file(ENHANCED_LOCK / name, provenance / "enhanced_lock_manifest" / name, rows, "manifest")

    # V3 package docs and manifests.
    component_rows = [
        {
            "package_section": "submission_ready/manuscript",
            "description": "Clean manuscript candidate with ED10 v1 and GSE310352 slide/FOV-level wording.",
            "upload_decision": "submission_candidate",
        },
        {
            "package_section": "submission_ready/figures",
            "description": "Main Figures 1-4 and Extended Data Figures 1-10 as PDFs.",
            "upload_decision": "submission_candidate",
        },
        {
            "package_section": "submission_ready/supplementary_tables",
            "description": "Base supplementary tables copied unchanged from the frozen submission candidate.",
            "upload_decision": "submission_candidate",
        },
        {
            "package_section": "submission_ready/source_data",
            "description": "Base source data package plus Source Data for Extended Data Fig. 10.",
            "upload_decision": "submission_candidate",
        },
        {
            "package_section": "provenance_internal",
            "description": "Internal review, rebuild, QC, and claim-boundary provenance; not intended as the journal upload set unless requested.",
            "upload_decision": "internal_review_only",
        },
    ]
    write_csv(
        PACKAGE / "manifest" / "submission_v3_component_map.csv",
        ["package_section", "description", "upload_decision"],
        component_rows,
    )

    write_csv(
        PACKAGE / "manifest" / "submission_v3_file_manifest.csv",
        ["component", "package_path", "source_lock", "source_relative_path", "status", "note"],
        rows,
    )

    source_rows = []
    for src_file in sorted((submission / "source_data").glob("*")):
        source_rows.append(
            {
                "file": str(src_file.relative_to(PACKAGE)).replace("\\", "/"),
                "bytes": str(src_file.stat().st_size),
                "ed10_added": "yes" if src_file.name == "Source_Data_Extended_Data_Figure_10_v1.csv" else "no",
            }
        )
    write_csv(PACKAGE / "manifest" / "submission_v3_source_data_manifest.csv", ["file", "bytes", "ed10_added"], source_rows)

    readme = """# Submission package v3

This package is a clean v3 submission candidate assembled from the frozen base lock and the enhanced ED10 v1 lock.

## Submission-ready contents

- `submission_ready/manuscript/Manuscript_submission_v3_with_ED10.docx`
- `submission_ready/figures/main_figures/`: Figures 1-4
- `submission_ready/figures/extended_data_figures/`: Extended Data Figures 1-10
- `submission_ready/supplementary_tables/`: Supplementary Tables copied unchanged from the base frozen candidate
- `submission_ready/source_data/`: base source data package plus Source Data for Extended Data Fig. 10

## Internal provenance

`provenance_internal/` contains ED10 claim-boundary notes, GSE310352 slide/FOV-level and cell-state transparency QC, Cho IMC source-only decision notes, strong-search summaries, and enhanced-lock manifests. These files are included for human review and traceability, not as the minimal journal upload set.

## Claim boundaries preserved

- GeoMx datasets support compartment-level CAF/matrix and immune/TME programs.
- GSE310352 CosMx supports slide/FOV-level CAF/matrix-associated TGF/EMT stromal-interface organization.
- GSE310352 is not interpreted as patient-level or specimen-level validation.
- GSE310352 cell states are rule-based because public processed files lacked author cell-type annotations.
- No causal signaling, direct SPP1-CD44 validation, tumor-intrinsic EMT, Visium gradient reconstruction, or LN immune-uncoupling validation is claimed from ED10.
- Cho IMC remains source-only and is not included in ED10 v1.

"""
    write_text(PACKAGE / "README_submission_package_v3.md", readme)

    qa = """# Submission v3 QA report

## Assembly

- Manuscript source: enhanced ED10 v1 submission-safe slide/FOV-level patch.
- Figures 1-4 and Extended Data Figures 1-9: copied unchanged from the base frozen submission candidate.
- Extended Data Figure 10: copied unchanged from the enhanced ED10 v1 lock.
- Source data: base source data package plus Source Data for Extended Data Fig. 10.

## Scientific scope

No new biological analysis was run while assembling this package. No figure panels, labels, source-data values, or manuscript scientific conclusions were changed.

## Remaining human checks

- Final Data availability repository URL/DOI must be filled by the authors if not already final.
- Final Code availability reviewer link must be filled by the authors if not already final.
- Journal-specific file naming, word limits, and upload categories should be checked before submission.

"""
    write_text(PACKAGE / "docs" / "submission_v3_QA_report.md", qa)

    # Local path audit for the package-facing manuscript and package docs.
    audit_targets = [
        submission / "manuscript" / "Manuscript_submission_v3_with_ED10.docx",
        PACKAGE / "README_submission_package_v3.md",
        PACKAGE / "docs" / "submission_v3_QA_report.md",
    ]
    audit_rows = audit_local_paths(audit_targets)
    write_csv(PACKAGE / "manifest" / "submission_v3_local_path_audit.csv", ["file", "pattern", "match_count", "status"], audit_rows)
    audit_md_lines = ["# Local path audit for submission package v3", ""]
    flagged = [r for r in audit_rows if r["status"] != "pass"]
    if flagged:
        audit_md_lines.append("Potential local-path strings were found and need review:")
        for row in flagged:
            audit_md_lines.append(f"- `{row['file']}` pattern `{row['pattern']}`: {row['match_count']}")
    else:
        audit_md_lines.append("No local-path strings were detected in the submission-facing manuscript or v3 package README/QA docs.")
    audit_md_lines.append("")
    audit_md_lines.append("Internal provenance files may contain relative provenance paths by design.")
    write_text(PACKAGE / "docs" / "local_path_audit_submission_v3.md", "\n".join(audit_md_lines) + "\n")

    # Checksums after all package files except the checksum itself are present.
    checksum_path = PACKAGE / "manifest" / "checksums_sha256_submission_v3.txt"
    checksum_lines: list[str] = []
    for path in sorted(p for p in PACKAGE.rglob("*") if p.is_file() and p != checksum_path):
        checksum_lines.append(f"{sha256(path)}  {path.relative_to(PACKAGE).as_posix()}")
    write_text(checksum_path, "\n".join(checksum_lines) + "\n")

    # Zip the package for transfer.
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
        for path in sorted(p for p in PACKAGE.rglob("*") if p.is_file()):
            zf.write(path, path.relative_to(PACKAGE).as_posix())

    print(f"PACKAGE={PACKAGE}")
    print(f"ZIP={ZIP_PATH}")
    print(f"FILES={sum(1 for _ in PACKAGE.rglob('*') if _.is_file())}")
    print(f"ZIP_BYTES={ZIP_PATH.stat().st_size}")
    print(f"LOCAL_PATH_FLAGS={len(flagged)}")


if __name__ == "__main__":
    main()
