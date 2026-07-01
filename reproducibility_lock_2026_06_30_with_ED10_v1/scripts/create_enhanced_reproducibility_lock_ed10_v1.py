from __future__ import annotations

import csv
import hashlib
import shutil
from datetime import datetime
from pathlib import Path


PROJECT = Path("pdac_spatial_ecology")
RESULTS = PROJECT / "results"
ORTHO = RESULTS / "orthogonal_validation_strong_search"
REV = RESULTS / "revision_2026_06_30_with_ED10_v1"
BASE_LOCK = PROJECT / "reproducibility_lock_2026_06_30"
LOCK = PROJECT / "reproducibility_lock_2026_06_30_with_ED10_v1"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def copy_file(src: Path, dest_rel: Path, category: str, description: str, rows: list[dict], required: bool = True):
    dest = LOCK / dest_rel
    if src.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        rows.append(
            {
                "category": category,
                "description": description,
                "status": "copied",
                "source_path": str(src.as_posix()),
                "lock_path": str(dest.as_posix()),
                "size_bytes": dest.stat().st_size,
                "sha256": sha256_file(dest),
            }
        )
    else:
        rows.append(
            {
                "category": category,
                "description": description,
                "status": "missing_required" if required else "not_available",
                "source_path": str(src.as_posix()),
                "lock_path": str(dest.as_posix()),
                "size_bytes": "",
                "sha256": "",
            }
        )


def write_csv(path: Path, rows: list[dict], fieldnames: list[str] | None = None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        keys: list[str] = []
        for row in rows:
            for key in row:
                if key not in keys:
                    keys.append(key)
        fieldnames = keys
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def copy_text_if_exists(src: Path, dest: Path):
    if src.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)


def main():
    if not BASE_LOCK.exists():
        raise SystemExit(f"Base lock not found: {BASE_LOCK}")
    if LOCK.exists():
        raise SystemExit(f"Enhanced lock already exists; refusing to overwrite: {LOCK}")
    LOCK.mkdir(parents=True)
    (LOCK / "docs").mkdir()
    (LOCK / "manifest").mkdir()

    rows: list[dict] = []

    requested_files = [
        # Manuscript candidate
        (REV / "Manuscript_with_ED10_v1_submission_safe_slidelevel_patch.docx", Path("manuscript/Manuscript_with_ED10_v1_submission_safe_slidelevel_patch.docx"), "manuscript", "Enhanced submission-safe manuscript candidate with GSE310352 slide/FOV-level patch", True),
        # ED10
        (ORTHO / "figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf", Path("ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf"), "ed10_figure", "Extended Data Fig. 10 v1 PDF", True),
        (ORTHO / "figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg", Path("ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg"), "ed10_figure", "Extended Data Fig. 10 v1 SVG", True),
        (ORTHO / "source_data/Source_Data_Extended_Data_Figure_10_v1.csv", Path("ed10/source_data/Source_Data_Extended_Data_Figure_10_v1.csv"), "ed10_source_data", "Source Data for Extended Data Fig. 10 v1", True),
        (ORTHO / "docs/ed10_v1_figure_legend_draft.md", Path("ed10/docs/ed10_v1_figure_legend_draft.md"), "ed10_doc", "ED10 v1 figure legend draft", True),
        (ORTHO / "docs/ed10_v1_methods_draft.md", Path("ed10/docs/ed10_v1_methods_draft.md"), "ed10_doc", "ED10 v1 methods draft", True),
        (ORTHO / "docs/ed10_v1_results_paragraph_draft.md", Path("ed10/docs/ed10_v1_results_paragraph_draft.md"), "ed10_doc", "ED10 v1 results paragraph draft", True),
        (ORTHO / "docs/ed10_v1_claim_boundary_notes.md", Path("ed10/docs/ed10_v1_claim_boundary_notes.md"), "ed10_doc", "ED10 v1 claim-boundary notes", True),
        (ORTHO / "manifest/ed10_v1_panel_map.csv", Path("ed10/manifest/ed10_v1_panel_map.csv"), "ed10_manifest", "ED10 v1 panel map", True),
        (ORTHO / "manifest/ed10_v1_source_data_manifest.csv", Path("ed10/manifest/ed10_v1_source_data_manifest.csv"), "ed10_manifest", "ED10 v1 source-data manifest", True),
        (ORTHO / "manifest/ed10_v1_script_manifest.csv", Path("ed10/manifest/ed10_v1_script_manifest.csv"), "ed10_manifest", "ED10 v1 script manifest", True),
        (ORTHO / "manifest/ed10_v1_parameter_manifest.csv", Path("ed10/manifest/ed10_v1_parameter_manifest.csv"), "ed10_manifest", "ED10 v1 parameter manifest", True),
        (ORTHO / "manifest/checksums_sha256_ed10_v1.txt", Path("ed10/manifest/checksums_sha256_ed10_v1.txt"), "ed10_manifest", "ED10 v1 original checksum file", True),
        # GSE310352 support docs/tables/figures
        (ORTHO / "docs/gse310352_patient_mapping_recovery_report.md", Path("gse310352/docs/gse310352_patient_mapping_recovery_report.md"), "gse310352_boundary", "GSE310352 patient/sample mapping recovery report", True),
        (ORTHO / "tables/gse310352_slide_patient_mapping_attempt.csv", Path("gse310352/tables/gse310352_slide_patient_mapping_attempt.csv"), "gse310352_boundary", "GSE310352 slide-patient mapping attempt table", True),
        (ORTHO / "docs/gse310352_cell_state_definition_transparency_report.md", Path("gse310352/docs/gse310352_cell_state_definition_transparency_report.md"), "gse310352_transparency", "GSE310352 cell-state definition transparency report", True),
        (ORTHO / "tables/gse310352_cell_state_definition_marker_enrichment.csv", Path("gse310352/tables/gse310352_cell_state_definition_marker_enrichment.csv"), "gse310352_transparency", "Marker enrichment table for rule-based cell states", True),
        (ORTHO / "tables/gse310352_cell_state_if_marker_qc.csv", Path("gse310352/tables/gse310352_cell_state_if_marker_qc.csv"), "gse310352_transparency", "IF marker QC table for rule-based cell states", True),
        (ORTHO / "tables/gse310352_cell_state_overlap_audit.csv", Path("gse310352/tables/gse310352_cell_state_overlap_audit.csv"), "gse310352_transparency", "Cell-state overlap audit table", True),
        (ORTHO / "tables/gse310352_tgfemt_spatial_identity_context.csv", Path("gse310352/tables/gse310352_tgfemt_spatial_identity_context.csv"), "gse310352_transparency", "TGF/EMT spatial identity context table", True),
        (ORTHO / "figures/gse310352_cell_state_definition_transparency_qc.pdf", Path("gse310352/figures/gse310352_cell_state_definition_transparency_qc.pdf"), "gse310352_transparency", "Compact cell-state transparency QC PDF", True),
        (ORTHO / "figures/gse310352_cell_state_definition_transparency_qc.svg", Path("gse310352/figures/gse310352_cell_state_definition_transparency_qc.svg"), "gse310352_transparency", "Compact cell-state transparency QC SVG", True),
        (ORTHO / "figures/gse310352_cell_state_marker_enrichment_heatmap.pdf", Path("gse310352/figures/gse310352_cell_state_marker_enrichment_heatmap.pdf"), "gse310352_transparency", "Marker enrichment heatmap", True),
        (ORTHO / "figures/gse310352_cell_state_if_marker_qc.pdf", Path("gse310352/figures/gse310352_cell_state_if_marker_qc.pdf"), "gse310352_transparency", "IF marker QC figure", True),
        (ORTHO / "figures/gse310352_cell_state_overlap_upset_or_barplot.pdf", Path("gse310352/figures/gse310352_cell_state_overlap_upset_or_barplot.pdf"), "gse310352_transparency", "Cell-state overlap barplot", True),
        (ORTHO / "figures/gse310352_tgfemt_spatial_identity_context.pdf", Path("gse310352/figures/gse310352_tgfemt_spatial_identity_context.pdf"), "gse310352_transparency", "Spatial identity context figure", True),
        # GSE310352 core ED10 source support
        (ORTHO / "source_data/Source_Data_GSE310352_CosMx.csv", Path("gse310352/source_data/Source_Data_GSE310352_CosMx.csv"), "gse310352_source_data", "GSE310352 CosMx source data used in ED10", True),
        (ORTHO / "tables/gse310352_nonoverlap_tgfemt_adjacency_results.csv", Path("gse310352/tables/gse310352_nonoverlap_tgfemt_adjacency_results.csv"), "gse310352_robustness", "Non-overlap TGF/EMT adjacency robustness", True),
        (ORTHO / "tables/gse310352_threshold_sensitivity.csv", Path("gse310352/tables/gse310352_threshold_sensitivity.csv"), "gse310352_robustness", "Threshold sensitivity", True),
        (ORTHO / "tables/gse310352_normalization_sensitivity.csv", Path("gse310352/tables/gse310352_normalization_sensitivity.csv"), "gse310352_robustness", "Normalization sensitivity", True),
        (ORTHO / "tables/gse310352_spatial_null_sensitivity.csv", Path("gse310352/tables/gse310352_spatial_null_sensitivity.csv"), "gse310352_robustness", "Spatial null sensitivity", True),
        (ORTHO / "tables/gse310352_leave_one_slide_out.csv", Path("gse310352/tables/gse310352_leave_one_slide_out.csv"), "gse310352_robustness", "Leave-one-slide-out support", True),
        (ORTHO / "tables/gse310352_caf_tgfemt_gene_overlap_audit.csv", Path("gse310352/tables/gse310352_caf_tgfemt_gene_overlap_audit.csv"), "gse310352_robustness", "CAF/TGF-EMT gene-overlap audit", True),
        # GeoMx validation docs/tables/source data
        (ORTHO / "source_data/Source_Data_GSE199102_GeoMx.csv", Path("geomx/gse199102/source_data/Source_Data_GSE199102_GeoMx.csv"), "geomx_source_data", "GSE199102 GeoMx source data", True),
        (ORTHO / "docs/gse199102_geomx_final_decision_report.md", Path("geomx/gse199102/docs/gse199102_geomx_final_decision_report.md"), "geomx_doc", "GSE199102 GeoMx final decision report", True),
        (ORTHO / "docs/gse199102_geomx_feasibility_report.md", Path("geomx/gse199102/docs/gse199102_geomx_feasibility_report.md"), "geomx_doc", "GSE199102 GeoMx feasibility report", True),
        (ORTHO / "docs/gse199102_geomx_methods_draft.md", Path("geomx/gse199102/docs/gse199102_geomx_methods_draft.md"), "geomx_doc", "GSE199102 GeoMx methods draft", True),
        (ORTHO / "docs/gse199102_geomx_results_paragraph_draft.md", Path("geomx/gse199102/docs/gse199102_geomx_results_paragraph_draft.md"), "geomx_doc", "GSE199102 GeoMx results draft", True),
        (ORTHO / "docs/gse199102_geomx_figure_legend_draft.md", Path("geomx/gse199102/docs/gse199102_geomx_figure_legend_draft.md"), "geomx_doc", "GSE199102 GeoMx figure legend draft", True),
        (ORTHO / "docs/gse199102_gse240078_concordance_interpretation.md", Path("geomx/concordance/docs/gse199102_gse240078_concordance_interpretation.md"), "geomx_concordance", "GSE199102-GSE240078 concordance interpretation", True),
        (ORTHO / "tables/gse199102_gse240078_direction_concordance.csv", Path("geomx/concordance/tables/gse199102_gse240078_direction_concordance.csv"), "geomx_concordance", "GSE199102-GSE240078 direction concordance table", True),
        (ORTHO / "figures/gse199102_gse240078_concordance_plot.pdf", Path("geomx/concordance/figures/gse199102_gse240078_concordance_plot.pdf"), "geomx_concordance", "GSE199102-GSE240078 concordance plot", True),
        (ORTHO / "figures/gse199102_geomx_candidate_panels.pdf", Path("geomx/gse199102/figures/gse199102_geomx_candidate_panels.pdf"), "geomx_figure", "GSE199102 candidate panels", True),
        (ORTHO / "figures/gse199102_geomx_positive_negative_controls.pdf", Path("geomx/gse199102/figures/gse199102_geomx_positive_negative_controls.pdf"), "geomx_figure", "GSE199102 positive/negative controls", True),
        (ORTHO / "tables/gse199102_geomx_module_scores.csv", Path("geomx/gse199102/tables/gse199102_geomx_module_scores.csv"), "geomx_table", "GSE199102 module score matrix", True),
        (ORTHO / "tables/gse199102_geomx_patient_level_deltas.csv", Path("geomx/gse199102/tables/gse199102_geomx_patient_level_deltas.csv"), "geomx_table", "GSE199102 patient-level deltas", True),
        (ORTHO / "tables/gse199102_geomx_paired_segment_results.csv", Path("geomx/gse199102/tables/gse199102_geomx_paired_segment_results.csv"), "geomx_table", "GSE199102 paired segment support", True),
        (ORTHO / "tables/gse199102_geomx_segment_comparison_results.csv", Path("geomx/gse199102/tables/gse199102_geomx_segment_comparison_results.csv"), "geomx_table", "GSE199102 segment comparison results", True),
        (ORTHO / "tables/gse199102_geomx_gene_coverage_by_module.csv", Path("geomx/gse199102/tables/gse199102_geomx_gene_coverage_by_module.csv"), "geomx_table", "GSE199102 gene coverage by module", True),
        (ORTHO / "manifest/gse199102_geomx_manifest.csv", Path("geomx/gse199102/manifest/gse199102_geomx_manifest.csv"), "geomx_manifest", "GSE199102 manifest", True),
        # GSE240078 is represented through ED10 source data and concordance outputs in this workspace.
        (ORTHO / "tables/gse240078_geomx_metadata_summary.csv", Path("geomx/gse240078/tables/gse240078_geomx_metadata_summary.csv"), "geomx_gse240078", "GSE240078 metadata summary if separately available", False),
        (ORTHO / "tables/gse240078_geomx_module_scores.csv", Path("geomx/gse240078/tables/gse240078_geomx_module_scores.csv"), "geomx_gse240078", "GSE240078 module score matrix if separately available", False),
        # Cho IMC source-only archive
        (ORTHO / "cho_imc/docs/cho_imc_final_decision_report.md", Path("cho_imc_source_only/docs/cho_imc_final_decision_report.md"), "cho_imc_source_only", "Cho IMC final decision report; source-only, not ED10", True),
        (ORTHO / "cho_imc/tables/cho_imc_marker_availability.csv", Path("cho_imc_source_only/tables/cho_imc_marker_availability.csv"), "cho_imc_source_only", "Cho IMC marker table", True),
        (ORTHO / "cho_imc/tables/cho_imc_caf_myeloid_adjacency_results.csv", Path("cho_imc_source_only/tables/cho_imc_caf_myeloid_adjacency_results.csv"), "cho_imc_source_only", "Cho IMC main adjacency table", True),
        (ORTHO / "cho_imc/tables/cho_imc_decision_summary.csv", Path("cho_imc_source_only/tables/cho_imc_decision_summary.csv"), "cho_imc_source_only", "Cho IMC decision summary", True),
        (ORTHO / "cho_imc/source_data/Source_Data_Cho_IMC.csv", Path("cho_imc_source_only/source_data/Source_Data_Cho_IMC.csv"), "cho_imc_source_only", "Cho IMC source data; source-only, not ED10", True),
        (ORTHO / "cho_imc/docs/negative_or_weak_results_log.md", Path("cho_imc_source_only/docs/negative_or_weak_results_log.md"), "cho_imc_source_only", "Cho IMC weak/negative/source-only log", True),
        (ORTHO / "cho_imc/manifest/cho_imc_output_manifest.csv", Path("cho_imc_source_only/manifest/cho_imc_output_manifest.csv"), "cho_imc_source_only", "Cho IMC output manifest", True),
        # Strong-search summary
        (ORTHO / "tables/signal_strength_ranking.csv", Path("strong_search/tables/signal_strength_ranking.csv"), "strong_search_summary", "Signal strength ranking table", True),
        (ORTHO / "manifest/signal_strength_ranking.csv", Path("strong_search/manifest/signal_strength_ranking.csv"), "strong_search_summary", "Signal strength ranking manifest copy", True),
        (ORTHO / "docs/candidate_results_for_ED_figure.md", Path("strong_search/docs/candidate_results_for_ED_figure.md"), "strong_search_summary", "Candidate results for ED figure", True),
        (ORTHO / "docs/negative_or_weak_results_log.md", Path("strong_search/docs/negative_or_weak_results_log.md"), "strong_search_summary", "Negative/weak results log", True),
        (ORTHO / "docs/final_orthogonal_validation_decision_summary.md", Path("strong_search/docs/final_orthogonal_validation_decision_summary.md"), "strong_search_summary", "Final orthogonal validation decision summary if available", False),
        (ORTHO / "docs/cross_platform_evidence_summary.md", Path("strong_search/docs/cross_platform_evidence_summary.md"), "strong_search_summary", "Cross-platform evidence summary if available", False),
        (ORTHO / "docs/signal_ranking_report.md", Path("strong_search/docs/signal_ranking_report.md"), "strong_search_summary", "Signal ranking report", True),
        # Scripts
        (ORTHO / "scripts/make_extended_data_fig10_v1.R", Path("scripts/make_extended_data_fig10_v1.R"), "script", "ED10 v1 figure/source generation script", True),
        (ORTHO / "scripts/analyze_gse199102_geomx.R", Path("scripts/analyze_gse199102_geomx.R"), "script", "GSE199102 GeoMx analysis script", True),
        (ORTHO / "scripts/analyze_gse310352_cosmx.R", Path("scripts/analyze_gse310352_cosmx.R"), "script", "GSE310352 CosMx initial analysis script", True),
        (ORTHO / "scripts/robustness_gse310352_cosmx.R", Path("scripts/robustness_gse310352_cosmx.R"), "script", "GSE310352 CosMx robustness script", True),
        (ORTHO / "scripts/gse310352_patient_mapping_recovery.R", Path("scripts/gse310352_patient_mapping_recovery.R"), "script", "GSE310352 patient/sample mapping recovery script", True),
        (ORTHO / "scripts/gse310352_cell_state_definition_transparency.R", Path("scripts/gse310352_cell_state_definition_transparency.R"), "script", "GSE310352 cell-state transparency QC script", True),
        (ORTHO / "cho_imc/scripts/01_inspect_cho_imc_object.R", Path("scripts/cho_imc/01_inspect_cho_imc_object.R"), "script", "Cho IMC object inspection script", True),
        (ORTHO / "cho_imc/scripts/02_analyze_cho_imc_spatial.R", Path("scripts/cho_imc/02_analyze_cho_imc_spatial.R"), "script", "Cho IMC spatial analysis script", True),
    ]

    for src, dest_rel, category, description, required in requested_files:
        copy_file(src, dest_rel, category, description, rows, required=required)

    # Convenience top-level derived manifests requested by the user.
    write_csv(LOCK / "enhanced_submission_file_manifest.csv", rows)

    # ED10 panel map copy with enhanced boundary note.
    panel_src = ORTHO / "manifest/ed10_v1_panel_map.csv"
    if panel_src.exists():
        panel_rows: list[dict] = []
        with panel_src.open("r", encoding="utf-8", newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                row["enhanced_lock_note"] = "ED10 v1 only; Cho IMC not included. GSE310352 panels F/G are slide/FOV-level and rule-based."
                panel_rows.append(row)
        write_csv(LOCK / "enhanced_ed10_panel_map.csv", panel_rows)

    source_rows = [
        r for r in rows
        if "source_data" in r["category"] or r["category"] in {"ed10_source_data", "geomx_source_data"}
    ]
    write_csv(LOCK / "enhanced_source_data_manifest.csv", source_rows)

    script_rows = [r for r in rows if r["category"] == "script"]
    write_csv(LOCK / "enhanced_script_manifest.csv", script_rows)

    dataset_rows = [
        {
            "dataset": "GSE240078",
            "platform": "GeoMx DSP",
            "role": "ED10 v1 compartment-level validation",
            "support_unit": "AOI compartment-level",
            "lock_evidence": "ED10 source data and GSE199102-GSE240078 concordance outputs",
            "boundary": "No causal signaling; no direct SPP1-CD44 stromal support",
        },
        {
            "dataset": "GSE199102",
            "platform": "GeoMx WTA",
            "role": "ED10 v1 independent compartment-level validation",
            "support_unit": "segment/ROI/patient-level summaries where public metadata allowed",
            "lock_evidence": "GSE199102 source data, final decision report, concordance table",
            "boundary": "Compartment-level support only",
        },
        {
            "dataset": "GSE310352",
            "platform": "CosMx SMI",
            "role": "ED10 v1 cell-level stromal-interface support",
            "support_unit": "slide/FOV-level only",
            "lock_evidence": "CosMx source data, robustness tables, patient-mapping boundary, cell-state transparency QC",
            "boundary": "Rule-based states; author annotations unavailable; not patient-level; not tumor-intrinsic EMT",
        },
        {
            "dataset": "Zenodo 15596960 Cho rapid-autopsy PDAC IMC",
            "platform": "IMC",
            "role": "source-only archive",
            "support_unit": "ROI/sample-level source-only",
            "lock_evidence": "Cho IMC final decision report and source tables",
            "boundary": "Not included in ED10 v1; not used to upgrade ED10; source-only",
        },
    ]
    write_csv(LOCK / "enhanced_dataset_manifest.csv", dataset_rows)

    parameter_rows = [
        {"analysis": "ED10 v1", "parameter": "version", "value": "v1", "note": "No ED10 v2 created"},
        {"analysis": "GSE310352 CosMx", "parameter": "support_unit", "value": "slide/FOV-level", "note": "Patient/specimen/tissue-block identifiers not recoverable from public processed-slide metadata"},
        {"analysis": "GSE310352 CosMx", "parameter": "cell_state_definition", "value": "rule-based", "note": "Processed CSV had no author cell-type annotation"},
        {"analysis": "GSE310352 CosMx", "parameter": "TGF/EMT label", "value": "stromal-interface", "note": "Not tumor-intrinsic EMT; not causal EMT induction"},
        {"analysis": "GSE310352 CosMx", "parameter": "gene-overlap boundary", "value": "MMP2 removed; MMP2+ITGA5 removed sensitivity retained 8/8 slide support", "note": "See non-overlap robustness table"},
        {"analysis": "GeoMx", "parameter": "claim_scope", "value": "compartment-level support", "note": "No Visium gradient reconstruction"},
        {"analysis": "Cho IMC", "parameter": "display_decision", "value": "source-only", "note": "Not included in ED10 v1"},
    ]
    write_csv(LOCK / "enhanced_parameter_manifest.csv", parameter_rows)

    unresolved = [
        {
            "item": "GSE310352 patient/specimen/tissue-block mapping",
            "status": "unresolved",
            "note": "Public metadata did not allow reliable recovery; keep slide/FOV-level wording.",
        },
        {
            "item": "GSE310352 author cell-type annotations",
            "status": "unavailable",
            "note": "Processed CSV lacked author annotations; labels remain rule-based.",
        },
        {
            "item": "GSE240078 standalone intermediate files",
            "status": "partially archived",
            "note": "No separately named GSE240078 intermediate summary files were present in this workspace; ED10 source data and concordance outputs preserve the displayed GSE240078 quantities.",
        },
        {
            "item": "Cho IMC ED10 inclusion",
            "status": "not included",
            "note": "Cho IMC remains source-only because support was modest/secondary.",
        },
        {
            "item": "Mechanism",
            "status": "not established",
            "note": "No causal signaling, tumor-intrinsic EMT, direct SPP1-CD44 validation, or Visium gradient reconstruction is claimed.",
        },
    ]
    unresolved_md = [
        "# Enhanced Lock Unresolved Items",
        "",
        *[f"- **{x['item']}**: {x['status']}. {x['note']}" for x in unresolved],
        "",
    ]
    (LOCK / "enhanced_unresolved_items.md").write_text("\n".join(unresolved_md), encoding="utf-8", newline="\n")

    claim_lines = [
        "# Enhanced Claim Boundary Checklist",
        "",
        "- [x] ED10 v1 is frozen as optional enhanced validation; no ED10 v2 created.",
        "- [x] GeoMx supports compartment-level CAF/matrix and immune/TME programs.",
        "- [x] GSE310352 CosMx supports CAF/matrix-associated TGF/EMT stromal-interface organization at slide/FOV level.",
        "- [x] GSE310352 is not interpreted as patient-level or specimen-level validation.",
        "- [x] GSE310352 cell states are rule-based because author cell-type annotations were unavailable.",
        "- [x] TGF/EMT stromal-interface is not claimed as tumor-intrinsic EMT.",
        "- [x] No causal signaling is claimed.",
        "- [x] No direct SPP1-CD44 validation is claimed.",
        "- [x] No Visium distance-gradient reconstruction is claimed.",
        "- [x] No lymph-node immune-uncoupling validation is claimed from ED10/GSE310352/Cho IMC.",
        "- [x] Cho IMC remains source-only and is not included in ED10 v1.",
        "",
    ]
    (LOCK / "enhanced_claim_boundary_checklist.md").write_text("\n".join(claim_lines), encoding="utf-8", newline="\n")

    rebuild_lines = [
        "# How To Rebuild ED10 v1",
        "",
        "This enhanced lock archives the ED10 v1 figure/source outputs and the scripts needed to regenerate or audit the components available in this workspace.",
        "",
        "## Inputs",
        "- GSE240078 GeoMx DSP quantities are preserved in `ed10/source_data/Source_Data_Extended_Data_Figure_10_v1.csv` and concordance outputs.",
        "- GSE199102 GeoMx WTA source tables and script are archived under `geomx/gse199102/` and `scripts/analyze_gse199102_geomx.R`.",
        "- GSE310352 CosMx source/robustness/transparency outputs are archived under `gse310352/` with scripts under `scripts/`.",
        "- Cho IMC is archived under `cho_imc_source_only/` but is not part of ED10 v1.",
        "",
        "## Rebuild Steps",
        "1. Re-run upstream public-dataset scripts only if the public input data are available locally.",
        "2. Re-run GSE199102: `Rscript scripts/analyze_gse199102_geomx.R` from the project root.",
        "3. Re-run GSE310352 initial/robustness/transparency scripts from the project root if needed.",
        "4. Re-run ED10 figure generation with `scripts/make_extended_data_fig10_v1.R` if present and inputs are available.",
        "5. Compare outputs to `checksums_sha256_enhanced_lock.txt`.",
        "",
        "## Claim Boundaries",
        "Keep GSE310352 slide/FOV-level, rule-based and stromal-interface. Do not convert it into patient-level validation, tumor-intrinsic EMT or causal signaling.",
        "",
    ]
    (LOCK / "docs/how_to_rebuild_ED10.md").write_text("\n".join(rebuild_lines), encoding="utf-8", newline="\n")

    summary_lines = [
        "# Enhanced Submission Summary",
        "",
        f"Created: {datetime.now().isoformat(timespec='seconds')}",
        "",
        "This lock freezes the enhanced submission candidate with ED10 v1 and the post-ED10 claim-boundary support files.",
        "",
        "## Included",
        "- Submission-safe manuscript candidate with GSE310352 slide/FOV-level wording.",
        "- Extended Data Fig. 10 v1 PDF/SVG and source data.",
        "- ED10 docs and manifests.",
        "- GSE310352 patient-mapping boundary and cell-state transparency QC.",
        "- GSE199102/GSE240078 GeoMx concordance support.",
        "- Cho IMC source-only archive, explicitly excluded from ED10 v1.",
        "- Strong-search summaries and negative/weak logs where available.",
        "",
        "## Status",
        "Ready for final human review. Remaining items are wording/provenance checks, not additional biological analysis.",
        "",
    ]
    (LOCK / "docs/enhanced_submission_summary.md").write_text("\n".join(summary_lines), encoding="utf-8", newline="\n")

    changed_lines = [
        "# What Changed From Base Lock",
        "",
        f"Base lock preserved unchanged at `{BASE_LOCK.as_posix()}`.",
        "",
        "## Added In Enhanced Lock",
        "- ED10 v1 optional orthogonal validation layer.",
        "- Submission-safe manuscript candidate with ED10 and GSE310352 slide/FOV-level boundary patch.",
        "- GeoMx/CosMx ED10 source data and docs.",
        "- GSE310352 patient/sample mapping recovery proving public metadata do not support patient/specimen-level interpretation.",
        "- GSE310352 cell-state transparency QC for rule-based CAF/matrix-like and TGF/EMT stromal-interface states.",
        "- Cho IMC source-only archive and final source-only decision.",
        "",
        "## Not Changed",
        "- Original `reproducibility_lock_2026_06_30/`.",
        "- Figure 1-4.",
        "- ED Fig. 1-10 v1 figure files.",
        "- Scientific conclusions and claim boundaries.",
        "",
    ]
    (LOCK / "docs/what_changed_from_base_lock.md").write_text("\n".join(changed_lines), encoding="utf-8", newline="\n")

    readme_lines = [
        "# Reproducibility Lock 2026-06-30 With ED10 v1",
        "",
        "Enhanced submission candidate archive. This directory is additive and does not overwrite the base lock.",
        "",
        "## Core Candidate",
        "- Manuscript: `manuscript/Manuscript_with_ED10_v1_submission_safe_slidelevel_patch.docx`",
        "- ED10 PDF: `ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.pdf`",
        "- ED10 SVG: `ed10/figures/Extended_Data_Figure_10_Orthogonal_GeoMx_CosMx_v1.svg`",
        "- ED10 source data: `ed10/source_data/Source_Data_Extended_Data_Figure_10_v1.csv`",
        "",
        "## Interpretation Boundaries",
        "- GSE310352 is slide/FOV-level support only.",
        "- GSE310352 cell states are rule-based because author annotations were unavailable.",
        "- TGF/EMT is interpreted as stromal-interface, not tumor-intrinsic EMT.",
        "- No causal signaling or direct SPP1-CD44 validation is claimed.",
        "- Cho IMC is source-only and not part of ED10 v1.",
        "",
        "## Manifests",
        "- `enhanced_submission_file_manifest.csv`",
        "- `enhanced_ed10_panel_map.csv`",
        "- `enhanced_source_data_manifest.csv`",
        "- `enhanced_script_manifest.csv`",
        "- `enhanced_dataset_manifest.csv`",
        "- `enhanced_parameter_manifest.csv`",
        "- `checksums_sha256_enhanced_lock.txt`",
        "",
    ]
    (LOCK / "README.md").write_text("\n".join(readme_lines), encoding="utf-8", newline="\n")

    # Checksums for every file in the enhanced lock except the checksum file itself.
    checksum_path = LOCK / "checksums_sha256_enhanced_lock.txt"
    checksum_entries = []
    for path in sorted(LOCK.rglob("*")):
        if path.is_file() and path != checksum_path:
            rel = path.relative_to(LOCK).as_posix()
            checksum_entries.append((sha256_file(path), rel))
    checksum_path.write_text(
        "".join(f"{sha}  {rel}\n" for sha, rel in checksum_entries),
        encoding="utf-8",
        newline="\n",
    )

    print(f"Created enhanced lock: {LOCK.resolve()}")
    print(f"Files manifest rows: {len(rows)}")
    print(f"Missing/not-available entries: {sum(1 for r in rows if r['status'] != 'copied')}")


if __name__ == "__main__":
    main()
