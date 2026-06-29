from __future__ import annotations

import csv
from pathlib import Path

from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures" / "submission"
OUT_TABLE = ROOT / "results" / "tables" / "submission_final_figure_manifest.csv"
OUT_REPORT = ROOT / "results" / "reports" / "top_journal_figure_package_audit_2026_06_28.md"


FIGURES = [
    ("Figure 1", "main", "figure1_submission_spatial_specificity_nc_style", "Spatial specificity", "CAF-myeloid cores are reproducible spatial organizing regions."),
    ("Figure 2", "main", "figure2_submission_metastatic_decoupling_nc_style", "Metastatic contrast", "Primary and liver metastases preserve CAF-core immune coupling, whereas lymph-node metastases decouple immune organization."),
    ("Figure 3", "main", "figure3_submission_ecotypes_mechanism_axes_nc_style", "Ecotype and interface biology", "CAF-core ecotypes nominate immune-decoupled invasive-interface states."),
    ("Figure 4", "main", "figure4_submission_multiresolution_validation_nc_style", "External multi-resolution validation", "Independent Visium and Xenium datasets support CAF-domain organization."),
    ("Extended Data Figure 1", "extended-data", "figure1_supplement_submission_post_nact_spatial_example", "Representative discovery section", "A post-NACT section shows inspectable CAF-core spatial organization."),
    ("Extended Data Figure 2", "extended-data", "figure2_supplement_submission_spatial_examples", "Representative validation sections", "Primary, liver-metastasis and lymph-node-metastasis maps show the site contrast at larger scale."),
    ("Extended Data Figure 3", "extended-data", "figure3_supplement_targeted_gene_axis_validation", "Targeted genes", "Candidate CAF-core axes are supported at targeted-gene level."),
    ("Extended Data Figure 4", "extended-data", "extended_data_figure4_he_morphology_bridge", "H&E bridge", "Simple H&E features provide an exploratory pathology bridge."),
    ("Extended Data Figure 5", "extended-data", "extended_data_figure5_external_anchor_robustness", "Robustness anchor", "Random-core, threshold and external-anchor tests reinforce spatial specificity."),
    ("Extended Data Figure 6", "extended-data", "extended_data_figure6_cell_state_reference_support", "Cell-state support", "Marker and reference-projection analyses support CAF/myeloid interpretation."),
    ("Extended Data Figure 7", "extended-data", "extended_data_tcga_paad_bulk_context", "TCGA bulk context", "Nominated axes map to broader non-spatial PAAD tumor-microenvironment variation."),
    ("Extended Data Figure 8", "extended-data", "extended_data_figure8_cxcl9_spp1_polarity", "SPP1/CXCL9 polarity", "CAF cores are SPP1/TAM-high, while immune-decoupled samples show attenuated CXCL9/IFN coupling."),
    ("Extended Data Figure 9", "extended-data", "extended_data_figure9_focused_interface_axes", "Focused interface axes", "Matrix-integrin, SPP1-CD44/integrin and TGF-beta/TGFBR are nominated communication axes."),
    ("Extended Data Figure 10", "extended-data", "extended_data_figure10_gse274557_external_validation", "External Visium validation", "GSE274557 validates broad CAF-core organization across metastatic contexts."),
    ("Extended Data Figure 11", "extended-data", "extended_data_figure11_gse274673_xenium_cell_resolution", "Xenium validation", "GSE274673 supports cell-resolution CAF-domain immune/myeloid centering."),
    ("Extended Data Figure 12", "extended-data", "extended_data_figure12_distance_to_caf_core_dynamics", "Distance dynamics", "Spatial programs decay or diverge with distance from CAF cores."),
    ("Extended Data Figure 13", "extended-data", "extended_data_figure13_xenium_cell_domain_maps", "Xenium maps", "Cell-level CAF-domain examples make the validation inspectable."),
    ("Extended Data Figure 14", "extended-data", "extended_data_figure14_spatial_atlas_overview", "Atlas overview", "The spatial atlas scale and cohort structure are transparent."),
    ("Extended Data Figure 15", "extended-data", "extended_data_figure15_local_spatial_program_maps", "Local program maps", "Representative local spatial programs show tissue-scale heterogeneity."),
    ("Extended Data Figure 16", "extended-data", "extended_data_figure16_interface_compartment_maps", "Interface compartments", "Interface regions and candidate programs are visible in selected sections."),
    ("Extended Data Figure 17", "extended-data", "extended_data_figure17_he_patch_examples", "H&E examples", "Patch examples ground the exploratory morphology bridge."),
    ("Extended Data Figure 18", "extended-data", "extended_data_figure18_xenium_program_neighborhoods", "Xenium neighborhoods", "Program-defined neighborhoods around CAF-SPP1/TAM domains support the cell-resolution layer."),
    ("Extended Data Figure 19", "extended-data", "extended_data_figure19_random_core_null_diagnostics", "Null diagnostics", "Representative random-core nulls expose the spatial-specificity baseline."),
    ("Extended Data Figure 20", "extended-data", "extended_data_figure20_ecotype_context_flow", "Ecotype flow", "CAF-core ecotype architecture varies across tissue contexts."),
    ("Extended Data Figure 21", "extended-data", "extended_data_figure21_mechanism_triangulation_priority", "Mechanism triangulation", "Matrix-integrin and SPP1-CD44/integrin are top perturbation-ready candidates."),
    ("Extended Data Figure 22", "extended-data", "extended_data_figure22_tcga_survival_context", "TCGA survival context", "Matrix/stromal-myeloid axes have exploratory adverse bulk survival-context associations."),
    ("Extended Data Figure 23", "extended-data", "extended_data_figure23_tls_maturity_stress_test", "TLS boundary", "Strict TLS-maturity support is limited; mature TLS is not the central claim."),
    ("Extended Data Figure 24", "extended-data", "extended_data_figure24_review_risk_resolution_nc_style", "Reviewer-risk synthesis", "Vulnerable claims are converted into tested or bounded claims."),
    ("Extended Data Figure 25", "extended-data", "extended_data_figure25_spatial_robustness_module_nc_style", "Spatial robustness synthesis", "Random-core, threshold, distance-gradient and external-validation evidence is consolidated."),
    ("Extended Data Figure 26", "extended-data", "extended_data_figure26_cell_state_reference_xenium_module_nc_style", "Cell-state synthesis", "Marker, reference-projection, NNLS-adjacent and Xenium views support cell-state interpretation."),
    ("Extended Data Figure 27", "extended-data", "extended_data_figure27_metastatic_immune_decoupling_module_nc_style", "Metastatic immune-decoupling synthesis", "Lymph-node metastases selectively weaken immune/IFN coupling while retaining stromal-tumor coupling."),
    ("Extended Data Figure 28", "extended-data", "extended_data_figure28_strict_nnls_deconvolution_sensitivity", "Strict NNLS sensitivity", "Per-spot NNLS preserves the key immune-decoupling and myCAF/matrix directions."),
    ("Extended Data Figure 29", "extended-data", "extended_data_figure29_mechanism_gene_interface_module_nc_style", "Mechanism gene/interface synthesis", "Targeted-gene, interface and priority evidence deepen the perturbation agenda."),
]

AUDIT_ONLY = [
    "figure1_submission_spatial_specificity",
    "figure2_submission_metastatic_decoupling",
    "figure3_submission_ecotypes_mechanism_axes",
    "figure4_submission_multiresolution_validation",
    "figure3_candidate_v2_ecotype_interface_story",
    "figure3_candidate_nc_style_ecotype_interface_story",
    "extended_data_gap1_cxcl9_spp1_polarity",
    "extended_data_gap2_cell_state_attribution",
    "extended_data_gap2_reference_projection_deconvolution",
    "extended_data_gap2_full_reference_projection_deconvolution",
    "extended_data_gap2_reference_projection_small_vs_full_comparison",
    "extended_data_gap3_focused_lr_interface",
]


def file_size(path: Path) -> int:
    return path.stat().st_size if path.exists() else 0


def image_dimensions(path: Path) -> tuple[int, int]:
    if not path.exists():
        return 0, 0
    with Image.open(path) as im:
        return im.size


def risk_for(label: str, base: str, width: int, height: int) -> str:
    dense_map_figures = {
        "figure2_submission_metastatic_decoupling_nc_style",
        "figure3_submission_ecotypes_mechanism_axes_nc_style",
        "figure4_submission_multiresolution_validation_nc_style",
        "figure2_supplement_submission_spatial_examples",
        "extended_data_figure13_xenium_cell_domain_maps",
        "extended_data_figure15_local_spatial_program_maps",
        "extended_data_figure16_interface_compartment_maps",
    }
    if base in dense_map_figures:
        return "readability-risk: map-heavy; inspect at final journal column size"
    if label in {"Extended Data Figure 24"}:
        return "synthesis-risk: text/summary content should stay supplement-only"
    if width < 1600 or height < 1000:
        return "production-risk: raster preview is small; rely on vector PDF/SVG"
    return "low"


def main() -> None:
    rows = []
    for label, role, base, figure_role, conclusion in FIGURES:
        pdf = FIG_DIR / f"{base}.pdf"
        svg = FIG_DIR / f"{base}.svg"
        png = FIG_DIR / f"{base}.png"
        width, height = image_dimensions(png)
        rows.append(
            {
                "label": label,
                "role": role,
                "figure_role": figure_role,
                "core_conclusion": conclusion,
                "pdf": str(pdf.relative_to(ROOT)),
                "svg": str(svg.relative_to(ROOT)),
                "png": str(png.relative_to(ROOT)),
                "pdf_exists": pdf.exists(),
                "svg_exists": svg.exists(),
                "png_exists": png.exists(),
                "png_width": width,
                "png_height": height,
                "pdf_bytes": file_size(pdf),
                "svg_bytes": file_size(svg),
                "png_bytes": file_size(png),
                "review_risk": risk_for(label, base, width, height),
            }
        )

    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TABLE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    missing = [
        row
        for row in rows
        if not (row["pdf_exists"] and row["svg_exists"] and row["png_exists"])
    ]
    risk_rows = [row for row in rows if row["review_risk"] != "low"]

    report = []
    report.append("# Top-Journal Figure Package Audit")
    report.append("")
    report.append("Last updated: 2026-06-28")
    report.append("")
    report.append("## Scope")
    report.append("")
    report.append("This audit defines the final intended figure set and separates it from retained audit-only or candidate graphics.")
    report.append("")
    report.append("## Final Figure Set")
    report.append("")
    report.append(f"- Main figures: {sum(row['role'] == 'main' for row in rows)}")
    report.append(f"- Extended Data figures: {sum(row['role'] == 'extended-data' for row in rows)}")
    report.append(f"- Figures with complete PDF/SVG/PNG exports: {len(rows) - len(missing)} / {len(rows)}")
    report.append(f"- Figures with readability or production risk flags: {len(risk_rows)}")
    report.append("")
    report.append("## Figure Contracts")
    report.append("")
    report.append("| figure | role | core conclusion | risk |")
    report.append("|---|---|---|---|")
    for row in rows:
        report.append(
            f"| {row['label']} | {row['figure_role']} | {row['core_conclusion']} | {row['review_risk']} |"
        )
    report.append("")
    report.append("## Audit-Only Figures")
    report.append("")
    report.append("The following files should be treated as provenance, candidate or earlier-version figures, not as the active submission set unless deliberately promoted:")
    report.append("")
    for base in AUDIT_ONLY:
        exists = (FIG_DIR / f"{base}.pdf").exists()
        status = "exists" if exists else "missing"
        report.append(f"- `{FIG_DIR.relative_to(ROOT) / (base + '.pdf')}` ({status})")
    report.append("")
    report.append("## Top-Journal Readability Priorities")
    report.append("")
    report.append("1. Keep Figures 1-4 as the main editorial path: discovery specificity, metastatic remodeling, ecotype/interface biology and external multi-resolution validation.")
    report.append("2. Keep ED24-ED29 as reviewer-facing synthesis modules, not as additional main-text claims.")
    report.append("3. Inspect map-heavy figures at final print size before submission, especially Figure 2, Figure 3, Figure 4 and ED13/ED15/ED16.")
    report.append("4. Do not submit earlier compact or gap figures as numbered figures unless the numbering and captions are deliberately rewritten.")
    report.append("5. Treat TCGA panels as context and boundary-setting material; spatial localization and clinical prediction should remain outside the claim.")
    report.append("")
    report.append("## Visual Inspection Decision")
    report.append("")
    report.append("The main NC-style figures were inspected at full raster-preview size after the automated manifest audit.")
    report.append("")
    report.append("- Figure 1 is ready as the opening discovery/specificity figure.")
    report.append("- Figure 2 is scientifically clear and should remain the main metastatic-remodeling figure. Minor polish could increase spacing around the atlas-count panel and reduce top whitespace, but no scientific redraw is required.")
    report.append("- Figure 3 is dense but coherent as an NC-style multi-panel main figure. If a journal requests a cleaner main figure, the contingency version should keep panels A-F in the main figure and move the representative spatial maps G-R to Extended Data.")
    report.append("- Figure 4 is readable and works as a multi-resolution validation figure. The Xenium maps are visually useful because they show cell-resolution support rather than merely repeating heatmaps.")
    report.append("- ED24 is intentionally a reviewer-risk synthesis module and should stay in Extended Data or supplement, not in the main figure set.")
    report.append("")
    report.append("## Output")
    report.append("")
    report.append(f"- Manifest table: `{OUT_TABLE.relative_to(ROOT)}`")
    report.append("- Visual contact sheet: `results/figures/qc/main_and_ed24_29_contact_sheet.png`")
    report.append("")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
