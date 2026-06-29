from __future__ import annotations

import csv
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures" / "submission"
QC_DIR = ROOT / "results" / "figures" / "qc" / "supplement_modules"
OUT_TABLE = ROOT / "results" / "tables" / "supplement_module_reorganization_manifest.csv"
OUT_REPORT = ROOT / "results" / "reports" / "supplement_module_reorganization_plan_2026_06_28.md"


MODULES = [
    {
        "module": "Supplementary Module 1",
        "short_name": "spatial_specificity_robustness",
        "title": "Spatial specificity and robustness of CAF-myeloid cores",
        "core_question": "Is the organizing-core signal spatially specific and robust to null models and threshold choices?",
        "recommended_final": "Use ED25 as the module backbone; retain ED1, ED5, ED12 and ED19 as component/detail plates if the journal allows.",
        "current_figures": [
            ("ED1", "figure1_supplement_submission_post_nact_spatial_example", "Representative discovery section"),
            ("ED5", "extended_data_figure5_external_anchor_robustness", "External paired-ST and robustness"),
            ("ED12", "extended_data_figure12_distance_to_caf_core_dynamics", "Distance dynamics"),
            ("ED19", "extended_data_figure19_random_core_null_diagnostics", "Random-core null diagnostics"),
            ("ED25", "extended_data_figure25_spatial_robustness_module_nc_style", "NC-style spatial robustness module"),
        ],
    },
    {
        "module": "Supplementary Module 2",
        "short_name": "metastatic_immune_decoupling",
        "title": "Metastatic-site remodeling and lymph-node immune decoupling",
        "core_question": "Does metastatic site preserve stromal-tumor coupling while changing immune coupling to CAF cores?",
        "recommended_final": "Use ED27 as the module backbone; retain ED2 and ED20 as large-scale visual/ecotype detail.",
        "current_figures": [
            ("ED2", "figure2_supplement_submission_spatial_examples", "Representative GSE272362 spatial examples"),
            ("ED20", "extended_data_figure20_ecotype_context_flow", "Ecotype context flow"),
            ("ED27", "extended_data_figure27_metastatic_immune_decoupling_module_nc_style", "NC-style metastatic immune-decoupling module"),
        ],
    },
    {
        "module": "Supplementary Module 3",
        "short_name": "cell_state_multiresolution_validation",
        "title": "Cell-state interpretation and multi-resolution validation",
        "core_question": "Do marker, reference-projection, NNLS and Xenium layers support the cellular interpretation of CAF-myeloid cores?",
        "recommended_final": "Use ED26 plus ED28 as the module backbone; retain ED10, ED11, ED13 and ED18 as validation/detail plates.",
        "current_figures": [
            ("ED6", "extended_data_figure6_cell_state_reference_support", "Marker and reference-projection support"),
            ("ED10", "extended_data_figure10_gse274557_external_validation", "Independent GSE274557 Visium validation"),
            ("ED11", "extended_data_figure11_gse274673_xenium_cell_resolution", "GSE274673 Xenium validation"),
            ("ED13", "extended_data_figure13_xenium_cell_domain_maps", "Representative Xenium cell maps"),
            ("ED18", "extended_data_figure18_xenium_program_neighborhoods", "Xenium program neighborhoods"),
            ("ED26", "extended_data_figure26_cell_state_reference_xenium_module_nc_style", "NC-style cell-state/reference/Xenium module"),
            ("ED28", "extended_data_figure28_strict_nnls_deconvolution_sensitivity", "Strict NNLS deconvolution sensitivity"),
        ],
    },
    {
        "module": "Supplementary Module 4",
        "short_name": "mechanism_interface_priority",
        "title": "Mechanism, interface biology and perturbation-priority axes",
        "core_question": "Which CAF-core/interface axes are strongest candidates for perturbation follow-up?",
        "recommended_final": "Use ED29 as the module backbone; retain ED3, ED8, ED9, ED16 and ED21 as gene/interface/triangulation detail.",
        "current_figures": [
            ("ED3", "figure3_supplement_targeted_gene_axis_validation", "Targeted gene-level support"),
            ("ED8", "extended_data_figure8_cxcl9_spp1_polarity", "SPP1/CXCL9 polarity"),
            ("ED9", "extended_data_figure9_focused_interface_axes", "Focused interface axes"),
            ("ED16", "extended_data_figure16_interface_compartment_maps", "Interface compartment maps"),
            ("ED21", "extended_data_figure21_mechanism_triangulation_priority", "Mechanism triangulation"),
            ("ED29", "extended_data_figure29_mechanism_gene_interface_module_nc_style", "NC-style mechanism gene/interface module"),
        ],
    },
    {
        "module": "Supplementary Module 5",
        "short_name": "pathology_clinical_tls_boundaries",
        "title": "Pathology bridge, TCGA context and TLS/claim boundaries",
        "core_question": "Which translational or alternative interpretations are supported, and which are explicitly bounded?",
        "recommended_final": "Use ED24 as the reviewer-risk backbone; retain ED4, ED7, ED17, ED22 and ED23 as pathology/TCGA/TLS detail.",
        "current_figures": [
            ("ED4", "extended_data_figure4_he_morphology_bridge", "Exploratory H&E morphology bridge"),
            ("ED7", "extended_data_tcga_paad_bulk_context", "TCGA PAAD bulk context"),
            ("ED17", "extended_data_figure17_he_patch_examples", "Representative H&E patches"),
            ("ED22", "extended_data_figure22_tcga_survival_context", "TCGA survival context"),
            ("ED23", "extended_data_figure23_tls_maturity_stress_test", "TLS maturity stress test"),
            ("ED24", "extended_data_figure24_review_risk_resolution_nc_style", "NC-style reviewer-risk resolution module"),
        ],
    },
    {
        "module": "Supplementary Module 6",
        "short_name": "atlas_map_plates",
        "title": "Atlas and inspectable spatial map plates",
        "core_question": "Can readers visually inspect the tissue-scale spatial ecology across representative contexts?",
        "recommended_final": "Use this as a map plate module only if the journal allows large supplementary figures; otherwise keep as source-detail figures.",
        "current_figures": [
            ("ED14", "extended_data_figure14_spatial_atlas_overview", "Spatial atlas overview"),
            ("ED15", "extended_data_figure15_local_spatial_program_maps", "Local spatial program maps"),
            ("ED16", "extended_data_figure16_interface_compartment_maps", "Interface compartment maps"),
            ("ED13", "extended_data_figure13_xenium_cell_domain_maps", "Xenium cell-domain maps"),
        ],
    },
]


def load_font(size: int):
    for name in ("arial.ttf", "DejaVuSans.ttf"):
        try:
            return ImageFont.truetype(name, size)
        except Exception:
            pass
    return ImageFont.load_default()


def make_contact_sheet(module: dict) -> Path:
    QC_DIR.mkdir(parents=True, exist_ok=True)
    entries = module["current_figures"]
    thumb_w, thumb_h = 520, 360
    label_h = 52
    pad = 26
    cols = 2
    rows = (len(entries) + cols - 1) // cols
    title_h = 92
    sheet = Image.new(
        "RGB",
        (cols * (thumb_w + pad) + pad, title_h + rows * (thumb_h + label_h + pad) + pad),
        "white",
    )
    draw = ImageDraw.Draw(sheet)
    title_font = load_font(24)
    label_font = load_font(15)
    draw.text((pad, 18), module["module"] + ": " + module["title"], fill=(0, 0, 0), font=title_font)
    draw.text((pad, 52), module["core_question"], fill=(60, 60, 60), font=label_font)
    for i, (label, base, caption) in enumerate(entries):
        png = FIG_DIR / f"{base}.png"
        col = i % cols
        row = i // cols
        x = pad + col * (thumb_w + pad)
        y = title_h + row * (thumb_h + label_h + pad)
        text = f"{label}: {caption}"
        draw.text((x, y), text, fill=(0, 0, 0), font=label_font)
        canvas = Image.new("RGB", (thumb_w, thumb_h), "white")
        if png.exists():
            im = Image.open(png).convert("RGB")
            im.thumbnail((thumb_w, thumb_h), Image.Resampling.LANCZOS)
            canvas.paste(im, ((thumb_w - im.width) // 2, (thumb_h - im.height) // 2))
        else:
            draw2 = ImageDraw.Draw(canvas)
            draw2.text((20, 20), f"Missing: {png.name}", fill=(160, 0, 0), font=label_font)
        sheet.paste(canvas, (x, y + label_h))
        draw.rectangle([x, y + label_h, x + thumb_w, y + label_h + thumb_h], outline=(190, 190, 190), width=1)
    out = QC_DIR / f"{module['short_name']}_contact_sheet.png"
    sheet.save(out, dpi=(180, 180))
    return out


def main() -> None:
    rows = []
    contact_sheets = []
    for module in MODULES:
        contact = make_contact_sheet(module)
        contact_sheets.append((module, contact))
        for order, (label, base, caption) in enumerate(module["current_figures"], start=1):
            rows.append(
                {
                    "module": module["module"],
                    "module_short_name": module["short_name"],
                    "module_title": module["title"],
                    "core_question": module["core_question"],
                    "recommended_final_structure": module["recommended_final"],
                    "panel_order": order,
                    "current_figure": label,
                    "current_base": base,
                    "current_pdf": str((FIG_DIR / f"{base}.pdf").relative_to(ROOT)),
                    "current_png": str((FIG_DIR / f"{base}.png").relative_to(ROOT)),
                    "component_role": caption,
                    "contact_sheet": str(contact.relative_to(ROOT)),
                }
            )
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TABLE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    report = []
    report.append("# Supplement Module Reorganization Plan")
    report.append("")
    report.append("Last updated: 2026-06-28")
    report.append("")
    report.append("## Short Answer")
    report.append("")
    report.append("Yes. The Extended Data package can be reorganized into a smaller number of module-style multi-panel figures. This is preferable to a linear ED1-ED29 list because it makes each supplement answer a reviewer-facing question.")
    report.append("")
    report.append("The recommended structure is not to stitch all existing ED figures into a few unreadable mega-panels. Instead, use module backbone figures, regenerate key panels from source data where needed, and retain map/detail plates only where visual inspection is necessary.")
    report.append("")
    report.append("## Recommended Modules")
    report.append("")
    for module, contact in contact_sheets:
        report.append(f"### {module['module']}: {module['title']}")
        report.append("")
        report.append(f"**Core question:** {module['core_question']}")
        report.append("")
        report.append(f"**Recommended final structure:** {module['recommended_final']}")
        report.append("")
        report.append(f"**Contact sheet:** `{contact.relative_to(ROOT)}`")
        report.append("")
        report.append("| current figure | role in module |")
        report.append("|---|---|")
        for label, _base, caption in module["current_figures"]:
            report.append(f"| {label} | {caption} |")
        report.append("")
    report.append("## Practical Recommendation")
    report.append("")
    report.append("For submission, the strongest configuration is:")
    report.append("")
    report.append("- Keep 4 main figures.")
    report.append("- Replace the long ED1-ED29 logic with 5-6 module-style Extended Data figures if the journal format permits.")
    report.append("- Use ED25-ED29 as the current backbone for robustness, cell-state, metastatic-decoupling, NNLS and mechanism modules.")
    report.append("- Keep atlas/map plates as high-resolution supplementary figures or source-detail figures because excessive downscaling would weaken readability.")
    report.append("- Do not submit Module 6 as a redundant numbered figure if its map panels are already included in Modules 3-4. Use it as a high-resolution map plate or source-detail appendix.")
    report.append("")
    report.append("## Recommended Renumbering If Fully Executed")
    report.append("")
    report.append("| new supplement figure | module | recommended content |")
    report.append("|---|---|---|")
    report.append("| Extended Data Fig. 1 | Spatial specificity and robustness | ED25 backbone plus selected ED1/ED5/ED12/ED19 panels |")
    report.append("| Extended Data Fig. 2 | Metastatic-site remodeling | ED27 backbone plus selected ED2/ED20 panels |")
    report.append("| Extended Data Fig. 3 | Cell-state and multi-resolution validation | ED26/ED28 backbone plus selected ED10/ED11/ED13/ED18 panels |")
    report.append("| Extended Data Fig. 4 | Mechanism and interface priority | ED29 backbone plus selected ED3/ED8/ED9/ED16/ED21 panels |")
    report.append("| Extended Data Fig. 5 | Pathology, TCGA and TLS boundaries | ED24 backbone plus selected ED4/ED7/ED17/ED22/ED23 panels |")
    report.append("| Supplementary Fig. or Source Data Plate | Atlas/map plates | ED14/ED15 and any non-duplicated high-resolution map panels |")
    report.append("")
    report.append("This renumbering should be executed only after final journal-format decisions, because it requires updating all in-text citations, legends, source-data filenames and package indices.")
    report.append("")
    report.append("## Next Execution Step")
    report.append("")
    report.append("The next code step should create true regenerated module figures from source data for Modules 1-5, rather than assembling thumbnails. The contact sheets generated here are planning/QC artifacts, not final publication figures.")
    report.append("")
    report.append("## Outputs")
    report.append("")
    report.append(f"- Manifest: `{OUT_TABLE.relative_to(ROOT)}`")
    for _module, contact in contact_sheets:
        report.append(f"- Contact sheet: `{contact.relative_to(ROOT)}`")
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
