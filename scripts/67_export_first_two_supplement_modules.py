from __future__ import annotations

import runpy
import shutil
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"
TABLE_DIR = ROOT / "results" / "tables"

MODULE_DEFS = [
    {
        "module_id": "Supplementary Module 1",
        "base": "supplementary_module1_spatial_specificity_robustness",
        "title": "Spatial specificity and robustness of CAF-myeloid cores",
        "script": ROOT / "scripts" / "60_make_nc_style_robustness_module.py",
        "backbone_base": "extended_data_figure25_spatial_robustness_module_nc_style",
        "source": "Source_Data_Extended_Data_Fig_25_spatial_robustness_module.csv",
        "report": "extended_data_figure25_spatial_robustness_module_notes.md",
        "component_figures": [
            "figure1_supplement_submission_post_nact_spatial_example",
            "extended_data_figure5_external_anchor_robustness",
            "extended_data_figure12_distance_to_caf_core_dynamics",
            "extended_data_figure19_random_core_null_diagnostics",
            "extended_data_figure25_spatial_robustness_module_nc_style",
        ],
        "claim": "CAF-myeloid cores are spatially specific organizing regions robust to same-size random-core controls, CAF-core threshold changes and independent validation layers.",
        "boundary": "This module does not establish perturbational causality, clinical prediction or histology-annotated compartment identity.",
    },
    {
        "module_id": "Supplementary Module 2",
        "base": "supplementary_module2_metastatic_immune_decoupling",
        "title": "Metastatic-site remodeling and lymph-node immune decoupling",
        "script": ROOT / "scripts" / "62_make_nc_style_metastatic_decoupling_module.py",
        "backbone_base": "extended_data_figure27_metastatic_immune_decoupling_module_nc_style",
        "source": "Source_Data_Extended_Data_Fig_27_metastatic_immune_decoupling_module.csv",
        "report": "extended_data_figure27_metastatic_immune_decoupling_module_notes.md",
        "component_figures": [
            "figure2_supplement_submission_spatial_examples",
            "extended_data_figure20_ecotype_context_flow",
            "extended_data_figure27_metastatic_immune_decoupling_module_nc_style",
        ],
        "claim": "Primary tumors and liver metastases preserve CAF-core immune/tumor organization, whereas lymph-node metastases retain stromal-tumor coupling while weakening immune/IFN coupling.",
        "boundary": "The lymph-node result is a spatial decoupling contrast from five lymph-node metastasis samples, not a clinical subtype or causal mechanism.",
    },
]


def copy_backbone(module: dict) -> None:
    src_base = FIG_DIR / module["backbone_base"]
    dst_base = FIG_DIR / module["base"]
    for ext in ("pdf", "svg", "png"):
        src = src_base.with_suffix(f".{ext}")
        dst = dst_base.with_suffix(f".{ext}")
        if not src.exists():
            raise FileNotFoundError(src)
        shutil.copy2(src, dst)

    src_source = SOURCE_DIR / module["source"]
    dst_source = SOURCE_DIR / f"Source_Data_{module['base']}.csv"
    if not src_source.exists():
        raise FileNotFoundError(src_source)
    shutil.copy2(src_source, dst_source)


def write_report(module: dict) -> Path:
    report_path = REPORT_DIR / f"{module['base']}_notes.md"
    components = "\n".join(
        f"- `{(FIG_DIR / (component + '.pdf')).relative_to(ROOT)}`"
        for component in module["component_figures"]
    )
    report_path.write_text(
        f"# {module['module_id']} Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        f"## Title\n\n{module['title']}\n\n"
        "## Role\n\n"
        "This is a module-style replacement candidate for the long linear Extended Data figure sequence. "
        "It uses a regenerated NC-style backbone figure and records the component/detail figures that can be retained as high-resolution supplementary plates if needed.\n\n"
        f"## Core claim\n\n{module['claim']}\n\n"
        f"## Boundary\n\n{module['boundary']}\n\n"
        "## Module backbone\n\n"
        f"- `{(FIG_DIR / (module['base'] + '.pdf')).relative_to(ROOT)}`\n"
        f"- `{(FIG_DIR / (module['base'] + '.svg')).relative_to(ROOT)}`\n"
        f"- `{(FIG_DIR / (module['base'] + '.png')).relative_to(ROOT)}`\n"
        f"- `{(SOURCE_DIR / ('Source_Data_' + module['base'] + '.csv')).relative_to(ROOT)}`\n\n"
        "## Component/detail figures retained for optional plates\n\n"
        f"{components}\n",
        encoding="utf-8",
    )
    return report_path


def main() -> None:
    rows = []
    for module in MODULE_DEFS:
        runpy.run_path(str(module["script"]), run_name="__main__")
        copy_backbone(module)
        report = write_report(module)
        for order, component in enumerate(module["component_figures"], start=1):
            rows.append(
                {
                    "module_id": module["module_id"],
                    "module_title": module["title"],
                    "module_base": module["base"],
                    "module_pdf": str((FIG_DIR / f"{module['base']}.pdf").relative_to(ROOT)),
                    "module_source_data": str((SOURCE_DIR / f"Source_Data_{module['base']}.csv").relative_to(ROOT)),
                    "module_notes": str(report.relative_to(ROOT)),
                    "component_order": order,
                    "component_base": component,
                    "component_pdf": str((FIG_DIR / f"{component}.pdf").relative_to(ROOT)),
                    "component_role": "backbone" if component == module["backbone_base"] else "optional_detail_plate",
                    "core_claim": module["claim"],
                    "claim_boundary": module["boundary"],
                }
            )
    out = TABLE_DIR / "supplementary_modules_1_2_export_manifest.csv"
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"Wrote {out}")
    for module in MODULE_DEFS:
        print(f"Wrote {FIG_DIR / (module['base'] + '.pdf')}")


if __name__ == "__main__":
    main()
