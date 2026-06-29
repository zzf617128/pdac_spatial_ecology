from __future__ import annotations

from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
MANUSCRIPT_IN = ROOT / "results" / "manuscript" / "pdac_caf_myeloid_spatial_niche_submission_v2_full.md"
MANUSCRIPT_OUT = ROOT / "results" / "manuscript" / "pdac_caf_myeloid_spatial_niche_nature_subjournal.md"
SI_OUT = ROOT / "results" / "manuscript" / "pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.md"
REPORT_OUT = ROOT / "results" / "reports" / "nature_subjournal_submission_readiness_2026_06_28.md"
RENUMBER_OUT = ROOT / "results" / "tables" / "nature_subjournal_display_item_renumbering.csv"


MAIN_FIGURE_LEGENDS = """## Figure legends

### Figure 1 | CAF-myeloid cores define reproducible spatial organizing regions in PDAC.

This figure tests whether CAF-myeloid-high regions behave as spatial organizing cores rather than arbitrary same-size tissue regions. **a,** Cohort scale and evidence roles across public PDAC ST-H&E, validation and cell-resolution datasets. **b,** Random-core specificity across GSE282302, GSE274103 and GSE235315. Bars show the median difference between the observed distance-to-core Spearman correlation and the random-core median. Labels indicate the number of samples in which the observed CAF-core gradient was more negative than the matched random-core median. **c,** CAF-core threshold sensitivity across top 15%, top 10% and top 5% within-sample CAF-core definitions. **d,** Representative same-size random-core null intervals showing observed CAF-core gradients against matched random-core distributions. **e,** Median program decay from CAF-core to far regions across discovery/support cohorts. **f-j,** Representative post-neoadjuvant GSE282302 section showing H&E, CAF-core rings, CAF-myeloid, IFN/MHC, tumor-aggressive and immune-core spatial programs. The figure supports spatial organization around CAF-myeloid cores but not causality, mature TLS formation or clinical outcome prediction.

### Figure 2 | Lymph-node metastases uncouple immune programs from CAF-myeloid cores.

This figure tests whether metastatic site changes the immune programs organized around CAF-myeloid cores. **a,** GSE272362 validation atlas composition by specimen site. **b,** Random-core validation across primary tumors, liver metastases and lymph-node metastases. Primary tumors and liver metastases validate CAF-core-centered IFN/MHC, immune-core and tumor-aggressive organization. Lymph-node metastases retain tumor-aggressive CAF association but show inconsistent CAF-core-centered IFN/MHC and immune-core organization. **c,** CAF-core subprogram decomposition across specimen sites. Negative values indicate enrichment nearer to the CAF core. **d,** Immune-decoupling index across tissue contexts, defined as stromal-tumor CAF-core enrichment minus immune CAF-core enrichment. **e-p,** Representative primary, liver-metastatic and lymph-node-metastatic sections showing H&E plus CAF-core rings, IFN/MHC, immune-core and tumor-aggressive programs. The figure supports a lymph-node immune-decoupling lead, not a clinical subtype or causal mechanism.

### Figure 3 | CAF-core ecotypes expose invasive-interface and immune-coupling axes.

This figure asks what biological states are encoded by the CAF-core architecture. **a,** CAF-core NMF ecotype loadings across stromal, myeloid, immune and tumor-state programs. **b,** Context-to-ecotype architecture summarizing the distribution of dominant CAF-core ecotypes across analyzed tissue contexts. **c,** Dominant CAF-core ecotype composition by cohort or specimen context. **d,** Immune-decoupling index by dominant CAF-core ecotype, defined as stromal-tumor CAF-core coupling minus immune CAF-core coupling. **e,** Candidate-axis enrichment in CAF cores and tumor-stroma interfaces. **f,** Spearman correlations between candidate-axis CAF-core enrichment and the immune-decoupling index across samples. **g-r,** Representative primary tumor, liver-metastasis and lymph-node-metastasis spatial maps showing compartment assignments, SPP1/TAM, TGF-beta and tumor-aggressive programs; cyan outlines mark tumor-stroma interface spots. The figure nominates pathway-level axes and representative spatial states for follow-up validation, not causal ligand-receptor signaling.

### Figure 4 | Independent Visium and Xenium data validate CAF-domain organization at complementary resolutions.

This figure tests whether the CAF-domain model generalizes across independent Visium sections and cell-resolution Xenium data. **a,** External validation scale across GSE274557 Visium primary and metastatic tissue contexts and GSE274673 Xenium treatment contexts. **b,** GSE274557 Visium CAF-core validation across primary PDAC, liver metastasis, lung metastasis and peritoneal metastasis. Tiles show median observed-minus-random distance-to-CAF-core Spearman rho after 1,000 same-size random cores per sample, with labels showing support counts. **c,** GSE274673 Xenium fixed-anchor cohort summary for CAF-APC and CAF-SPP1/TAM anchors. **d,** CAF-SPP1/TAM anchor deltas across all six Xenium sections. Negative values indicate stronger target-program centering around CAF-domain anchors than around matched random cell anchors. **e-l,** Representative treatment-naive and chemoradiotherapy-treated Xenium sections showing CAF-SPP1/TAM anchor, SPP1/TAM, IFN/APC and tumor epithelial programs at cell resolution. The figure supports independent Visium and cell-resolution validation of CAF-domain immune/myeloid organization but does not establish lymph-node-specific validation in GSE274557, direct CAF-to-tumor epithelial proximity in GSE274673 or causal signaling.
"""


EXTENDED_DATA_LEGENDS = """## Extended Data figure legends

### Extended Data Figure 1 | Spatial specificity and robustness of CAF-myeloid cores.

Regenerated module-style synthesis of the spatial-specificity evidence. **a,** Cross-cohort random-core specificity for CAF-core-centered target programs. **b,** CAF-core threshold sensitivity across top 15%, top 10% and top 5% definitions. **c,d,** Distance-gradient and core-to-far summaries from CAF-core-proximal to distal regions. **e,f,** Independent GSE274557 Visium and GSE274673 Xenium support layers. **g,** Cross-context support fraction. This module supports spatial specificity and robustness to null-core and threshold perturbations; it remains observational and program-defined.

### Extended Data Figure 2 | Metastatic-site remodeling and lymph-node immune decoupling.

Regenerated module-style synthesis of the GSE272362 metastatic-site analysis. **a,** Site-level scale across primary tumors, liver metastases and lymph-node metastases. **b,** Same-size random-core support for IFN/MHC, immune-core and tumor-aggressive programs by site. **c,** CAF-core subprogram distance gradients. **d,** Stromal-tumor coupling, immune-core coupling and immune-decoupling index. **e,f,** Patient-matched primary-to-metastasis deltas where available. **g,** Interpretation boundary for the lymph-node result. This module supports a metastatic-site immune-remodeling lead, not a clinical subtype, prognosis marker or causal mechanism.

### Extended Data Figure 3 | Cell-state interpretation and multi-resolution validation.

Regenerated module-style synthesis of marker, reference-projection, strict NNLS and Xenium evidence. **a,b,** Marker-state CAF-core enrichment and spot-level association with CAF-myeloid score. **c,d,** Marker-level, full-reference and strict NNLS associations with immune decoupling. **e,** Small-reference versus full-reference projection stability. **f,g,** Xenium CAF-domain support and signature/sample coverage. **h,** Claim ladder distinguishing supported cell-state interpretation from unsupported single-cell-resolved ground truth. This module strengthens cellular interpretation but does not replace immunostaining, image segmentation or orthogonal cell-abundance validation.

### Extended Data Figure 4 | Mechanism, interface biology and perturbation-priority axes.

Regenerated module-style synthesis of mechanism nomination. **a,** Evidence-layer matrix across matrix-integrin, SPP1-CD44/integrin, TGF-beta/TGFBR and IL6-OSM/LIF-JAKSTAT. **b,** Perturbation-priority rank. **c,d,** Targeted-gene CAF-core and interface enrichment across contexts. **e,** Focused ligand-core, receptor-interface, response-interface and directional metrics. **f,** Directional-score associations with immune decoupling and stromal-tumor coupling. **g,** External-context anchors. **h,** Claim ladder for mechanism follow-up. This module nominates matrix-integrin and SPP1-CD44/integrin as leading perturbation-ready axes, but it does not establish causal ligand-receptor signaling.

### Extended Data Figure 5 | Pathology bridge, TCGA context and TLS/claim boundaries.

Regenerated module-style synthesis of translational context and claim boundaries. **a,b,** H&E patch-model performance and feature-direction map. **c,d,** TCGA PAAD bulk-axis correlation and exploratory univariable survival context. **e,f,** TLS-compatibility and immune-compartment stress testing. **g,** Boundary map for reviewer-risk claims. This module supports translational plausibility and explicit claim control; it does not establish mature TLS biology, clinical-grade pathology prediction, prognosis, therapy response, spatial localization from TCGA or causal signaling.
"""


STATEMENTS = """## Acknowledgements

[Author to complete: funding bodies, grant numbers, institutional support and data/resource acknowledgements.]

## Author contributions

[Author to complete using CRediT roles, for example: Conceptualization, Data curation, Formal analysis, Funding acquisition, Investigation, Methodology, Project administration, Resources, Software, Supervision, Validation, Visualization, Writing - original draft, Writing - review and editing.]

## Competing interests

[Author to complete. If none: The authors declare no competing interests.]
"""


def section(text: str, start: str, end: str | None = None) -> str:
    start_pat = re.escape(start)
    if end is None:
        m = re.search(start_pat + r"\n(.*)", text, re.S)
    else:
        m = re.search(start_pat + r"\n(.*?)(?=\n" + re.escape(end) + r"\n)", text, re.S)
    return m.group(0 if end is None else 0).strip() if m else ""


def extract_until(text: str, end_heading: str) -> str:
    idx = text.index("\n" + end_heading + "\n")
    return text[:idx].rstrip()


def build_manuscript() -> str:
    text = MANUSCRIPT_IN.read_text(encoding="utf-8")
    main = extract_until(text, "## Main figure legends")
    main = main.replace("## Data and code availability", "## Data availability")
    split_token = "\nPrimary scripts include "
    if split_token in main:
        data_part, code_part = main.split(split_token, 1)
        code_part = "Primary scripts include " + code_part
        main = (
            data_part.rstrip()
            + "\n\n## Code availability\n\n"
            + code_part.strip()
            + "\n\nCustom analysis scripts are currently indexed locally in `results/reports/submission_package_index.md`. "
            "For submission, the analysis code should be deposited in a persistent public repository with a DOI or made available to editors and reviewers through the journal system."
        )
    out = "\n\n".join([main, STATEMENTS, MAIN_FIGURE_LEGENDS, EXTENDED_DATA_LEGENDS])
    return out + "\n"


def build_si() -> str:
    return """# Supplementary Information

## Supplementary Figure Guide

For a Nature Communications-style submission, the regenerated module figures can be submitted as Supplementary Figures 1-6 instead of Extended Data Figures 1-6. The scientific content and source-data mapping should remain unchanged.

| Nature-family item | Nature Communications equivalent | file |
|---|---|---|
| Extended Data Fig. 1 | Supplementary Fig. 1 | `results/figures/submission/supplementary_module1_spatial_specificity_robustness.pdf` |
| Extended Data Fig. 2 | Supplementary Fig. 2 | `results/figures/submission/supplementary_module2_metastatic_immune_decoupling.pdf` |
| Extended Data Fig. 3 | Supplementary Fig. 3 | `results/figures/submission/supplementary_module3_cell_state_multiresolution_validation.pdf` |
| Extended Data Fig. 4 | Supplementary Fig. 4 | `results/figures/submission/supplementary_module4_mechanism_interface_priority.pdf` |
| Extended Data Fig. 5 | Supplementary Fig. 5 | `results/figures/submission/supplementary_module5_pathology_tcga_tls_boundaries.pdf` |
| Extended Data Fig. 6 | Supplementary Fig. 6 | `results/figures/submission/supplementary_module6_spatial_architecture_mechanism_deepening.pdf` |

## Source-Detail Figure Archive

The previous ED1-ED33 figure suite remains a provenance and inspection archive. It should not be submitted as 33 active Extended Data figures for a Nature-family journal unless the editor specifically requests additional source-detail plates.

Primary source mapping is maintained in `results/reports/submission_figure_captions_and_source_map.md`.
"""


def build_renumbering() -> str:
    rows = [
        ("Main Figure 1", "Figure 1", "results/figures/submission/figure1_submission_spatial_specificity_nc_style.pdf", "Core discovery and robustness"),
        ("Main Figure 2", "Figure 2", "results/figures/submission/figure2_submission_metastatic_decoupling_nc_style.pdf", "Metastatic-site contrast"),
        ("Main Figure 3", "Figure 3", "results/figures/submission/figure3_submission_ecotypes_mechanism_axes_nc_style.pdf", "Ecotypes and interface biology"),
        ("Main Figure 4", "Figure 4", "results/figures/submission/figure4_submission_multiresolution_validation_nc_style.pdf", "Independent Visium/Xenium validation"),
        ("Supplementary Module 1", "Extended Data Fig. 1 / Supplementary Fig. 1", "results/figures/submission/supplementary_module1_spatial_specificity_robustness.pdf", "Spatial specificity and robustness"),
        ("Supplementary Module 2", "Extended Data Fig. 2 / Supplementary Fig. 2", "results/figures/submission/supplementary_module2_metastatic_immune_decoupling.pdf", "Metastatic immune decoupling"),
        ("Supplementary Module 3", "Extended Data Fig. 3 / Supplementary Fig. 3", "results/figures/submission/supplementary_module3_cell_state_multiresolution_validation.pdf", "Cell-state validation"),
        ("Supplementary Module 4", "Extended Data Fig. 4 / Supplementary Fig. 4", "results/figures/submission/supplementary_module4_mechanism_interface_priority.pdf", "Mechanism priority"),
        ("Supplementary Module 5", "Extended Data Fig. 5 / Supplementary Fig. 5", "results/figures/submission/supplementary_module5_pathology_tcga_tls_boundaries.pdf", "Pathology/TCGA/TLS boundaries"),
    ]
    lines = ["current_item,nature_family_item,file,role"]
    lines += [",".join(f'"{x}"' for x in row) for row in rows]
    return "\n".join(lines) + "\n"


def build_report() -> str:
    return """# Nature-Subjournal Submission Readiness Plan

Last updated: 2026-06-28

## Target Interpretation

Default target: Nature Communications / Nature-family oncology or biology journal initial submission.

The active display set should be constrained to:

- 4 main figures.
- 6 regenerated multi-panel module figures, submitted either as Extended Data Figures 1-6 for Nature-family research journals or as Supplementary Figures 1-6 for Nature Communications.
- Source Data files for every main and module figure.

## Actions Completed

- Created `results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal.md`.
- Created `results/manuscript/pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.md`.
- Created `results/tables/nature_subjournal_display_item_renumbering.csv`.
- Reduced the active supplementary display set from a 33-figure Extended Data archive to 6 regenerated module figures.
- Preserved the ED1-ED33 archive as provenance/source-detail material, not the active submission structure.

## Requirements Mapped To Local Package

| Nature-family requirement | Local status | Next action |
|---|---|---|
| Cover letter | Not yet drafted in final form | Draft concise cover letter with novelty, fit and boundaries |
| Manuscript text file | Prepared as Markdown and Word export candidate | Export Nature-subjournal Word file |
| Individual figure files | Available as PDF/SVG/PNG | Submit PDF/SVG or journal-requested production formats |
| Supplementary Information | Prepared as SI guide | Convert to Word if targeting Nature Communications |
| Source Data | Available for main and module figures | Package into a Source Data folder/zip |
| Data availability | Present but should cite persistent public repositories and dataset identifiers | Add formal dataset citations and repository accession links |
| Code availability | Present as local script index, but not DOI-backed | Deposit code or prepare private reviewer link |
| Reporting Summary | Required for life-science submission | Fill Bio/life sciences reporting summary; add sample-size, randomization, blinding and statistics answers |
| Author contributions | Placeholder only | Author input required |
| Competing interests | Placeholder only | Author input required |
| Acknowledgements/funding | Placeholder only | Author input required |

## Main Claim Boundary

Safe central claim: CAF-myeloid stromal cores organize inflammatory, immune-core and tumor-aggressive spatial programs in PDAC, with metastatic-site remodeling and prioritized interface candidates for perturbation follow-up.

Do not claim mature TLS biology, clinical-grade prediction, survival/prognosis validation, treatment response, spatial localization from TCGA, direct ligand-receptor signaling or experimentally proven causality.
"""


def main() -> None:
    MANUSCRIPT_OUT.parent.mkdir(parents=True, exist_ok=True)
    REPORT_OUT.parent.mkdir(parents=True, exist_ok=True)
    RENUMBER_OUT.parent.mkdir(parents=True, exist_ok=True)
    MANUSCRIPT_OUT.write_text(build_manuscript(), encoding="utf-8")
    SI_OUT.write_text(build_si(), encoding="utf-8")
    RENUMBER_OUT.write_text(build_renumbering(), encoding="utf-8")
    REPORT_OUT.write_text(build_report(), encoding="utf-8")
    print(f"Wrote {MANUSCRIPT_OUT}")
    print(f"Wrote {SI_OUT}")
    print(f"Wrote {RENUMBER_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
