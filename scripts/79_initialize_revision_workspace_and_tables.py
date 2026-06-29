from __future__ import annotations

import ast
import csv
import shutil
from datetime import date
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT / "results"
REVISION = RESULTS / "revision_2026_06_29"
TABLES = RESULTS / "tables"
MANUSCRIPT = RESULTS / "manuscript"
STAGED = RESULTS / "nature_subjournal_submission"


def ensure_dirs() -> None:
    for rel in [
        "manuscript",
        "figures",
        "supplementary_tables",
        "analysis_outputs",
        "docs",
    ]:
        (REVISION / rel).mkdir(parents=True, exist_ok=True)


def copy_if_exists(src: Path, dst: Path) -> None:
    if src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


def copy_current_package_files() -> None:
    copy_if_exists(
        MANUSCRIPT / "pdac_caf_myeloid_spatial_niche_nature_subjournal.md",
        REVISION / "manuscript" / "Manuscript_NatureSubjournal_current.md",
    )
    copy_if_exists(
        MANUSCRIPT / "pdac_caf_myeloid_spatial_niche_nature_subjournal.docx",
        REVISION / "manuscript" / "Manuscript_NatureSubjournal_current.docx",
    )
    copy_if_exists(
        MANUSCRIPT / "pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.md",
        REVISION / "manuscript" / "Supplementary_Information_current.md",
    )
    copy_if_exists(
        MANUSCRIPT / "pdac_caf_myeloid_spatial_niche_nature_subjournal_supplementary_information.docx",
        REVISION / "manuscript" / "Supplementary_Information_current.docx",
    )
    fig_dir = STAGED / "figures"
    if fig_dir.exists():
        for fig in fig_dir.glob("*.pdf"):
            copy_if_exists(fig, REVISION / "figures" / fig.name.replace(".pdf", "_current.pdf"))


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def parse_signatures() -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    for line in (PROJECT / "config" / "signatures.yaml").read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        name, genes = line.split(":", 1)
        genes = genes.strip()
        if genes.startswith("[") and genes.endswith("]"):
            signatures[name.strip()] = [
                item.strip().strip("'\"")
                for item in genes.strip("[]").split(",")
                if item.strip()
            ]
    return signatures


def build_dataset_table() -> None:
    rows = [
        {
            "dataset": "GSE282302",
            "accession": "GSE282302",
            "platform": "10x Visium ST-H&E",
            "sample_type": "post-NACT PDAC sections",
            "n_sections_or_samples": 108,
            "treatment_context": "post-NACT / FOLFIRINOX-context metadata",
            "tissue_sites": "primary PDAC sections",
            "role_in_study": "discovery cohort for CAF-myeloid core, same-size random-core null, threshold sensitivity and H&E bridge",
            "notes": "Used as section-level spatial evidence; patient-level clinical outcome not inferred.",
        },
        {
            "dataset": "GSE274103",
            "accession": "GSE274103",
            "platform": "10x Visium ST-H&E",
            "sample_type": "treatment-naive PDAC sections",
            "n_sections_or_samples": 5,
            "treatment_context": "treatment-naive",
            "tissue_sites": "primary PDAC",
            "role_in_study": "treatment-naive support for CAF-core spatial organization and H&E bridge",
            "notes": "Small support cohort; interpreted directionally.",
        },
        {
            "dataset": "GSE235315",
            "accession": "GSE235315",
            "platform": "10x Visium paired ST-H&E",
            "sample_type": "PDAC paired spatial samples",
            "n_sections_or_samples": 7,
            "treatment_context": "metadata pending",
            "tissue_sites": "PDAC",
            "role_in_study": "external paired-ST anchor for CAF-core random-core support",
            "notes": "Used as spatial-state support, not patient-level or treatment-context evidence.",
        },
        {
            "dataset": "GSE272362",
            "accession": "GSE272362 / Zenodo PDAC_Updated.rds",
            "platform": "spatial transcriptomics RDS-derived atlas",
            "sample_type": "primary, liver metastasis, lymph-node metastasis and normal pancreas specimens",
            "n_sections_or_samples": 30,
            "treatment_context": "mixed public atlas context",
            "tissue_sites": "primary tumor; liver metastasis; lymph-node metastasis; normal pancreas",
            "role_in_study": "independent primary/metastatic validation and lymph-node immune-decoupling lead",
            "notes": "LN metastasis subset has n=5 and must be framed as hypothesis-generating.",
        },
        {
            "dataset": "GSE274557",
            "accession": "GSE274557",
            "platform": "10x Visium",
            "sample_type": "primary and metastatic PDAC sections",
            "n_sections_or_samples": 55,
            "treatment_context": "public metastatic PDAC atlas",
            "tissue_sites": "primary PDAC; liver metastasis; lung metastasis; peritoneal metastasis",
            "role_in_study": "external multi-site Visium validation of broad CAF-core organization",
            "notes": "Does not test lymph-node-specific immune decoupling because lymph-node metastases were not included in the analyzed Visium set.",
        },
        {
            "dataset": "GSE274673",
            "accession": "GSE274673",
            "platform": "Xenium targeted in situ transcriptomics",
            "sample_type": "PDAC cell-resolution sections",
            "n_sections_or_samples": 6,
            "treatment_context": "3 treatment-naive; 3 chemoradiotherapy-treated",
            "tissue_sites": "primary PDAC sections",
            "role_in_study": "cell-resolution CAF-domain validation and cell-neighborhood support",
            "notes": "Targeted gene panel; supports cell-resolution spatial organization, not whole-transcriptome cell-state abundance.",
        },
        {
            "dataset": "GSE202051",
            "accession": "GSE202051",
            "platform": "single-cell/single-nucleus RNA-seq h5ad references",
            "sample_type": "PDAC reference cells/nuclei",
            "n_sections_or_samples": "reference profiles; selected marker-high cells per state",
            "treatment_context": "public reference",
            "tissue_sites": "PDAC reference atlas",
            "role_in_study": "marker-constrained reference projection and strict NNLS sensitivity",
            "notes": "Computational projection only; not spatial single-cell ground truth.",
        },
        {
            "dataset": "TCGA-PAAD",
            "accession": "TCGA-PAAD",
            "platform": "bulk RNA-seq with clinical annotations",
            "sample_type": "primary pancreatic cancer bulk tumors",
            "n_sections_or_samples": 178,
            "treatment_context": "TCGA primary tumor cohort",
            "tissue_sites": "primary tumor",
            "role_in_study": "non-spatial external biological and exploratory survival context",
            "notes": "Does not validate spatial localization or provide a clinical prediction model.",
        },
    ]
    write_csv(
        REVISION / "supplementary_tables" / "Supplementary_Table_1_Datasets.csv",
        rows,
        [
            "dataset",
            "accession",
            "platform",
            "sample_type",
            "n_sections_or_samples",
            "treatment_context",
            "tissue_sites",
            "role_in_study",
            "notes",
        ],
    )


def build_gene_module_table() -> None:
    signatures = parse_signatures()
    aliases = {
        "ifn_antigen_presentation": "IFN/MHC antigen-presentation",
        "immune_core": "immune-core",
        "spp1_tam": "SPP1/TAM",
        "tgfb_pathway": "TGF-beta",
        "emt_invasion": "EMT/invasion",
        "mycaf": "myCAF",
        "icaf": "iCAF",
        "apcaf": "apCAF",
        "pan_caf": "panCAF",
        "tumor_epithelial": "tumor epithelial",
        "pdac_basal_like": "basal-like",
        "t_cell": "T cell",
        "b_cell": "B cell",
        "dc_apc": "DC/APC",
        "plasma_cell": "plasma cell",
        "tls_chemokine": "TLS chemokine",
        "fdc_gc_like": "FDC/GC-like",
    }
    module_defs: dict[str, list[str]] = {}
    for key, genes in signatures.items():
        module_defs[aliases.get(key, key)] = genes

    composite = {
        "CAF-myeloid": signatures["mycaf"] + signatures["pan_caf"] + signatures["myeloid"] + signatures["spp1_tam"],
        "immune-core": signatures["ifn_antigen_presentation"] + signatures["t_cell"] + signatures["b_cell"] + signatures["dc_apc"] + signatures["plasma_cell"],
        "tumor-aggressive": signatures["pdac_basal_like"] + signatures["emt_invasion"] + signatures["hypoxia"],
        "immune-maturity-like": signatures["tls_chemokine"] + signatures["fdc_gc_like"] + signatures["b_cell"] + signatures["plasma_cell"] + signatures["dc_apc"],
        "myCAF/matrix": signatures["mycaf"] + signatures["pan_caf"],
        "matrix-integrin": ["COL1A1", "COL1A2", "COL3A1", "FN1", "DCN", "LUM", "ITGA5", "ITGAV", "ITGB1", "ITGA3"],
        "SPP1-CD44/integrin": ["SPP1", "CD44", "ITGAV", "ITGB1", "ITGA5", "TREM2", "APOE", "LGALS3"],
        "TGF-beta/TGFBR": ["TGFB1", "TGFBI", "TGFBR1", "TGFBR2", "SERPINE1", "CTGF", "INHBA", "SMAD3"],
        "IL6-OSM/LIF-JAKSTAT": ["IL6", "OSM", "LIF", "IL6ST", "JAK1", "JAK2", "STAT3", "SOCS3"],
    }
    for name, genes in composite.items():
        module_defs[name] = list(dict.fromkeys(genes))

    rows: list[dict[str, object]] = []
    for module_name in sorted(module_defs):
        source_type = "curated composite marker set" if module_name in composite else "curated marker set"
        used = infer_used_in_figures(module_name)
        for gene in sorted(set(module_defs[module_name])):
            rows.append(
                {
                    "module_name": module_name,
                    "gene_symbol": gene,
                    "source_reference": "To be finalized in formatted references; current set derived from project curated signatures.yaml and manuscript-defined composite axes.",
                    "source_type": source_type,
                    "used_in_figures": used,
                    "notes": "Composite modules are explicitly listed so overlap sensitivity can be audited.",
                }
            )
    write_csv(
        REVISION / "supplementary_tables" / "Supplementary_Table_2_Gene_Modules.csv",
        rows,
        ["module_name", "gene_symbol", "source_reference", "source_type", "used_in_figures", "notes"],
    )


def infer_used_in_figures(module_name: str) -> str:
    if module_name in {"CAF-myeloid", "IFN/MHC antigen-presentation", "immune-core", "tumor-aggressive", "immune-maturity-like"}:
        return "Figure 1; Figure 2; Extended Data Figure 1; Extended Data Figure 2"
    if module_name in {"myCAF/matrix", "SPP1/TAM", "TGF-beta/EMT", "matrix-integrin", "SPP1-CD44/integrin", "TGF-beta/TGFBR", "IL6-OSM/LIF-JAKSTAT"}:
        return "Figure 3; Figure 4; Extended Data Figure 3; Extended Data Figure 4; Extended Data Figure 6"
    if module_name in {"TLS chemokine", "FDC/GC-like"}:
        return "Extended Data Figure 5"
    return "Methods; supplementary analyses"


def normalize_random_rows(rows: list[dict[str, str]], source: str) -> list[dict[str, object]]:
    out: list[dict[str, object]] = []
    for r in rows:
        dataset = r.get("dataset_id", "")
        target = r.get("target") or r.get("target_label") or r.get("target_program", "")
        sample_id = r.get("sample_id") or r.get("geo_accession") or r.get("title", "")
        tissue_site = r.get("specimen_type") or r.get("tissue") or r.get("treatment", "")
        observed = r.get("observed_rho", "")
        random_median = r.get("null_median_rho") or r.get("random_median_rho") or r.get("null_median", "")
        delta = r.get("delta_vs_null_median") or r.get("delta_vs_random_median") or r.get("delta_vs_random", "")
        p = r.get("empirical_p_more_negative", "")
        support = r.get("observed_more_negative_than_null") or r.get("observed_more_negative_than_random_median", "")
        if support == "" and delta not in {"", None}:
            try:
                support = str(float(delta) < 0)
            except ValueError:
                support = ""
        out.append(
            {
                "dataset": dataset,
                "sample_id": sample_id,
                "tissue_site": tissue_site,
                "target_program": target,
                "observed_rho": observed,
                "random_median_rho": random_median,
                "delta_vs_random": delta,
                "empirical_p": p,
                "support": support,
                "source_file": source,
            }
        )
    return out


def build_random_core_table() -> None:
    rows: list[dict[str, object]] = []
    rows.extend(normalize_random_rows(read_csv(TABLES / "mvp_random_core_permutation_sample_stats.csv"), "mvp_random_core_permutation_sample_stats.csv"))
    rows.extend(normalize_random_rows(read_csv(TABLES / "gse272362_rds_random_core_permutation_sample_stats.csv"), "gse272362_rds_random_core_permutation_sample_stats.csv"))
    rows.extend(normalize_random_rows(read_csv(TABLES / "gse235315_random_core_anchor_summary.csv"), "gse235315_random_core_anchor_summary.csv"))
    rows.extend(normalize_random_rows(read_csv(TABLES / "gse274557_full_caf_core_gradients.csv"), "gse274557_full_caf_core_gradients.csv"))
    rows.extend(normalize_random_rows(read_csv(TABLES / "gse274673_xenium_fixed_anchor_gradients.csv"), "gse274673_xenium_fixed_anchor_gradients.csv"))
    write_csv(
        REVISION / "supplementary_tables" / "Supplementary_Table_3_Random_Core_Results.csv",
        rows,
        [
            "dataset",
            "sample_id",
            "tissue_site",
            "target_program",
            "observed_rho",
            "random_median_rho",
            "delta_vs_random",
            "empirical_p",
            "support",
            "source_file",
        ],
    )


def write_docs() -> None:
    today = date.today().isoformat()
    (REVISION / "docs" / "revision_log.md").write_text(
        "\n".join(
            [
                "# Revision log",
                "",
                f"Initialized: {today}",
                "",
                "## Completed",
                "",
                "- Created revision workspace.",
                "- Copied current manuscript, supplementary information and staged PDF display items.",
                "- Generated Supplementary Table 1: datasets.",
                "- Generated Supplementary Table 2: gene modules.",
                "- Generated Supplementary Table 3: random-core results.",
                "",
                "## Pending analyses",
                "",
                "- Stronger spatial null sensitivity.",
                "- CAF-only / myeloid-only / CAF-myeloid anchor comparison.",
                "- Gene-module overlap sensitivity.",
                "- LN metastasis leave-one-out.",
                "- NMF rank stability.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (REVISION / "docs" / "analysis_todo.md").write_text(
        "\n".join(
            [
                "# Analysis todo",
                "",
                "## P1 required analyses",
                "",
                "1. Spatially contiguous random-core null.",
                "2. CAF-only, myeloid-only and CAF-myeloid combined anchor comparison.",
                "3. Gene-module overlap matrix and overlap-removed gradient sensitivity.",
                "4. LN metastasis leave-one-out and individual-sample summary.",
                "5. NMF rank stability and rank-4 justification.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (REVISION / "docs" / "claim_boundary_checklist.md").write_text(
        "\n".join(
            [
                "# Claim boundary checklist",
                "",
                "## Supported after current evidence",
                "",
                "- CAF-myeloid-high regions mark reproducible spatial cores in public PDAC spatial-transcriptomic cohorts.",
                "- Primary tumors and liver metastases show immune-coupled CAF-core architecture.",
                "- A five-sample lymph-node metastasis subset suggests selective immune/IFN decoupling.",
                "- Matrix-integrin and SPP1-CD44/integrin are candidate axes for future perturbation.",
                "",
                "## Not supported",
                "",
                "- Causal CAF-core organization of immune or tumor programs.",
                "- Definitive lymph-node clinical subtype.",
                "- Mature TLS biology without histologic/FDC/GC validation.",
                "- Clinical-grade H&E prediction.",
                "- TCGA spatial validation.",
                "",
            ]
        ),
        encoding="utf-8",
    )
    (REVISION / "docs" / "figure_panel_map.md").write_text(
        "\n".join(
            [
                "# Figure panel map",
                "",
                "- Figure 1: discovery, random-core specificity, threshold sensitivity and representative maps.",
                "- Figure 2: primary/liver/LN metastatic-site remodeling and LN boundary.",
                "- Figure 3: CAF-core ecotypes, candidate axes and cell-state interpretation.",
                "- Figure 4: external Visium and Xenium validation.",
                "- Extended Data Figure 1: robustness.",
                "- Extended Data Figure 2: metastatic immune decoupling.",
                "- Extended Data Figure 3: marker/reference/NNLS/Xenium cell-state support.",
                "- Extended Data Figure 4: mechanism/interface prioritization.",
                "- Extended Data Figure 5: H&E, TCGA and TLS boundaries.",
                "- Extended Data Figure 6: deep spatial architecture.",
                "",
            ]
        ),
        encoding="utf-8",
    )


def write_summary() -> None:
    inventory_rows = []
    for path in sorted(REVISION.rglob("*")):
        if path.is_file():
            inventory_rows.append(
                {
                    "relative_path": path.relative_to(REVISION).as_posix(),
                    "size_bytes": path.stat().st_size,
                }
            )
    write_csv(REVISION / "docs" / "revision_workspace_inventory.csv", inventory_rows, ["relative_path", "size_bytes"])


def main() -> None:
    ensure_dirs()
    copy_current_package_files()
    build_dataset_table()
    build_gene_module_table()
    build_random_core_table()
    write_docs()
    write_summary()
    print(f"Revision workspace: {REVISION}")
    print(f"Files: {sum(1 for p in REVISION.rglob('*') if p.is_file())}")
    print("Generated Supplementary Tables 1-3.")


if __name__ == "__main__":
    main()
