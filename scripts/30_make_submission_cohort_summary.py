from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
REPORT_DIR = PROJECT / "results" / "reports"


def _fmt_int(value: int) -> str:
    return f"{int(value):,}"


def main() -> None:
    mvp = pd.read_csv(TABLE_DIR / "mvp_sample_level_scores.csv")
    gse272362 = pd.read_csv(TABLE_DIR / "gse272362_rds_sample_level_scores.csv")
    gse235315 = pd.read_csv(TABLE_DIR / "gse235315_sample_level_scores.csv")

    rows = []

    for dataset_id, label, role, limitation in [
        (
            "GSE282302",
            "post-neoadjuvant PDAC ST-H&E discovery cohort",
            "Discovery cohort for CAF-myeloid cores, random-core specificity, threshold sensitivity and H&E morphology bridge.",
            "Post-neoadjuvant/FOLFIRINOX context; patient-level clinical outcome and patient heterogeneity are not audited.",
        ),
        (
            "GSE274103",
            "treatment-naive PDAC ST-H&E support cohort",
            "Treatment-naive support for CAF-core spatial organization and H&E morphology bridge.",
            "Small five-section support cohort; used for directional support rather than patient-level inference.",
        ),
    ]:
        sub = mvp.loc[mvp["dataset_id"].eq(dataset_id)].copy()
        rows.append(
            {
                "cohort": dataset_id,
                "context": label,
                "specimen_group": "PDAC sections",
                "n_samples": len(sub),
                "n_spots": int(sub["n_spots_qc"].sum()),
                "analysis_role": role,
                "main_figures": "Figure 1; Extended Data Figures 1 and 4",
                "key_limitation": limitation,
            }
        )

    site_labels = {
        "primary_tumor": "primary tumors",
        "liver_metastasis": "liver metastases",
        "lymph_node_metastasis": "lymph-node metastases",
        "normal_pancreas": "normal pancreas",
    }
    site_roles = {
        "primary_tumor": "Independent validation of CAF-core-centered inflammatory and tumor-aggressive organization.",
        "liver_metastasis": "Independent metastatic validation of CAF-core-centered organization.",
        "lymph_node_metastasis": "Biological contrast showing preserved tumor-aggressive CAF coupling but immune/IFN decoupling.",
        "normal_pancreas": "Contextual normal-pancreas reference; not used as the main tumor-organization claim.",
    }
    site_limits = {
        "primary_tumor": "RDS-derived analysis; edge/background image QC is not identical to directly processed cohorts.",
        "liver_metastasis": "RDS-derived analysis; edge/background image QC is not identical to directly processed cohorts.",
        "lymph_node_metastasis": "Five specimens; best framed as a strong spatial biology lead rather than a definitive clinical subtype.",
        "normal_pancreas": "Small contextual group; not used for clinical or tumor-subtype conclusions.",
    }
    for specimen_type, sub in gse272362.groupby("specimen_type", sort=True):
        rows.append(
            {
                "cohort": "GSE272362",
                "context": "independent primary/metastatic PDAC atlas",
                "specimen_group": site_labels.get(specimen_type, specimen_type),
                "n_samples": len(sub),
                "n_spots": int(sub["n_spots_qc"].sum()),
                "analysis_role": site_roles.get(specimen_type, "Independent spatial validation context."),
                "main_figures": "Figure 2; Extended Data Figures 2, 3 and 5",
                "key_limitation": site_limits.get(specimen_type, "RDS-derived spatial analysis."),
            }
        )

    rows.append(
        {
            "cohort": "GSE235315",
            "context": "external paired-ST anchor cohort",
            "specimen_group": "PDAC paired ST-H&E samples",
            "n_samples": len(gse235315),
            "n_spots": int(gse235315["n_spots_qc"].sum()),
            "analysis_role": "External paired-ST anchor for CAF-core random-core support.",
            "main_figures": "Figure 1; Extended Data Figure 5",
            "key_limitation": "Sample metadata are not fully curated; use as spatial-state support, not patient-level or treatment-context evidence.",
        }
    )

    out = pd.DataFrame(rows)
    out_path = TABLE_DIR / "submission_cohort_summary.csv"
    out.to_csv(out_path, index=False)

    md_lines = [
        "# Submission Cohort Summary",
        "",
        "Last updated: 2026-06-25",
        "",
        "This table is intended for the manuscript, supplement, cover letter and reviewer-facing data provenance summary.",
        "",
        "| cohort | context | specimen group | n samples | n spots | analysis role | main figures | key limitation |",
        "|---|---|---|---:|---:|---|---|---|",
    ]
    for row in out.to_dict("records"):
        md_lines.append(
            "| {cohort} | {context} | {specimen_group} | {n_samples} | {n_spots} | {analysis_role} | {main_figures} | {key_limitation} |".format(
                cohort=row["cohort"],
                context=row["context"],
                specimen_group=row["specimen_group"],
                n_samples=_fmt_int(row["n_samples"]),
                n_spots=_fmt_int(row["n_spots"]),
                analysis_role=row["analysis_role"],
                main_figures=row["main_figures"],
                key_limitation=row["key_limitation"],
            )
        )

    md_lines.extend(
        [
            "",
            "## Manuscript Use",
            "",
            "- Use this as a compact Supplementary Table for cohort provenance.",
            "- Keep the limitation column visible; it prevents overclaiming patient-level, clinical-outcome or causal conclusions.",
            "- The active quantitative source tables remain under `results/source_data` and are mapped in `results/reports/submission_figure_captions_and_source_map.md`.",
        ]
    )
    report_path = REPORT_DIR / "submission_cohort_summary.md"
    report_path.write_text("\n".join(md_lines) + "\n", encoding="utf-8")

    print(f"Wrote {out_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
