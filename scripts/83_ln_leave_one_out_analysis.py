from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "results" / "tables"
REVISION = PROJECT / "results" / "revision_2026_06_29"
ANALYSIS_OUT = REVISION / "analysis_outputs"
SUPP_TABLES = REVISION / "supplementary_tables"
DOCS = REVISION / "docs"


TARGET_MAP = {
    "ifn_mhc": "IFN_MHC",
    "immune_core": "immune_core",
    "tumor_aggressive": "tumor_aggressive",
}


def load_ln_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    random_stats = pd.read_csv(TABLES / "gse272362_rds_random_core_permutation_sample_stats.csv", encoding="utf-8-sig")
    random_stats = random_stats[random_stats["specimen_type"].eq("lymph_node_metastasis")].copy()
    random_stats = random_stats[random_stats["target"].isin(TARGET_MAP)].copy()
    random_stats["target_clean"] = random_stats["target"].map(TARGET_MAP)

    pivot = random_stats.pivot_table(
        index=["dataset_id", "sample_id", "specimen_type"],
        columns="target_clean",
        values=["observed_rho", "null_median_rho", "delta_vs_null_median", "empirical_p_more_negative"],
        aggfunc="first",
    )
    pivot.columns = [f"{metric}__{target}" for metric, target in pivot.columns]
    pivot = pivot.reset_index()

    ecotype = pd.read_csv(TABLES / "spatial_ecotype_sample_summary.csv", encoding="utf-8-sig")
    ecotype_ln = ecotype[
        ecotype["dataset_id"].eq("GSE272362") & ecotype["cohort_context"].eq("lymph_node_metastasis")
    ][
        [
            "sample_id",
            "stromal_tumor_core_coupling",
            "immune_core_coupling",
            "immune_decoupling_index",
            "core_enrichment__IFN/MHC",
            "core_enrichment__immune core",
            "core_enrichment__tumor aggressive",
            "dominant_nmf_ecotype",
            "dominant_nmf_ecotype_label",
        ]
    ].copy()
    individual = pivot.merge(ecotype_ln, on="sample_id", how="left")
    return random_stats, individual


def build_leave_one_out(individual: pd.DataFrame) -> pd.DataFrame:
    rows = []
    samples = list(individual["sample_id"])
    metrics = {
        "IFN_MHC_median_delta": "delta_vs_null_median__IFN_MHC",
        "immune_core_median_delta": "delta_vs_null_median__immune_core",
        "tumor_aggressive_median_delta": "delta_vs_null_median__tumor_aggressive",
        "immune_decoupling_index": "immune_decoupling_index",
        "stromal_tumor_core_coupling": "stromal_tumor_core_coupling",
        "immune_core_coupling": "immune_core_coupling",
    }

    def summarize(sub: pd.DataFrame, left_out: str) -> dict[str, object]:
        row: dict[str, object] = {
            "analysis": "leave_one_out" if left_out != "none" else "all_LN_samples",
            "left_out_sample": left_out,
            "n_samples_included": len(sub),
        }
        for out_name, col in metrics.items():
            row[out_name] = float(np.nanmedian(sub[col])) if col in sub.columns else np.nan
        row["IFN_MHC_support_n"] = int(np.nansum(sub["delta_vs_null_median__IFN_MHC"] < 0))
        row["immune_core_support_n"] = int(np.nansum(sub["delta_vs_null_median__immune_core"] < 0))
        row["tumor_aggressive_support_n"] = int(np.nansum(sub["delta_vs_null_median__tumor_aggressive"] < 0))
        return row

    rows.append(summarize(individual, "none"))
    for sample in samples:
        sub = individual[individual["sample_id"].ne(sample)].copy()
        rows.append(summarize(sub, sample))
    return pd.DataFrame(rows)


def write_report(individual: pd.DataFrame, loo: pd.DataFrame) -> None:
    all_row = loo[loo["left_out_sample"].eq("none")].iloc[0]
    loo_rows = loo[~loo["left_out_sample"].eq("none")].copy()
    lines = [
        "# LN metastasis leave-one-out analysis",
        "",
        "This analysis tests whether the five-sample lymph-node metastasis immune-decoupling pattern is driven by a single specimen.",
        "",
        "## All LN samples",
        "",
        f"- IFN/MHC median delta: {all_row['IFN_MHC_median_delta']:.3f}; support {int(all_row['IFN_MHC_support_n'])}/5.",
        f"- Immune-core median delta: {all_row['immune_core_median_delta']:.3f}; support {int(all_row['immune_core_support_n'])}/5.",
        f"- Tumor-aggressive median delta: {all_row['tumor_aggressive_median_delta']:.3f}; support {int(all_row['tumor_aggressive_support_n'])}/5.",
        f"- Median immune-decoupling index: {all_row['immune_decoupling_index']:.3f}.",
        "",
        "## Leave-one-out range",
        "",
        f"- IFN/MHC median delta range: {loo_rows['IFN_MHC_median_delta'].min():.3f} to {loo_rows['IFN_MHC_median_delta'].max():.3f}.",
        f"- Immune-core median delta range: {loo_rows['immune_core_median_delta'].min():.3f} to {loo_rows['immune_core_median_delta'].max():.3f}.",
        f"- Tumor-aggressive median delta range: {loo_rows['tumor_aggressive_median_delta'].min():.3f} to {loo_rows['tumor_aggressive_median_delta'].max():.3f}.",
        f"- Immune-decoupling index range: {loo_rows['immune_decoupling_index'].min():.3f} to {loo_rows['immune_decoupling_index'].max():.3f}.",
        "",
        "## Interpretation",
        "",
        "The LN subset consistently retains tumor-aggressive CAF-core coupling while IFN/MHC and immune-core deltas remain weak or unstable. Because this subset contains only five samples, the manuscript should describe the finding as a hypothesis-generating lymph-node metastasis lead rather than a definitive clinical subtype.",
        "",
        "## Individual LN samples",
        "",
    ]
    for _, row in individual.iterrows():
        lines.append(
            f"- {row['sample_id']}: IFN/MHC delta {row['delta_vs_null_median__IFN_MHC']:.3f}; "
            f"immune-core delta {row['delta_vs_null_median__immune_core']:.3f}; "
            f"tumor-aggressive delta {row['delta_vs_null_median__tumor_aggressive']:.3f}; "
            f"decoupling index {row['immune_decoupling_index']:.3f}."
        )
    (DOCS / "ln_metastasis_leave_one_out_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    ANALYSIS_OUT.mkdir(parents=True, exist_ok=True)
    SUPP_TABLES.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    _, individual = load_ln_tables()
    loo = build_leave_one_out(individual)
    individual.to_csv(ANALYSIS_OUT / "ln_metastasis_individual_sample_summary.csv", index=False, encoding="utf-8")
    loo.to_csv(ANALYSIS_OUT / "ln_metastasis_leave_one_out_summary.csv", index=False, encoding="utf-8")
    loo.to_csv(SUPP_TABLES / "Supplementary_Table_6_LN_Leave_One_Out.csv", index=False, encoding="utf-8")
    write_report(individual, loo)
    print("Wrote LN leave-one-out outputs.")


if __name__ == "__main__":
    main()
