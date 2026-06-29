from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_27 = PROJECT_ROOT / "scripts" / "27_targeted_gene_axis_validation.py"
OUT_SAMPLE = PROJECT_ROOT / "results/tables/gap1_cxcl9_spp1_polarity_per_sample.csv"
OUT_SUMMARY = PROJECT_ROOT / "results/tables/gap1_cxcl9_spp1_polarity_context_summary.csv"
OUT_CORR = PROJECT_ROOT / "results/tables/gap1_cxcl9_spp1_polarity_decoupling_correlations.csv"
OUT_FIG = PROJECT_ROOT / "results/figures/submission/extended_data_gap1_cxcl9_spp1_polarity"
OUT_REPORT = PROJECT_ROOT / "results/reports/gap1_cxcl9_spp1_polarity_report.md"
STATUS = PROJECT_ROOT / "results/logs/stage_32_gap1_cxcl9_spp1_polarity_status.json"

GENES = [
    "SPP1",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "APOE",
    "TREM2",
    "LGALS3",
    "CD68",
    "C1QA",
    "C1QB",
    "C1QC",
    "LST1",
    "FCGR3A",
    "MARCO",
    "MRC1",
    "STAT1",
    "IRF1",
]

SPP1_TAM_GENES = ["SPP1", "APOE", "TREM2", "LGALS3", "CD68", "C1QA", "C1QB", "C1QC", "LST1", "FCGR3A"]
CXCL9_IFN_GENES = ["CXCL9", "CXCL10", "CXCL11", "STAT1", "IRF1"]

FEATURE_ORDER = [
    "SPP1",
    "CXCL9",
    "CXCL10",
    "CXCL11",
    "APOE",
    "TREM2",
    "CD68",
    "SPP1 TAM program",
    "CXCL9 IFN program",
    "SPP1-high/CXCL9-low polarity",
]

CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
    "external_paired_st_anchor",
]

CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "external_paired_st_anchor": "GSE235315",
    "normal_pancreas": "normal",
}

COLORS = {
    "post_neoadjuvant_sections": "#8C6D31",
    "treatment_naive_primary": "#B55A30",
    "primary_tumor": "#4E79A7",
    "liver_metastasis": "#59A14F",
    "lymph_node_metastasis": "#9C755F",
    "external_paired_st_anchor": "#7B6BB1",
}


def load_stage27_module():
    spec = importlib.util.spec_from_file_location("stage27_targeted", SCRIPT_27)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {SCRIPT_27}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def zscore(values: pd.Series) -> pd.Series:
    arr = values.to_numpy(float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(arr)), index=values.index)
    return pd.Series((arr - np.nanmean(arr)) / sd, index=values.index)


def load_expression() -> pd.DataFrame:
    stage27 = load_stage27_module()
    stage27.TARGET_GENES = sorted(set(GENES))
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    mvp = stage27.load_mvp_spots()
    anchor = stage27.load_gse235315_spots()
    combined_spots = pd.concat([mvp, anchor], ignore_index=True)
    h5_expr = stage27.extract_h5_gene_values(manifest, combined_spots)
    rds_expr = stage27.extract_gse272362_gene_values()
    frames = [frame for frame in [h5_expr, rds_expr] if not frame.empty]
    if not frames:
        raise RuntimeError("No expression frames available for Gap 1 polarity analysis.")
    expr = pd.concat(frames, ignore_index=True)
    if "cohort_context" not in expr.columns:
        expr["cohort_context"] = expr.apply(stage27.context_for_row, axis=1)
    else:
        missing = expr["cohort_context"].isna() | expr["cohort_context"].astype(str).eq("")
        if missing.any():
            expr.loc[missing, "cohort_context"] = expr.loc[missing].apply(stage27.context_for_row, axis=1)
    expr = stage27.add_regions(expr)
    return expr


def add_polarity_features(expr: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    present_genes = [gene for gene in GENES if gene in out.columns and out[gene].notna().any()]
    for gene in present_genes:
        out[f"z_{gene}"] = out.groupby("sample_id", group_keys=False)[gene].apply(zscore)

    spp1_present = [f"z_{gene}" for gene in SPP1_TAM_GENES if f"z_{gene}" in out.columns]
    cxcl9_present = [f"z_{gene}" for gene in CXCL9_IFN_GENES if f"z_{gene}" in out.columns]
    out["SPP1 TAM program"] = out[spp1_present].mean(axis=1) if spp1_present else np.nan
    out["CXCL9 IFN program"] = out[cxcl9_present].mean(axis=1) if cxcl9_present else np.nan
    out["SPP1-high/CXCL9-low polarity"] = out["SPP1 TAM program"] - out["CXCL9 IFN program"]
    for gene in present_genes:
        out[gene] = out[f"z_{gene}"]
    return out


def summarize_samples(expr: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    features = [feature for feature in FEATURE_ORDER if feature in expr.columns]
    for sample_id, sample in expr.groupby("sample_id", sort=False):
        if "is_caf_core_top10" not in sample or sample["is_caf_core_top10"].sum() < 3:
            continue
        core = sample["is_caf_core_top10"].astype(bool)
        context = sample["cohort_context"].iloc[0]
        dataset_id = sample["dataset_id"].iloc[0] if "dataset_id" in sample else ""
        specimen_type = sample["specimen_type"].iloc[0] if "specimen_type" in sample else ""
        row = {
            "dataset_id": dataset_id,
            "sample_id": sample_id,
            "specimen_type": specimen_type,
            "cohort_context": context,
            "n_spots": int(len(sample)),
            "n_caf_core_spots": int(core.sum()),
            "n_interface_spots": int(sample["is_interface"].sum()) if "is_interface" in sample else 0,
        }
        for feature in features:
            values = sample[feature]
            row[f"core_mean__{feature}"] = float(values[core].mean())
            row[f"noncore_mean__{feature}"] = float(values[~core].mean())
            row[f"core_enrichment__{feature}"] = float(values[core].mean() - values[~core].mean())
            if "is_interface" in sample and sample["is_interface"].sum() >= 3:
                interface = sample["is_interface"].astype(bool)
                row[f"interface_enrichment__{feature}"] = float(values[interface].mean() - values[~interface].mean())
            else:
                row[f"interface_enrichment__{feature}"] = np.nan
        rows.append(row)

    per_sample = pd.DataFrame(rows)
    decoupling = pd.read_csv(PROJECT_ROOT / "results/tables/mechanism_candidate_axis_sample_summary.csv")
    keep = [
        "sample_id",
        "stromal_tumor_core_coupling",
        "immune_core_coupling",
        "immune_decoupling_index",
    ]
    per_sample = per_sample.merge(decoupling[keep], on="sample_id", how="left")
    return per_sample


def summarize_context(per_sample: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for context, sub in per_sample.groupby("cohort_context", sort=False):
        for feature in FEATURE_ORDER:
            col = f"core_enrichment__{feature}"
            if col not in sub:
                continue
            values = sub[col].dropna()
            if values.empty:
                continue
            rows.append(
                {
                    "cohort_context": context,
                    "feature": feature,
                    "n_samples": int(values.shape[0]),
                    "median_core_enrichment": float(values.median()),
                    "n_core_positive": int((values > 0).sum()),
                    "mean_core_enrichment": float(values.mean()),
                }
            )
    return pd.DataFrame(rows)


def compute_correlations(per_sample: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    feature_cols = [
        "core_enrichment__SPP1-high/CXCL9-low polarity",
        "core_mean__SPP1-high/CXCL9-low polarity",
        "core_enrichment__SPP1",
        "core_enrichment__CXCL9",
        "core_enrichment__SPP1 TAM program",
        "core_enrichment__CXCL9 IFN program",
    ]
    groups = {
        "all_non_normal": per_sample[~per_sample["cohort_context"].eq("normal_pancreas")],
        "tumor_contexts_no_external": per_sample[
            per_sample["cohort_context"].isin(
                [
                    "post_neoadjuvant_sections",
                    "treatment_naive_primary",
                    "primary_tumor",
                    "liver_metastasis",
                    "lymph_node_metastasis",
                ]
            )
        ],
        "gse272362_tumor_sites": per_sample[
            per_sample["cohort_context"].isin(["primary_tumor", "liver_metastasis", "lymph_node_metastasis"])
        ],
    }
    for group_name, sub in groups.items():
        for col in feature_cols:
            if col not in sub:
                continue
            tmp = sub[[col, "immune_decoupling_index"]].dropna()
            if tmp.shape[0] < 4:
                continue
            rho, p = spearmanr(tmp[col], tmp["immune_decoupling_index"])
            rows.append(
                {
                    "analysis_group": group_name,
                    "feature": col.replace("core_enrichment__", "").replace("core_mean__", "core_mean "),
                    "n_samples": int(tmp.shape[0]),
                    "spearman_rho_vs_immune_decoupling": float(rho),
                    "p_value": float(p),
                }
            )
    return pd.DataFrame(rows)


def plot_gap1(summary: pd.DataFrame, per_sample: pd.DataFrame, corr: pd.DataFrame) -> None:
    contexts = [c for c in CONTEXT_ORDER if c in set(summary["cohort_context"])]
    features = [f for f in FEATURE_ORDER if f in set(summary["feature"])]
    mat = (
        summary.pivot_table(index="feature", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
        .reindex(index=features, columns=contexts)
    )

    fig = plt.figure(figsize=(12.5, 6.6), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1.0, 1.12], height_ratios=[1.0, 0.08])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])

    fig.suptitle("Gap 1: CXCL9:SPP1 TAM polarity around CAF-myeloid cores", fontsize=15, fontweight="bold")

    im = ax_a.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-1.1, vmax=1.1, aspect="auto")
    ax_a.set_title("A  CAF-core gene/program enrichment", loc="left", fontsize=10.5, fontweight="bold")
    ax_a.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax_a.set_yticks(np.arange(len(features)), features)
    ax_a.tick_params(axis="both", labelsize=8)
    for i, feature in enumerate(features):
        for j, context in enumerate(contexts):
            val = mat.loc[feature, context]
            if np.isfinite(val):
                row = summary[summary["feature"].eq(feature) & summary["cohort_context"].eq(context)].iloc[0]
                ax_a.text(j, i, f"{val:.2f}\n{int(row.n_core_positive)}/{int(row.n_samples)}", ha="center", va="center", fontsize=5.7)
    cb = fig.colorbar(im, ax=ax_a, fraction=0.045, pad=0.02)
    cb.set_label("median core enrichment", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    ycol = "core_enrichment__SPP1-high/CXCL9-low polarity"
    plot_df = per_sample[per_sample["cohort_context"].isin(contexts)].copy()
    positions = np.arange(len(contexts))
    medians = [plot_df.loc[plot_df["cohort_context"].eq(c), ycol].median() for c in contexts]
    ax_b.bar(positions, medians, color=[COLORS.get(c, "#777777") for c in contexts], alpha=0.82, width=0.62)
    rng = np.random.default_rng(13)
    for x, context in zip(positions, contexts):
        vals = plot_df.loc[plot_df["cohort_context"].eq(context), ycol].dropna().to_numpy(float)
        jitter = rng.uniform(-0.13, 0.13, size=len(vals))
        ax_b.scatter(np.full(len(vals), x) + jitter, vals, s=14, color="black", alpha=0.38, linewidths=0)
    ax_b.axhline(0, color="#333333", lw=0.8)
    ax_b.set_title("B  SPP1-high / CXCL9-low polarity", loc="left", fontsize=10.5, fontweight="bold")
    ax_b.set_xticks(positions, [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax_b.set_ylabel("CAF-core enrichment", fontsize=8.5)
    ax_b.tick_params(axis="both", labelsize=8)
    ax_b.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax_b.spines[["top", "right"]].set_visible(False)

    scatter_df = per_sample[
        per_sample["cohort_context"].isin(
            ["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
        )
    ][[ycol, "immune_decoupling_index", "cohort_context"]].dropna()
    for context, sub in scatter_df.groupby("cohort_context"):
        ax_c.scatter(
            sub[ycol],
            sub["immune_decoupling_index"],
            s=28 if len(sub) < 20 else 18,
            color=COLORS.get(context, "#777777"),
            alpha=0.72,
            label=CONTEXT_LABELS.get(context, context),
            edgecolor="white",
            linewidth=0.3,
        )
    if len(scatter_df) >= 4:
        x = scatter_df[ycol].to_numpy(float)
        y = scatter_df["immune_decoupling_index"].to_numpy(float)
        coeff = np.polyfit(x, y, deg=1)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ax_c.plot(xx, coeff[0] * xx + coeff[1], color="#222222", lw=1.1)
        rho, p = spearmanr(x, y)
        ax_c.text(
            0.03,
            0.97,
            f"Spearman rho = {rho:.2f}\np = {p:.2g}",
            transform=ax_c.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#DDDDDD", lw=0.6),
        )
    ax_c.set_title("C  Polarity versus immune decoupling", loc="left", fontsize=10.5, fontweight="bold")
    ax_c.set_xlabel("CAF-core polarity enrichment", fontsize=8.5)
    ax_c.set_ylabel("immune-decoupling index", fontsize=8.5)
    ax_c.tick_params(axis="both", labelsize=8)
    ax_c.grid(color="#E5E5E5", lw=0.7)
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.legend(frameon=False, fontsize=6.7, loc="lower right")

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = OUT_FIG.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(summary: pd.DataFrame, per_sample: pd.DataFrame, corr: pd.DataFrame) -> None:
    polarity = "SPP1-high/CXCL9-low polarity"
    pol_col = f"core_enrichment__{polarity}"
    context_rows = []
    for context in CONTEXT_ORDER:
        sub = per_sample[per_sample["cohort_context"].eq(context)]
        if sub.empty or pol_col not in sub:
            continue
        context_rows.append(
            f"- {CONTEXT_LABELS.get(context, context)}: median CAF-core polarity enrichment {sub[pol_col].median():.3f}; positive in {(sub[pol_col] > 0).sum()}/{sub[pol_col].notna().sum()} samples."
        )

    main_corr = corr[
        corr["analysis_group"].eq("tumor_contexts_no_external")
        & corr["feature"].eq("SPP1-high/CXCL9-low polarity")
    ]
    corr_text = "not available"
    if not main_corr.empty:
        row = main_corr.iloc[0]
        corr_text = (
            f"Spearman rho {row['spearman_rho_vs_immune_decoupling']:.3f}, "
            f"p = {row['p_value']:.3g}, n = {int(row['n_samples'])}"
        )

    report = [
        "# Gap 1 CXCL9:SPP1 TAM Polarity Report",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Gap 1 Definition",
        "",
        "The current manuscript nominates SPP1-TAM/matrix as a CAF-core candidate axis, but had not tested whether CAF-core immune decoupling aligns with a CXCL9:SPP1 macrophage-polarity framework.",
        "",
        "## Result Summary",
        "",
        *context_rows,
        "",
        f"- Across tumor contexts excluding the external anchor, polarity enrichment versus immune-decoupling index: {corr_text}.",
        "",
        "## Generated Outputs",
        "",
        f"- `{OUT_SAMPLE.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_SUMMARY.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_CORR.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.pdf').relative_to(PROJECT_ROOT)}`",
        "",
        "## Claim Boundary",
        "",
        "This analysis supports or refutes alignment between CAF-core immune decoupling and an SPP1-high/CXCL9-low TAM polarity axis. It remains targeted gene-level spatial evidence, not causal macrophage-state perturbation or ligand-receptor proof.",
    ]
    OUT_REPORT.write_text("\n".join(report) + "\n", encoding="utf-8")


def main() -> int:
    expr = load_expression()
    expr = add_polarity_features(expr)
    per_sample = summarize_samples(expr)
    summary = summarize_context(per_sample)
    corr = compute_correlations(per_sample)

    OUT_SAMPLE.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(OUT_SAMPLE, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    corr.to_csv(OUT_CORR, index=False)
    plot_gap1(summary, per_sample, corr)
    write_report(summary, per_sample, corr)

    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "32_gap1_cxcl9_spp1_polarity",
                "status": "success",
                "n_samples": int(per_sample["sample_id"].nunique()),
                "n_contexts": int(per_sample["cohort_context"].nunique()),
                "outputs": [
                    str(OUT_SAMPLE),
                    str(OUT_SUMMARY),
                    str(OUT_CORR),
                    str(OUT_FIG.with_suffix(".pdf")),
                    str(OUT_REPORT),
                ],
                "claim_boundary": "Targeted gene-level CXCL9:SPP1 polarity analysis; not causal macrophage perturbation.",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote Gap 1 CXCL9:SPP1 polarity analysis for {per_sample['sample_id'].nunique()} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
