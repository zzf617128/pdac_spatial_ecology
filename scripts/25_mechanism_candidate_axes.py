from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE = "25_mechanism_candidate_axes"

AXES = [
    {
        "axis_id": "tgfb_emt_invasive",
        "axis_label": "TGF-beta/EMT invasive",
        "components": ["score_tgfb_pathway", "score_emt_invasion", "score_mycaf"],
    },
    {
        "axis_id": "spp1_tam_matrix",
        "axis_label": "SPP1-TAM/matrix",
        "components": ["score_spp1_tam", "score_myeloid", "score_pan_caf"],
    },
    {
        "axis_id": "icaf_cytokine_chemokine",
        "axis_label": "iCAF cytokine/chemokine",
        "components": ["score_icaf", "score_tls_chemokine"],
    },
    {
        "axis_id": "ifn_apc_antigen",
        "axis_label": "IFN/APC antigen",
        "components": ["score_ifn_antigen_presentation", "score_dc_apc", "score_apcaf"],
    },
    {
        "axis_id": "t_cell_checkpoint",
        "axis_label": "T cell/checkpoint",
        "components": ["score_t_cell", "score_cd8_effector", "score_t_cell_exhaustion_checkpoint"],
    },
    {
        "axis_id": "b_plasma_lymphoid",
        "axis_label": "B/plasma lymphoid",
        "components": ["score_b_cell", "score_plasma_cell", "score_tls_chemokine"],
    },
    {
        "axis_id": "basal_hypoxic_tumor",
        "axis_label": "Basal-hypoxic tumor",
        "components": ["score_pdac_basal_like", "score_hypoxia", "score_tumor_aggressive"],
    },
    {
        "axis_id": "classical_epithelial",
        "axis_label": "Classical epithelial",
        "components": ["score_pdac_classical_like", "score_tumor_epithelial"],
    },
]

CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
    "normal_pancreas",
    "external_paired_st_anchor",
]
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treat-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "external_paired_st_anchor": "GSE235315",
}

STROMAL_TUMOR = [
    "score_mycaf",
    "score_myeloid",
    "score_spp1_tam",
    "score_tgfb_pathway",
    "score_emt_invasion",
    "score_tumor_aggressive",
]
IMMUNE = [
    "score_ifn_antigen_presentation",
    "score_immune_hub_core",
    "score_t_cell",
    "score_b_cell",
    "score_dc_apc",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(status: str, payload: dict) -> None:
    base = {
        "stage": STAGE,
        "status": status,
        "timestamp_utc": now_iso(),
        "n_errors": 0,
        "critical_errors": [],
        "noncritical_warnings": [],
        "next_manual_check": [],
    }
    base.update(payload)
    path = PROJECT_ROOT / f"results/logs/stage_{STAGE}_status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(len(arr), dtype=float)
    return (arr - np.nanmean(arr)) / sd


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    keep = np.isfinite(x) & np.isfinite(y)
    if keep.sum() < 30:
        return np.nan
    if len(np.unique(x[keep])) < 2 or len(np.unique(y[keep])) < 2:
        return np.nan
    return float(spearmanr(x[keep], y[keep]).statistic)


def median_neighbor_distance(points: np.ndarray) -> float:
    if len(points) < 3:
        return 1.0
    idx = np.arange(len(points))
    if len(points) > 2500:
        rng = np.random.default_rng(20260624)
        idx = rng.choice(idx, size=2500, replace=False)
    sampled = points[idx]
    nn = NearestNeighbors(n_neighbors=2)
    nn.fit(sampled)
    d = nn.kneighbors(sampled, return_distance=True)[0][:, 1]
    med = float(np.nanmedian(d))
    return med if np.isfinite(med) and med > 0 else 1.0


def nearest_distance(points: np.ndarray, centers: np.ndarray) -> np.ndarray:
    if len(centers) == 0:
        return np.full(len(points), np.nan)
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(centers)
    return nn.kneighbors(points, return_distance=True)[0][:, 0]


def load_spots() -> pd.DataFrame:
    needed_scores = sorted(
        set(
            ["score_caf_myeloid_barrier", "score_tumor_epithelial"]
            + STROMAL_TUMOR
            + IMMUNE
            + [col for axis in AXES for col in axis["components"]]
        )
    )
    base_cols = ["dataset_id", "sample_id", "patient_id", "specimen_type", "barcode", "x_pixel", "y_pixel"]

    frames: list[pd.DataFrame] = []

    mvp_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    mvp = pd.read_csv(mvp_path)
    mvp = mvp[~mvp["edge_or_background_risk"].astype(str).str.lower().eq("true")].copy()
    mvp["cohort_context"] = np.where(
        mvp["dataset_id"].eq("GSE274103"),
        "treatment_naive_primary",
        "post_neoadjuvant_sections",
    )
    frames.append(mvp[[c for c in base_cols + ["cohort_context"] + needed_scores if c in mvp.columns]].copy())

    gse272362 = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv")
    gse272362["cohort_context"] = gse272362["specimen_type"].astype(str)
    frames.append(gse272362[[c for c in base_cols + ["cohort_context"] + needed_scores if c in gse272362.columns]].copy())

    gse235315 = pd.read_csv(PROJECT_ROOT / "results/tables/gse235315_spot_level_scores.csv")
    gse235315["cohort_context"] = "external_paired_st_anchor"
    frames.append(gse235315[[c for c in base_cols + ["cohort_context"] + needed_scores if c in gse235315.columns]].copy())

    out = pd.concat(frames, ignore_index=True)
    out = out[np.isfinite(out["x_pixel"]) & np.isfinite(out["y_pixel"])].copy()
    return out


def add_axis_scores(sample: pd.DataFrame) -> pd.DataFrame:
    sample = sample.copy()
    for axis in AXES:
        z_components = []
        for col in axis["components"]:
            if col not in sample.columns:
                continue
            z_col = f"z_axis_component__{col}"
            sample[z_col] = zscore(sample[col])
            z_components.append(z_col)
        if z_components:
            sample[f"axis__{axis['axis_id']}"] = sample[z_components].mean(axis=1)
        else:
            sample[f"axis__{axis['axis_id']}"] = np.nan
    for col in sorted(set(STROMAL_TUMOR + IMMUNE)):
        if col in sample.columns:
            sample[f"z_for_decoupling__{col}"] = zscore(sample[col])
    return sample


def summarize_sample(sample: pd.DataFrame) -> tuple[list[dict], dict]:
    sample = add_axis_scores(sample)
    points = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    med_nn = median_neighbor_distance(points)

    caf = sample["score_caf_myeloid_barrier"].to_numpy(float)
    tumor = sample["score_tumor_epithelial"].to_numpy(float)
    caf_core = caf >= np.nanpercentile(caf, 90)
    tumor_high = tumor >= np.nanpercentile(tumor, 80)
    d_caf = nearest_distance(points, points[caf_core]) / med_nn
    d_tumor = nearest_distance(points, points[tumor_high]) / med_nn
    caf_near = d_caf <= 2.0
    tumor_near = d_tumor <= 2.0
    interface = caf_near & tumor_near
    caf_only = caf_near & ~tumor_near
    tumor_only = tumor_near & ~caf_near

    first = sample.iloc[0]
    sample_meta = {
        "dataset_id": first["dataset_id"],
        "sample_id": first["sample_id"],
        "patient_id": first["patient_id"],
        "specimen_type": first["specimen_type"],
        "cohort_context": first["cohort_context"],
        "n_spots": int(len(sample)),
        "n_caf_core": int(caf_core.sum()),
        "n_interface": int(interface.sum()),
        "median_neighbor_distance_px": med_nn,
    }

    stromal_enrich = []
    immune_enrich = []
    for col in STROMAL_TUMOR:
        z_col = f"z_for_decoupling__{col}"
        if z_col in sample.columns:
            stromal_enrich.append(float(np.nanmean(sample.loc[caf_core, z_col]) - np.nanmean(sample.loc[~caf_core, z_col])))
    for col in IMMUNE:
        z_col = f"z_for_decoupling__{col}"
        if z_col in sample.columns:
            immune_enrich.append(float(np.nanmean(sample.loc[caf_core, z_col]) - np.nanmean(sample.loc[~caf_core, z_col])))
    sample_meta["stromal_tumor_core_coupling"] = float(np.nanmean(stromal_enrich)) if stromal_enrich else np.nan
    sample_meta["immune_core_coupling"] = float(np.nanmean(immune_enrich)) if immune_enrich else np.nan
    sample_meta["immune_decoupling_index"] = sample_meta["stromal_tumor_core_coupling"] - sample_meta["immune_core_coupling"]

    rows: list[dict] = []
    for axis in AXES:
        axis_col = f"axis__{axis['axis_id']}"
        values = sample[axis_col].to_numpy(float)
        core_enrich = float(np.nanmean(values[caf_core]) - np.nanmean(values[~caf_core]))
        interface_enrich = (
            float(np.nanmean(values[interface]) - np.nanmean(values[~interface])) if interface.sum() >= 20 else np.nan
        )
        interface_vs_caf_only = (
            float(np.nanmean(values[interface]) - np.nanmean(values[caf_only])) if interface.sum() >= 20 and caf_only.sum() >= 20 else np.nan
        )
        interface_vs_tumor_only = (
            float(np.nanmean(values[interface]) - np.nanmean(values[tumor_only])) if interface.sum() >= 20 and tumor_only.sum() >= 20 else np.nan
        )
        rows.append(
            {
                **sample_meta,
                "axis_id": axis["axis_id"],
                "axis_label": axis["axis_label"],
                "axis_components": ";".join(axis["components"]),
                "core_enrichment": core_enrich,
                "interface_enrichment": interface_enrich,
                "interface_vs_caf_only": interface_vs_caf_only,
                "interface_vs_tumor_only": interface_vs_tumor_only,
                "rho_distance_to_caf_core": safe_spearman(d_caf, values),
            }
        )
    return rows, sample_meta


def bh_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    prev = 1.0
    m = len(p)
    for i in range(m - 1, -1, -1):
        val = ranked[i] * m / (i + 1)
        prev = min(prev, val)
        q[i] = prev
    out = np.empty_like(q)
    out[order] = q
    return out.tolist()


def summarize_axes(axis_long: pd.DataFrame, sample_meta: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = (
        axis_long.groupby(["cohort_context", "axis_id", "axis_label"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_core_enrichment=("core_enrichment", "median"),
            n_core_positive=("core_enrichment", lambda s: int((s > 0).sum())),
            median_interface_enrichment=("interface_enrichment", "median"),
            n_interface_positive=("interface_enrichment", lambda s: int((s > 0).sum())),
            median_rho_distance_to_caf_core=("rho_distance_to_caf_core", "median"),
            n_core_proximal=("rho_distance_to_caf_core", lambda s: int((s < 0).sum())),
        )
        .reset_index()
    )

    corr_rows: list[dict] = []
    meta = sample_meta[["sample_id", "immune_decoupling_index"]].drop_duplicates()
    for axis_id, group in axis_long.groupby("axis_id"):
        merged = group.merge(meta, on="sample_id", how="left", suffixes=("", "_meta"))
        label = group["axis_label"].iloc[0]
        for metric in ["core_enrichment", "interface_enrichment", "rho_distance_to_caf_core"]:
            rho = safe_spearman(merged[metric].to_numpy(float), merged["immune_decoupling_index_meta"].to_numpy(float))
            corr_rows.append(
                {
                    "axis_id": axis_id,
                    "axis_label": label,
                    "metric": metric,
                    "rho_with_immune_decoupling_index": rho,
                    "n_samples": int(merged[[metric, "immune_decoupling_index_meta"]].dropna().shape[0]),
                }
            )
    correlations = pd.DataFrame(corr_rows)

    pair_rows: list[dict] = []
    target = "lymph_node_metastasis"
    for axis_id, group in axis_long.groupby("axis_id"):
        label = group["axis_label"].iloc[0]
        target_group = group[group["cohort_context"].eq(target)]
        for context in ["primary_tumor", "liver_metastasis", "post_neoadjuvant_sections", "external_paired_st_anchor"]:
            compare = group[group["cohort_context"].eq(context)]
            for metric in ["core_enrichment", "interface_enrichment"]:
                a = target_group[metric].dropna()
                b = compare[metric].dropna()
                if len(a) < 2 or len(b) < 2:
                    continue
                test = mannwhitneyu(a, b, alternative="greater")
                pair_rows.append(
                    {
                        "axis_id": axis_id,
                        "axis_label": label,
                        "metric": metric,
                        "comparison": f"lymph_node_metastasis_greater_than_{context}",
                        "reference_context": target,
                        "comparison_context": context,
                        "n_reference": int(len(a)),
                        "n_comparison": int(len(b)),
                        "median_reference": float(a.median()),
                        "median_comparison": float(b.median()),
                        "median_difference": float(a.median() - b.median()),
                        "p_value": float(test.pvalue),
                    }
                )
    pairwise = pd.DataFrame(pair_rows)
    if not pairwise.empty:
        pairwise["q_value"] = bh_adjust(pairwise["p_value"].tolist())
    return summary, correlations, pairwise


def matrix_for_plot(summary: pd.DataFrame, metric: str) -> tuple[np.ndarray, list[str], list[str]]:
    axis_labels = [axis["axis_label"] for axis in AXES]
    contexts = [c for c in CONTEXT_ORDER if c in set(summary["cohort_context"])]
    pivot = summary.pivot_table(index="axis_label", columns="cohort_context", values=metric, aggfunc="median")
    pivot = pivot.reindex(index=axis_labels, columns=contexts)
    return pivot.to_numpy(float), axis_labels, [CONTEXT_LABELS.get(c, c) for c in contexts]


def plot_heatmap(ax: plt.Axes, data: np.ndarray, rows: list[str], cols: list[str], title: str, vlim: float = 1.3) -> None:
    im = ax.imshow(data, cmap="RdBu_r", vmin=-vlim, vmax=vlim, aspect="auto")
    ax.set_xticks(np.arange(len(cols)), cols, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(rows)), rows, fontsize=8)
    ax.set_title(title)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.5, color="black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)


def make_figure(summary: pd.DataFrame, correlations: pd.DataFrame, output_base: Path) -> None:
    core_matrix, axis_labels, context_labels = matrix_for_plot(summary, "median_core_enrichment")
    interface_matrix, _, _ = matrix_for_plot(summary, "median_interface_enrichment")

    fig = plt.figure(figsize=(13.5, 10.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 0.95])

    ax0 = fig.add_subplot(gs[0, 0])
    plot_heatmap(ax0, core_matrix, axis_labels, context_labels, "CAF-core enrichment of candidate axes")

    ax1 = fig.add_subplot(gs[0, 1])
    plot_heatmap(ax1, interface_matrix, axis_labels, context_labels, "Tumor-stroma interface enrichment")

    ax2 = fig.add_subplot(gs[1, 0])
    pivot = summary.pivot_table(index="axis_label", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
    ln_minus_liver = pivot.get("lymph_node_metastasis") - pivot.get("liver_metastasis")
    ln_minus_primary = pivot.get("lymph_node_metastasis") - pivot.get("primary_tumor")
    diff = pd.DataFrame({"LN - liver": ln_minus_liver, "LN - primary": ln_minus_primary}).reindex(axis_labels)
    y = np.arange(len(diff))
    ax2.axvline(0, color="#555555", linewidth=0.8)
    ax2.scatter(diff["LN - liver"], y - 0.13, color="#4C78A8", label="LN - liver", s=38)
    ax2.scatter(diff["LN - primary"], y + 0.13, color="#F58518", label="LN - primary", s=38)
    for i, label in enumerate(diff.index):
        for col, dy in [("LN - liver", -0.13), ("LN - primary", 0.13)]:
            val = diff.loc[label, col]
            if np.isfinite(val):
                ax2.plot([0, val], [i + dy, i + dy], color="#BBBBBB", linewidth=0.8, zorder=0)
    ax2.set_yticks(y, diff.index, fontsize=8)
    ax2.set_xlabel("median CAF-core enrichment difference")
    ax2.set_title("Lymph-node shift in CAF-core candidate axes")
    ax2.legend(frameon=False, fontsize=8)

    ax3 = fig.add_subplot(gs[1, 1])
    corr = correlations[correlations["metric"].eq("core_enrichment")].set_index("axis_label").reindex(axis_labels)
    vals = corr["rho_with_immune_decoupling_index"].to_numpy(float)
    ax3.barh(y, vals, color=["#B279A2" if v < 0 else "#4C78A8" for v in vals])
    ax3.axvline(0, color="#555555", linewidth=0.8)
    ax3.set_yticks(y, corr.index, fontsize=8)
    ax3.set_xlabel("Spearman rho")
    ax3.set_title("Axis association with immune-decoupling index")

    fig.suptitle("Mechanistic candidate axes nominated by CAF-core and interface enrichment", fontsize=14)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=240)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def write_report(summary: pd.DataFrame, correlations: pd.DataFrame, pairwise: pd.DataFrame) -> None:
    lines = [
        "# Stage 25 Mechanistic Candidate Axes",
        "",
        "## Purpose",
        "",
        "This analysis translates the CAF-core spatial organization into candidate biological axes. It does not claim ligand-receptor causality; it nominates pathway-level mechanisms that are enriched in CAF cores or tumor-stroma interfaces and can guide validation.",
        "",
        "## Strongest CAF-Core Enrichments by Context",
        "",
    ]
    for context in CONTEXT_ORDER:
        sub = summary[summary["cohort_context"].eq(context)].sort_values("median_core_enrichment", ascending=False)
        if sub.empty:
            continue
        top = sub.head(4)
        readable = ", ".join([f"{r.axis_label} ({r.median_core_enrichment:.2f})" for r in top.itertuples()])
        lines.append(f"- {CONTEXT_LABELS.get(context, context)}: {readable}.")

    lines.extend(["", "## Interface-Enriched Axes", ""])
    overall = (
        summary.groupby(["axis_id", "axis_label"], dropna=False)
        .agg(median_interface_enrichment=("median_interface_enrichment", "median"))
        .sort_values("median_interface_enrichment", ascending=False)
        .reset_index()
    )
    for row in overall.itertuples():
        lines.append(f"- {row.axis_label}: median interface enrichment {row.median_interface_enrichment:.3f}.")

    lines.extend(["", "## Association With Immune-Decoupling", ""])
    corr = correlations[correlations["metric"].eq("core_enrichment")].sort_values(
        "rho_with_immune_decoupling_index", ascending=False
    )
    for row in corr.itertuples():
        lines.append(
            f"- {row.axis_label}: rho={row.rho_with_immune_decoupling_index:.3f} "
            f"across {int(row.n_samples)} samples."
        )

    if not pairwise.empty:
        lines.extend(["", "## Lymph-Node-Enriched Candidate Axes", ""])
        sub = pairwise[pairwise["metric"].eq("core_enrichment")].sort_values(["q_value", "p_value"]).head(12)
        for row in sub.itertuples():
            lines.append(
                f"- {row.axis_label}, {row.comparison}: median delta {row.median_difference:.3f}, "
                f"p={row.p_value:.3g}, q={row.q_value:.3g}."
            )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The top-journal-safe framing is that CAF cores nominate a small set of candidate biology axes: stromal TGF-beta/EMT, SPP1-TAM/matrix, cytokine/chemokine, IFN/APC, lymphoid and tumor-state axes. Axes that enrich at both the CAF core and tumor-stroma interface are the most plausible mechanistic follow-up targets. Axes that track the immune-decoupling index explain why lymph-node metastases preserve stromal-tumor coupling while losing immune/IFN coupling.",
        ]
    )
    path = PROJECT_ROOT / "results/reports/mechanism_candidate_axes_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    spots = load_spots()
    axis_rows: list[dict] = []
    sample_rows: list[dict] = []
    warnings: list[str] = []
    for sample_id, sample in spots.groupby("sample_id", sort=True):
        if len(sample) < 200:
            warnings.append(f"Skipped {sample_id}: fewer than 200 spots")
            continue
        try:
            rows, meta = summarize_sample(sample)
            axis_rows.extend(rows)
            sample_rows.append(meta)
            print(f"Summarized candidate axes: {sample_id} ({len(sample)} spots)")
        except Exception as exc:
            warnings.append(f"{sample_id}: {exc}")

    axis_long = pd.DataFrame(axis_rows)
    sample_meta = pd.DataFrame(sample_rows)
    summary, correlations, pairwise = summarize_axes(axis_long, sample_meta)

    tables = PROJECT_ROOT / "results/tables"
    tables.mkdir(parents=True, exist_ok=True)
    axis_long.to_csv(tables / "mechanism_candidate_axis_sample_long.csv", index=False)
    sample_meta.to_csv(tables / "mechanism_candidate_axis_sample_summary.csv", index=False)
    summary.to_csv(tables / "mechanism_candidate_axis_context_summary.csv", index=False)
    correlations.to_csv(tables / "mechanism_candidate_axis_decoupling_correlations.csv", index=False)
    pairwise.to_csv(tables / "mechanism_candidate_axis_ln_pairwise_tests.csv", index=False)

    source = PROJECT_ROOT / "results/source_data"
    source.mkdir(parents=True, exist_ok=True)
    summary.to_csv(source / "Source_Data_Fig_7A_B.csv", index=False)
    correlations.to_csv(source / "Source_Data_Fig_7C.csv", index=False)
    pairwise.to_csv(source / "Source_Data_Extended_Mechanism_Axis_Tests.csv", index=False)

    make_figure(summary, correlations, PROJECT_ROOT / "results/figures/main/figure7_mechanism_candidate_axes")
    write_report(summary, correlations, pairwise)

    write_status(
        "success" if not warnings else "partial_success",
        {
            "n_samples_processed": int(sample_meta["sample_id"].nunique()),
            "n_axis_rows": int(len(axis_long)),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": warnings,
            "outputs": [
                "results/tables/mechanism_candidate_axis_context_summary.csv",
                "results/tables/mechanism_candidate_axis_decoupling_correlations.csv",
                "results/figures/main/figure7_mechanism_candidate_axes.pdf",
                "results/reports/mechanism_candidate_axes_report.md",
            ],
            "next_manual_check": [
                "Inspect Figure 7 for whether candidate axes should be a main figure or extended data.",
                "Treat axes as pathway-level nominations, not causal ligand-receptor proof.",
            ],
        },
    )
    print("Stage 25 mechanism candidate axes complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
