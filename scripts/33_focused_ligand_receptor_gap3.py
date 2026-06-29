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
OUT_SAMPLE = PROJECT_ROOT / "results/tables/gap3_focused_lr_interface_per_sample.csv"
OUT_SUMMARY = PROJECT_ROOT / "results/tables/gap3_focused_lr_interface_context_summary.csv"
OUT_CORR = PROJECT_ROOT / "results/tables/gap3_focused_lr_interface_correlations.csv"
OUT_FIG = PROJECT_ROOT / "results/figures/submission/extended_data_gap3_focused_lr_interface"
OUT_REPORT = PROJECT_ROOT / "results/reports/gap3_focused_ligand_receptor_interface_report.md"
STATUS = PROJECT_ROOT / "results/logs/stage_33_gap3_focused_lr_interface_status.json"

AXES = {
    "SPP1-CD44/integrin": {
        "ligand": ["SPP1"],
        "receptor": ["CD44", "ITGAV", "ITGB1", "ITGB5", "ITGA5"],
        "response": ["FN1", "SPARC", "MMP14", "LGALS3"],
    },
    "TGF-beta/TGFBR": {
        "ligand": ["TGFB1", "TGFB2", "TGFB3", "INHBA"],
        "receptor": ["TGFBR1", "TGFBR2", "TGFBR3"],
        "response": ["TGFBI", "CTGF", "SERPINE1", "VIM", "MMP14"],
    },
    "matrix-integrin": {
        "ligand": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "FN1", "SPARC", "POSTN"],
        "receptor": ["ITGA5", "ITGAV", "ITGB1", "ITGB5"],
        "response": ["VIM", "MMP14", "SERPINE1", "CTGF"],
    },
    "IL6-OSM/LIF-JAKSTAT": {
        "ligand": ["IL6", "OSM", "LIF"],
        "receptor": ["IL6R", "IL6ST", "OSMR", "LIFR"],
        "response": ["JAK1", "STAT3", "SOCS3"],
    },
}

GENES = sorted({gene for axis in AXES.values() for genes in axis.values() for gene in genes})
CONTEXT_ORDER = ["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
}
CONTEXT_COLORS = {
    "post_neoadjuvant_sections": "#8C6D31",
    "treatment_naive_primary": "#B55A30",
    "primary_tumor": "#4E79A7",
    "liver_metastasis": "#59A14F",
    "lymph_node_metastasis": "#9C755F",
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
        raise RuntimeError("No expression frames available for Gap 3 LR/interface analysis.")
    expr = pd.concat(frames, ignore_index=True)
    if "cohort_context" not in expr.columns:
        expr["cohort_context"] = expr.apply(stage27.context_for_row, axis=1)
    else:
        missing = expr["cohort_context"].isna() | expr["cohort_context"].astype(str).eq("")
        if missing.any():
            expr.loc[missing, "cohort_context"] = expr.loc[missing].apply(stage27.context_for_row, axis=1)
    expr = stage27.add_regions(expr)
    return expr


def add_axis_scores(expr: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    for gene in GENES:
        if gene in out.columns and out[gene].notna().any():
            out[f"z_{gene}"] = out.groupby("sample_id", group_keys=False)[gene].apply(zscore)
    for axis_name, parts in AXES.items():
        for part_name, genes in parts.items():
            cols = [f"z_{gene}" for gene in genes if f"z_{gene}" in out.columns]
            out[f"{axis_name}::{part_name}"] = out[cols].mean(axis=1) if cols else np.nan
        out[f"{axis_name}::directional_score"] = (
            out[f"{axis_name}::ligand"] + out[f"{axis_name}::receptor"] + out[f"{axis_name}::response"]
        ) / 3.0
    return out


def summarize_samples(expr: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    decoupling = pd.read_csv(PROJECT_ROOT / "results/tables/mechanism_candidate_axis_sample_summary.csv")
    dec_cols = [
        "sample_id",
        "stromal_tumor_core_coupling",
        "immune_core_coupling",
        "immune_decoupling_index",
    ]
    for sample_id, sample in expr.groupby("sample_id", sort=False):
        if "is_caf_core_top10" not in sample or sample["is_caf_core_top10"].sum() < 3:
            continue
        core = sample["is_caf_core_top10"].astype(bool)
        interface = sample["is_interface"].astype(bool) if "is_interface" in sample else pd.Series(False, index=sample.index)
        if interface.sum() < 3:
            continue
        row_base = {
            "dataset_id": sample["dataset_id"].iloc[0] if "dataset_id" in sample else "",
            "sample_id": sample_id,
            "specimen_type": sample["specimen_type"].iloc[0] if "specimen_type" in sample else "",
            "cohort_context": sample["cohort_context"].iloc[0],
            "n_spots": int(len(sample)),
            "n_caf_core_spots": int(core.sum()),
            "n_interface_spots": int(interface.sum()),
        }
        for axis_name in AXES:
            ligand = sample[f"{axis_name}::ligand"]
            receptor = sample[f"{axis_name}::receptor"]
            response = sample[f"{axis_name}::response"]
            directional = sample[f"{axis_name}::directional_score"]
            rows.append(
                {
                    **row_base,
                    "axis": axis_name,
                    "ligand_genes": ";".join(AXES[axis_name]["ligand"]),
                    "receptor_genes": ";".join(AXES[axis_name]["receptor"]),
                    "response_genes": ";".join(AXES[axis_name]["response"]),
                    "ligand_core_enrichment": float(ligand[core].mean() - ligand[~core].mean()),
                    "receptor_interface_enrichment": float(receptor[interface].mean() - receptor[~interface].mean()),
                    "response_interface_enrichment": float(response[interface].mean() - response[~interface].mean()),
                    "directional_interface_score": float(directional[interface].mean() - directional[~interface].mean()),
                    "directional_core_to_interface_score": float(
                        (ligand[core].mean() - ligand[~core].mean())
                        + (receptor[interface].mean() - receptor[~interface].mean())
                        + (response[interface].mean() - response[~interface].mean())
                    ),
                }
            )
    per_sample = pd.DataFrame(rows)
    per_sample = per_sample.merge(decoupling[dec_cols], on="sample_id", how="left")
    return per_sample


def summarize_context(per_sample: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "ligand_core_enrichment",
        "receptor_interface_enrichment",
        "response_interface_enrichment",
        "directional_core_to_interface_score",
    ]
    rows: list[dict] = []
    for (context, axis), sub in per_sample.groupby(["cohort_context", "axis"], sort=False):
        for metric in metrics:
            values = sub[metric].dropna()
            if values.empty:
                continue
            rows.append(
                {
                    "cohort_context": context,
                    "axis": axis,
                    "metric": metric,
                    "n_samples": int(values.shape[0]),
                    "median_value": float(values.median()),
                    "n_positive": int((values > 0).sum()),
                    "mean_value": float(values.mean()),
                }
            )
    return pd.DataFrame(rows)


def compute_correlations(per_sample: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    tumor = per_sample[per_sample["cohort_context"].isin(CONTEXT_ORDER)].copy()
    for axis, sub_axis in tumor.groupby("axis"):
        for metric in ["ligand_core_enrichment", "receptor_interface_enrichment", "response_interface_enrichment", "directional_core_to_interface_score"]:
            for target in ["immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]:
                tmp = sub_axis[[metric, target]].dropna()
                if tmp.shape[0] < 5:
                    continue
                rho, p = spearmanr(tmp[metric], tmp[target])
                rows.append(
                    {
                        "axis": axis,
                        "metric": metric,
                        "target": target,
                        "n_samples": int(tmp.shape[0]),
                        "spearman_rho": float(rho),
                        "p_value": float(p),
                    }
                )
    return pd.DataFrame(rows)


def plot_gap3(summary: pd.DataFrame, per_sample: pd.DataFrame, corr: pd.DataFrame) -> None:
    contexts = [c for c in CONTEXT_ORDER if c in set(summary["cohort_context"])]
    axes = list(AXES)
    metrics = ["ligand_core_enrichment", "receptor_interface_enrichment", "response_interface_enrichment"]
    metric_labels = ["ligand\nCAF core", "receptor\ninterface", "response\ninterface"]

    fig = plt.figure(figsize=(13.2, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.15, 1.15, 1.0])
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[:, 1])
    ax_c = fig.add_subplot(gs[:, 2])
    fig.suptitle("Gap 3: focused ligand-receptor/interface nomination", fontsize=15, fontweight="bold")

    rows = []
    row_labels = []
    for axis in axes:
        for metric, label in zip(metrics, metric_labels):
            vals = []
            for context in contexts:
                hit = summary[
                    summary["axis"].eq(axis)
                    & summary["metric"].eq(metric)
                    & summary["cohort_context"].eq(context)
                ]
                vals.append(hit["median_value"].iloc[0] if not hit.empty else np.nan)
            rows.append(vals)
            row_labels.append(f"{axis} | {label}")
    mat = np.asarray(rows, float)
    im = ax_a.imshow(mat, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax_a.set_title("A  Directional component enrichment", loc="left", fontsize=10.5, fontweight="bold")
    ax_a.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax_a.set_yticks(np.arange(len(row_labels)), row_labels)
    ax_a.tick_params(axis="both", labelsize=7.2)
    for i, axis_metric in enumerate(row_labels):
        axis = axis_metric.split(" | ")[0]
        metric_idx = i % len(metrics)
        metric = metrics[metric_idx]
        for j, context in enumerate(contexts):
            hit = summary[
                summary["axis"].eq(axis)
                & summary["metric"].eq(metric)
                & summary["cohort_context"].eq(context)
            ]
            if hit.empty:
                continue
            val = hit["median_value"].iloc[0]
            ax_a.text(j, i, f"{val:.2f}\n{int(hit['n_positive'].iloc[0])}/{int(hit['n_samples'].iloc[0])}", ha="center", va="center", fontsize=5.2)
    cb = fig.colorbar(im, ax=ax_a, fraction=0.045, pad=0.02)
    cb.set_label("median enrichment", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    y_metric = "directional_core_to_interface_score"
    positions = np.arange(len(axes))
    width = 0.13
    for k, context in enumerate(contexts):
        vals = [
            per_sample.loc[per_sample["axis"].eq(axis) & per_sample["cohort_context"].eq(context), y_metric].median()
            for axis in axes
        ]
        ax_b.bar(positions + (k - (len(contexts) - 1) / 2) * width, vals, width=width, color=CONTEXT_COLORS.get(context, "#777777"), label=CONTEXT_LABELS.get(context, context), alpha=0.82)
    ax_b.axhline(0, color="#333333", lw=0.8)
    ax_b.set_title("B  Combined directional nomination", loc="left", fontsize=10.5, fontweight="bold")
    ax_b.set_xticks(positions, axes, rotation=35, ha="right")
    ax_b.set_ylabel("ligand core + receptor/response interface", fontsize=8.5)
    ax_b.tick_params(axis="both", labelsize=8)
    ax_b.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.legend(frameon=False, fontsize=6.4, ncols=1, loc="upper left")

    score_col = y_metric
    scatter_axis = "matrix-integrin"
    scatter = per_sample[per_sample["axis"].eq(scatter_axis) & per_sample["cohort_context"].isin(contexts)][
        [score_col, "immune_decoupling_index", "cohort_context"]
    ].dropna()
    for context, sub in scatter.groupby("cohort_context"):
        ax_c.scatter(
            sub[score_col],
            sub["immune_decoupling_index"],
            s=24 if len(sub) < 20 else 16,
            color=CONTEXT_COLORS.get(context, "#777777"),
            alpha=0.72,
            edgecolor="white",
            linewidth=0.3,
            label=CONTEXT_LABELS.get(context, context),
        )
    if len(scatter) >= 5:
        x = scatter[score_col].to_numpy(float)
        y = scatter["immune_decoupling_index"].to_numpy(float)
        coeff = np.polyfit(x, y, deg=1)
        xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
        ax_c.plot(xx, coeff[0] * xx + coeff[1], color="#222222", lw=1.1)
        rho, p = spearmanr(x, y)
        ax_c.text(
            0.03,
            0.97,
            f"matrix-integrin\nrho = {rho:.2f}\np = {p:.2g}",
            transform=ax_c.transAxes,
            va="top",
            ha="left",
            fontsize=8,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="#DDDDDD", lw=0.6),
        )
    ax_c.set_title("C  Candidate axis vs immune decoupling", loc="left", fontsize=10.5, fontweight="bold")
    ax_c.set_xlabel("matrix-integrin directional nomination score", fontsize=8.5)
    ax_c.set_ylabel("immune-decoupling index", fontsize=8.5)
    ax_c.tick_params(axis="both", labelsize=8)
    ax_c.grid(color="#E5E5E5", lw=0.7)
    ax_c.spines[["top", "right"]].set_visible(False)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = OUT_FIG.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(summary: pd.DataFrame, corr: pd.DataFrame) -> None:
    lines = [
        "# Gap 3 Focused Ligand-Receptor / Interface Report",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Gap 3 Definition",
        "",
        "SPP1-TAM/matrix and TGF-beta/EMT axes are spatially nominated but not yet shown as directional interaction candidates.",
        "",
        "## Main Results",
        "",
    ]
    for axis in AXES:
        hit = summary[
            summary["axis"].eq(axis)
            & summary["metric"].eq("directional_core_to_interface_score")
            & summary["cohort_context"].isin(CONTEXT_ORDER)
        ]
        if hit.empty:
            continue
        best = hit.sort_values("median_value", ascending=False).iloc[0]
        lines.append(
            f"- {axis}: strongest median directional score in {CONTEXT_LABELS.get(best['cohort_context'], best['cohort_context'])} "
            f"({best['median_value']:.3f}; positive in {int(best['n_positive'])}/{int(best['n_samples'])} samples)."
        )

    lines.extend(["", "## Correlation With Immune Decoupling", ""])
    dec = corr[corr["target"].eq("immune_decoupling_index") & corr["metric"].eq("directional_core_to_interface_score")].copy()
    for _, row in dec.sort_values("spearman_rho", ascending=False).iterrows():
        lines.append(
            f"- {row['axis']}: rho {row['spearman_rho']:.3f}, p = {row['p_value']:.3g}, n = {int(row['n_samples'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "This is a spatially constrained nomination of ligand, receptor and response programs across CAF cores and tumor-stroma interfaces. It is not causal ligand-receptor proof and should be described as candidate communication biology.",
            "",
            "## Generated Outputs",
            "",
            f"- `{OUT_SAMPLE.relative_to(PROJECT_ROOT)}`",
            f"- `{OUT_SUMMARY.relative_to(PROJECT_ROOT)}`",
            f"- `{OUT_CORR.relative_to(PROJECT_ROOT)}`",
            f"- `{OUT_FIG.with_suffix('.pdf').relative_to(PROJECT_ROOT)}`",
        ]
    )
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    expr = load_expression()
    expr = add_axis_scores(expr)
    per_sample = summarize_samples(expr)
    summary = summarize_context(per_sample)
    corr = compute_correlations(per_sample)

    OUT_SAMPLE.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(OUT_SAMPLE, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    corr.to_csv(OUT_CORR, index=False)
    plot_gap3(summary, per_sample, corr)
    write_report(summary, corr)

    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "33_gap3_focused_lr_interface",
                "status": "success",
                "n_samples": int(per_sample["sample_id"].nunique()),
                "n_axes": int(per_sample["axis"].nunique()),
                "outputs": [
                    str(OUT_SAMPLE),
                    str(OUT_SUMMARY),
                    str(OUT_CORR),
                    str(OUT_FIG.with_suffix(".pdf")),
                    str(OUT_REPORT),
                ],
                "claim_boundary": "Spatial ligand-receptor/interface nomination only; not causal proof.",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote Gap 3 focused LR/interface analysis for {per_sample['sample_id'].nunique()} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
