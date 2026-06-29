from __future__ import annotations

import importlib.util
import json
import time
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.65,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from scipy.optimize import nnls
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
SCRIPT_36 = ROOT / "scripts" / "36_reference_projection_deconvolution_gap2.py"
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
REPORT_DIR = ROOT / "results" / "reports"
SOURCE_DIR = ROOT / "results" / "source_data"
LOG_DIR = ROOT / "results" / "logs"

for directory in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR, LOG_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

SIGNATURE_FILE = TABLE_DIR / "gse202051_full_reference_projection_signature_matrix.csv"
OLD_PER_SAMPLE = TABLE_DIR / "gap2_full_reference_projection_deconvolution_per_sample.csv"
DECOUPLING_TABLE = TABLE_DIR / "mechanism_candidate_axis_sample_summary.csv"

OUT_PER_SAMPLE = TABLE_DIR / "strict_nnls_reference_deconvolution_per_sample.csv"
OUT_CONTEXT = TABLE_DIR / "strict_nnls_reference_deconvolution_context_summary.csv"
OUT_CORR = TABLE_DIR / "strict_nnls_reference_deconvolution_correlations.csv"
OUT_COMPARE = TABLE_DIR / "strict_nnls_vs_projection_comparison.csv"
OUT_SOURCE = SOURCE_DIR / "Source_Data_Extended_Data_Fig_28_strict_nnls_deconvolution.csv"
OUT_FIG = FIG_DIR / "extended_data_figure28_strict_nnls_deconvolution_sensitivity"
OUT_REPORT = REPORT_DIR / "strict_nnls_reference_deconvolution_report.md"
STATUS = LOG_DIR / "stage_63_strict_nnls_reference_deconvolution_status.json"

STATE_ORDER = ["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "B_plasma", "epithelial_tumor"]
STATE_LABELS = {
    "myCAF_matrix": "myCAF/matrix",
    "iCAF_inflammatory": "iCAF",
    "SPP1_TAM": "SPP1/TAM",
    "DC_APC": "DC/APC",
    "T_NK": "T/NK",
    "B_plasma": "B/plasma",
    "epithelial_tumor": "epithelial",
    "endothelial": "endothelial",
    "neural_schwann": "neural",
    "acinar_normal": "acinar",
}
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


def load_stage36():
    spec = importlib.util.spec_from_file_location("stage36_reference_projection", SCRIPT_36)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {SCRIPT_36}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_signature() -> pd.DataFrame:
    if not SIGNATURE_FILE.exists():
        stage36 = load_stage36()
        stage36.configure_outputs("full")
        stage36.build_reference_signature()
    long = pd.read_csv(SIGNATURE_FILE)
    sig = long.pivot_table(index="cell_state", columns="gene", values="reference_expression", aggfunc="mean").fillna(0.0)
    sig = sig.loc[(sig.sum(axis=1) > 0), (sig.max(axis=0) > 0)]
    return sig


def load_spatial_expression(genes: list[str]) -> pd.DataFrame:
    stage36 = load_stage36()
    return stage36.load_spatial_expression(genes)


def prepare_design(sig: pd.DataFrame, available: list[str]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    a = np.nan_to_num(sig[available].T.to_numpy(float), nan=0.0, posinf=0.0, neginf=0.0)
    gene_scale = np.nanmax(a, axis=1)
    gene_scale[~np.isfinite(gene_scale) | (gene_scale <= 0)] = 1.0
    a = a / gene_scale[:, None]
    col_norm = np.linalg.norm(a, axis=0)
    col_norm[~np.isfinite(col_norm) | (col_norm <= 0)] = 1.0
    a = a / col_norm[None, :]
    return a, gene_scale, col_norm


def nnls_coefficients(sample: pd.DataFrame, sig: pd.DataFrame, genes: list[str]) -> tuple[np.ndarray, list[str]]:
    available = [gene for gene in genes if gene in sample.columns and sample[gene].notna().any()]
    if len(available) < 15:
        return np.zeros((0, sig.shape[0])), []
    a, gene_scale, col_norm = prepare_design(sig, available)
    y = np.nan_to_num(sample[available].fillna(0).to_numpy(float), nan=0.0, posinf=0.0, neginf=0.0)
    y = y / gene_scale[None, :]
    coef = np.zeros((y.shape[0], a.shape[1]), dtype=float)
    for i in range(y.shape[0]):
        coef_i, _ = nnls(a, y[i, :], maxiter=500)
        coef[i, :] = coef_i / col_norm
    denom = coef.sum(axis=1)
    frac = np.divide(coef, denom[:, None], out=np.zeros_like(coef), where=denom[:, None] > 0)
    return frac, available


def summarize_sample(sample: pd.DataFrame, sig: pd.DataFrame, genes: list[str]) -> tuple[pd.DataFrame, dict[str, object]]:
    states = list(sig.index)
    frac, available = nnls_coefficients(sample, sig, genes)
    if frac.shape[0] == 0:
        return pd.DataFrame(), {"status": "skipped", "reason": "insufficient_genes", "n_genes": len(available)}
    core = sample["is_caf_core_top10"].astype(bool).to_numpy()
    if core.sum() < 3 or (~core).sum() < 3:
        return pd.DataFrame(), {"status": "skipped", "reason": "insufficient_core_spots", "n_genes": len(available)}
    rows = []
    caf_score = sample["score_caf_myeloid_barrier"].to_numpy(float)
    for idx, state in enumerate(states):
        values = frac[:, idx]
        rho, p = spearmanr(caf_score, values, nan_policy="omit")
        rows.append(
            {
                "dataset_id": sample["dataset_id"].iloc[0],
                "sample_id": sample["sample_id"].iloc[0],
                "specimen_type": sample["specimen_type"].iloc[0] if "specimen_type" in sample.columns else "",
                "cohort_context": sample["cohort_context"].iloc[0],
                "n_spots": int(len(sample)),
                "n_caf_core_spots": int(core.sum()),
                "n_genes_used": int(len(available)),
                "cell_state": state,
                "nnls_core_enrichment": float(values[core].mean() - values[~core].mean()),
                "nnls_core_mean_fraction": float(values[core].mean()),
                "nnls_noncore_mean_fraction": float(values[~core].mean()),
                "nnls_spot_spearman_with_caf_myeloid": float(rho) if np.isfinite(rho) else np.nan,
                "nnls_spot_spearman_p": float(p) if np.isfinite(p) else np.nan,
            }
        )
    return pd.DataFrame(rows), {"status": "done", "n_genes": len(available), "n_spots": len(sample)}


def append_rows(path: Path, rows: pd.DataFrame) -> None:
    if rows.empty:
        return
    header = not path.exists()
    rows.to_csv(path, mode="a", header=header, index=False)


def run_nnls(sig: pd.DataFrame, expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    completed: set[str] = set()
    if OUT_PER_SAMPLE.exists():
        try:
            completed = set(pd.read_csv(OUT_PER_SAMPLE, usecols=["sample_id"])["sample_id"].astype(str).unique())
        except Exception:
            completed = set()

    total = expr["sample_id"].nunique()
    start = time.time()
    processed = 0
    status_rows = []
    for sample_id, sample in expr.groupby("sample_id", sort=False):
        sample_id = str(sample_id)
        if sample_id in completed:
            processed += 1
            continue
        t0 = time.time()
        rows, status = summarize_sample(sample.copy(), sig, genes)
        if not rows.empty:
            append_rows(OUT_PER_SAMPLE, rows)
        status_rows.append({"sample_id": sample_id, **status, "elapsed_sec": round(time.time() - t0, 3)})
        processed += 1
        STATUS.write_text(
            json.dumps(
                {
                    "stage": "63_strict_nnls_reference_deconvolution",
                    "status": "running",
                    "processed_samples": processed,
                    "total_samples": int(total),
                    "last_sample": sample_id,
                    "elapsed_min": round((time.time() - start) / 60, 2),
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        if processed % 10 == 0:
            print(f"Processed {processed}/{total} samples; last={sample_id}; elapsed={(time.time()-start)/60:.1f} min", flush=True)
    if status_rows:
        status_path = LOG_DIR / "stage_63_strict_nnls_sample_status.csv"
        pd.DataFrame(status_rows).to_csv(status_path, mode="a", header=not status_path.exists(), index=False)
    return pd.read_csv(OUT_PER_SAMPLE)


def summarize_outputs(per_sample: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    context_rows = []
    for (context, state), sub in per_sample.groupby(["cohort_context", "cell_state"], sort=False):
        vals = sub["nnls_core_enrichment"].dropna()
        rhos = sub["nnls_spot_spearman_with_caf_myeloid"].dropna()
        context_rows.append(
            {
                "cohort_context": context,
                "cell_state": state,
                "n_samples": int(vals.shape[0]),
                "median_nnls_core_enrichment": float(vals.median()) if len(vals) else np.nan,
                "n_nnls_core_positive": int((vals > 0).sum()),
                "median_nnls_spot_spearman_with_caf_myeloid": float(rhos.median()) if len(rhos) else np.nan,
                "n_nnls_spearman_positive": int((rhos > 0).sum()),
                "median_n_genes_used": float(sub["n_genes_used"].median()),
            }
        )
    context = pd.DataFrame(context_rows)
    context.to_csv(OUT_CONTEXT, index=False)

    corr_rows = []
    if DECOUPLING_TABLE.exists():
        dec = pd.read_csv(DECOUPLING_TABLE)[["sample_id", "immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]]
        merged = per_sample.merge(dec, on="sample_id", how="left")
        tumor = merged[~merged["cohort_context"].eq("normal_pancreas")].copy()
        for state, sub in tumor.groupby("cell_state"):
            for target in ["immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]:
                tmp = sub[["nnls_core_enrichment", target]].dropna()
                if len(tmp) < 5:
                    continue
                rho, p = spearmanr(tmp["nnls_core_enrichment"], tmp[target])
                corr_rows.append(
                    {
                        "cell_state": state,
                        "metric": "nnls_core_enrichment",
                        "target": target,
                        "n_samples": int(len(tmp)),
                        "spearman_rho": float(rho),
                        "p_value": float(p),
                    }
                )
    corr = pd.DataFrame(corr_rows)
    corr.to_csv(OUT_CORR, index=False)

    compare_rows = []
    if OLD_PER_SAMPLE.exists():
        old = pd.read_csv(OLD_PER_SAMPLE)
        merged = per_sample.merge(old, on=["sample_id", "cell_state"], how="inner")
        for state, sub in merged.groupby("cell_state"):
            tmp = sub[["nnls_core_enrichment", "core_enrichment"]].dropna()
            if len(tmp) < 5:
                continue
            rho, p = spearmanr(tmp["nnls_core_enrichment"], tmp["core_enrichment"])
            compare_rows.append(
                {
                    "cell_state": state,
                    "n_pairs": int(len(tmp)),
                    "nnls_vs_projection_spearman_rho": float(rho),
                    "p_value": float(p),
                    "same_direction_fraction": float((np.sign(tmp["nnls_core_enrichment"]) == np.sign(tmp["core_enrichment"])).mean()),
                    "median_nnls_core_enrichment": float(tmp["nnls_core_enrichment"].median()),
                    "median_projection_core_enrichment": float(tmp["core_enrichment"].median()),
                }
            )
    compare = pd.DataFrame(compare_rows)
    compare.to_csv(OUT_COMPARE, index=False)
    return context, corr, compare


def plot_results(context: pd.DataFrame, corr: pd.DataFrame, compare: pd.DataFrame) -> None:
    states = [s for s in STATE_ORDER if s in set(context["cell_state"])]
    contexts = [c for c in CONTEXT_ORDER if c in set(context["cohort_context"])]
    mat = (
        context.pivot_table(index="cell_state", columns="cohort_context", values="median_nnls_core_enrichment", aggfunc="median")
        .reindex(index=states, columns=contexts)
    )
    rho_mat = (
        context.pivot_table(index="cell_state", columns="cohort_context", values="median_nnls_spot_spearman_with_caf_myeloid", aggfunc="median")
        .reindex(index=states, columns=contexts)
    )
    corr_plot = corr[corr["target"].eq("immune_decoupling_index") & corr["cell_state"].isin(states)].copy()
    corr_plot["cell_state"] = pd.Categorical(corr_plot["cell_state"], categories=states, ordered=True)
    corr_plot = corr_plot.sort_values("cell_state")
    compare_plot = compare[compare["cell_state"].isin(states)].copy()
    compare_plot["cell_state"] = pd.Categorical(compare_plot["cell_state"], categories=states, ordered=True)
    compare_plot = compare_plot.sort_values("cell_state")

    fig = plt.figure(figsize=(14.8, 10.8), constrained_layout=False)
    gs = GridSpec(3, 6, figure=fig, height_ratios=[1.0, 0.95, 0.46], hspace=0.72, wspace=0.82)
    fig.suptitle("Strict NNLS reference-deconvolution sensitivity", fontsize=15, fontweight="bold", y=0.985)
    ax_a = fig.add_subplot(gs[0, 0:3])
    ax_b = fig.add_subplot(gs[0, 3:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:6])
    ax_note = fig.add_subplot(gs[2, :])

    im = ax_a.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.08, vmax=0.08, aspect="auto")
    ax_a.set_title("A  NNLS CAF-core enrichment", loc="left", fontsize=10.5, fontweight="bold")
    ax_a.set_xticks(np.arange(len(contexts)))
    ax_a.set_xticklabels([CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right", fontsize=7)
    ax_a.set_yticks(np.arange(len(states)))
    ax_a.set_yticklabels([STATE_LABELS.get(s, s) for s in states], fontsize=7)
    for i, state in enumerate(states):
        for j, ctx in enumerate(contexts):
            val = mat.loc[state, ctx]
            if np.isfinite(val):
                row = context[(context["cell_state"].eq(state)) & (context["cohort_context"].eq(ctx))]
                support = ""
                if not row.empty:
                    support = f"\n{int(row['n_nnls_core_positive'].iloc[0])}/{int(row['n_samples'].iloc[0])}"
                ax_a.text(j, i, f"{val:.2f}{support}", ha="center", va="center", fontsize=5.7, color="#FFFFFF" if abs(val) > 0.05 else "#222222")
    cb = fig.colorbar(im, ax=ax_a, fraction=0.04, pad=0.02)
    cb.set_label("core - noncore fraction", fontsize=7)
    cb.ax.tick_params(labelsize=6)

    im2 = ax_b.imshow(rho_mat.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax_b.set_title("B  NNLS fraction vs CAF-myeloid score", loc="left", fontsize=10.5, fontweight="bold")
    ax_b.set_xticks(np.arange(len(contexts)))
    ax_b.set_xticklabels([CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right", fontsize=7)
    ax_b.set_yticks(np.arange(len(states)))
    ax_b.set_yticklabels([STATE_LABELS.get(s, s) for s in states], fontsize=7)
    for i, state in enumerate(states):
        for j, ctx in enumerate(contexts):
            val = rho_mat.loc[state, ctx]
            if np.isfinite(val):
                ax_b.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > 0.30 else "#222222")
    cb2 = fig.colorbar(im2, ax=ax_b, fraction=0.04, pad=0.02)
    cb2.set_label("median Spearman rho", fontsize=7)
    cb2.ax.tick_params(labelsize=6)

    y = np.arange(len(corr_plot))
    colors = np.where(corr_plot["spearman_rho"] > 0, "#B23A48", "#4C78A8")
    ax_c.axvline(0, color="#333333", lw=0.7)
    ax_c.hlines(y, 0, corr_plot["spearman_rho"], color=colors, lw=1.8)
    ax_c.scatter(corr_plot["spearman_rho"], y, color=colors, s=42, edgecolor="white", linewidth=0.5)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([STATE_LABELS.get(str(s), str(s)) for s in corr_plot["cell_state"]], fontsize=7)
    ax_c.invert_yaxis()
    ax_c.set_xlim(-0.82, 0.82)
    ax_c.set_xlabel("rho with immune-decoupling index", fontsize=7)
    ax_c.set_title("C  NNLS links to immune decoupling", loc="left", fontsize=10.5, fontweight="bold")
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.grid(axis="x", color="#E7E7E7", linewidth=0.55)

    y2 = np.arange(len(compare_plot))
    ax_d.barh(y2, compare_plot["nnls_vs_projection_spearman_rho"], color="#2C7A51", height=0.52)
    ax_d.set_yticks(y2)
    ax_d.set_yticklabels([STATE_LABELS.get(str(s), str(s)) for s in compare_plot["cell_state"]], fontsize=7)
    ax_d.invert_yaxis()
    ax_d.set_xlim(-0.1, 1.0)
    ax_d.set_xlabel("NNLS vs projection per-sample rho", fontsize=7)
    ax_d.set_title("D  Agreement with previous projection", loc="left", fontsize=10.5, fontweight="bold")
    ax_d.spines[["top", "right"]].set_visible(False)
    ax_d.grid(axis="x", color="#E7E7E7", linewidth=0.55)
    for i, row in compare_plot.reset_index(drop=True).iterrows():
        ax_d.text(row["nnls_vs_projection_spearman_rho"] + 0.025, i, f"{row['same_direction_fraction']:.2f}", va="center", fontsize=6.5)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.72,
        "Interpretation: strict per-spot NNLS deconvolution uses the full GSE202051-derived state signatures and non-negative coefficients. "
        "It is used as an orthogonal sensitivity analysis for ED26, asking whether the same CAF-core cell-state directions survive a stricter fitting constraint.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.28,
        "Boundary: this is still reference-dependent computational deconvolution, not immunostaining, segmentation or single-cell-resolved ground truth.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT_FIG.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_source_data(context: pd.DataFrame, corr: pd.DataFrame, compare: pd.DataFrame) -> None:
    rows: list[dict[str, object]] = []
    for _, row in context.iterrows():
        rows.append({"panel": "A-B", "item": f"{row['cohort_context']}|{row['cell_state']}", "metric": "median_nnls_core_enrichment", "value": row["median_nnls_core_enrichment"]})
        rows.append({"panel": "A-B", "item": f"{row['cohort_context']}|{row['cell_state']}", "metric": "median_nnls_spot_spearman_with_caf_myeloid", "value": row["median_nnls_spot_spearman_with_caf_myeloid"]})
    for _, row in corr[corr["target"].eq("immune_decoupling_index")].iterrows():
        rows.append({"panel": "C", "item": row["cell_state"], "metric": "rho_with_immune_decoupling", "value": row["spearman_rho"]})
    for _, row in compare.iterrows():
        rows.append({"panel": "D", "item": row["cell_state"], "metric": "nnls_vs_projection_spearman_rho", "value": row["nnls_vs_projection_spearman_rho"]})
        rows.append({"panel": "D", "item": row["cell_state"], "metric": "same_direction_fraction", "value": row["same_direction_fraction"]})
    pd.DataFrame(rows).to_csv(OUT_SOURCE, index=False)


def write_report(per_sample: pd.DataFrame, context: pd.DataFrame, corr: pd.DataFrame, compare: pd.DataFrame) -> None:
    key_states = [s for s in STATE_ORDER if s in set(context["cell_state"])]
    lines = [
        "# Strict NNLS Reference-Deconvolution Sensitivity Report",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Purpose",
        "",
        "This analysis reruns the GSE202051 full-reference cell-state support layer using per-spot non-negative least-squares coefficients. It is stricter than the earlier clipped least-squares projection and is used as an orthogonal sensitivity analysis for the cell-state interpretation.",
        "",
        "## Scale",
        "",
        f"- Samples analyzed: {per_sample['sample_id'].nunique()}",
        f"- Cell states: {per_sample['cell_state'].nunique()}",
        f"- Median matched genes per sample: {per_sample['n_genes_used'].median():.0f}",
        "",
        "## Context Summary",
        "",
        "| context | state | n samples | median NNLS core enrichment | positive samples | median spot rho |",
        "|---|---|---:|---:|---:|---:|",
    ]
    show = context[context["cell_state"].isin(key_states)].copy()
    show["context_label"] = show["cohort_context"].map(CONTEXT_LABELS).fillna(show["cohort_context"])
    show["state_label"] = show["cell_state"].map(STATE_LABELS).fillna(show["cell_state"])
    for row in show.sort_values(["cell_state", "cohort_context"]).itertuples(index=False):
        lines.append(
            f"| {row.context_label} | {row.state_label} | {int(row.n_samples)} | {row.median_nnls_core_enrichment:.4f} | {int(row.n_nnls_core_positive)}/{int(row.n_samples)} | {row.median_nnls_spot_spearman_with_caf_myeloid:.3f} |"
        )
    lines.extend(["", "## Immune-Decoupling Correlations", "", "| state | n | rho | p |", "|---|---:|---:|---:|"])
    show_corr = corr[corr["target"].eq("immune_decoupling_index") & corr["cell_state"].isin(key_states)].copy()
    show_corr["state_label"] = show_corr["cell_state"].map(STATE_LABELS).fillna(show_corr["cell_state"])
    for row in show_corr.sort_values("spearman_rho").itertuples(index=False):
        lines.append(f"| {row.state_label} | {int(row.n_samples)} | {row.spearman_rho:.3f} | {row.p_value:.3g} |")
    lines.extend(["", "## Agreement With Earlier Full-Reference Projection", "", "| state | n pairs | rho | same-direction fraction |", "|---|---:|---:|---:|"])
    comp = compare[compare["cell_state"].isin(key_states)].copy()
    comp["state_label"] = comp["cell_state"].map(STATE_LABELS).fillna(comp["cell_state"])
    for row in comp.sort_values("nnls_vs_projection_spearman_rho", ascending=False).itertuples(index=False):
        lines.append(f"| {row.state_label} | {int(row.n_pairs)} | {row.nnls_vs_projection_spearman_rho:.3f} | {row.same_direction_fraction:.3f} |")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "This is stricter than the prior projection because coefficients are constrained to be non-negative per spot. It remains a reference-dependent computational deconvolution sensitivity analysis, not immunostaining, segmentation or single-cell-resolved ground truth.",
            "",
            "## Generated Outputs",
            "",
            f"- `{OUT_PER_SAMPLE.relative_to(ROOT).as_posix()}`",
            f"- `{OUT_CONTEXT.relative_to(ROOT).as_posix()}`",
            f"- `{OUT_CORR.relative_to(ROOT).as_posix()}`",
            f"- `{OUT_COMPARE.relative_to(ROOT).as_posix()}`",
            f"- `{OUT_FIG.relative_to(ROOT).as_posix()}.pdf`",
            f"- `{OUT_SOURCE.relative_to(ROOT).as_posix()}`",
        ]
    )
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    sig = load_signature()
    genes = list(sig.columns)
    expr = load_spatial_expression(genes)
    per_sample = run_nnls(sig, expr, genes)
    context, corr, compare = summarize_outputs(per_sample)
    plot_results(context, corr, compare)
    write_source_data(context, corr, compare)
    write_report(per_sample, context, corr, compare)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "63_strict_nnls_reference_deconvolution",
                "status": "success",
                "n_samples": int(per_sample["sample_id"].nunique()),
                "n_states": int(per_sample["cell_state"].nunique()),
                "outputs": [
                    str(OUT_PER_SAMPLE),
                    str(OUT_CONTEXT),
                    str(OUT_CORR),
                    str(OUT_COMPARE),
                    str(OUT_FIG.with_suffix(".pdf")),
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote strict NNLS deconvolution sensitivity for {per_sample['sample_id'].nunique()} samples")


if __name__ == "__main__":
    main()
