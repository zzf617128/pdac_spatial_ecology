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
OUT_SAMPLE = PROJECT_ROOT / "results/tables/gap2_cell_state_marker_attribution_per_sample.csv"
OUT_SUMMARY = PROJECT_ROOT / "results/tables/gap2_cell_state_marker_attribution_context_summary.csv"
OUT_CORR = PROJECT_ROOT / "results/tables/gap2_cell_state_marker_attribution_correlations.csv"
OUT_FIG = PROJECT_ROOT / "results/figures/submission/extended_data_gap2_cell_state_attribution"
OUT_REPORT = PROJECT_ROOT / "results/reports/gap2_cell_state_marker_attribution_report.md"
STATUS = PROJECT_ROOT / "results/logs/stage_34_gap2_cell_state_marker_attribution_status.json"

CELL_STATES = {
    "myCAF/matrix": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "POSTN", "PDGFRB"],
    "iCAF/inflammatory CAF": ["IL6", "CXCL12", "CXCL14", "CFD", "DPT", "HAS1", "PDGFRA", "LIF"],
    "SPP1/TAM": ["SPP1", "TREM2", "APOE", "LGALS3", "CD68", "C1QA", "C1QB", "C1QC", "LST1", "FCGR3A"],
    "DC/APC": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "LAMP3", "CCR7", "CLEC10A", "FCER1A"],
    "T/NK cell": ["CD3D", "CD3E", "TRAC", "CD8A", "CD4", "IL7R", "NKG7", "GZMB"],
    "B/plasma cell": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "CXCL13"],
    "epithelial/tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "MSLN", "CEACAM6"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "neural/Schwann": ["SOX10", "S100B", "PLP1", "MPZ", "NRXN1"],
    "acinar/normal pancreas": ["PRSS1", "CPA1", "REG1A"],
}

GENES = sorted({gene for genes in CELL_STATES.values() for gene in genes})
CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
    "normal_pancreas",
]
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
}
COLORS = {
    "post_neoadjuvant_sections": "#8C6D31",
    "treatment_naive_primary": "#B55A30",
    "primary_tumor": "#4E79A7",
    "liver_metastasis": "#59A14F",
    "lymph_node_metastasis": "#9C755F",
    "normal_pancreas": "#8A8A8A",
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
        raise RuntimeError("No expression frames available for Gap 2 marker attribution analysis.")
    expr = pd.concat(frames, ignore_index=True)
    if "cohort_context" not in expr.columns:
        expr["cohort_context"] = expr.apply(stage27.context_for_row, axis=1)
    else:
        missing = expr["cohort_context"].isna() | expr["cohort_context"].astype(str).eq("")
        if missing.any():
            expr.loc[missing, "cohort_context"] = expr.loc[missing].apply(stage27.context_for_row, axis=1)
    expr = stage27.add_regions(expr)
    return expr


def add_cell_state_scores(expr: pd.DataFrame) -> pd.DataFrame:
    out = expr.copy()
    for gene in GENES:
        if gene in out.columns and out[gene].notna().any():
            out[f"z_{gene}"] = out.groupby("sample_id", group_keys=False)[gene].apply(zscore)
    for state, genes in CELL_STATES.items():
        cols = [f"z_{gene}" for gene in genes if f"z_{gene}" in out.columns]
        out[state] = out[cols].mean(axis=1) if cols else np.nan
    return out


def summarize_samples(expr: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for sample_id, sample in expr.groupby("sample_id", sort=False):
        if "is_caf_core_top10" not in sample or sample["is_caf_core_top10"].sum() < 3:
            continue
        core = sample["is_caf_core_top10"].astype(bool)
        context = sample["cohort_context"].iloc[0]
        row_base = {
            "dataset_id": sample["dataset_id"].iloc[0] if "dataset_id" in sample else "",
            "sample_id": sample_id,
            "specimen_type": sample["specimen_type"].iloc[0] if "specimen_type" in sample else "",
            "cohort_context": context,
            "n_spots": int(len(sample)),
            "n_caf_core_spots": int(core.sum()),
        }
        for state, genes in CELL_STATES.items():
            values = sample[state]
            present = [gene for gene in genes if f"z_{gene}" in sample.columns and sample[f"z_{gene}"].notna().any()]
            if values.notna().sum() < 5:
                continue
            rho, p = spearmanr(sample["score_caf_myeloid_barrier"], values, nan_policy="omit")
            rows.append(
                {
                    **row_base,
                    "cell_state": state,
                    "n_genes_present": len(present),
                    "present_genes": ";".join(present),
                    "core_enrichment": float(values[core].mean() - values[~core].mean()),
                    "core_mean": float(values[core].mean()),
                    "noncore_mean": float(values[~core].mean()),
                    "spot_spearman_with_caf_myeloid": float(rho) if np.isfinite(rho) else np.nan,
                    "spot_spearman_p": float(p) if np.isfinite(p) else np.nan,
                }
            )
    return pd.DataFrame(rows)


def summarize_context(per_sample: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (context, state), sub in per_sample.groupby(["cohort_context", "cell_state"], sort=False):
        values = sub["core_enrichment"].dropna()
        corr = sub["spot_spearman_with_caf_myeloid"].dropna()
        rows.append(
            {
                "cohort_context": context,
                "cell_state": state,
                "n_samples": int(values.shape[0]),
                "median_core_enrichment": float(values.median()) if len(values) else np.nan,
                "n_core_positive": int((values > 0).sum()),
                "median_spot_spearman_with_caf_myeloid": float(corr.median()) if len(corr) else np.nan,
                "n_spearman_positive": int((corr > 0).sum()),
            }
        )
    return pd.DataFrame(rows)


def compute_correlations(per_sample: pd.DataFrame) -> pd.DataFrame:
    dec = pd.read_csv(PROJECT_ROOT / "results/tables/mechanism_candidate_axis_sample_summary.csv")
    dec = dec[["sample_id", "immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]]
    merged = per_sample.merge(dec, on="sample_id", how="left")
    rows: list[dict] = []
    tumor = merged[merged["cohort_context"].isin(CONTEXT_ORDER[:-1])].copy()
    for state, sub in tumor.groupby("cell_state"):
        for target in ["immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]:
            tmp = sub[["core_enrichment", target]].dropna()
            if len(tmp) < 5:
                continue
            rho, p = spearmanr(tmp["core_enrichment"], tmp[target])
            rows.append(
                {
                    "cell_state": state,
                    "metric": "core_enrichment",
                    "target": target,
                    "n_samples": int(len(tmp)),
                    "spearman_rho": float(rho),
                    "p_value": float(p),
                }
            )
    return pd.DataFrame(rows)


def plot_gap2(summary: pd.DataFrame, per_sample: pd.DataFrame, corr: pd.DataFrame) -> None:
    contexts = [c for c in CONTEXT_ORDER if c in set(summary["cohort_context"])]
    states = list(CELL_STATES)
    mat = (
        summary.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
        .reindex(index=states, columns=contexts)
    )
    corr_mat = (
        summary.pivot_table(index="cell_state", columns="cohort_context", values="median_spot_spearman_with_caf_myeloid", aggfunc="median")
        .reindex(index=states, columns=contexts)
    )

    fig = plt.figure(figsize=(13.2, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.38, 1.02, 1.02])
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    fig.suptitle("Gap 2: marker-level cell-state attribution of CAF-myeloid cores", fontsize=15, fontweight="bold")

    im = ax_a.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax_a.set_title("A  CAF-core marker enrichment", loc="left", fontsize=10.5, fontweight="bold")
    ax_a.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax_a.set_yticks(np.arange(len(states)), states)
    ax_a.tick_params(axis="both", labelsize=7.5)
    for i, state in enumerate(states):
        for j, context in enumerate(contexts):
            val = mat.loc[state, context]
            if np.isfinite(val):
                row = summary[summary["cell_state"].eq(state) & summary["cohort_context"].eq(context)].iloc[0]
                ax_a.text(j, i, f"{val:.2f}\n{int(row.n_core_positive)}/{int(row.n_samples)}", ha="center", va="center", fontsize=5.2)
    cb = fig.colorbar(im, ax=ax_a, fraction=0.045, pad=0.02)
    cb.set_label("median core enrichment", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    im2 = ax_b.imshow(corr_mat.to_numpy(float), cmap="RdBu_r", vmin=-0.65, vmax=0.65, aspect="auto")
    ax_b.set_title("B  Spot-level correlation with CAF-myeloid score", loc="left", fontsize=10.5, fontweight="bold")
    ax_b.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax_b.set_yticks(np.arange(len(states)), states)
    ax_b.tick_params(axis="both", labelsize=7.3)
    for i, state in enumerate(states):
        for j, context in enumerate(contexts):
            val = corr_mat.loc[state, context]
            if np.isfinite(val):
                ax_b.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.4)
    cb2 = fig.colorbar(im2, ax=ax_b, fraction=0.045, pad=0.02)
    cb2.set_label("median Spearman rho", fontsize=8)
    cb2.ax.tick_params(labelsize=7)

    plot_states = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK cell", "B/plasma cell", "epithelial/tumor"]
    positions = np.arange(len(plot_states))
    width = 0.14
    plot_contexts = [c for c in ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"] if c in contexts]
    for k, context in enumerate(plot_contexts):
        vals = [
            per_sample.loc[
                per_sample["cell_state"].eq(state) & per_sample["cohort_context"].eq(context),
                "core_enrichment",
            ].median()
            for state in plot_states
        ]
        ax_c.bar(
            positions + (k - (len(plot_contexts) - 1) / 2) * width,
            vals,
            width=width,
            color=COLORS.get(context, "#777777"),
            alpha=0.82,
            label=CONTEXT_LABELS.get(context, context),
        )
    ax_c.axhline(0, color="#333333", lw=0.8)
    ax_c.set_title("C  Primary/metastatic cell-state contrast", loc="left", fontsize=10.5, fontweight="bold")
    ax_c.set_xticks(positions, plot_states, rotation=35, ha="right")
    ax_c.set_ylabel("CAF-core enrichment", fontsize=8.5)
    ax_c.tick_params(axis="both", labelsize=8)
    ax_c.grid(axis="y", color="#E5E5E5", lw=0.7)
    ax_c.spines[["top", "right"]].set_visible(False)
    ax_c.legend(frameon=False, fontsize=7, loc="upper right")

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = OUT_FIG.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(summary: pd.DataFrame, corr: pd.DataFrame) -> None:
    def get_line(context: str, state: str) -> str:
        hit = summary[summary["cohort_context"].eq(context) & summary["cell_state"].eq(state)]
        if hit.empty:
            return "NA"
        row = hit.iloc[0]
        return f"{row['median_core_enrichment']:.3f}; {int(row['n_core_positive'])}/{int(row['n_samples'])}"

    lines = [
        "# Gap 2 Cell-State Marker Attribution Report",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Gap 2 Definition",
        "",
        "The CAF-myeloid core is currently module-defined. This analysis tests whether CAF cores are supported by marker-level enrichment of CAF and myeloid/TAM states, and whether immune marker states decouple in lymph-node metastases.",
        "",
        "## Main Results",
        "",
        f"- myCAF/matrix marker enrichment in CAF cores: post-NACT {get_line('post_neoadjuvant_sections', 'myCAF/matrix')}; primary {get_line('primary_tumor', 'myCAF/matrix')}; liver met {get_line('liver_metastasis', 'myCAF/matrix')}; LN met {get_line('lymph_node_metastasis', 'myCAF/matrix')}.",
        f"- SPP1/TAM marker enrichment in CAF cores: post-NACT {get_line('post_neoadjuvant_sections', 'SPP1/TAM')}; primary {get_line('primary_tumor', 'SPP1/TAM')}; liver met {get_line('liver_metastasis', 'SPP1/TAM')}; LN met {get_line('lymph_node_metastasis', 'SPP1/TAM')}.",
        f"- DC/APC marker enrichment in CAF cores: primary {get_line('primary_tumor', 'DC/APC')}; liver met {get_line('liver_metastasis', 'DC/APC')}; LN met {get_line('lymph_node_metastasis', 'DC/APC')}.",
        f"- T/NK marker enrichment in CAF cores: primary {get_line('primary_tumor', 'T/NK cell')}; liver met {get_line('liver_metastasis', 'T/NK cell')}; LN met {get_line('lymph_node_metastasis', 'T/NK cell')}.",
        "",
        "## Interpretation Boundary",
        "",
        "This is marker-level cell-state attribution, not formal deconvolution, segmentation or orthogonal immunostaining. It supports whether the module-defined CAF-myeloid core is aligned with expected marker programs.",
        "",
        "## Generated Outputs",
        "",
        f"- `{OUT_SAMPLE.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_SUMMARY.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_CORR.relative_to(PROJECT_ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.pdf').relative_to(PROJECT_ROOT)}`",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    expr = load_expression()
    expr = add_cell_state_scores(expr)
    per_sample = summarize_samples(expr)
    summary = summarize_context(per_sample)
    corr = compute_correlations(per_sample)

    OUT_SAMPLE.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(OUT_SAMPLE, index=False)
    summary.to_csv(OUT_SUMMARY, index=False)
    corr.to_csv(OUT_CORR, index=False)
    plot_gap2(summary, per_sample, corr)
    write_report(summary, corr)

    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "34_gap2_cell_state_marker_attribution",
                "status": "success",
                "n_samples": int(per_sample["sample_id"].nunique()),
                "n_cell_states": int(per_sample["cell_state"].nunique()),
                "outputs": [
                    str(OUT_SAMPLE),
                    str(OUT_SUMMARY),
                    str(OUT_CORR),
                    str(OUT_FIG.with_suffix(".pdf")),
                    str(OUT_REPORT),
                ],
                "claim_boundary": "Marker-level attribution only; not formal deconvolution or orthogonal validation.",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote Gap 2 marker attribution analysis for {per_sample['sample_id'].nunique()} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
