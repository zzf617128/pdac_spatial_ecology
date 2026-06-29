from __future__ import annotations

import argparse
import importlib.util
import json
from pathlib import Path

import h5py
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
from scipy import sparse
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT_27 = PROJECT / "scripts" / "27_targeted_gene_axis_validation.py"
H5AD_FILES = [
    PROJECT / "data/external/GSE202051/GSE202051_adata_010nuc_10x.h5ad",
    PROJECT / "data/external/GSE202051/GSE202051_adata_010orgCRT_10x.h5ad",
]
FULL_H5AD_FILES = [
    PROJECT / "data/external/GSE202051/GSE202051_totaldata-final-toshare.h5ad",
]
REFERENCE_LABEL = "GSE202051 h5ad subsets"

OUT_SIGNATURE = PROJECT / "results/tables/gse202051_reference_projection_signature_matrix.csv"
OUT_PER_SAMPLE = PROJECT / "results/tables/gap2_reference_projection_deconvolution_per_sample.csv"
OUT_CONTEXT = PROJECT / "results/tables/gap2_reference_projection_deconvolution_context_summary.csv"
OUT_CORR = PROJECT / "results/tables/gap2_reference_projection_deconvolution_correlations.csv"
OUT_REPORT = PROJECT / "results/reports/gap2_reference_projection_deconvolution_report.md"
OUT_FIG = PROJECT / "results/figures/submission/extended_data_gap2_reference_projection_deconvolution"
STATUS = PROJECT / "results/logs/stage_36_reference_projection_deconvolution_status.json"

STATE_GENES = {
    "myCAF_matrix": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "POSTN", "PDGFRB", "SPARC", "FN1"],
    "iCAF_inflammatory": ["IL6", "CXCL12", "CXCL14", "CFD", "DPT", "HAS1", "PDGFRA", "LIF", "CTGF"],
    "SPP1_TAM": ["SPP1", "TREM2", "APOE", "LGALS3", "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LST1", "FCGR3A", "ITGAM", "MRC1"],
    "DC_APC": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "LAMP3", "CCR7", "CLEC10A", "FCER1A", "CXCL9", "CXCL10", "CXCL11"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "CD8A", "CD4", "IL7R", "NKG7", "GZMB", "PDCD1", "CTLA4", "LAG3", "TIGIT"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "CXCL13", "BANK1"],
    "epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "MSLN", "CEACAM6"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "neural_schwann": ["SOX10", "S100B", "PLP1", "MPZ", "NRXN1"],
    "acinar_normal": ["PRSS1", "CPA1", "REG1A"],
}
GENES = sorted({gene for genes in STATE_GENES.values() for gene in genes})

REFERENCE_SCORE_COLUMNS = {
    "myCAF_matrix": ["PanCAF", "myCAF", "Tuveson_mCAF", "Pan_Stellate", "Activated_Stellate", "FIBROBLASTS"],
    "iCAF_inflammatory": ["iCAF", "Tuveson_iCAF", "Neuzillet_CAFb"],
    "SPP1_TAM": ["Macrophage", "M2", "Monocytes_1", "Monocytes_2"],
    "DC_APC": ["APC", "cDC1", "cDC2", "DC_activated", "pDC"],
    "T_NK": ["CD8_Tcells", "CD4_Tcells", "NK", "CD8_exhausted", "CD4_memory_activated"],
    "B_plasma": ["B_cell", "Plasma", "Bcell_naive", "Bcell_memory"],
    "epithelial_tumor": ["MALIGNANT CELLS", "DUCTAL", "ductal14", "ductal2", "ductal3", "ductal4"],
    "endothelial": ["ENDOTHELIAL"],
    "neural_schwann": ["PANCREATIC_SCHWANN_CELLS"],
    "acinar_normal": ["ACINAR"],
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
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "external_paired_st_anchor": "GSE235315",
}


def load_stage27():
    spec = importlib.util.spec_from_file_location("stage27_targeted", SCRIPT_27)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {SCRIPT_27}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decode(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def configure_outputs(reference: str) -> None:
    global H5AD_FILES
    global REFERENCE_LABEL
    global OUT_SIGNATURE
    global OUT_PER_SAMPLE
    global OUT_CONTEXT
    global OUT_CORR
    global OUT_REPORT
    global OUT_FIG
    global STATUS

    if reference == "full":
        H5AD_FILES = FULL_H5AD_FILES
        REFERENCE_LABEL = "GSE202051 total h5ad reference"
        OUT_SIGNATURE = PROJECT / "results/tables/gse202051_full_reference_projection_signature_matrix.csv"
        OUT_PER_SAMPLE = PROJECT / "results/tables/gap2_full_reference_projection_deconvolution_per_sample.csv"
        OUT_CONTEXT = PROJECT / "results/tables/gap2_full_reference_projection_deconvolution_context_summary.csv"
        OUT_CORR = PROJECT / "results/tables/gap2_full_reference_projection_deconvolution_correlations.csv"
        OUT_REPORT = PROJECT / "results/reports/gap2_full_reference_projection_deconvolution_report.md"
        OUT_FIG = PROJECT / "results/figures/submission/extended_data_gap2_full_reference_projection_deconvolution"
        STATUS = PROJECT / "results/logs/stage_36_full_reference_projection_deconvolution_status.json"
    else:
        H5AD_FILES = [
            PROJECT / "data/external/GSE202051/GSE202051_adata_010nuc_10x.h5ad",
            PROJECT / "data/external/GSE202051/GSE202051_adata_010orgCRT_10x.h5ad",
        ]
        REFERENCE_LABEL = "GSE202051 h5ad subsets"
        OUT_SIGNATURE = PROJECT / "results/tables/gse202051_reference_projection_signature_matrix.csv"
        OUT_PER_SAMPLE = PROJECT / "results/tables/gap2_reference_projection_deconvolution_per_sample.csv"
        OUT_CONTEXT = PROJECT / "results/tables/gap2_reference_projection_deconvolution_context_summary.csv"
        OUT_CORR = PROJECT / "results/tables/gap2_reference_projection_deconvolution_correlations.csv"
        OUT_REPORT = PROJECT / "results/reports/gap2_reference_projection_deconvolution_report.md"
        OUT_FIG = PROJECT / "results/figures/submission/extended_data_gap2_reference_projection_deconvolution"
        STATUS = PROJECT / "results/logs/stage_36_reference_projection_deconvolution_status.json"


def read_h5ad_var_names(handle: h5py.File) -> list[str]:
    var_obj = handle["var"]
    if isinstance(var_obj, h5py.Dataset):
        var = var_obj[:]
        key = "index" if "index" in var.dtype.names else "_index"
        return decode(var[key])
    if "_index" in var_obj:
        return decode(var_obj["_index"][:])
    if "index" in var_obj:
        return decode(var_obj["index"][:])
    raise RuntimeError("Could not find h5ad var index.")


def read_h5ad_obs_columns(handle: h5py.File, columns: set[str]) -> pd.DataFrame:
    obs_obj = handle["obs"]
    if isinstance(obs_obj, h5py.Dataset):
        obs = pd.DataFrame.from_records(obs_obj[:])
        for col in obs.columns:
            if obs[col].dtype == object:
                obs[col] = obs[col].map(lambda value: value.decode("utf-8") if isinstance(value, bytes) else value)
        return obs

    data: dict[str, object] = {}
    if "_index" in obs_obj:
        data["_index"] = decode(obs_obj["_index"][:])
    elif "index" in obs_obj:
        data["index"] = decode(obs_obj["index"][:])
    for col in sorted(columns):
        if col in obs_obj:
            values = obs_obj[col][:]
            if values.dtype == object:
                data[col] = decode(values)
            else:
                data[col] = values
    return pd.DataFrame(data)


def read_h5ad_gene_matrix(path: Path, genes: list[str]) -> tuple[pd.DataFrame, np.ndarray, list[str]]:
    with h5py.File(path, "r") as handle:
        gene_names = read_h5ad_var_names(handle)
        gene_to_idx = {gene.upper(): idx for idx, gene in enumerate(gene_names)}
        indices = [gene_to_idx[gene.upper()] for gene in genes if gene.upper() in gene_to_idx]
        present = [gene for gene in genes if gene.upper() in gene_to_idx]
        x_group = handle["raw.X"] if "raw.X" in handle else handle["X"]
        data = x_group["data"][()]
        idx = x_group["indices"][()]
        indptr = x_group["indptr"][()]
        n_obs = len(indptr) - 1
        n_vars = len(gene_names)
        x = sparse.csr_matrix((data, idx, indptr), shape=(n_obs, n_vars), dtype=float)[:, indices].tocsr()
        wanted_obs = {col for cols in REFERENCE_SCORE_COLUMNS.values() for col in cols}
        obs = read_h5ad_obs_columns(handle, wanted_obs)
    # If the matrix looks count-like, library-normalize the selected genes.
    if x.data.size and np.nanmax(x.data) > 50:
        totals = np.asarray(x.sum(axis=1)).ravel()
        scale = np.divide(10000.0, totals, out=np.zeros_like(totals, dtype=float), where=totals > 0)
        x = x.multiply(scale[:, None]).tocsr()
        x.data = np.log1p(x.data)
    return obs, pd.DataFrame.sparse.from_spmatrix(x, columns=present).sparse.to_dense().to_numpy(float), present


def build_reference_signature() -> tuple[pd.DataFrame, list[str], list[str]]:
    signatures = []
    selection_rows = []
    matrices = []
    obs_frames = []
    present_genes: list[str] | None = None
    for path in H5AD_FILES:
        obs, x, present = read_h5ad_gene_matrix(path, GENES)
        obs["reference_file"] = path.name
        if present_genes is None:
            present_genes = present
        elif present != present_genes:
            common = [gene for gene in present_genes if gene in set(present)]
            left_idx = [present_genes.index(gene) for gene in common]
            right_idx = [present.index(gene) for gene in common]
            matrices = [mat[:, left_idx] for mat in matrices]
            x = x[:, right_idx]
            present_genes = common
        matrices.append(x)
        obs_frames.append(obs)
    obs_all = pd.concat(obs_frames, ignore_index=True)
    x_all = np.vstack(matrices)
    if present_genes is None:
        raise RuntimeError("No reference genes were found in GSE202051 h5ad files.")

    for state, columns in REFERENCE_SCORE_COLUMNS.items():
        usable = [col for col in columns if col in obs_all.columns]
        if not usable:
            continue
        scores = obs_all[usable].astype(float).max(axis=1).to_numpy()
        order = np.argsort(scores)[::-1]
        positive = np.where(scores > 0)[0]
        n_select = min(250, max(25, int(0.06 * len(scores))))
        if len(positive) >= 25:
            selected = positive[np.argsort(scores[positive])[::-1][: min(n_select, len(positive))]]
        else:
            selected = order[:n_select]
        profile = np.nan_to_num(np.nanmean(x_all[selected, :], axis=0), nan=0.0, posinf=0.0, neginf=0.0)
        signatures.append(profile)
        selection_rows.append(
            {
                "cell_state": state,
                "n_reference_cells_selected": int(len(selected)),
                "median_reference_score": float(np.nanmedian(scores[selected])),
                "score_columns_used": ";".join(usable),
            }
        )

    states = [row["cell_state"] for row in selection_rows]
    sig = pd.DataFrame(np.vstack(signatures), index=states, columns=present_genes)
    for state in sig.index:
        allowed = {gene for gene in STATE_GENES.get(state, []) if gene in sig.columns}
        blocked = [gene for gene in sig.columns if gene not in allowed]
        sig.loc[state, blocked] = 0.0
    # Keep genes with nonzero reference signal in at least one state.
    sig = sig.replace([np.inf, -np.inf], np.nan).fillna(0.0)
    nonzero = sig.max(axis=0) > 0
    sig = sig.loc[:, nonzero]
    selection = pd.DataFrame(selection_rows)
    long = sig.reset_index(names="cell_state").melt(id_vars="cell_state", var_name="gene", value_name="reference_expression")
    long = long.merge(selection, on="cell_state", how="left")
    OUT_SIGNATURE.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_SIGNATURE, index=False)
    return sig, list(sig.columns), list(sig.index)


def load_spatial_expression(genes: list[str]) -> pd.DataFrame:
    stage27 = load_stage27()
    stage27.TARGET_GENES = sorted(set(genes))
    manifest = pd.read_csv(PROJECT / "metadata/dataset_manifest_curated.csv")
    mvp = stage27.load_mvp_spots()
    anchor = stage27.load_gse235315_spots()
    combined = pd.concat([mvp, anchor], ignore_index=True)
    h5_expr = stage27.extract_h5_gene_values(manifest, combined)
    rds_expr = stage27.extract_gse272362_gene_values()
    frames = [frame for frame in [h5_expr, rds_expr] if not frame.empty]
    expr = pd.concat(frames, ignore_index=True)
    if "cohort_context" not in expr.columns:
        expr["cohort_context"] = expr.apply(stage27.context_for_row, axis=1)
    else:
        missing = expr["cohort_context"].isna() | expr["cohort_context"].astype(str).eq("")
        expr.loc[missing, "cohort_context"] = expr.loc[missing].apply(stage27.context_for_row, axis=1)
    expr = stage27.add_regions(expr)
    return expr


def project_sample(sample: pd.DataFrame, sig: pd.DataFrame, genes: list[str], states: list[str]) -> pd.DataFrame:
    available = [gene for gene in genes if gene in sample.columns and sample[gene].notna().any()]
    if len(available) < 15:
        return pd.DataFrame()
    a = np.nan_to_num(sig[available].T.to_numpy(float), nan=0.0, posinf=0.0, neginf=0.0)
    gene_scale = np.nanmax(a, axis=1)
    gene_scale[~np.isfinite(gene_scale) | (gene_scale <= 0)] = 1.0
    a = a / gene_scale[:, None]
    # Unit-normalize state columns so high-expression states do not dominate the fit.
    col_norm = np.linalg.norm(a, axis=0)
    col_norm[~np.isfinite(col_norm) | (col_norm <= 0)] = 1.0
    a = a / col_norm[None, :]
    y = np.nan_to_num(sample[available].fillna(0).to_numpy(float).T, nan=0.0, posinf=0.0, neginf=0.0) / gene_scale[:, None]
    coef = np.linalg.lstsq(a, y, rcond=None)[0].T
    coef = np.clip(coef / col_norm[None, :], 0, None)
    denom = coef.sum(axis=1)
    frac = np.divide(coef, denom[:, None], out=np.zeros_like(coef), where=denom[:, None] > 0)
    out = sample[["dataset_id", "sample_id", "specimen_type", "cohort_context", "barcode", "score_caf_myeloid_barrier", "is_caf_core_top10", "is_interface"]].copy()
    for i, state in enumerate(states):
        out[f"proj_{state}"] = frac[:, i]
    return out


def project_spatial(expr: pd.DataFrame, sig: pd.DataFrame, genes: list[str], states: list[str]) -> pd.DataFrame:
    frames = []
    for _, sample in expr.groupby("sample_id", sort=False):
        projected = project_sample(sample, sig, genes, states)
        if not projected.empty:
            frames.append(projected)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize_projection(projected: pd.DataFrame, states: list[str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for sample_id, sample in projected.groupby("sample_id", sort=False):
        core = sample["is_caf_core_top10"].astype(bool)
        if core.sum() < 3:
            continue
        base = {
            "dataset_id": sample["dataset_id"].iloc[0],
            "sample_id": sample_id,
            "specimen_type": sample["specimen_type"].iloc[0],
            "cohort_context": sample["cohort_context"].iloc[0],
            "n_spots": int(len(sample)),
            "n_caf_core_spots": int(core.sum()),
        }
        for state in states:
            col = f"proj_{state}"
            values = sample[col].astype(float)
            rho, p = spearmanr(sample["score_caf_myeloid_barrier"], values, nan_policy="omit")
            rows.append(
                {
                    **base,
                    "cell_state": state,
                    "core_enrichment": float(values[core].mean() - values[~core].mean()),
                    "core_mean_projection": float(values[core].mean()),
                    "noncore_mean_projection": float(values[~core].mean()),
                    "spot_spearman_with_caf_myeloid": float(rho) if np.isfinite(rho) else np.nan,
                    "spot_spearman_p": float(p) if np.isfinite(p) else np.nan,
                }
            )
    per_sample = pd.DataFrame(rows)

    context_rows = []
    for (context, state), sub in per_sample.groupby(["cohort_context", "cell_state"], sort=False):
        vals = sub["core_enrichment"].dropna()
        rhos = sub["spot_spearman_with_caf_myeloid"].dropna()
        context_rows.append(
            {
                "cohort_context": context,
                "cell_state": state,
                "n_samples": int(len(vals)),
                "median_core_enrichment": float(vals.median()) if len(vals) else np.nan,
                "n_core_positive": int((vals > 0).sum()),
                "median_spot_spearman_with_caf_myeloid": float(rhos.median()) if len(rhos) else np.nan,
                "n_spearman_positive": int((rhos > 0).sum()),
            }
        )
    context = pd.DataFrame(context_rows)

    dec_path = PROJECT / "results/tables/mechanism_candidate_axis_sample_summary.csv"
    corr_rows = []
    if dec_path.exists():
        dec = pd.read_csv(dec_path)[["sample_id", "immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]]
        merged = per_sample.merge(dec, on="sample_id", how="left")
        tumor = merged[~merged["cohort_context"].eq("normal_pancreas")].copy()
        for state, sub in tumor.groupby("cell_state"):
            for target in ["immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]:
                tmp = sub[["core_enrichment", target]].dropna()
                if len(tmp) < 5:
                    continue
                rho, p = spearmanr(tmp["core_enrichment"], tmp[target])
                corr_rows.append(
                    {
                        "cell_state": state,
                        "metric": "core_enrichment",
                        "target": target,
                        "n_samples": int(len(tmp)),
                        "spearman_rho": float(rho),
                        "p_value": float(p),
                    }
                )
    corr = pd.DataFrame(corr_rows)
    return per_sample, context, corr


def plot_results(context: pd.DataFrame, per_sample: pd.DataFrame, corr: pd.DataFrame, states: list[str]) -> None:
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    contexts = [c for c in CONTEXT_ORDER if c in set(context["cohort_context"])]
    mat = context.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median").reindex(index=states, columns=contexts)
    rho_mat = context.pivot_table(index="cell_state", columns="cohort_context", values="median_spot_spearman_with_caf_myeloid", aggfunc="median").reindex(index=states, columns=contexts)

    fig = plt.figure(figsize=(13.5, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.3, 1.0, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[0, 2])
    fig.suptitle("Gap 2: GSE202051 reference-projection deconvolution prototype", fontsize=15, fontweight="bold")

    im1 = ax1.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.08, vmax=0.08, aspect="auto")
    ax1.set_title("A  CAF-core projection enrichment", loc="left", fontsize=10.5, fontweight="bold")
    ax1.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax1.set_yticks(np.arange(len(states)), states)
    ax1.tick_params(labelsize=7)
    for i, state in enumerate(states):
        for j, ctx in enumerate(contexts):
            val = mat.loc[state, ctx]
            if np.isfinite(val):
                row = context[context["cell_state"].eq(state) & context["cohort_context"].eq(ctx)].iloc[0]
                ax1.text(j, i, f"{val:.2f}\n{int(row.n_core_positive)}/{int(row.n_samples)}", ha="center", va="center", fontsize=5.0)
    cb1 = fig.colorbar(im1, ax=ax1, fraction=0.045, pad=0.02)
    cb1.set_label("median projection enrichment", fontsize=8)
    cb1.ax.tick_params(labelsize=7)

    im2 = ax2.imshow(rho_mat.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax2.set_title("B  Projection vs CAF-myeloid score", loc="left", fontsize=10.5, fontweight="bold")
    ax2.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax2.set_yticks(np.arange(len(states)), states)
    ax2.tick_params(labelsize=7)
    for i, state in enumerate(states):
        for j, ctx in enumerate(contexts):
            val = rho_mat.loc[state, ctx]
            if np.isfinite(val):
                ax2.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.2)
    cb2 = fig.colorbar(im2, ax=ax2, fraction=0.045, pad=0.02)
    cb2.set_label("median Spearman rho", fontsize=8)
    cb2.ax.tick_params(labelsize=7)

    plot_states = ["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "epithelial_tumor"]
    x = np.arange(len(plot_states))
    width = 0.22
    plot_contexts = [c for c in ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"] if c in contexts]
    colors = {"primary_tumor": "#4E79A7", "liver_metastasis": "#59A14F", "lymph_node_metastasis": "#9C755F"}
    for k, ctx in enumerate(plot_contexts):
        vals = [
            per_sample.loc[per_sample["cohort_context"].eq(ctx) & per_sample["cell_state"].eq(state), "core_enrichment"].median()
            for state in plot_states
        ]
        ax3.bar(x + (k - 1) * width, vals, width=width, color=colors.get(ctx, "#999999"), label=CONTEXT_LABELS.get(ctx, ctx))
    ax3.axhline(0, color="#333333", lw=0.8)
    ax3.set_title("C  Metastatic-site contrast", loc="left", fontsize=10.5, fontweight="bold")
    ax3.set_xticks(x, plot_states, rotation=35, ha="right")
    ax3.set_ylabel("median CAF-core projection enrichment", fontsize=8)
    ax3.tick_params(labelsize=7)
    ax3.legend(frameon=False, fontsize=7)

    for ext in [".pdf", ".png", ".svg"]:
        fig.savefig(f"{OUT_FIG}{ext}", dpi=320 if ext == ".png" else None)
    plt.close(fig)


def write_report(per_sample: pd.DataFrame, context: pd.DataFrame, corr: pd.DataFrame, states: list[str]) -> None:
    key = context[context["cell_state"].isin(["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "epithelial_tumor"])].copy()
    key["context_label"] = key["cohort_context"].map(CONTEXT_LABELS).fillna(key["cohort_context"])
    lines = [
        "# Gap 2 Reference-Projection Deconvolution Report",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Purpose",
        "",
        f"This analysis uses the downloaded {REFERENCE_LABEL} as a reference-projection prototype. It estimates broad cell-state projections in spatial spots using curated marker genes and a clipped least-squares projection against reference signatures.",
        "",
        "## Main Context Summary",
        "",
        "| context | cell state | n samples | median CAF-core enrichment | positive samples |",
        "|---|---|---:|---:|---:|",
    ]
    for row in key.sort_values(["cell_state", "cohort_context"]).itertuples(index=False):
        lines.append(
            f"| {row.context_label} | {row.cell_state} | {int(row.n_samples)} | {row.median_core_enrichment:.4f} | {int(row.n_core_positive)}/{int(row.n_samples)} |"
        )

    if not corr.empty:
        lines.extend(["", "## Immune-Decoupling Correlations", "", "| cell state | target | n | rho | p |", "|---|---|---:|---:|---:|"])
        show = corr[corr["target"].eq("immune_decoupling_index")].sort_values("spearman_rho")
        for row in show.itertuples(index=False):
            lines.append(f"| {row.cell_state} | {row.target} | {int(row.n_samples)} | {row.spearman_rho:.3f} | {row.p_value:.3g} |")

    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "This is a reference-projection prototype, not a final validated deconvolution model. The projection uses selected marker genes shared with the spatial expression tables and clipped least-squares coefficients normalized per spot. It is useful for checking whether the marker-level Gap 2 conclusion is directionally supported by an external PDAC snRNA reference.",
            "",
            "## Generated Outputs",
            "",
            f"- `{OUT_SIGNATURE.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_PER_SAMPLE.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_CONTEXT.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_CORR.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_FIG.relative_to(PROJECT).as_posix()}.pdf`",
        ]
    )
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description="Run Gap 2 GSE202051 marker-constrained reference projection.")
    parser.add_argument("--reference", choices=["small", "full"], default="small", help="Use small h5ad subsets or the total GSE202051 h5ad.")
    args = parser.parse_args()
    configure_outputs(args.reference)

    sig, genes, states = build_reference_signature()
    expr = load_spatial_expression(genes)
    projected = project_spatial(expr, sig, genes, states)
    per_sample, context, corr = summarize_projection(projected, states)

    OUT_PER_SAMPLE.parent.mkdir(parents=True, exist_ok=True)
    per_sample.to_csv(OUT_PER_SAMPLE, index=False)
    context.to_csv(OUT_CONTEXT, index=False)
    corr.to_csv(OUT_CORR, index=False)
    plot_results(context, per_sample, corr, states)
    write_report(per_sample, context, corr, states)
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "36_reference_projection_deconvolution",
                "reference": args.reference,
                "status": "success",
                "n_samples": int(per_sample["sample_id"].nunique()) if not per_sample.empty else 0,
                "n_states": len(states),
                "n_genes": len(genes),
                "outputs": [
                    str(OUT_SIGNATURE),
                    str(OUT_PER_SAMPLE),
                    str(OUT_CONTEXT),
                    str(OUT_CORR),
                    f"{OUT_FIG}.pdf",
                ],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote reference projection deconvolution prototype for {per_sample['sample_id'].nunique()} samples")


if __name__ == "__main__":
    main()
