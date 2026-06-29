from __future__ import annotations

import gzip
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
from scipy.stats import mannwhitneyu, spearmanr


PROJECT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT / "data/external/TCGA_PAAD"
EXPR_PATH = DATA_DIR / "HiSeqV2.gz"
CLIN_PATH = DATA_DIR / "PAAD_clinicalMatrix.tsv"
SURV_PATH = DATA_DIR / "TCGA-PAAD.survival.tsv.gz"

OUT_SCORES = PROJECT / "results/tables/tcga_paad_bulk_context_scores.csv"
OUT_CORR = PROJECT / "results/tables/tcga_paad_bulk_context_axis_correlations.csv"
OUT_CLIN = PROJECT / "results/tables/tcga_paad_bulk_context_clinical_exploratory.csv"
OUT_COVERAGE = PROJECT / "results/tables/tcga_paad_bulk_context_gene_coverage.csv"
OUT_REPORT = PROJECT / "results/reports/tcga_paad_bulk_context_report.md"
OUT_FIG = PROJECT / "results/figures/submission/extended_data_tcga_paad_bulk_context"
STATUS = PROJECT / "results/logs/stage_38_tcga_paad_bulk_context_status.json"

GENE_SETS = {
    "myCAF_matrix": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "POSTN", "PDGFRB", "SPARC", "FN1"],
    "iCAF_inflammatory": ["IL6", "CXCL12", "CXCL14", "CFD", "DPT", "HAS1", "PDGFRA", "LIF", "CTGF"],
    "SPP1_TAM": ["SPP1", "TREM2", "APOE", "LGALS3", "CD68", "CD163", "C1QA", "C1QB", "C1QC", "LST1", "FCGR3A", "ITGAM", "MRC1"],
    "DC_APC": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "LAMP3", "CCR7", "CLEC10A", "FCER1A", "CXCL9", "CXCL10", "CXCL11"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "CD8A", "CD4", "IL7R", "NKG7", "GZMB", "PDCD1", "CTLA4", "LAG3", "TIGIT"],
    "B_plasma": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "CXCL13", "BANK1"],
    "epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "MSLN", "CEACAM6"],
    "TGFb_EMT": ["TGFB1", "TGFBI", "CTGF", "INHBA", "VIM", "ITGA5", "SERPINE1", "MMP14"],
    "matrix_integrin": ["COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3", "FN1", "SPARC", "POSTN", "ITGA5", "ITGAV", "ITGB1", "ITGB5", "CD44"],
}

PLOT_AXES = ["myCAF_matrix", "SPP1_TAM", "TGFb_EMT", "matrix_integrin", "DC_APC", "T_NK", "B_plasma", "epithelial_tumor"]


def zscore_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    values = matrix.to_numpy(float)
    mean = np.nanmean(values, axis=1, keepdims=True)
    sd = np.nanstd(values, axis=1, keepdims=True)
    sd[~np.isfinite(sd) | (sd == 0)] = 1.0
    return pd.DataFrame((values - mean) / sd, index=matrix.index, columns=matrix.columns)


def load_expression() -> pd.DataFrame:
    expr = pd.read_csv(EXPR_PATH, sep="\t", compression="gzip")
    gene_col = expr.columns[0]
    expr[gene_col] = expr[gene_col].astype(str).str.upper()
    expr = expr.groupby(gene_col, as_index=True).mean(numeric_only=True)
    # Retain primary tumor aliquots; legacy Xena sample ids are mostly TCGA-..-01.
    expr = expr[[col for col in expr.columns if "-01" in col]]
    return expr


def load_clinical() -> tuple[pd.DataFrame, pd.DataFrame]:
    clin = pd.read_csv(CLIN_PATH, sep="\t", dtype=str)
    with gzip.open(SURV_PATH, "rt", encoding="utf-8", errors="replace") as handle:
        surv = pd.read_csv(handle, sep="\t", dtype=str)
    return clin, surv


def score_axes(expr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_genes = sorted({gene for genes in GENE_SETS.values() for gene in genes})
    present = [gene for gene in all_genes if gene in expr.index]
    z = zscore_rows(expr.loc[present])
    rows = []
    scores = pd.DataFrame(index=expr.columns)
    for axis, genes in GENE_SETS.items():
        available = [gene for gene in genes if gene in z.index]
        scores[axis] = z.loc[available].mean(axis=0) if available else np.nan
        rows.append(
            {
                "axis": axis,
                "n_genes_defined": len(genes),
                "n_genes_present": len(available),
                "present_genes": ";".join(available),
                "missing_genes": ";".join([gene for gene in genes if gene not in z.index]),
            }
        )
    scores.index.name = "sample_id"
    scores = scores.reset_index()
    scores["patient_id"] = scores["sample_id"].str.slice(0, 12)
    scores["bulk_stromal_myeloid_index"] = scores[["myCAF_matrix", "SPP1_TAM", "TGFb_EMT", "matrix_integrin"]].mean(axis=1)
    scores["bulk_immune_projection_index"] = scores[["DC_APC", "T_NK", "B_plasma"]].mean(axis=1)
    scores["bulk_decoupling_like_index"] = scores["bulk_stromal_myeloid_index"] - scores["bulk_immune_projection_index"]
    return scores, pd.DataFrame(rows)


def correlate_axes(scores: pd.DataFrame) -> pd.DataFrame:
    axes = PLOT_AXES + ["bulk_stromal_myeloid_index", "bulk_immune_projection_index", "bulk_decoupling_like_index"]
    rows = []
    for a in axes:
        for b in axes:
            if a == b:
                tmp = scores[[a]].dropna()
                rho, p = (1.0, 0.0) if len(tmp) >= 1 else (np.nan, np.nan)
            else:
                tmp = scores[[a, b]].dropna()
                rho, p = spearmanr(tmp[a], tmp[b]) if len(tmp) >= 5 else (np.nan, np.nan)
            rows.append({"axis_x": a, "axis_y": b, "n_samples": int(len(tmp)), "spearman_rho": float(rho), "p_value": float(p)})
    return pd.DataFrame(rows)


def clinical_exploratory(scores: pd.DataFrame, clin: pd.DataFrame, surv: pd.DataFrame) -> pd.DataFrame:
    clinical_cols = [
        "sampleID",
        "_PATIENT",
        "pathologic_stage",
        "neoplasm_histologic_grade",
        "histological_type",
        "residual_tumor",
        "primary_therapy_outcome_success",
        "vital_status",
    ]
    merged = scores.merge(clin[[c for c in clinical_cols if c in clin.columns]], left_on="sample_id", right_on="sampleID", how="left")
    surv = surv.rename(columns={"sample": "survival_sample_id"})
    surv["patient_id"] = surv["survival_sample_id"].str.slice(0, 12)
    merged = merged.merge(surv[["patient_id", "OS.time", "OS"]], on="patient_id", how="left")

    rows = []
    for axis in ["bulk_decoupling_like_index", "bulk_stromal_myeloid_index", "bulk_immune_projection_index", "SPP1_TAM", "myCAF_matrix", "DC_APC", "T_NK"]:
        for cat_col in ["pathologic_stage", "neoplasm_histologic_grade", "residual_tumor", "primary_therapy_outcome_success"]:
            if cat_col not in merged.columns:
                continue
            sub = merged[[axis, cat_col]].dropna()
            counts = sub[cat_col].value_counts()
            if len(counts) < 2:
                continue
            groups = counts.index[:2].tolist()
            a = sub.loc[sub[cat_col].eq(groups[0]), axis]
            b = sub.loc[sub[cat_col].eq(groups[1]), axis]
            if len(a) < 5 or len(b) < 5:
                continue
            stat = mannwhitneyu(a, b, alternative="two-sided")
            rows.append(
                {
                    "axis": axis,
                    "clinical_variable": cat_col,
                    "group_1": groups[0],
                    "n_group_1": int(len(a)),
                    "median_group_1": float(a.median()),
                    "group_2": groups[1],
                    "n_group_2": int(len(b)),
                    "median_group_2": float(b.median()),
                    "mannwhitney_p": float(stat.pvalue),
                }
            )

        tmp = merged[[axis, "OS.time", "OS"]].dropna()
        if len(tmp) >= 30:
            tmp["OS.time"] = pd.to_numeric(tmp["OS.time"], errors="coerce")
            tmp["OS"] = pd.to_numeric(tmp["OS"], errors="coerce")
            tmp = tmp.dropna()
            if len(tmp) >= 30:
                rho, p = spearmanr(tmp[axis], tmp["OS.time"])
                rows.append(
                    {
                        "axis": axis,
                        "clinical_variable": "OS.time_exploratory_spearman",
                        "group_1": "all",
                        "n_group_1": int(len(tmp)),
                        "median_group_1": float(tmp[axis].median()),
                        "group_2": "",
                        "n_group_2": int(tmp["OS"].sum()),
                        "median_group_2": float(tmp["OS.time"].median()),
                        "mannwhitney_p": float(p),
                        "spearman_rho": float(rho),
                    }
                )
    return pd.DataFrame(rows), merged


def plot_context(scores: pd.DataFrame, corr: pd.DataFrame) -> None:
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    axes = PLOT_AXES
    mat = corr[corr["axis_x"].isin(axes) & corr["axis_y"].isin(axes)].pivot(index="axis_y", columns="axis_x", values="spearman_rho").reindex(index=axes, columns=axes)

    fig = plt.figure(figsize=(13.4, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.18, 1.0], height_ratios=[1.0, 0.9])
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    fig.suptitle("External bulk context in TCGA PAAD", fontsize=15, fontweight="bold")

    im = ax1.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1)
    ax1.set_title("A  Bulk gene-set correlation structure", loc="left", fontsize=10.5, fontweight="bold")
    ax1.set_xticks(np.arange(len(axes)), axes, rotation=35, ha="right")
    ax1.set_yticks(np.arange(len(axes)), axes)
    ax1.tick_params(labelsize=8)
    for i, y in enumerate(axes):
        for j, x in enumerate(axes):
            val = mat.loc[y, x]
            if np.isfinite(val):
                ax1.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6)
    cb = fig.colorbar(im, ax=ax1, fraction=0.045, pad=0.02)
    cb.set_label("Spearman rho", fontsize=8)
    cb.ax.tick_params(labelsize=7)

    sc = ax2.scatter(scores["bulk_stromal_myeloid_index"], scores["bulk_immune_projection_index"], c=scores["bulk_decoupling_like_index"], cmap="RdBu_r", s=36, edgecolor="white", linewidth=0.4)
    rho, p = spearmanr(scores["bulk_stromal_myeloid_index"], scores["bulk_immune_projection_index"])
    ax2.set_title("B  Stromal-myeloid vs immune bulk axes", loc="left", fontsize=10.5, fontweight="bold")
    ax2.set_xlabel("bulk stromal-myeloid index")
    ax2.set_ylabel("bulk immune projection index")
    ax2.text(0.02, 0.95, f"rho = {rho:.2f}", transform=ax2.transAxes, ha="left", va="top", fontsize=8)
    ax2.tick_params(labelsize=8)
    cb2 = fig.colorbar(sc, ax=ax2, fraction=0.045, pad=0.02)
    cb2.set_label("bulk decoupling-like index", fontsize=8)
    cb2.ax.tick_params(labelsize=7)

    ordered = scores.sort_values("bulk_decoupling_like_index").reset_index(drop=True)
    x = np.arange(len(ordered))
    for axis, color in [("myCAF_matrix", "#B07AA1"), ("SPP1_TAM", "#9C755F"), ("DC_APC", "#4E79A7"), ("T_NK", "#59A14F")]:
        rolling = ordered[axis].rolling(window=15, min_periods=5, center=True).mean()
        ax3.plot(x, rolling, label=axis, lw=1.8, color=color)
    ax3.set_title("C  Samples ranked by bulk decoupling-like index", loc="left", fontsize=10.5, fontweight="bold")
    ax3.set_xlabel("TCGA PAAD samples, low to high decoupling-like index")
    ax3.set_ylabel("rolling mean score")
    ax3.legend(frameon=False, fontsize=7, ncol=2)
    ax3.tick_params(labelsize=8)

    for ext in [".pdf", ".png", ".svg"]:
        fig.savefig(f"{OUT_FIG}{ext}", dpi=320 if ext == ".png" else None)
    plt.close(fig)


def write_report(scores: pd.DataFrame, coverage: pd.DataFrame, corr: pd.DataFrame, clinical: pd.DataFrame) -> None:
    pairs = [
        ("myCAF_matrix", "SPP1_TAM"),
        ("myCAF_matrix", "DC_APC"),
        ("SPP1_TAM", "DC_APC"),
        ("bulk_decoupling_like_index", "DC_APC"),
        ("bulk_decoupling_like_index", "T_NK"),
        ("bulk_decoupling_like_index", "B_plasma"),
        ("bulk_decoupling_like_index", "myCAF_matrix"),
    ]
    corr_lookup = corr.set_index(["axis_x", "axis_y"])
    lines = [
        "# TCGA PAAD Bulk Context Report",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Purpose",
        "",
        "This analysis places the spatially nominated CAF/TAM, matrix, immune and epithelial axes into an external TCGA PAAD bulk RNA-seq context. It is used as orthogonal biological context, not as spatial or clinical validation.",
        "",
        "## Dataset",
        "",
        f"- Expression matrix: `{EXPR_PATH.relative_to(PROJECT).as_posix()}`",
        f"- Samples scored: {scores.shape[0]} primary tumor RNA-seq samples.",
        "",
        "## Selected Bulk Correlations",
        "",
        "| axis 1 | axis 2 | rho | p |",
        "|---|---|---:|---:|",
    ]
    for a, b in pairs:
        row = corr_lookup.loc[(a, b)]
        lines.append(f"| {a} | {b} | {row.spearman_rho:.3f} | {row.p_value:.3g} |")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "TCGA PAAD supports a broad stromal-myeloid/matrix expression continuum that can be compared with immune and epithelial bulk axes. Because bulk RNA-seq lacks spatial information, these results should be framed as external context for the nominated axes rather than validation of CAF-core localization.",
            "",
            "## Generated Outputs",
            "",
            f"- `{OUT_SCORES.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_COVERAGE.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_CORR.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_CLIN.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_FIG.relative_to(PROJECT).as_posix()}.pdf`",
        ]
    )
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    expr = load_expression()
    clin, surv = load_clinical()
    scores, coverage = score_axes(expr)
    corr = correlate_axes(scores)
    clinical, merged = clinical_exploratory(scores, clin, surv)

    OUT_SCORES.parent.mkdir(parents=True, exist_ok=True)
    scores.to_csv(OUT_SCORES, index=False)
    corr.to_csv(OUT_CORR, index=False)
    clinical.to_csv(OUT_CLIN, index=False)
    coverage.to_csv(OUT_COVERAGE, index=False)
    plot_context(scores, corr)
    write_report(scores, coverage, corr, clinical)

    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "38_tcga_paad_bulk_context",
                "status": "success",
                "n_samples": int(scores.shape[0]),
                "outputs": [str(OUT_SCORES), str(OUT_CORR), str(OUT_CLIN), f"{OUT_FIG}.pdf"],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote TCGA PAAD bulk context for {scores.shape[0]} samples")


if __name__ == "__main__":
    main()
