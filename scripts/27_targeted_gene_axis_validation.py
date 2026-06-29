from __future__ import annotations

import json
import math
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
from scipy.spatial import cKDTree


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_TABLE = PROJECT_ROOT / "results/tables/targeted_gene_axis_validation_summary.csv"
OUT_LONG = PROJECT_ROOT / "results/tables/targeted_gene_axis_validation_per_sample.csv"
OUT_FIG = PROJECT_ROOT / "results/figures/submission/figure3_supplement_targeted_gene_axis_validation"
STATUS = PROJECT_ROOT / "results/logs/stage_27_targeted_gene_axis_validation_status.json"

GENE_SETS = {
    "SPP1-TAM/matrix": ["SPP1", "TREM2", "APOE", "LGALS3", "COL1A1", "COL1A2", "FN1", "SPARC"],
    "TGF-beta/EMT invasive": ["TGFB1", "TGFBI", "CTGF", "INHBA", "VIM", "ITGA5", "SERPINE1", "MMP14"],
    "IFN/APC antigen": ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "CXCL9", "CXCL10", "STAT1", "IRF1"],
    "B/plasma lymphoid": ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "IGHG1", "CXCL13", "BANK1"],
    "T cell/checkpoint": ["CD3D", "CD3E", "CD8A", "TRAC", "PDCD1", "CTLA4", "LAG3", "TIGIT"],
}
TARGET_GENES = sorted({gene for genes in GENE_SETS.values() for gene in genes})

CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "external_paired_st_anchor": "GSE235315",
}


def decode_array(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def read_10x_h5(path: Path) -> tuple[sparse.csc_matrix, list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        data = matrix["data"][()]
        indices = matrix["indices"][()]
        indptr = matrix["indptr"][()]
        shape = tuple(matrix["shape"][()])
        barcodes = decode_array(matrix["barcodes"][()])
        gene_names = decode_array(matrix["features"]["name"][()])
    return sparse.csc_matrix((data, indices, indptr), shape=shape), gene_names, barcodes


def zscore_by_sample(df: pd.DataFrame, value_col: str, group_col: str = "sample_id") -> pd.Series:
    def _z(values: pd.Series) -> pd.Series:
        arr = values.to_numpy(float)
        sd = np.nanstd(arr)
        if not np.isfinite(sd) or sd == 0:
            return pd.Series(np.zeros(len(arr)), index=values.index)
        return pd.Series((arr - np.nanmean(arr)) / sd, index=values.index)

    return df.groupby(group_col, group_keys=False)[value_col].apply(_z)


def context_for_row(row: pd.Series) -> str:
    dataset = str(row.get("dataset_id", ""))
    specimen = str(row.get("specimen_type", ""))
    treatment = str(row.get("treatment_context", ""))
    if dataset == "GSE282302":
        return "post_neoadjuvant_sections"
    if dataset == "GSE274103":
        return "treatment_naive_primary"
    if dataset == "GSE235315":
        return "external_paired_st_anchor"
    if specimen in {"primary_tumor", "liver_metastasis", "lymph_node_metastasis", "normal_pancreas"}:
        return specimen
    if "naive" in treatment.lower():
        return "treatment_naive_primary"
    return dataset or "unknown"


def load_mvp_spots() -> pd.DataFrame:
    usecols = [
        "dataset_id",
        "sample_id",
        "barcode",
        "x_pixel",
        "y_pixel",
        "score_caf_myeloid_barrier",
        "score_tumor_epithelial",
    ]
    path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    if not path.exists():
        path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv"
    optional = ["edge_or_background_risk"]
    header = pd.read_csv(path, nrows=0).columns
    cols = [col for col in usecols + optional if col in header]
    df = pd.read_csv(path, usecols=cols)
    if "edge_or_background_risk" in df.columns:
        df = df[~df["edge_or_background_risk"].fillna(False).astype(bool)].copy()
    return df


def load_gse235315_spots() -> pd.DataFrame:
    path = PROJECT_ROOT / "results/tables/gse235315_spot_level_scores.csv"
    if not path.exists():
        return pd.DataFrame()
    header = pd.read_csv(path, nrows=0).columns
    cols = [
        col
        for col in [
            "dataset_id",
            "sample_id",
            "barcode",
            "x_pixel",
            "y_pixel",
            "score_caf_myeloid_barrier",
            "score_tumor_epithelial",
        ]
        if col in header
    ]
    return pd.read_csv(path, usecols=cols)


def load_gse272362_spots() -> pd.DataFrame:
    path = PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv"
    header = pd.read_csv(path, nrows=0).columns
    cols = [
        col
        for col in [
            "dataset_id",
            "sample_id",
            "specimen_type",
            "barcode",
            "x_pixel",
            "y_pixel",
            "score_caf_myeloid_barrier",
            "score_tumor_epithelial",
        ]
        if col in header
    ]
    return pd.read_csv(path, usecols=cols)


def extract_h5_gene_values(manifest: pd.DataFrame, spot: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    targets = set(spot[["dataset_id", "sample_id"]].drop_duplicates().itertuples(index=False, name=None))
    manifest = manifest[manifest[["dataset_id", "sample_id"]].apply(tuple, axis=1).isin(targets)].copy()
    for _, meta in manifest.iterrows():
        expression_path = Path(str(meta["expression_path"]))
        if not expression_path.exists() or expression_path.suffix.lower() != ".h5":
            continue
        sample_spots = spot[(spot["dataset_id"].eq(meta["dataset_id"])) & (spot["sample_id"].eq(meta["sample_id"]))].copy()
        if sample_spots.empty:
            continue
        counts, genes, barcodes = read_10x_h5(expression_path)
        gene_to_idx: dict[str, list[int]] = {}
        for idx, gene in enumerate(genes):
            gene_to_idx.setdefault(gene.upper(), []).append(idx)
        barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
        sample_spots["matrix_idx"] = sample_spots["barcode"].map(barcode_to_idx)
        sample_spots = sample_spots.dropna(subset=["matrix_idx"]).copy()
        if sample_spots.empty:
            continue
        spot_idx = sample_spots["matrix_idx"].astype(int).to_numpy()
        n_counts = np.asarray(counts[:, spot_idx].sum(axis=0)).ravel()
        scale = np.divide(10000.0, n_counts, out=np.zeros_like(n_counts, dtype=float), where=n_counts > 0)
        out = sample_spots.drop(columns=["matrix_idx"]).copy()
        for gene in TARGET_GENES:
            idx = gene_to_idx.get(gene.upper(), [])
            if not idx:
                out[gene] = np.nan
                continue
            values = np.asarray(counts[idx, :][:, spot_idx].sum(axis=0)).ravel().astype(float)
            out[gene] = np.log1p(values * scale)
        out["specimen_type"] = str(meta.get("specimen_type", ""))
        out["cohort_context"] = out.apply(context_for_row, axis=1)
        rows.append(out)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def extract_gse272362_gene_values() -> pd.DataFrame:
    # The GSE272362 RDS expression matrix has already been reduced to module scores.
    # To keep this script single-language and reproducible for quick iteration, use a
    # companion CSV if it has been exported; otherwise return an empty frame and log it.
    path = PROJECT_ROOT / "results/tables/gse272362_rds_target_gene_expression.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path)


def add_regions(df: pd.DataFrame) -> pd.DataFrame:
    out: list[pd.DataFrame] = []
    required = {"x_pixel", "y_pixel", "score_caf_myeloid_barrier", "score_tumor_epithelial"}
    for sample_id, sample in df.groupby("sample_id", sort=False):
        sample = sample.copy()
        if not required.issubset(sample.columns) or len(sample) < 30:
            continue
        coords = sample[["x_pixel", "y_pixel"]].to_numpy(float)
        finite = np.isfinite(coords).all(axis=1)
        sample = sample.loc[finite].copy()
        coords = coords[finite]
        if len(sample) < 30:
            continue
        caf = sample["score_caf_myeloid_barrier"].to_numpy(float)
        tumor = sample["score_tumor_epithelial"].to_numpy(float)
        caf_core = caf >= np.nanpercentile(caf, 90)
        tumor_high = tumor >= np.nanpercentile(tumor, 75)
        sample["is_caf_core_top10"] = caf_core
        sample["is_tumor_high_top25"] = tumor_high
        if caf_core.sum() >= 3 and tumor_high.sum() >= 3:
            med_nn = np.nanmedian(cKDTree(coords).query(coords, k=2)[0][:, 1])
            radius = 2.0 * med_nn if np.isfinite(med_nn) and med_nn > 0 else np.inf
            caf_d = cKDTree(coords[caf_core]).query(coords, k=1)[0]
            tumor_d = cKDTree(coords[tumor_high]).query(coords, k=1)[0]
            sample["is_interface"] = (caf_d <= radius) & (tumor_d <= radius)
        else:
            sample["is_interface"] = False
        out.append(sample)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


def summarize_gene_enrichment(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    per_sample: list[dict] = []
    genes_present = [gene for gene in TARGET_GENES if gene in df.columns]
    for sample_id, sample in df.groupby("sample_id", sort=False):
        if sample["is_caf_core_top10"].sum() < 3:
            continue
        context = sample["cohort_context"].iloc[0]
        specimen = sample["specimen_type"].iloc[0] if "specimen_type" in sample else ""
        for axis, genes in GENE_SETS.items():
            present = [gene for gene in genes if gene in genes_present and sample[gene].notna().any()]
            if not present:
                continue
            axis_values = sample[present].mean(axis=1)
            z = (axis_values - axis_values.mean()) / axis_values.std(ddof=0) if axis_values.std(ddof=0) else axis_values * 0
            sample_axis = sample.assign(axis_z=z)
            core_enrichment = float(sample_axis.loc[sample_axis["is_caf_core_top10"], "axis_z"].mean() - sample_axis.loc[~sample_axis["is_caf_core_top10"], "axis_z"].mean())
            if sample_axis["is_interface"].sum() >= 3:
                interface_enrichment = float(sample_axis.loc[sample_axis["is_interface"], "axis_z"].mean() - sample_axis.loc[~sample_axis["is_interface"], "axis_z"].mean())
            else:
                interface_enrichment = math.nan
            per_sample.append(
                {
                    "sample_id": sample_id,
                    "cohort_context": context,
                    "specimen_type": specimen,
                    "axis_label": axis,
                    "n_genes_present": len(present),
                    "present_genes": ";".join(present),
                    "core_enrichment": core_enrichment,
                    "interface_enrichment": interface_enrichment,
                    "n_spots": int(len(sample_axis)),
                    "n_caf_core_spots": int(sample_axis["is_caf_core_top10"].sum()),
                    "n_interface_spots": int(sample_axis["is_interface"].sum()),
                }
            )
            for gene in present:
                values = sample[gene]
                sd = values.std(ddof=0)
                gene_z = (values - values.mean()) / sd if sd else values * 0
                gene_df = sample.assign(gene_z=gene_z)
                per_sample.append(
                    {
                        "sample_id": sample_id,
                        "cohort_context": context,
                        "specimen_type": specimen,
                        "axis_label": f"{axis}::{gene}",
                        "n_genes_present": 1,
                        "present_genes": gene,
                        "core_enrichment": float(gene_df.loc[gene_df["is_caf_core_top10"], "gene_z"].mean() - gene_df.loc[~gene_df["is_caf_core_top10"], "gene_z"].mean()),
                        "interface_enrichment": float(gene_df.loc[gene_df["is_interface"], "gene_z"].mean() - gene_df.loc[~gene_df["is_interface"], "gene_z"].mean()) if gene_df["is_interface"].sum() >= 3 else math.nan,
                        "n_spots": int(len(gene_df)),
                        "n_caf_core_spots": int(gene_df["is_caf_core_top10"].sum()),
                        "n_interface_spots": int(gene_df["is_interface"].sum()),
                    }
                )
    long = pd.DataFrame(per_sample)
    axis_long = long[~long["axis_label"].str.contains("::", regex=False)].copy()
    summary = (
        axis_long.groupby(["cohort_context", "axis_label"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_core_enrichment=("core_enrichment", "median"),
            n_core_positive=("core_enrichment", lambda s: int((s > 0).sum())),
            median_interface_enrichment=("interface_enrichment", "median"),
            n_interface_positive=("interface_enrichment", lambda s: int((s > 0).sum())),
            median_n_genes_present=("n_genes_present", "median"),
        )
        .sort_values(["cohort_context", "axis_label"])
    )
    return long, summary


def plot_summary(summary: pd.DataFrame) -> None:
    contexts = [
        "post_neoadjuvant_sections",
        "treatment_naive_primary",
        "primary_tumor",
        "liver_metastasis",
        "lymph_node_metastasis",
        "external_paired_st_anchor",
    ]
    contexts = [context for context in contexts if context in set(summary["cohort_context"])]
    axes = list(GENE_SETS)
    mat = (
        summary.pivot_table(index="axis_label", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
        .reindex(index=axes, columns=contexts)
    )
    fig, ax = plt.subplots(figsize=(8.3, 4.8), constrained_layout=True)
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-1.1, vmax=1.1, aspect="auto")
    ax.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS.get(c, c) for c in contexts], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(axes)), axes)
    ax.set_title("Targeted gene-level support for candidate CAF-core axes", fontsize=11, fontweight="bold", loc="left")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                row = summary[(summary["axis_label"].eq(mat.index[i])) & (summary["cohort_context"].eq(mat.columns[j]))]
                label = f"{val:.2f}\n{int(row['n_core_positive'].iloc[0])}/{int(row['n_samples'].iloc[0])}"
                ax.text(j, i, label, ha="center", va="center", fontsize=6.4)
    cb = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("median CAF-core enrichment", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    for ext in ["pdf", "svg", "png"]:
        path = OUT_FIG.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    mvp = load_mvp_spots()
    anchor = load_gse235315_spots()
    combined_spot = pd.concat([mvp, anchor], ignore_index=True)
    h5_expr = extract_h5_gene_values(manifest, combined_spot)
    rds_expr = extract_gse272362_gene_values()
    frames = [frame for frame in [h5_expr, rds_expr] if not frame.empty]
    if not frames:
        raise RuntimeError("No target-gene expression tables could be generated.")
    expr = pd.concat(frames, ignore_index=True)
    if "cohort_context" not in expr.columns:
        expr["cohort_context"] = expr.apply(context_for_row, axis=1)
    else:
        missing_context = expr["cohort_context"].isna() | expr["cohort_context"].astype(str).eq("")
        if missing_context.any():
            expr.loc[missing_context, "cohort_context"] = expr.loc[missing_context].apply(context_for_row, axis=1)
    expr = add_regions(expr)
    long, summary = summarize_gene_enrichment(expr)
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    long.to_csv(OUT_LONG, index=False)
    summary.to_csv(OUT_TABLE, index=False)
    plot_summary(summary)
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "27_targeted_gene_axis_validation",
                "status": "success",
                "n_spots_with_target_genes": int(len(expr)),
                "n_samples": int(expr["sample_id"].nunique()),
                "n_contexts": int(expr["cohort_context"].nunique()),
                "used_gse272362_gene_export": bool(not rds_expr.empty),
                "outputs": [str(OUT_LONG), str(OUT_TABLE), str(OUT_FIG.with_suffix(".pdf"))],
                "claim_boundary": "Targeted gene-level support for candidate axes; not ligand-receptor or causal validation.",
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Wrote targeted gene-axis validation for {expr['sample_id'].nunique()} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
