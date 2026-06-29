from __future__ import annotations

from pathlib import Path

import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT / "data" / "external" / "GSE274673"
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "internal"
REPORT_DIR = PROJECT / "results" / "reports"
for path in [TABLE_DIR, FIG_DIR, REPORT_DIR]:
    path.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260627)

SIGNATURES = {
    "CAF_matrix": ["FAP", "LUM", "DCN", "ACTA2", "PDGFRA", "PDGFRB", "COL1A2", "COL6A1", "COL6A2", "COL6A3", "FN1", "THBS2", "MMP2", "ADAM12"],
    "SPP1_TAM": ["SPP1", "CD68", "CD163", "APOE", "C1QA", "C1QB", "C1QC", "LYZ", "AIF1", "TYROBP", "CSF1R", "LGALS3", "MARCO"],
    "Tumor_epithelial": ["EPCAM", "KRT19", "KRT8", "KRT18", "KRT7", "TACSTD2", "CEACAM6", "S100P", "CLDN4", "MUC1"],
    "SPP1_tumor_like": ["SPP1", "EPCAM", "KRT19", "KRT8", "KRT18", "CEACAM6", "S100P"],
    "IFN_APC": ["CXCL9", "CXCL10", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "B2M", "STAT1", "IFIT1", "ISG15", "TAP1"],
    "T_NK": ["CD3D", "CD3E", "CD3G", "CD2", "CD4", "CD8A", "CD8B", "NKG7", "GNLY", "GZMB", "GZMK"],
    "TGFb_EMT": ["TGFB1", "TGFBR1", "TGFBR2", "VIM", "SNAI1", "SNAI2", "ITGA5", "ITGB1", "MMP14", "TGFBI", "SERPINH1"],
}


def find_pilot_base() -> Path:
    bases = list((DATA_DIR / "pilot_GSM8454448").glob("output-*"))
    if not bases:
        raise FileNotFoundError("GSE274673 pilot output directory not found")
    return bases[0]


def read_10x_h5(path: Path) -> tuple[csc_matrix, list[str], list[str]]:
    with h5py.File(path, "r") as f:
        group = f["matrix"]
        data = group["data"][:]
        indices = group["indices"][:]
        indptr = group["indptr"][:]
        shape = tuple(group["shape"][:])
        genes = [x.decode("utf-8") for x in group["features"]["name"][:]]
        barcodes = [x.decode("utf-8") for x in group["barcodes"][:]]
    matrix = csc_matrix((data, indices, indptr), shape=shape)
    return matrix, genes, barcodes


def zscore(values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    mean = np.nanmean(values)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros_like(values, dtype=float)
    return (values - mean) / sd


def score_modules(matrix: csc_matrix, genes: list[str], cells: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    gene_to_idx = {g.upper(): i for i, g in enumerate(genes)}
    counts = matrix.T.tocsr()
    lib = np.asarray(counts.sum(axis=1)).ravel()
    lib[lib <= 0] = 1
    norm = counts.multiply(1000 / lib[:, None]).tocsr()
    log_norm = norm.copy()
    log_norm.data = np.log1p(log_norm.data)

    coverage_rows = []
    out = cells.copy()
    for name, gene_list in SIGNATURES.items():
        present = [g for g in gene_list if g.upper() in gene_to_idx]
        coverage_rows.append(
            {
                "signature": name,
                "n_genes": len(gene_list),
                "n_present": len(present),
                "present_genes": ";".join(present),
                "missing_genes": ";".join([g for g in gene_list if g.upper() not in gene_to_idx]),
            }
        )
        if not present:
            out[f"score_{name}"] = np.nan
            continue
        idx = [gene_to_idx[g.upper()] for g in present]
        gene_expr = log_norm[:, idx].toarray()
        gene_z = np.apply_along_axis(zscore, 0, gene_expr)
        out[f"score_{name}"] = np.nanmean(gene_z, axis=1)
    return out, pd.DataFrame(coverage_rows)


def median_nn_scale(xy: np.ndarray) -> float:
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    val = float(np.nanmedian(dists[:, 1]))
    return val if np.isfinite(val) and val > 0 else 1.0


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def nearest_dist(xy: np.ndarray, anchor_mask: np.ndarray) -> np.ndarray:
    tree = cKDTree(xy[anchor_mask])
    dists, _ = tree.query(xy, k=1)
    return dists


def run_gradient_test(scores: pd.DataFrame, n_random: int = 1000) -> tuple[pd.DataFrame, pd.DataFrame]:
    xy = scores[["x_centroid", "y_centroid"]].to_numpy(float)
    scale = median_nn_scale(xy)
    anchor_score = scores["score_CAF_matrix"].to_numpy(float)
    n_anchor = max(50, int(np.ceil(0.10 * len(scores))))
    threshold = np.nanquantile(anchor_score, 1 - n_anchor / len(scores))
    anchor_mask = anchor_score >= threshold
    if anchor_mask.sum() < 10:
        anchor_mask[np.argsort(anchor_score)[-n_anchor:]] = True

    dist = nearest_dist(xy, anchor_mask) / scale
    rows = []
    for sig in ["SPP1_TAM", "Tumor_epithelial", "SPP1_tumor_like", "IFN_APC", "T_NK", "TGFb_EMT"]:
        col = f"score_{sig}"
        vals = scores[col].to_numpy(float)
        obs = safe_spearman(dist, vals)
        random_rhos = []
        for _ in range(n_random):
            rand_mask = np.zeros(len(scores), dtype=bool)
            rand_mask[RNG.choice(len(scores), size=anchor_mask.sum(), replace=False)] = True
            rand_dist = nearest_dist(xy, rand_mask) / scale
            random_rhos.append(safe_spearman(rand_dist, vals))
        random_rhos = np.array(random_rhos, dtype=float)
        rand_med = float(np.nanmedian(random_rhos))
        rows.append(
            {
                "dataset_id": "GSE274673",
                "geo_accession": "GSM8454448",
                "title": "Patient 3 PDAC_3",
                "treatment": "chemoradiotherapy-treated",
                "anchor": "CAF_matrix_top10pct",
                "target_program": sig,
                "n_cells": int(len(scores)),
                "n_anchor_cells": int(anchor_mask.sum()),
                "observed_rho": obs,
                "random_median_rho": rand_med,
                "delta_vs_random_median": obs - rand_med,
                "observed_more_negative_than_random_median": bool(obs < rand_med),
            }
        )
    cell_table = scores.copy()
    cell_table["is_CAF_matrix_anchor"] = anchor_mask
    cell_table["dist_to_CAF_matrix_anchor_norm"] = dist
    return pd.DataFrame(rows), cell_table


def make_figure(summary: pd.DataFrame, cells: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1.2, 1.0]})
    fig.suptitle("GSE274673 Xenium pilot | expression-defined CAF domain test", fontsize=15, fontweight="bold")

    ax = axes[0]
    plot = cells.sample(min(len(cells), 12000), random_state=7)
    sc = ax.scatter(plot["x_centroid"], plot["y_centroid"], c=plot["score_CAF_matrix"], s=1, cmap="viridis", rasterized=True)
    ax.set_title("A  CAF/matrix score", loc="left", fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)

    ax = axes[1]
    order = ["SPP1_TAM", "Tumor_epithelial", "SPP1_tumor_like", "IFN_APC", "T_NK", "TGFb_EMT"]
    data = summary.set_index("target_program").reindex(order)
    colors = ["#2f78b7" if v else "#d55e5e" for v in data["observed_more_negative_than_random_median"]]
    ax.barh(range(len(order)), data["delta_vs_random_median"], color=colors)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels(order)
    ax.invert_yaxis()
    ax.set_xlabel("observed-minus-random rho")
    ax.set_title("B  Distance-to-CAF domain vs random", loc="left", fontweight="bold")
    for i, val in enumerate(data["delta_vs_random_median"]):
        ax.text(val, i, f" {val:.2f}", va="center", ha="left" if val >= 0 else "right", fontsize=8)

    axes[2].axis("off")
    axes[2].set_title("C  Interpretation boundary", loc="left", fontweight="bold")
    axes[2].text(
        0,
        0.92,
        "This pilot uses expression-defined\n"
        "CAF/matrix anchor cells rather than\n"
        "coarse annotations.\n\n"
        "Negative values indicate that target\n"
        "programs are more centered on the\n"
        "CAF/matrix domain than expected from\n"
        "same-size random cell anchors.\n\n"
        "Use as feasibility evidence only until\n"
        "a treatment-naive sample is added.",
        va="top",
        fontsize=10,
        linespacing=1.25,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIG_DIR / "gse274673_xenium_pilot_expression_domain"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame, coverage: pd.DataFrame) -> None:
    lines = [
        "# GSE274673 Xenium Pilot Expression-Domain Validation",
        "",
        "Last updated: 2026-06-27",
        "",
        "## Scope",
        "",
        "Downloaded and parsed one GSE274673 Xenium chemoradiotherapy-treated PDAC sample, GSM8454448 / Patient 3 PDAC_3.",
        "",
        "## Outputs",
        "",
        "- `results/tables/gse274673_xenium_pilot_signature_coverage.csv`",
        "- `results/tables/gse274673_xenium_pilot_caf_domain_gradients.csv`",
        "- `results/tables/gse274673_xenium_pilot_cell_scores.csv`",
        "- `results/figures/internal/gse274673_xenium_pilot_expression_domain.pdf`",
        "",
        "## Initial Interpretation",
        "",
        "Negative delta values indicate that the target program is more CAF/matrix-domain centered than expected from same-size random cell anchors.",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- {row.target_program}: observed rho {row.observed_rho:.3f}, delta vs random {row.delta_vs_random_median:.3f}, support={row.observed_more_negative_than_random_median}."
        )
    lines.extend(
        [
            "",
            "## Coverage",
            "",
        ]
    )
    for _, row in coverage.iterrows():
        lines.append(f"- {row.signature}: {int(row.n_present)}/{int(row.n_genes)} genes present.")
    lines.extend(
        [
            "",
            "## Next Step",
            "",
            "Download one treatment-naive GSE274673 sample and rerun the same expression-domain analysis. Promote this route only if the CAF/matrix-domain centering is reproducible across at least one naive and one treated sample.",
        ]
    )
    (REPORT_DIR / "gse274673_xenium_pilot_expression_domain.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    base = find_pilot_base()
    cells = pd.read_csv(base / "cells.csv.gz")
    matrix, genes, barcodes = read_10x_h5(base / "cell_feature_matrix.h5")
    if list(cells["cell_id"]) != barcodes:
        cells = cells.set_index("cell_id").reindex(barcodes).reset_index()
    scores, coverage = score_modules(matrix, genes, cells)
    summary, cell_table = run_gradient_test(scores)
    coverage.to_csv(TABLE_DIR / "gse274673_xenium_pilot_signature_coverage.csv", index=False)
    summary.to_csv(TABLE_DIR / "gse274673_xenium_pilot_caf_domain_gradients.csv", index=False)
    cell_table.to_csv(TABLE_DIR / "gse274673_xenium_pilot_cell_scores.csv", index=False)
    make_figure(summary, cell_table)
    write_report(summary, coverage)
    print("Wrote GSE274673 Xenium pilot outputs.")


if __name__ == "__main__":
    main()
