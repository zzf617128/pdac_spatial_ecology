from __future__ import annotations

import importlib.util
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT / "data" / "external" / "GSE274673"
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
REPORT_DIR = PROJECT / "results" / "reports"
SOURCE_DIR = PROJECT / "results" / "source_data"
for path in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR]:
    path.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260627)
N_RANDOM = 1000


def load_stage46():
    path = PROJECT / "scripts" / "46_gse274673_xenium_pilot_expression_domain.py"
    spec = importlib.util.spec_from_file_location("stage46", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def load_sample_meta() -> dict[str, dict]:
    meta = pd.read_csv(TABLE_DIR / "gse274673_xenium_sample_metadata.csv")
    out = {}
    for _, row in meta.iterrows():
        treatment = str(row["treatment"]).replace("treatment-naïve", "treatment-naive")
        out[row["geo_accession"]] = {
            "title": row["title"],
            "treatment": treatment,
            "patient_id": row["title"].replace("Patient ", "P").replace(" PDAC_", "_PDAC"),
            "pilot_dir": f"pilot_{row['geo_accession']}",
        }
    return out


def sample_base(pilot_dir: str) -> Path:
    bases = list((DATA_DIR / pilot_dir).glob("output-*"))
    if not bases:
        raise FileNotFoundError(f"No output directory found under {pilot_dir}")
    return bases[0]


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 10 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def nearest_dist(xy: np.ndarray, anchor_mask: np.ndarray) -> np.ndarray:
    tree = cKDTree(xy[anchor_mask])
    dists, _ = tree.query(xy, k=1)
    return dists


def median_nn_scale(xy: np.ndarray) -> float:
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    val = float(np.nanmedian(dists[:, 1]))
    return val if np.isfinite(val) and val > 0 else 1.0


def score_sample(stage46, accession: str, meta: dict) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base = sample_base(meta["pilot_dir"])
    cells = pd.read_csv(base / "cells.csv.gz")
    matrix, genes, barcodes = stage46.read_10x_h5(base / "cell_feature_matrix.h5")
    if list(cells["cell_id"]) != barcodes:
        cells = cells.set_index("cell_id").reindex(barcodes).reset_index()
    scores, coverage = stage46.score_modules(matrix, genes, cells)
    scores["anchor_CAF_APC"] = (scores["score_CAF_matrix"] + scores["score_IFN_APC"]) / 2
    scores["anchor_CAF_SPP1TAM"] = (scores["score_CAF_matrix"] + scores["score_SPP1_TAM"]) / 2

    xy = scores[["x_centroid", "y_centroid"]].to_numpy(float)
    scale = median_nn_scale(xy)
    rows = []
    anchors = {
        "CAF_APC": "anchor_CAF_APC",
        "CAF_SPP1TAM": "anchor_CAF_SPP1TAM",
    }
    targets = {
        "SPP1_TAM": "score_SPP1_TAM",
        "IFN_APC": "score_IFN_APC",
        "T_NK": "score_T_NK",
        "TGFb_EMT": "score_TGFb_EMT",
        "Tumor_epithelial": "score_Tumor_epithelial",
        "SPP1_tumor_like": "score_SPP1_tumor_like",
    }
    for anchor_label, anchor_col in anchors.items():
        anchor_values = scores[anchor_col].to_numpy(float)
        n_anchor = max(50, int(np.ceil(0.10 * len(scores))))
        anchor_mask = anchor_values >= np.nanquantile(anchor_values, 1 - n_anchor / len(scores))
        observed_dist = nearest_dist(xy, anchor_mask) / scale
        for target_label, target_col in targets.items():
            target_values = scores[target_col].to_numpy(float)
            observed = safe_spearman(observed_dist, target_values)
            random_rhos = []
            for _ in range(N_RANDOM):
                random_mask = np.zeros(len(scores), dtype=bool)
                random_mask[RNG.choice(len(scores), size=anchor_mask.sum(), replace=False)] = True
                random_rhos.append(safe_spearman(nearest_dist(xy, random_mask) / scale, target_values))
            random_rhos = np.array(random_rhos, dtype=float)
            random_median = float(np.nanmedian(random_rhos))
            rows.append(
                {
                    "dataset_id": "GSE274673_Xenium",
                    "geo_accession": accession,
                    "title": meta["title"],
                    "patient_id": meta["patient_id"],
                    "treatment": meta["treatment"],
                    "anchor": anchor_label,
                    "target_program": target_label,
                    "n_cells": int(len(scores)),
                    "n_anchor_cells": int(anchor_mask.sum()),
                    "observed_rho": observed,
                    "random_median_rho": random_median,
                    "delta_vs_random_median": observed - random_median,
                    "observed_more_negative_than_random_median": bool(observed < random_median),
                    "random_p05_rho": float(np.nanpercentile(random_rhos, 5)),
                    "random_p95_rho": float(np.nanpercentile(random_rhos, 95)),
                }
            )

    cell_score_cols = [
        "cell_id",
        "x_centroid",
        "y_centroid",
        "total_counts",
        "score_CAF_matrix",
        "score_IFN_APC",
        "score_SPP1_TAM",
        "score_T_NK",
        "score_TGFb_EMT",
        "score_Tumor_epithelial",
        "score_SPP1_tumor_like",
        "anchor_CAF_APC",
        "anchor_CAF_SPP1TAM",
    ]
    cell_scores = scores[cell_score_cols].copy()
    cell_scores["geo_accession"] = accession
    cell_scores["title"] = meta["title"]
    cell_scores["treatment"] = meta["treatment"]

    coverage["geo_accession"] = accession
    coverage["title"] = meta["title"]
    coverage["treatment"] = meta["treatment"]
    return pd.DataFrame(rows), coverage, cell_scores


def make_figure(results: pd.DataFrame, composition: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    sample_order = ["GSM8454446", "GSM8454449", "GSM8454450", "GSM8454447", "GSM8454448", "GSM8454451"]
    sample_labels = ["Naive\nP1", "Naive\nP4", "Naive\nP5", "CRT\nP2", "CRT\nP3", "CRT\nP6"]
    target_order = ["SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT", "Tumor_epithelial", "SPP1_tumor_like"]

    fig = plt.figure(figsize=(12.8, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.0], height_ratios=[1.15, 0.78], wspace=0.52, hspace=0.55)
    fig.suptitle("Extended Data Fig. 11 | GSE274673 Xenium cell-resolution CAF-domain validation", fontsize=16, fontweight="bold")

    for ax, anchor, title in [
        (fig.add_subplot(gs[0, 0]), "CAF_APC", "A  CAF-APC anchor"),
        (fig.add_subplot(gs[0, 1]), "CAF_SPP1TAM", "B  CAF-SPP1/TAM anchor"),
    ]:
        sub = results[results["anchor"].eq(anchor)]
        mat = sub.pivot(index="target_program", columns="geo_accession", values="delta_vs_random_median").reindex(
            index=target_order, columns=sample_order
        )
        labels = mat.map(lambda v: f"{v:.2f}" if pd.notna(v) else "NA")
        im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xticks(range(len(sample_order)))
        ax.set_xticklabels(sample_labels, fontsize=8)
        ax.set_yticks(range(len(target_order)))
        ax.set_yticklabels(target_order)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                value = mat.iat[i, j]
                text_color = "white" if pd.notna(value) and abs(value) >= 0.34 else "#111111"
                ax.text(j, i, labels.iat[i, j], ha="center", va="center", fontsize=8, color=text_color)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("observed-minus-random rho")

    ax_c = fig.add_subplot(gs[1, 0])
    focus = results[
        results["target_program"].isin(["SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT"])
        & results["anchor"].isin(["CAF_APC", "CAF_SPP1TAM"])
    ]
    summary = (
        focus.groupby(["anchor", "target_program"], as_index=False)
        .agg(median_delta=("delta_vs_random_median", "median"), n_support=("observed_more_negative_than_random_median", "sum"))
    )
    summary["label"] = summary["anchor"] + " -> " + summary["target_program"]
    ax_c.barh(range(len(summary)), summary["median_delta"], color="#3B78A8")
    ax_c.axvline(0, color="#333333", lw=0.8)
    ax_c.set_yticks(range(len(summary)))
    ax_c.set_yticklabels(summary["label"], fontsize=8)
    ax_c.invert_yaxis()
    ax_c.set_xlabel("median delta vs random")
    ax_c.set_title("C  Cohort-level support", loc="left", fontweight="bold")
    for i, row in summary.iterrows():
        ax_c.text(row["median_delta"], i, f" {int(row['n_support'])}/6", va="center", fontsize=8)

    ax_d = fig.add_subplot(gs[1, 1])
    comp = composition.set_index("geo_accession").reindex(sample_order)
    x = np.arange(len(sample_order))
    ax_d.bar(x, comp["n_cells"], color="#7A9E9F")
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(sample_labels, fontsize=8)
    ax_d.set_ylabel("cells")
    ax_d.set_title("D  Xenium sample scale", loc="left", fontweight="bold")
    ax_d.set_ylim(0, float(comp["n_cells"].max()) * 1.14)
    for i, val in enumerate(comp["n_cells"]):
        ax_d.text(i, val, f"{int(val/1000)}k", ha="center", va="bottom", fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    out = FIG_DIR / "extended_data_figure11_gse274673_xenium_cell_resolution"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(results: pd.DataFrame, context: pd.DataFrame, composition: pd.DataFrame) -> None:
    lines = [
        "# GSE274673 Xenium Fixed-Anchor Validation",
        "",
        "Last updated: 2026-06-27",
        "",
        "## Scope",
        "",
        "Analyzed all six GSE274673 Xenium PDAC sections with fixed CAF-APC and CAF-SPP1/TAM anchor definitions and 1,000 same-size random cell anchors per sample.",
        "",
        "## Outputs",
        "",
        "- `results/tables/gse274673_xenium_fixed_anchor_gradients.csv`",
        "- `results/tables/gse274673_xenium_fixed_anchor_context_summary.csv`",
        "- `results/tables/gse274673_xenium_fixed_anchor_sample_composition.csv`",
        "- `results/tables/gse274673_xenium_fixed_anchor_signature_coverage.csv`",
        "- `results/source_data/Source_Data_Extended_Data_Fig_11A_B.csv`",
        "- `results/source_data/Source_Data_Extended_Data_Fig_11C.csv`",
        "- `results/source_data/Source_Data_Extended_Data_Fig_11D.csv`",
        "- `results/figures/submission/extended_data_figure11_gse274673_xenium_cell_resolution.pdf`",
        "",
        "## Main Results",
        "",
        "Negative delta values indicate target-program centering around expression-defined anchors beyond same-size random cell anchors.",
    ]
    for _, row in context.sort_values(["anchor", "target_program"]).iterrows():
        lines.append(
            f"- {row.anchor} -> {row.target_program}: median delta {row.median_delta_vs_random:.3f}, support {int(row.n_support)}/{int(row.n_samples)} samples."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "GSE274673 supports a cell-resolution immune/myeloid CAF-domain layer. CAF-APC and CAF-SPP1/TAM domains organize SPP1/TAM, IFN/APC and T/NK programs across treatment-naive and chemoradiotherapy-treated sections. Tumor epithelial and SPP1-tumor-like programs are not centered on these anchors, so this analysis should not be used to claim direct CAF-to-tumor epithelial proximity.",
        ]
    )
    (REPORT_DIR / "gse274673_xenium_fixed_anchor_validation.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    stage46 = load_stage46()
    sample_meta = load_sample_meta()
    all_results = []
    all_coverage = []
    all_cells = []
    for accession, meta in sample_meta.items():
        results, coverage, cell_scores = score_sample(stage46, accession, meta)
        all_results.append(results)
        all_coverage.append(coverage)
        all_cells.append(cell_scores)
    results = pd.concat(all_results, ignore_index=True)
    coverage = pd.concat(all_coverage, ignore_index=True)
    cells = pd.concat(all_cells, ignore_index=True)
    composition = (
        cells.groupby(["geo_accession", "title", "treatment"], as_index=False)
        .agg(n_cells=("cell_id", "count"), median_total_counts=("total_counts", "median"))
    )
    context = (
        results.groupby(["anchor", "target_program"], as_index=False)
        .agg(
            n_samples=("geo_accession", "nunique"),
            median_delta_vs_random=("delta_vs_random_median", "median"),
            n_support=("observed_more_negative_than_random_median", "sum"),
            median_observed_rho=("observed_rho", "median"),
        )
    )

    results.to_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_gradients.csv", index=False)
    context.to_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_context_summary.csv", index=False)
    composition.to_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv", index=False)
    coverage.to_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_signature_coverage.csv", index=False)
    cells.to_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_cell_scores.csv", index=False)

    results.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_11A_B.csv", index=False)
    context.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_11C.csv", index=False)
    composition.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_11D.csv", index=False)

    make_figure(results, composition)
    write_report(results, context, composition)
    print("Wrote GSE274673 Xenium fixed-anchor validation outputs.")


if __name__ == "__main__":
    main()
