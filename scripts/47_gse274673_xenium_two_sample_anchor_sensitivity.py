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
for path in [TABLE_DIR, FIG_DIR, REPORT_DIR]:
    path.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260627)
N_RANDOM = 100


def load_stage46():
    path = PROJECT / "scripts" / "46_gse274673_xenium_pilot_expression_domain.py"
    spec = importlib.util.spec_from_file_location("stage46", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def sample_base(pilot_dir: str) -> Path:
    bases = list((DATA_DIR / pilot_dir).glob("output-*"))
    if not bases:
        raise FileNotFoundError(f"No output directory found under {pilot_dir}")
    return bases[0]


def load_sample_meta() -> dict[str, dict]:
    meta = pd.read_csv(TABLE_DIR / "gse274673_xenium_sample_metadata.csv")
    out = {}
    for _, row in meta.iterrows():
        out[row["geo_accession"]] = {
            "title": row["title"],
            "treatment": str(row["treatment"]).replace("treatment-naïve", "treatment-naive"),
            "pilot_dir": f"pilot_{row['geo_accession']}",
        }
    return out


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
    scores["anchor_SPP1TAM"] = scores["score_SPP1_TAM"]
    scores["anchor_TGFb_EMT"] = scores["score_TGFb_EMT"]

    xy = scores[["x_centroid", "y_centroid"]].to_numpy(float)
    scale = median_nn_scale(xy)
    rows = []
    anchors = ["score_CAF_matrix", "anchor_CAF_APC", "anchor_CAF_SPP1TAM"]
    targets = ["score_SPP1_TAM", "score_Tumor_epithelial", "score_SPP1_tumor_like", "score_IFN_APC", "score_T_NK", "score_TGFb_EMT"]
    for anchor in anchors:
        anchor_values = scores[anchor].to_numpy(float)
        n_anchor = max(50, int(np.ceil(0.10 * len(scores))))
        anchor_mask = anchor_values >= np.nanquantile(anchor_values, 1 - n_anchor / len(scores))
        dist = nearest_dist(xy, anchor_mask) / scale
        for target in targets:
            if target == anchor:
                continue
            target_values = scores[target].to_numpy(float)
            observed = safe_spearman(dist, target_values)
            random_rhos = []
            for _ in range(N_RANDOM):
                random_mask = np.zeros(len(scores), dtype=bool)
                random_mask[RNG.choice(len(scores), size=anchor_mask.sum(), replace=False)] = True
                random_rhos.append(safe_spearman(nearest_dist(xy, random_mask) / scale, target_values))
            random_median = float(np.nanmedian(random_rhos))
            rows.append(
                {
                    "dataset_id": "GSE274673",
                    "geo_accession": accession,
                    "title": meta["title"],
                    "treatment": meta["treatment"],
                    "anchor": anchor,
                    "target_program": target.replace("score_", ""),
                    "n_cells": int(len(scores)),
                    "n_anchor_cells": int(anchor_mask.sum()),
                    "observed_rho": observed,
                    "random_median_rho": random_median,
                    "delta_vs_random_median": observed - random_median,
                    "observed_more_negative_than_random_median": bool(observed < random_median),
                }
            )
    coverage["geo_accession"] = accession
    coverage["title"] = meta["title"]
    coverage["treatment"] = meta["treatment"]
    return pd.DataFrame(rows), coverage, scores


def make_figure(results: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(15.4, 4.8), gridspec_kw={"width_ratios": [1.5, 1.5, 1.0]})
    fig.suptitle("GSE274673 Xenium | six-sample expression-domain feasibility", fontsize=15, fontweight="bold")

    target_order = ["SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT", "Tumor_epithelial", "SPP1_tumor_like"]
    sample_order = ["GSM8454446", "GSM8454449", "GSM8454450", "GSM8454447", "GSM8454448", "GSM8454451"]
    for ax, anchor, title in [
        (axes[0], "anchor_CAF_APC", "A  CAF-APC anchor"),
        (axes[1], "anchor_CAF_SPP1TAM", "B  CAF-SPP1/TAM anchor"),
    ]:
        sub = results[results["anchor"].eq(anchor)].copy()
        mat = sub.pivot(index="target_program", columns="geo_accession", values="delta_vs_random_median").reindex(
            index=target_order, columns=sample_order
        )
        labels = (
            sub.assign(label=lambda x: x["delta_vs_random_median"].map(lambda v: f"{v:.2f}"))
            .pivot(index="target_program", columns="geo_accession", values="label")
            .reindex(index=target_order, columns=sample_order)
        )
        im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.35, vmax=0.35, aspect="auto")
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_xticks(range(len(sample_order)))
        ax.set_xticklabels(["Naive\nP1", "Naive\nP4", "Naive\nP5", "CRT\nP2", "CRT\nP3", "CRT\nP6"], fontsize=8)
        ax.set_yticks(range(len(target_order)))
        ax.set_yticklabels(target_order)
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                ax.text(j, i, labels.iat[i, j], ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04).set_label("delta vs random rho")

    axes[2].axis("off")
    axes[2].set_title("C  Decision rule", loc="left", fontweight="bold")
    axes[2].text(
        0,
        0.95,
        "Negative values support expression-domain\n"
        "centering beyond random cell anchors.\n\n"
        "This test distinguishes two questions:\n"
        "1. pure matrix CAF domains;\n"
        "2. antigen-presenting or SPP1/TAM-linked\n"
        "CAF domains.\n\n"
        "Promote this dataset only if the same\n"
        "anchor-target direction reproduces across\n"
        "the six-sample cohort.",
        va="top",
        fontsize=10,
        linespacing=1.25,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    out = FIG_DIR / "gse274673_xenium_anchor_sensitivity"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(results: pd.DataFrame, coverage: pd.DataFrame) -> None:
    lines = [
        "# GSE274673 Xenium Anchor Sensitivity",
        "",
        "Last updated: 2026-06-27",
        "",
        "## Scope",
        "",
        "Analyzed all six available GSE274673 Xenium PDAC samples: three treatment-naive and three chemoradiotherapy-treated sections.",
        "",
        "## Outputs",
        "",
        "- `results/tables/gse274673_xenium_anchor_sensitivity.csv`",
        "- `results/tables/gse274673_xenium_anchor_context_summary.csv`",
        "- `results/tables/gse274673_xenium_signature_coverage.csv`",
        "- `results/figures/submission/gse274673_xenium_anchor_sensitivity.pdf`",
        "",
        "## Initial Interpretation",
        "",
        "Negative delta values indicate target-program centering around expression-defined anchors beyond same-size random cell anchors.",
    ]
    focus = results[results["anchor"].isin(["anchor_CAF_APC", "anchor_CAF_SPP1TAM"])]
    for _, row in focus.sort_values(["anchor", "target_program", "geo_accession"]).iterrows():
        lines.append(
            f"- {row.geo_accession} / {row.anchor} -> {row.target_program}: delta {row.delta_vs_random_median:.3f}, support={row.observed_more_negative_than_random_median}."
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            f"This is a six-sample feasibility analysis using {N_RANDOM} random anchors per sample. If promoted to the manuscript, rerun the fixed anchor definitions with 1,000 random anchors per sample and package panel-level source data.",
        ]
    )
    (REPORT_DIR / "gse274673_xenium_anchor_sensitivity.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    stage46 = load_stage46()
    sample_meta = load_sample_meta()
    all_results = []
    all_coverage = []
    for accession, meta in sample_meta.items():
        results, coverage, _ = score_sample(stage46, accession, meta)
        all_results.append(results)
        all_coverage.append(coverage)
    results = pd.concat(all_results, ignore_index=True)
    coverage = pd.concat(all_coverage, ignore_index=True)
    context = (
        results.groupby(["anchor", "target_program"], as_index=False)
        .agg(
            n_samples=("geo_accession", "nunique"),
            median_delta_vs_random=("delta_vs_random_median", "median"),
            n_support=("observed_more_negative_than_random_median", "sum"),
        )
    )
    results.to_csv(TABLE_DIR / "gse274673_xenium_anchor_sensitivity.csv", index=False)
    context.to_csv(TABLE_DIR / "gse274673_xenium_anchor_context_summary.csv", index=False)
    coverage.to_csv(TABLE_DIR / "gse274673_xenium_signature_coverage.csv", index=False)
    make_figure(results)
    write_report(results, coverage)
    print("Wrote GSE274673 Xenium anchor sensitivity outputs.")


if __name__ == "__main__":
    main()
