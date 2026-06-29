from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_SCORES = {
    "immune_core": "score_immune_hub_core",
    "ifn_mhc": "z_ifn_antigen_presentation",
    "tumor_aggressive": "score_tumor_aggressive",
    "immune_maturity": "score_immune_hub_maturity",
}


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 30:
        return np.nan
    if len(np.unique(x[mask])) < 2 or len(np.unique(y[mask])) < 2:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).correlation)


def median_neighbor_distance(coords: np.ndarray) -> float:
    if len(coords) < 3:
        return np.nan
    tree = cKDTree(coords)
    distances, _ = tree.query(coords, k=2)
    return float(np.nanmedian(distances[:, 1]))


def analyze_sample(group: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    safe = group[~group["edge_or_background_risk"].fillna(False).astype(bool)].copy()
    if len(safe) < 100:
        return [], []
    coords = safe[["x_pixel", "y_pixel"]].to_numpy(float)
    nn = median_neighbor_distance(coords)
    if not np.isfinite(nn) or nn <= 0:
        return [], []

    caf_threshold = np.nanpercentile(safe["score_caf_myeloid_barrier"], 90)
    core = safe[safe["score_caf_myeloid_barrier"] >= caf_threshold]
    if len(core) < 10:
        return [], []
    core_tree = cKDTree(core[["x_pixel", "y_pixel"]].to_numpy(float))
    distances_px, _ = core_tree.query(coords, k=1)
    safe["distance_to_caf_core_ring"] = distances_px / nn
    safe["caf_core_bin"] = pd.cut(
        safe["distance_to_caf_core_ring"],
        bins=[-0.01, 1.5, 3, 6, np.inf],
        labels=["core_0_1.5", "near_1.5_3", "mid_3_6", "far_gt6"],
    )

    gradient_rows: list[dict] = []
    bin_rows: list[dict] = []
    base = {
        "dataset_id": safe["dataset_id"].iloc[0],
        "sample_id": safe["sample_id"].iloc[0],
        "n_safe_spots": len(safe),
        "n_caf_core_spots": len(core),
        "caf_core_threshold": float(caf_threshold),
        "median_neighbor_distance_px": float(nn),
    }
    for label, score_col in TARGET_SCORES.items():
        rho = safe_spearman(safe["distance_to_caf_core_ring"].to_numpy(float), safe[score_col].to_numpy(float))
        gradient_rows.append({**base, "target": label, "rho_distance_to_caf_core": rho})
        for bin_name, bin_df in safe.groupby("caf_core_bin", observed=True):
            bin_rows.append(
                {
                    **base,
                    "target": label,
                    "distance_bin": str(bin_name),
                    "n_spots_bin": len(bin_df),
                    "mean_score": float(np.nanmean(bin_df[score_col])),
                    "median_score": float(np.nanmedian(bin_df[score_col])),
                }
            )
    return gradient_rows, bin_rows


def write_status(stage: str, status: str, payload: dict) -> None:
    base = {
        "stage": stage,
        "status": status,
        "timestamp_utc": now_iso(),
        "n_datasets_processed": 0,
        "n_samples_processed": 0,
        "n_errors": 0,
        "critical_errors": [],
        "noncritical_warnings": [],
        "next_manual_check": [],
    }
    base.update(payload)
    path = PROJECT_ROOT / f"results/logs/stage_{stage}_status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def make_figure(gradients: pd.DataFrame, bins: pd.DataFrame) -> None:
    out_dir = PROJECT_ROOT / "results/figures/mvp/caf_myeloid_niche"
    out_dir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)

    ax = axes[0]
    targets = list(TARGET_SCORES.keys())
    data = [gradients[gradients["target"] == target]["rho_distance_to_caf_core"].dropna() for target in targets]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(targets) + 1), [target.replace("_", "\n") for target in targets])
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Spearman rho: distance to CAF-myeloid core vs score")
    ax.set_title("Negative values indicate co-localization with CAF-myeloid core")

    ax = axes[1]
    order = ["core_0_1.5", "near_1.5_3", "mid_3_6", "far_gt6"]
    for target in ["immune_core", "ifn_mhc", "tumor_aggressive"]:
        sub = bins[bins["target"] == target]
        means = []
        errs = []
        for bin_name in order:
            values = sub[sub["distance_bin"] == bin_name]["mean_score"].dropna().to_numpy()
            means.append(float(np.nanmedian(values)) if len(values) else np.nan)
            errs.append(float(np.nanpercentile(values, 75) - np.nanpercentile(values, 25)) / 2 if len(values) else np.nan)
        ax.errorbar(range(len(order)), means, yerr=errs, marker="o", capsize=3, label=target.replace("_", " "))
    ax.set_xticks(range(len(order)), order, rotation=25, ha="right")
    ax.set_ylabel("Median sample-level mean score")
    ax.set_title("Spatial programs around high CAF-myeloid cores")
    ax.legend(frameon=False)

    fig.suptitle("CAF-myeloid Niche Distance Gradients")
    fig.savefig(out_dir / "caf_myeloid_niche_gradients.png", dpi=180)
    fig.savefig(out_dir / "caf_myeloid_niche_gradients.pdf")
    plt.close(fig)


def main() -> int:
    spot_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    spot = pd.read_csv(spot_path)
    gradient_rows: list[dict] = []
    bin_rows: list[dict] = []
    for _, group in spot.groupby(["dataset_id", "sample_id"], sort=False):
        g_rows, b_rows = analyze_sample(group)
        gradient_rows.extend(g_rows)
        bin_rows.extend(b_rows)
    gradients = pd.DataFrame(gradient_rows)
    bins = pd.DataFrame(bin_rows)
    gradients.to_csv(PROJECT_ROOT / "results/tables/caf_myeloid_niche_gradient_stats.csv", index=False)
    bins.to_csv(PROJECT_ROOT / "results/tables/caf_myeloid_niche_distance_bins.csv", index=False)
    make_figure(gradients, bins)

    write_status(
        "09_caf_myeloid_niche_gradients",
        "success",
        {
            "n_datasets_processed": int(gradients["dataset_id"].nunique()) if not gradients.empty else 0,
            "n_samples_processed": int(gradients["sample_id"].nunique()) if not gradients.empty else 0,
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": [
                "This is a distance-to-high-score-region analysis, not a causal barrier model."
            ],
            "next_manual_check": [
                "Inspect caf_myeloid_niche_gradients.png.",
                "If gradients are consistent, use this as the core spatial ecology result.",
            ],
        },
    )
    print(f"Computed CAF-myeloid niche gradients for {gradients['sample_id'].nunique()} samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

