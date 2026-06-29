from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OVERLAY_SCORES = [
    ("score_caf_myeloid_barrier", "CAF-myeloid"),
    ("score_immune_hub_core", "immune core"),
    ("score_tumor_aggressive", "tumor aggressive"),
    ("z_ifn_antigen_presentation", "IFN/MHC z"),
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0))


def choose_samples(sample: pd.DataFrame, top_n: int) -> pd.DataFrame:
    edge_summary_path = PROJECT_ROOT / "results/tables/mvp_edge_qc_sample_summary.csv"
    if edge_summary_path.exists():
        edge = pd.read_csv(edge_summary_path)
        sample = sample.merge(edge, on=["dataset_id", "sample_id"], how="left")
    candidates: list[pd.Series] = []
    ranking_cols = [
        ("p90_caf_myeloid_barrier", "high_caf_myeloid", "score_caf_myeloid_barrier_high_fraction_risk"),
        ("p90_immune_hub_core", "high_immune_core", "score_immune_hub_core_high_fraction_risk"),
        ("p90_tumor_aggressive", "high_tumor_aggressive", "score_tumor_aggressive_high_fraction_risk"),
        ("p90_ifn_antigen_presentation", "high_ifn_mhc", "z_ifn_antigen_presentation_high_fraction_risk"),
    ]
    for col, reason, risk_col in ranking_cols:
        if col not in sample.columns:
            continue
        eligible = sample.copy()
        if "fraction_edge_or_background_risk" in eligible.columns and risk_col in eligible.columns:
            filtered = eligible[
                (eligible["fraction_edge_or_background_risk"].fillna(1.0) <= 0.40)
                & (eligible[risk_col].fillna(1.0) <= 0.35)
            ].copy()
            if len(filtered) >= top_n:
                eligible = filtered
        ranked = eligible.sort_values(col, ascending=False).head(top_n)
        for _, row in ranked.iterrows():
            item = row.copy()
            item["selection_reason"] = reason
            candidates.append(item)
    selected = pd.DataFrame(candidates)
    if selected.empty:
        return selected
    selected = selected.drop_duplicates(subset=["dataset_id", "sample_id"]).reset_index(drop=True)
    return selected


def render_overlay(
    sample_id: str,
    image_path: Path,
    scalefactors_path: Path,
    spots: pd.DataFrame,
    output_path: Path,
) -> None:
    image = Image.open(image_path).convert("RGB")
    scale = read_scale(scalefactors_path)
    x = spots["x_pixel"].to_numpy(dtype=float) * scale
    y = spots["y_pixel"].to_numpy(dtype=float) * scale

    fig, axes = plt.subplots(1, len(OVERLAY_SCORES) + 1, figsize=(5 * (len(OVERLAY_SCORES) + 1), 5))
    axes[0].imshow(image)
    axes[0].set_title("H&E")
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    for ax, (score, title) in zip(axes[1:], OVERLAY_SCORES):
        ax.imshow(image)
        values = spots[score].to_numpy(dtype=float)
        lo = np.nanpercentile(values, 2)
        hi = np.nanpercentile(values, 98)
        scatter = ax.scatter(
            x,
            y,
            c=values,
            s=7,
            cmap="viridis",
            alpha=0.85,
            linewidths=0,
            vmin=lo,
            vmax=hi,
        )
        if "edge_or_background_risk" in spots.columns:
            risk = spots["edge_or_background_risk"].fillna(False).astype(bool).to_numpy()
            if risk.any():
                ax.scatter(
                    x[risk],
                    y[risk],
                    facecolors="none",
                    edgecolors="#d62728",
                    s=18,
                    linewidths=0.35,
                    alpha=0.65,
                )
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(scatter, ax=ax, fraction=0.046, pad=0.02)

    fig.suptitle(sample_id)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Render H&E overlays for top MVP samples.")
    parser.add_argument("--top-n", type=int, default=3)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sample = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_sample_level_scores.csv")
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    selected = choose_samples(sample, args.top_n)

    usecols = ["dataset_id", "sample_id", "x_pixel", "y_pixel"] + [score for score, _ in OVERLAY_SCORES]
    edge_qc_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    if edge_qc_path.exists():
        usecols_with_qc = usecols + ["edge_or_background_risk", "local_white_background_fraction", "nearest_out_tissue_distance_px"]
        spot = pd.read_csv(edge_qc_path, usecols=usecols_with_qc)
    else:
        spot = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv", usecols=usecols)

    rows: list[dict] = []
    errors: list[str] = []
    for _, selected_row in selected.iterrows():
        dataset_id = selected_row["dataset_id"]
        sample_id = selected_row["sample_id"]
        manifest_row = manifest[(manifest["dataset_id"] == dataset_id) & (manifest["sample_id"] == sample_id)]
        if manifest_row.empty:
            errors.append(f"Manifest row missing: {dataset_id} {sample_id}")
            continue
        manifest_row = manifest_row.iloc[0]
        image_path = Path(manifest_row["image_path"])
        scalefactors_path = Path(manifest_row["scalefactors_path"])
        sample_spots = spot[(spot["dataset_id"] == dataset_id) & (spot["sample_id"] == sample_id)]
        if not image_path.exists() or sample_spots.empty:
            errors.append(f"Missing image or spots: {dataset_id} {sample_id}")
            continue
        safe_sample = sample_id.replace("/", "_").replace("\\", "_")
        output_path = PROJECT_ROOT / "results/figures/mvp/he_overlays" / dataset_id / f"{safe_sample}_he_overlay.png"
        render_overlay(sample_id, image_path, scalefactors_path, sample_spots, output_path)
        rows.append(
            {
                "dataset_id": dataset_id,
                "sample_id": sample_id,
                "selection_reason": selected_row["selection_reason"],
                "overlay_path": str(output_path),
                "image_path": str(image_path),
            }
        )

    pd.DataFrame(rows).to_csv(PROJECT_ROOT / "results/tables/mvp_overlay_manifest.csv", index=False)
    status = "success" if not errors else "partial_success"
    write_status(
        "05_he_overlays",
        status,
        {
            "n_datasets_processed": int(selected["dataset_id"].nunique()) if not selected.empty else 0,
            "n_samples_processed": len(rows),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Open representative overlays and record whether high-score regions match plausible histology.",
                "Do not claim H&E predictability from these overlays alone.",
            ],
        },
    )
    print(f"Wrote {len(rows)} H&E overlays")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
