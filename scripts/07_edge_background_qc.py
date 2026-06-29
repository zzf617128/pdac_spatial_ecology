from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
from scipy.spatial import cKDTree


PROJECT_ROOT = Path(__file__).resolve().parents[1]
SCORE_COLS = [
    "score_caf_myeloid_barrier",
    "score_immune_hub_core",
    "score_tumor_aggressive",
    "z_ifn_antigen_presentation",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_positions(path: Path) -> pd.DataFrame:
    positions = pd.read_csv(path)
    if "barcode" not in positions.columns and positions.shape[1] == 6:
        positions = pd.read_csv(
            path,
            header=None,
            names=[
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ],
        )
    positions = positions.rename(
        columns={"pxl_col_in_fullres": "x_pixel", "pxl_row_in_fullres": "y_pixel"}
    )
    return positions


def read_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0))


def estimate_neighbor_distance(coords: np.ndarray) -> float:
    if len(coords) < 3:
        return float("nan")
    tree = cKDTree(coords)
    dist, _ = tree.query(coords, k=2)
    return float(np.nanmedian(dist[:, 1]))


def edge_distances(positions: pd.DataFrame) -> pd.DataFrame:
    in_tissue = positions[positions["in_tissue"].astype(int) == 1].copy()
    out_tissue = positions[positions["in_tissue"].astype(int) == 0].copy()
    in_coords = in_tissue[["x_pixel", "y_pixel"]].to_numpy(float)
    out_coords = out_tissue[["x_pixel", "y_pixel"]].to_numpy(float)
    nn_dist = estimate_neighbor_distance(in_coords)
    if len(out_coords) == 0:
        distances = np.full(len(in_tissue), np.nan)
    else:
        tree = cKDTree(out_coords)
        distances, _ = tree.query(in_coords, k=1)
    in_tissue["nearest_out_tissue_distance_px"] = distances
    in_tissue["median_neighbor_distance_px"] = nn_dist
    in_tissue["edge_proximal"] = distances <= (2.5 * nn_dist)
    return in_tissue[["barcode", "nearest_out_tissue_distance_px", "median_neighbor_distance_px", "edge_proximal"]]


def integral_image(mask: np.ndarray) -> np.ndarray:
    return np.pad(mask.astype(np.int32).cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)))


def window_fraction(integral: np.ndarray, x: np.ndarray, y: np.ndarray, radius: int) -> np.ndarray:
    height = integral.shape[0] - 1
    width = integral.shape[1] - 1
    x0 = np.clip(np.floor(x - radius).astype(int), 0, width)
    x1 = np.clip(np.ceil(x + radius).astype(int), 0, width)
    y0 = np.clip(np.floor(y - radius).astype(int), 0, height)
    y1 = np.clip(np.ceil(y + radius).astype(int), 0, height)
    area = np.maximum((x1 - x0) * (y1 - y0), 1)
    values = integral[y1, x1] - integral[y0, x1] - integral[y1, x0] + integral[y0, x0]
    return values / area


def background_fractions(image_path: Path, scale: float, positions: pd.DataFrame) -> pd.DataFrame:
    image = np.asarray(Image.open(image_path).convert("RGB"))
    brightness = image.mean(axis=2)
    channel_range = image.max(axis=2) - image.min(axis=2)
    white_mask = (brightness > 235) & (channel_range < 25)
    integ = integral_image(white_mask)
    x = positions["x_pixel"].to_numpy(float) * scale
    y = positions["y_pixel"].to_numpy(float) * scale
    radius = 16
    frac = window_fraction(integ, x, y, radius=radius)
    return pd.DataFrame({"barcode": positions["barcode"].to_numpy(), "local_white_background_fraction": frac})


def per_sample_qc(
    sample_scores: pd.DataFrame,
    positions_path: Path,
    image_path: Path,
    scalefactors_path: Path,
) -> tuple[pd.DataFrame, dict]:
    positions = read_positions(positions_path)
    edge = edge_distances(positions)
    bg = background_fractions(image_path, read_scale(scalefactors_path), positions)
    annotated = sample_scores.merge(edge, on="barcode", how="left").merge(bg, on="barcode", how="left")
    annotated["background_risk"] = annotated["local_white_background_fraction"] > 0.50
    annotated["edge_or_background_risk"] = annotated["edge_proximal"].fillna(False) | annotated[
        "background_risk"
    ].fillna(False)

    summary = {
        "dataset_id": sample_scores["dataset_id"].iloc[0],
        "sample_id": sample_scores["sample_id"].iloc[0],
        "n_spots": len(annotated),
        "fraction_edge_proximal": float(np.nanmean(annotated["edge_proximal"])),
        "fraction_background_risk": float(np.nanmean(annotated["background_risk"])),
        "fraction_edge_or_background_risk": float(np.nanmean(annotated["edge_or_background_risk"])),
    }
    safe = annotated[~annotated["edge_or_background_risk"]]
    risk = annotated[annotated["edge_or_background_risk"]]
    for score in SCORE_COLS:
        high_thr = np.nanpercentile(annotated[score], 90)
        high = annotated[annotated[score] >= high_thr]
        summary[f"{score}_high_fraction_risk"] = float(np.nanmean(high["edge_or_background_risk"]))
        summary[f"{score}_safe_p90"] = float(np.nanpercentile(safe[score], 90)) if len(safe) else np.nan
        summary[f"{score}_risk_p90"] = float(np.nanpercentile(risk[score], 90)) if len(risk) else np.nan
        summary[f"{score}_risk_minus_safe_p90"] = (
            summary[f"{score}_risk_p90"] - summary[f"{score}_safe_p90"]
            if np.isfinite(summary[f"{score}_safe_p90"]) and np.isfinite(summary[f"{score}_risk_p90"])
            else np.nan
        )
    return annotated, summary


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


def main() -> int:
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    spot = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv")
    scored_keys = spot[["dataset_id", "sample_id"]].drop_duplicates()
    manifest = manifest.merge(scored_keys, on=["dataset_id", "sample_id"], how="inner")

    annotated_rows: list[pd.DataFrame] = []
    summary_rows: list[dict] = []
    errors: list[str] = []
    for _, row in manifest.iterrows():
        sample_scores = spot[(spot["dataset_id"] == row["dataset_id"]) & (spot["sample_id"] == row["sample_id"])]
        try:
            annotated, summary = per_sample_qc(
                sample_scores,
                Path(row["coordinates_path"]),
                Path(row["image_path"]),
                Path(row["scalefactors_path"]),
            )
            annotated_rows.append(annotated)
            summary_rows.append(summary)
        except Exception as exc:
            errors.append(f"{row['dataset_id']} {row['sample_id']}: {exc}")

    if annotated_rows:
        pd.concat(annotated_rows, ignore_index=True).to_csv(
            PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv", index=False
        )
    pd.DataFrame(summary_rows).to_csv(
        PROJECT_ROOT / "results/tables/mvp_edge_qc_sample_summary.csv", index=False
    )

    status = "success" if not errors else "partial_success"
    write_status(
        "07_edge_background_qc",
        status,
        {
            "n_datasets_processed": int(manifest["dataset_id"].nunique()),
            "n_samples_processed": len(summary_rows),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Inspect results/tables/mvp_edge_qc_sample_summary.csv.",
                "If high-score regions are mostly risk spots, revise story toward tissue-interface/edge-aware analysis.",
            ],
        },
    )
    print(f"Computed edge/background QC for {len(summary_rows)} samples")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
