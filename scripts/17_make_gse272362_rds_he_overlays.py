from __future__ import annotations

import json
import re
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data/raw/GSE272362"
SCORES = [
    ("score_caf_myeloid_barrier", "CAF-myeloid"),
    ("z_ifn_antigen_presentation", "IFN/MHC"),
    ("score_tumor_aggressive", "tumor aggressive"),
    ("score_immune_hub_core", "immune core"),
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(stage: str, status: str, payload: dict) -> None:
    base = {
        "stage": stage,
        "status": status,
        "timestamp_utc": now_iso(),
        "n_errors": 0,
        "critical_errors": [],
        "noncritical_warnings": [],
        "next_manual_check": [],
    }
    base.update(payload)
    path = PROJECT_ROOT / f"results/logs/stage_{stage}_status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def bare_barcode(value: str) -> str:
    return str(value).split("_")[-1]


def read_positions(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    if "barcode" not in df.columns and df.shape[1] == 6:
        df = pd.read_csv(
            path,
            header=None,
            names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
        )
    return df[["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"]].copy()


def read_hires_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0))


def base_from_position_path(path: Path) -> str:
    return path.name[: -len("_tissue_positions.csv")]


def build_position_cache() -> dict[str, dict]:
    cache: dict[str, dict] = {}
    for pos_path in sorted(RAW_DIR.glob("*_tissue_positions.csv")):
        base = base_from_position_path(pos_path)
        positions = read_positions(pos_path)
        positions["bare_barcode"] = positions["barcode"].map(bare_barcode)
        cache[base] = {
            "positions_path": pos_path,
            "positions": positions,
            "image_path": RAW_DIR / f"{base}_tissue_hires_image.png",
            "scalefactors_path": RAW_DIR / f"{base}_scalefactors_json.json",
        }
    return cache


def match_sample_to_geo(sample_spots: pd.DataFrame, position_cache: dict[str, dict]) -> dict:
    query = sample_spots[["barcode", "x_pixel", "y_pixel"]].copy()
    query["bare_barcode"] = query["barcode"].map(bare_barcode)
    if len(query) > 400:
        query = query.sample(n=400, random_state=20260624)
    best: dict | None = None
    for base, payload in position_cache.items():
        pos = payload["positions"]
        merged = query.merge(pos, on="bare_barcode", how="inner")
        if merged.empty:
            continue
        exact = (
            (merged["x_pixel"].astype(float) - merged["pxl_col_in_fullres"].astype(float)).abs() <= 1
        ) & (
            (merged["y_pixel"].astype(float) - merged["pxl_row_in_fullres"].astype(float)).abs() <= 1
        )
        n_exact = int(exact.sum())
        n_in_tissue_exact = int((exact & (merged["in_tissue"].astype(int) == 1)).sum())
        candidate = {
            **payload,
            "geo_base": base,
            "n_barcode_matches": int(len(merged)),
            "n_exact_coordinate_matches": n_exact,
            "n_in_tissue_exact_matches": n_in_tissue_exact,
            "match_fraction": n_exact / max(1, len(query)),
        }
        if best is None or candidate["n_exact_coordinate_matches"] > best["n_exact_coordinate_matches"]:
            best = candidate
    if best is None or best["n_exact_coordinate_matches"] < 20:
        raise ValueError(f"Could not confidently match sample {sample_spots['sample_id'].iloc[0]} to GEO image.")
    return best


def select_representative_samples(stats: pd.DataFrame) -> pd.DataFrame:
    rows: list[pd.Series] = []
    rules = [
        ("primary_tumor", "tumor_aggressive", True, "primary_strong_tumor_caf_core"),
        ("primary_tumor", "ifn_mhc", True, "primary_strong_ifn_caf_core"),
        ("liver_metastasis", "ifn_mhc", True, "liver_strong_ifn_caf_core"),
        ("liver_metastasis", "tumor_aggressive", True, "liver_strong_tumor_caf_core"),
        ("lymph_node_metastasis", "tumor_aggressive", True, "lnm_tumor_caf_core"),
        ("lymph_node_metastasis", "immune_core", False, "lnm_immune_divergent"),
    ]
    used: set[str] = set()
    for specimen, target, ascending_delta, reason in rules:
        sub = stats[(stats["specimen_type"] == specimen) & (stats["target"] == target)].copy()
        if sub.empty:
            continue
        if ascending_delta:
            sub = sub.sort_values(["delta_vs_null_median", "empirical_p_more_negative"], ascending=[True, True])
        else:
            sub = sub.sort_values(["delta_vs_null_median", "observed_rho"], ascending=[False, False])
        chosen = None
        for _, row in sub.iterrows():
            if row["sample_id"] not in used:
                chosen = row.copy()
                break
        if chosen is None:
            chosen = sub.iloc[0].copy()
        chosen["selection_reason"] = reason
        used.add(str(chosen["sample_id"]))
        rows.append(chosen)
    return pd.DataFrame(rows)


def add_caf_core_and_distance(spots: pd.DataFrame) -> pd.DataFrame:
    spots = spots.copy()
    caf = spots["score_caf_myeloid_barrier"].to_numpy(float)
    threshold = np.nanpercentile(caf, 90)
    spots["is_caf_core_top10"] = caf >= threshold
    return spots


def render_overlay(sample_id: str, specimen_type: str, reason: str, spots: pd.DataFrame, match: dict, output_path: Path) -> None:
    image_path = Path(match["image_path"])
    image = Image.open(image_path).convert("RGB")
    scale = read_hires_scale(Path(match["scalefactors_path"]))
    spots = add_caf_core_and_distance(spots)
    x = spots["x_pixel"].to_numpy(float) * scale
    y = spots["y_pixel"].to_numpy(float) * scale
    core = spots["is_caf_core_top10"].to_numpy(bool)

    fig, axes = plt.subplots(1, len(SCORES) + 1, figsize=(5 * (len(SCORES) + 1), 5.2), constrained_layout=True)
    axes[0].imshow(image)
    axes[0].scatter(x, y, s=3, c="#222222", alpha=0.15, linewidths=0)
    axes[0].scatter(x[core], y[core], s=14, facecolors="none", edgecolors="#e31a1c", linewidths=0.55, alpha=0.85)
    axes[0].set_title("H&E + CAF core")
    axes[0].set_xticks([])
    axes[0].set_yticks([])

    for ax, (score_col, title) in zip(axes[1:], SCORES):
        ax.imshow(image)
        values = spots[score_col].to_numpy(float)
        lo = np.nanpercentile(values, 2)
        hi = np.nanpercentile(values, 98)
        sc = ax.scatter(
            x,
            y,
            c=values,
            s=6,
            cmap="viridis",
            alpha=0.82,
            linewidths=0,
            vmin=lo,
            vmax=hi,
        )
        ax.scatter(x[core], y[core], s=15, facecolors="none", edgecolors="#e31a1c", linewidths=0.45, alpha=0.75)
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)

    fig.suptitle(
        f"{sample_id} | {specimen_type} | {reason} | {match['geo_base']} | exact matches={match['n_exact_coordinate_matches']}",
        fontsize=13,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> int:
    spot = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv")
    stats = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_random_core_permutation_sample_stats.csv")
    selected = select_representative_samples(stats)
    position_cache = build_position_cache()

    rows: list[dict] = []
    errors: list[str] = []
    for _, row in selected.iterrows():
        sample_id = row["sample_id"]
        sample_spots = spot[spot["sample_id"] == sample_id].copy()
        try:
            match = match_sample_to_geo(sample_spots, position_cache)
            safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(sample_id))
            output_path = (
                PROJECT_ROOT
                / "results/figures/mvp/gse272362_rds/he_overlays"
                / f"{safe_id}_{row['selection_reason']}_he_overlay.png"
            )
            render_overlay(sample_id, row["specimen_type"], row["selection_reason"], sample_spots, match, output_path)
            rows.append(
                {
                    "sample_id": sample_id,
                    "specimen_type": row["specimen_type"],
                    "selection_reason": row["selection_reason"],
                    "target": row["target"],
                    "observed_rho": row["observed_rho"],
                    "delta_vs_null_median": row["delta_vs_null_median"],
                    "empirical_p_more_negative": row["empirical_p_more_negative"],
                    "geo_base": match["geo_base"],
                    "n_exact_coordinate_matches": match["n_exact_coordinate_matches"],
                    "image_path": str(match["image_path"]),
                    "overlay_path": str(output_path),
                }
            )
        except Exception as exc:
            errors.append(f"{sample_id}: {exc}")

    out_table = PROJECT_ROOT / "results/tables/gse272362_rds_overlay_manifest.csv"
    pd.DataFrame(rows).to_csv(out_table, index=False)
    status = "success" if not errors else "partial_success"
    write_status(
        "17_gse272362_rds_he_overlays",
        status,
        {
            "n_samples_processed": len(rows),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Inspect GSE272362 primary, liver, and lymph node overlays.",
                "Confirm red-ring CAF cores correspond to plausible tissue compartments.",
                "For lymph node metastasis overlays, check whether immune-rich areas follow lymphoid architecture rather than CAF cores.",
            ],
        },
    )
    print(f"Wrote {len(rows)} GSE272362 H&E overlays")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

