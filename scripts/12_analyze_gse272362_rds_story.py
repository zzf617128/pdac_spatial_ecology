from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
TARGET_SCORES = {
    "immune_core": "score_immune_hub_core",
    "ifn_mhc": "z_ifn_antigen_presentation",
    "tumor_aggressive": "score_tumor_aggressive",
    "immune_maturity": "score_immune_hub_maturity",
}


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


def bh_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    q = np.full(len(p), np.nan)
    mask = np.isfinite(p)
    if not mask.any():
        return q.tolist()
    idx = np.where(mask)[0]
    order = idx[np.argsort(p[idx])]
    ranked = p[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    q[order] = np.clip(adjusted, 0, 1)
    return q.tolist()


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 30:
        return np.nan
    if len(np.unique(x[mask])) < 2 or len(np.unique(y[mask])) < 2:
        return np.nan
    xr = pd.Series(x[mask]).rank(method="average").to_numpy(float)
    yr = pd.Series(y[mask]).rank(method="average").to_numpy(float)
    return float(np.corrcoef(xr, yr)[0, 1])


def median_neighbor_distance(coords: np.ndarray) -> float:
    if len(coords) < 3:
        return np.nan
    sample = coords
    if len(coords) > 1200:
        rng = np.random.default_rng(20260624)
        sample = coords[rng.choice(len(coords), size=1200, replace=False)]
    nearest = nearest_distance(sample, sample, chunk_size=256, exclude_self=True)
    return float(np.nanmedian(nearest))


def nearest_distance(query: np.ndarray, reference: np.ndarray, chunk_size: int = 512, exclude_self: bool = False) -> np.ndarray:
    out = np.empty(len(query), dtype=float)
    for start in range(0, len(query), chunk_size):
        chunk = query[start : start + chunk_size]
        diff = chunk[:, None, :] - reference[None, :, :]
        d2 = np.sum(diff * diff, axis=2)
        if exclude_self and len(query) == len(reference):
            rows = np.arange(start, min(start + chunk_size, len(query)))
            d2[np.arange(len(rows)), rows] = np.inf
        out[start : start + chunk_size] = np.sqrt(np.nanmin(d2, axis=1))
    return out


def normal_two_sided_p_from_z(z: float) -> float:
    if not np.isfinite(z):
        return np.nan
    return float(math.erfc(abs(z) / math.sqrt(2)))


def mann_whitney_approx_p(x: np.ndarray, y: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return np.nan
    ranks = pd.Series(np.concatenate([x, y])).rank(method="average").to_numpy(float)
    r1 = np.sum(ranks[:n1])
    u1 = r1 - n1 * (n1 + 1) / 2
    mean_u = n1 * n2 / 2
    sd_u = math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    if sd_u == 0:
        return np.nan
    return normal_two_sided_p_from_z((u1 - mean_u) / sd_u)


def exact_sign_test_p(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    values = values[values != 0]
    n = len(values)
    if n == 0:
        return np.nan
    k = int(min(np.sum(values > 0), np.sum(values < 0)))
    prob = sum(math.comb(n, i) for i in range(k + 1)) / (2**n)
    return float(min(1.0, 2 * prob))


def normalize_specimen_type(values: pd.Series) -> pd.Series:
    text = values.fillna("metadata_required").astype(str).str.lower()
    out = pd.Series("metadata_required", index=values.index, dtype=object)
    out[text.str.contains("normal", regex=False)] = "normal_pancreas"
    out[text.str.contains("primary|tumou?r|pdac", regex=True)] = "primary_tumor"
    out[text.str.contains("liver|hepatic", regex=True)] = "liver_metastasis"
    out[text.str.contains("lymph|node|ln", regex=True)] = "lymph_node_metastasis"
    return out


def summarize_by_sample_and_specimen(spot: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    score_cols = [
        "score_caf_myeloid_barrier",
        "score_immune_hub_core",
        "score_immune_hub_maturity",
        "score_tumor_aggressive",
        "z_ifn_antigen_presentation",
    ]
    for (sample_id, specimen_type), group in spot.groupby(["sample_id", "specimen_type"], sort=False):
        row = {
            "dataset_id": "GSE272362",
            "sample_id": sample_id,
            "specimen_type": specimen_type,
            "n_spots": len(group),
            "median_counts": float(np.nanmedian(group["n_counts"])),
            "median_genes": float(np.nanmedian(group["n_genes"])),
        }
        for col in score_cols:
            if col in group:
                label = col.replace("score_", "").replace("z_", "")
                row[f"mean_{label}"] = float(np.nanmean(group[col]))
                row[f"p90_{label}"] = float(np.nanpercentile(group[col], 90))
                row[f"fraction_{label}_gt1"] = float(np.nanmean(group[col] > 1.0))
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_specimen_contrasts(sample_summary: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "p90_caf_myeloid_barrier",
        "fraction_caf_myeloid_barrier_gt1",
        "p90_immune_hub_core",
        "fraction_immune_hub_core_gt1",
        "p90_immune_hub_maturity",
        "fraction_immune_hub_maturity_gt1",
        "p90_tumor_aggressive",
        "fraction_tumor_aggressive_gt1",
        "p90_ifn_antigen_presentation",
        "fraction_ifn_antigen_presentation_gt1",
    ]
    rows: list[dict] = []
    for metric in metrics:
        if metric not in sample_summary:
            continue
        for specimen, group in sample_summary.groupby("specimen_type", sort=False):
            values = group[metric].dropna().to_numpy(float)
            rows.append(
                {
                    "contrast_type": "specimen_summary",
                    "metric": metric,
                    "specimen_type": specimen,
                    "comparison": "",
                    "n_samples": len(values),
                    "median": float(np.nanmedian(values)) if len(values) else np.nan,
                    "iqr": float(np.nanpercentile(values, 75) - np.nanpercentile(values, 25)) if len(values) else np.nan,
                    "p_value": np.nan,
                }
            )
        primary = sample_summary[sample_summary["specimen_type"] == "primary_tumor"][metric].dropna().to_numpy(float)
        for specimen in ["liver_metastasis", "lymph_node_metastasis", "normal_pancreas"]:
            other = sample_summary[sample_summary["specimen_type"] == specimen][metric].dropna().to_numpy(float)
            p_value = np.nan
            if len(primary) >= 2 and len(other) >= 2:
                p_value = mann_whitney_approx_p(primary, other)
            rows.append(
                {
                    "contrast_type": "primary_vs_other",
                    "metric": metric,
                    "specimen_type": specimen,
                    "comparison": f"primary_tumor_vs_{specimen}",
                    "n_samples": len(other),
                    "median": float(np.nanmedian(other)) if len(other) else np.nan,
                    "iqr": float(np.nanpercentile(other, 75) - np.nanpercentile(other, 25)) if len(other) else np.nan,
                    "p_value": p_value,
                }
            )
    contrast = pd.DataFrame(rows)
    contrast["q_value"] = bh_adjust(contrast["p_value"].tolist()) if not contrast.empty else []
    return contrast


def analyze_sample_gradient(group: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    group = group.copy()
    group = group[np.isfinite(group["x_pixel"]) & np.isfinite(group["y_pixel"])]
    if len(group) < 100:
        return [], []
    coords = group[["x_pixel", "y_pixel"]].to_numpy(float)
    nn = median_neighbor_distance(coords)
    if not np.isfinite(nn) or nn <= 0:
        return [], []
    caf = group["score_caf_myeloid_barrier"].to_numpy(float)
    threshold = np.nanpercentile(caf, 90)
    core = group[caf >= threshold]
    if len(core) < 10:
        return [], []
    distances = nearest_distance(coords, core[["x_pixel", "y_pixel"]].to_numpy(float))
    group["distance_to_caf_core_ring"] = distances / nn
    group["caf_core_bin"] = pd.cut(
        group["distance_to_caf_core_ring"],
        bins=[-0.01, 1.5, 3, 6, np.inf],
        labels=["core_0_1.5", "near_1.5_3", "mid_3_6", "far_gt6"],
    )
    base = {
        "dataset_id": "GSE272362",
        "sample_id": group["sample_id"].iloc[0],
        "specimen_type": group["specimen_type"].iloc[0],
        "n_spots": len(group),
        "n_caf_core_spots": len(core),
        "caf_core_threshold": float(threshold),
        "median_neighbor_distance_px": float(nn),
    }
    gradient_rows: list[dict] = []
    bin_rows: list[dict] = []
    for label, score_col in TARGET_SCORES.items():
        rho = safe_spearman(group["distance_to_caf_core_ring"].to_numpy(float), group[score_col].to_numpy(float))
        gradient_rows.append({**base, "target": label, "rho_distance_to_caf_core": rho})
        for bin_name, bin_df in group.groupby("caf_core_bin", observed=True):
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


def summarize_gradient_tests(gradients: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    for (specimen_type, target), group in gradients.groupby(["specimen_type", "target"], sort=False):
        values = group["rho_distance_to_caf_core"].dropna().to_numpy(float)
        p_value = np.nan
        if len(values) >= 5 and np.any(values != 0):
            p_value = exact_sign_test_p(values)
        rows.append(
            {
                "dataset_id": "GSE272362",
                "specimen_type": specimen_type,
                "target": target,
                "n_samples": len(values),
                "median_rho_distance_to_caf_core": float(np.nanmedian(values)) if len(values) else np.nan,
                "n_negative": int(np.sum(values < 0)) if len(values) else 0,
                "n_positive": int(np.sum(values > 0)) if len(values) else 0,
                "p_value_sign_test_vs_zero": p_value,
            }
        )
    out = pd.DataFrame(rows)
    out["q_value"] = bh_adjust(out["p_value_sign_test_vs_zero"].tolist()) if not out.empty else []
    return out


def make_figure(sample_summary: pd.DataFrame, gradients: pd.DataFrame, bins: pd.DataFrame) -> None:
    # Figure generation is intentionally handled elsewhere when matplotlib is available.
    return None


def main() -> int:
    spot_path = PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv"
    if not spot_path.exists():
        write_status(
            "12_gse272362_rds_story",
            "blocked",
            {
                "critical_errors": [f"Missing {spot_path}. Run scripts/11_score_gse272362_pdac_updated_rds.R first."],
                "n_errors": 1,
            },
        )
        print(f"Missing {spot_path}")
        print("Run scripts/11_score_gse272362_pdac_updated_rds.R first.")
        return 1

    spot = pd.read_csv(spot_path)
    spot["specimen_type"] = normalize_specimen_type(spot.get("specimen_type", pd.Series(index=spot.index)))
    sample_summary = summarize_by_sample_and_specimen(spot)
    specimen_contrasts = summarize_specimen_contrasts(sample_summary)

    gradient_rows: list[dict] = []
    bin_rows: list[dict] = []
    for _, group in spot.groupby(["sample_id", "specimen_type"], sort=False):
        gradients, bins = analyze_sample_gradient(group)
        gradient_rows.extend(gradients)
        bin_rows.extend(bins)
    gradients = pd.DataFrame(gradient_rows)
    bins = pd.DataFrame(bin_rows)
    gradient_tests = summarize_gradient_tests(gradients) if not gradients.empty else pd.DataFrame()

    out_dir = PROJECT_ROOT / "results/tables"
    sample_summary.to_csv(out_dir / "gse272362_rds_sample_specimen_summary.csv", index=False)
    specimen_contrasts.to_csv(out_dir / "gse272362_rds_specimen_contrasts.csv", index=False)
    gradients.to_csv(out_dir / "gse272362_rds_caf_myeloid_gradient_stats.csv", index=False)
    bins.to_csv(out_dir / "gse272362_rds_caf_myeloid_distance_bins.csv", index=False)
    gradient_tests.to_csv(out_dir / "gse272362_rds_caf_myeloid_gradient_model_stats.csv", index=False)
    if not sample_summary.empty and not gradients.empty and not bins.empty:
        make_figure(sample_summary, gradients, bins)

    write_status(
        "12_gse272362_rds_story",
        "success",
        {
            "n_samples": int(sample_summary["sample_id"].nunique()) if not sample_summary.empty else 0,
            "n_specimen_types": int(sample_summary["specimen_type"].nunique()) if not sample_summary.empty else 0,
            "n_gradient_samples": int(gradients["sample_id"].nunique()) if not gradients.empty else 0,
            "noncritical_warnings": [
                "Specimen labels are inferred from RDS metadata columns and must be manually audited before manuscript claims.",
                "No edge/background image QC is applied to GSE272362 RDS-derived scores unless exported coordinates/images support it.",
            ],
            "next_manual_check": [
                "Open metadata/gse272362_rds_metadata_columns.csv and confirm which column encodes tissue site.",
                "Inspect results/figures/mvp/gse272362_rds/gse272362_specimen_ecology.png.",
            ],
        },
    )
    print(f"Analyzed {sample_summary['sample_id'].nunique()} GSE272362 samples")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
