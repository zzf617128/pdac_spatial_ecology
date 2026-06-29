from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "results" / "tables"
REVISION = PROJECT / "results" / "revision_2026_06_29"
ANALYSIS_OUT = REVISION / "analysis_outputs"
SUPP_TABLES = REVISION / "supplementary_tables"
REPORTS = REVISION / "docs"
RNG = np.random.default_rng(20260629)

DATASETS = {
    "mvp": (TABLES / "mvp_spot_level_scores_with_edge_qc.csv", None, None, ["edge_or_background_risk"]),
    "gse272362": (TABLES / "gse272362_rds_spot_level_scores.csv", "GSE272362", "specimen_type", ["specimen_type"]),
    "gse235315": (TABLES / "gse235315_spot_level_scores.csv", "GSE235315", "specimen_type", ["specimen_type"]),
    "gse274557": (TABLES / "gse274557_full_spot_scores.csv", "GSE274557", "tissue", ["tissue", "treatment", "geo_accession"]),
}

USECOLS = [
    "dataset_id",
    "sample_id",
    "patient_id",
    "x_pixel",
    "y_pixel",
    "score_caf_myeloid_barrier",
    "score_immune_hub_core",
    "score_tumor_aggressive",
    "score_tumor_epithelial",
    "z_ifn_antigen_presentation",
    "z_mycaf",
    "z_pan_caf",
    "z_myeloid",
    "z_spp1_tam",
    "z_tgfb_pathway",
    "z_emt_invasion",
]


def available_usecols(path: Path, requested: list[str]) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        header = next(csv.reader(f))
    return [c for c in requested if c in header]


def normalize_specimen(value: object) -> str:
    lower = str(value).strip().lower()
    if "lymph" in lower or "ln" in lower:
        return "lymph_node_metastasis"
    if "liver" in lower or "hepatic" in lower:
        return "liver_metastasis"
    if "lung" in lower:
        return "lung_metastasis"
    if "peritoneal" in lower:
        return "peritoneal_metastasis"
    if "normal" in lower:
        return "normal_pancreas"
    if "primary" in lower or "pdac" in lower or "tumor" in lower or "tumour" in lower:
        return "primary_tumor"
    return str(value).strip() or "metadata_required"


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 30:
        return np.nan
    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def median_nn_scale(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return 1.0
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    val = float(np.nanmedian(dists[:, 1]))
    return val if np.isfinite(val) and val > 0 else 1.0


def nearest_dist(xy: np.ndarray, idx: np.ndarray, scale: float) -> np.ndarray:
    tree = cKDTree(xy[idx])
    dists, _ = tree.query(xy, k=1)
    return dists / scale


def anchor_scores(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        "CAF-myeloid combined": df["score_caf_myeloid_barrier"].to_numpy(float),
        "CAF-only": df[["z_mycaf", "z_pan_caf"]].mean(axis=1, skipna=True).to_numpy(float),
        "myeloid-only": df[["z_myeloid", "z_spp1_tam"]].mean(axis=1, skipna=True).to_numpy(float),
        "panCAF": df["z_pan_caf"].to_numpy(float),
        "SPP1/TAM": df["z_spp1_tam"].to_numpy(float),
        "immune-high": df["score_immune_hub_core"].to_numpy(float),
        "tumor-high": df["score_tumor_epithelial"].to_numpy(float),
    }


def target_scores(df: pd.DataFrame) -> dict[str, np.ndarray]:
    return {
        "IFN/MHC": df["z_ifn_antigen_presentation"].to_numpy(float),
        "immune-core": df["score_immune_hub_core"].to_numpy(float),
        "tumor-aggressive": df["score_tumor_aggressive"].to_numpy(float),
        "SPP1/TAM": df["z_spp1_tam"].to_numpy(float),
        "TGF-beta/EMT": df[["z_tgfb_pathway", "z_emt_invasion"]].mean(axis=1, skipna=True).to_numpy(float),
        "tumor epithelial": df["score_tumor_epithelial"].to_numpy(float),
    }


def analyze_sample(df: pd.DataFrame, dataset: str, tissue_site: str, n_perm: int) -> list[dict[str, object]]:
    if "edge_or_background_risk" in df.columns:
        risk = df["edge_or_background_risk"].astype(str).str.lower().isin(["true", "1", "yes"])
        df = df[~risk].copy()
    df = df[np.isfinite(df["x_pixel"]) & np.isfinite(df["y_pixel"])].copy()
    if len(df) < 100:
        return []
    xy = df[["x_pixel", "y_pixel"]].to_numpy(float)
    scale = median_nn_scale(xy)
    n_anchor = max(10, int(math.ceil(0.10 * len(df))))
    anchors = anchor_scores(df)
    targets = target_scores(df)
    rows: list[dict[str, object]] = []
    random_indices = [RNG.choice(len(df), size=n_anchor, replace=False) for _ in range(n_perm)]
    random_null_by_target: dict[str, np.ndarray] = {name: [] for name in targets}
    for ridx in random_indices:
        rd = nearest_dist(xy, ridx, scale)
        for target_name, tv in targets.items():
            random_null_by_target[target_name].append(safe_spearman(rd, tv))
    random_null_by_target = {
        target_name: np.array(values, dtype=float)
        for target_name, values in random_null_by_target.items()
    }
    sample_id = str(df["sample_id"].iloc[0])
    patient_id = str(df["patient_id"].iloc[0]) if "patient_id" in df.columns else ""

    for anchor_name, av in anchors.items():
        if np.isfinite(av).sum() < n_anchor:
            continue
        idx = np.argsort(av)[-n_anchor:]
        dist = nearest_dist(xy, idx, scale)
        for target_name, tv in targets.items():
            observed = safe_spearman(dist, tv)
            null = random_null_by_target[target_name]
            n_finite = int(np.isfinite(null).sum())
            median = float(np.nanmedian(null)) if n_finite else np.nan
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": sample_id,
                    "patient_id": patient_id,
                    "tissue_site": tissue_site,
                    "anchor_type": anchor_name,
                    "target_program": target_name,
                    "n_spots": len(df),
                    "n_anchor_spots": n_anchor,
                    "n_perm": n_perm,
                    "observed_rho": observed,
                    "random_median_rho": median,
                    "delta_vs_random": observed - median if np.isfinite(observed) and np.isfinite(median) else np.nan,
                    "random_p05": float(np.nanquantile(null, 0.05)) if n_finite else np.nan,
                    "random_p95": float(np.nanquantile(null, 0.95)) if n_finite else np.nan,
                    "empirical_p": (1 + int(np.nansum(null <= observed))) / (1 + n_finite) if np.isfinite(observed) and n_finite else np.nan,
                    "support": bool(observed < median) if np.isfinite(observed) and np.isfinite(median) else "",
                }
            )
    return rows


def load_dataset(key: str) -> pd.DataFrame:
    path, _, _, extras = DATASETS[key]
    cols = available_usecols(path, USECOLS + extras)
    return pd.read_csv(path, usecols=cols, encoding="utf-8-sig")


def summarize(detail: pd.DataFrame) -> pd.DataFrame:
    return (
        detail.groupby(["anchor_type", "target_program"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_delta=("delta_vs_random", "median"),
            median_observed_rho=("observed_rho", "median"),
            support_n=("support", lambda x: int(np.sum([v is True or str(v).lower() == "true" for v in x]))),
            support_fraction=("support", lambda x: float(np.mean([v is True or str(v).lower() == "true" for v in x]))),
            median_empirical_p=("empirical_p", "median"),
        )
        .reset_index()
    )


def write_report(summary: pd.DataFrame, n_perm: int) -> None:
    lines = [
        "# CAF-only, myeloid-only and CAF-myeloid anchor comparison",
        "",
        f"n_perm per sample and anchor: {n_perm}",
        "",
        "This analysis asks whether the combined CAF-myeloid core adds value beyond CAF-only or myeloid-only anchors.",
        "Negative delta indicates stronger target-program centering around the biological anchor than around same-size random anchors.",
        "",
        "## Key anchor summary",
        "",
    ]
    key = summary[
        summary["anchor_type"].isin(["CAF-myeloid combined", "CAF-only", "myeloid-only", "tumor-high", "immune-high"])
        & summary["target_program"].isin(["IFN/MHC", "immune-core", "tumor-aggressive", "SPP1/TAM", "TGF-beta/EMT", "tumor epithelial"])
    ].copy()
    for _, row in key.iterrows():
        lines.append(
            f"- {row['anchor_type']} -> {row['target_program']}: "
            f"median delta {row['median_delta']:.3f}, support {int(row['support_n'])}/{int(row['n_samples'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation rule",
            "",
            "If CAF-only and CAF-myeloid combined anchors perform similarly, the manuscript should frame the architecture as a CAF-dominant stromal core with myeloid enrichment rather than claiming that the combined score is uniquely superior.",
            "",
        ]
    )
    (REPORTS / "caf_myeloid_component_anchor_comparison_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--max-samples", type=int, default=0)
    args = parser.parse_args()

    ANALYSIS_OUT.mkdir(parents=True, exist_ok=True)
    SUPP_TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    counter = 0
    for key in args.datasets:
        path, fixed_label, specimen_col, _ = DATASETS[key]
        print(f"Reading {key}: {path}", flush=True)
        df = load_dataset(key)
        for sample_id, sdf in df.groupby("sample_id", sort=True):
            if args.max_samples and counter >= args.max_samples:
                break
            dataset = fixed_label or str(sdf["dataset_id"].iloc[0])
            tissue_site = normalize_specimen(sdf[specimen_col].iloc[0]) if specimen_col and specimen_col in sdf.columns else "metadata_required"
            counter += 1
            print(f"[{counter}] {dataset} {sample_id} {tissue_site} n={len(sdf)}", flush=True)
            rows.extend(analyze_sample(sdf, dataset, tissue_site, args.n_perm))
        if args.max_samples and counter >= args.max_samples:
            break

    detail = pd.DataFrame(rows)
    summary = summarize(detail)
    detail_path = ANALYSIS_OUT / "caf_myeloid_component_anchor_comparison_per_sample.csv"
    summary_path = ANALYSIS_OUT / "caf_myeloid_component_anchor_comparison_summary.csv"
    supp_path = SUPP_TABLES / "Supplementary_Table_4B_Anchor_Component_Comparison.csv"
    detail.to_csv(detail_path, index=False, encoding="utf-8")
    summary.to_csv(summary_path, index=False, encoding="utf-8")
    detail.to_csv(supp_path, index=False, encoding="utf-8")
    write_report(summary, args.n_perm)
    print(f"Wrote {detail_path}", flush=True)
    print(f"Wrote {summary_path}", flush=True)
    print(f"Wrote {supp_path}", flush=True)


if __name__ == "__main__":
    main()
