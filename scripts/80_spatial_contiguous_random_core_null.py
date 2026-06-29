from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import rankdata


PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "results" / "tables"
REVISION = PROJECT / "results" / "revision_2026_06_29"
ANALYSIS_OUT = REVISION / "analysis_outputs"
SUPP_TABLES = REVISION / "supplementary_tables"
REPORTS = REVISION / "docs"

RNG = np.random.default_rng(20260629)

BASE_USECOLS = [
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
    "z_spp1_tam",
    "z_tgfb_pathway",
    "z_emt_invasion",
    "z_mycaf",
    "z_pan_caf",
]


DATASETS = {
    "mvp": {
        "path": TABLES / "mvp_spot_level_scores_with_edge_qc.csv",
        "extra_usecols": ["edge_or_background_risk"],
        "specimen_col": None,
        "dataset_label": None,
    },
    "gse272362": {
        "path": TABLES / "gse272362_rds_spot_level_scores.csv",
        "extra_usecols": ["specimen_type"],
        "specimen_col": "specimen_type",
        "dataset_label": "GSE272362",
    },
    "gse235315": {
        "path": TABLES / "gse235315_spot_level_scores.csv",
        "extra_usecols": ["specimen_type"],
        "specimen_col": "specimen_type",
        "dataset_label": "GSE235315",
    },
    "gse274557": {
        "path": TABLES / "gse274557_full_spot_scores.csv",
        "extra_usecols": ["tissue", "treatment", "geo_accession"],
        "specimen_col": "tissue",
        "dataset_label": "GSE274557",
    },
}


def available_usecols(path: Path, requested: list[str]) -> list[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as f:
        header = next(csv.reader(f))
    return [c for c in requested if c in header]


def normalize_specimen(value: object) -> str:
    text = str(value).strip()
    lower = text.lower()
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
    return text or "metadata_required"


def finite_rank(values: np.ndarray) -> np.ndarray:
    ranks = np.full(values.shape, np.nan, dtype=float)
    mask = np.isfinite(values)
    if mask.sum() >= 3:
        ranks[mask] = rankdata(values[mask], method="average")
    return ranks


def safe_rank_corr(rank_x: np.ndarray, rank_y: np.ndarray) -> float:
    mask = np.isfinite(rank_x) & np.isfinite(rank_y)
    if mask.sum() < 30:
        return np.nan
    x = rank_x[mask]
    y = rank_y[mask]
    x = x - x.mean()
    y = y - y.mean()
    denom = math.sqrt(float(np.dot(x, x) * np.dot(y, y)))
    if denom <= 0:
        return np.nan
    return float(np.dot(x, y) / denom)


def nearest_distance_rank(xy: np.ndarray, core_idx: np.ndarray) -> np.ndarray:
    tree = cKDTree(xy[core_idx])
    dists, _ = tree.query(xy, k=1)
    return finite_rank(dists.astype(float))


def build_knn(xy: np.ndarray, k: int = 7) -> list[np.ndarray]:
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=min(k + 1, len(xy)))
    neighbors: list[np.ndarray] = []
    for row_i, row in enumerate(np.atleast_2d(idx)):
        neighbors.append(np.array([j for j in row if j != row_i], dtype=int))
    return neighbors


def contiguous_region(neighbors: list[np.ndarray], n_core: int) -> np.ndarray:
    n = len(neighbors)
    seed = int(RNG.integers(0, n))
    selected = {seed}
    frontier = set(int(x) for x in neighbors[seed])
    while len(selected) < n_core:
        if frontier:
            options = np.fromiter(frontier, dtype=int)
            chosen = int(options[int(RNG.integers(0, len(options)))])
            frontier.discard(chosen)
        else:
            remaining = np.array(list(set(range(n)) - selected), dtype=int)
            if len(remaining) == 0:
                break
            chosen = int(remaining[int(RNG.integers(0, len(remaining)))])
        if chosen in selected:
            continue
        selected.add(chosen)
        for nb in neighbors[chosen]:
            nb_int = int(nb)
            if nb_int not in selected:
                frontier.add(nb_int)
    return np.fromiter(selected, dtype=int)


def target_arrays(df: pd.DataFrame) -> dict[str, np.ndarray]:
    tgfb_emt = df[["z_tgfb_pathway", "z_emt_invasion"]].mean(axis=1, skipna=True).to_numpy(float)
    mycaf_matrix = df[["z_mycaf", "z_pan_caf"]].mean(axis=1, skipna=True).to_numpy(float)
    return {
        "IFN/MHC": df["z_ifn_antigen_presentation"].to_numpy(float),
        "immune-core": df["score_immune_hub_core"].to_numpy(float),
        "tumor-aggressive": df["score_tumor_aggressive"].to_numpy(float),
        "SPP1/TAM": df["z_spp1_tam"].to_numpy(float),
        "TGF-beta/EMT": tgfb_emt,
        "tumor epithelial": df["score_tumor_epithelial"].to_numpy(float),
        "myCAF/matrix": mycaf_matrix,
    }


def analyze_sample(df: pd.DataFrame, dataset_label: str, specimen_type: str, n_perm: int) -> list[dict[str, object]]:
    df = df[np.isfinite(df["x_pixel"]) & np.isfinite(df["y_pixel"]) & np.isfinite(df["score_caf_myeloid_barrier"])].copy()
    if "edge_or_background_risk" in df.columns:
        risk = df["edge_or_background_risk"].astype(str).str.lower().isin(["true", "1", "yes"])
        df = df[~risk].copy()
    if len(df) < 100:
        return []

    xy = df[["x_pixel", "y_pixel"]].to_numpy(float)
    caf = df["score_caf_myeloid_barrier"].to_numpy(float)
    n_core = max(10, int(math.ceil(0.10 * len(df))))
    core_idx = np.argsort(caf)[-n_core:]
    observed_dist_rank = nearest_distance_rank(xy, core_idx)

    targets = target_arrays(df)
    target_ranks = {name: finite_rank(vals) for name, vals in targets.items()}
    observed_rhos = {name: safe_rank_corr(observed_dist_rank, ranks) for name, ranks in target_ranks.items()}

    neighbors = build_knn(xy, k=7)
    null = {name: [] for name in targets}
    for _ in range(n_perm):
        random_idx = contiguous_region(neighbors, n_core)
        dist_rank = nearest_distance_rank(xy, random_idx)
        for name, ranks in target_ranks.items():
            null[name].append(safe_rank_corr(dist_rank, ranks))

    rows: list[dict[str, object]] = []
    sample_id = str(df["sample_id"].iloc[0])
    patient_id = str(df["patient_id"].iloc[0]) if "patient_id" in df.columns else ""
    for name, values in null.items():
        arr = np.array(values, dtype=float)
        observed = observed_rhos[name]
        n_finite = int(np.isfinite(arr).sum())
        null_median = float(np.nanmedian(arr)) if n_finite else np.nan
        empirical_p = (1 + int(np.nansum(arr <= observed))) / (1 + n_finite) if np.isfinite(observed) and n_finite else np.nan
        rows.append(
            {
                "dataset": dataset_label,
                "sample_id": sample_id,
                "patient_id": patient_id,
                "tissue_site": specimen_type,
                "target_program": name,
                "null_type": "spatially_contiguous_random_core",
                "n_spots": len(df),
                "n_core_spots": n_core,
                "n_perm": n_perm,
                "observed_rho": observed,
                "null_median": null_median,
                "null_p05": float(np.nanquantile(arr, 0.05)) if n_finite else np.nan,
                "null_p95": float(np.nanquantile(arr, 0.95)) if n_finite else np.nan,
                "delta": observed - null_median if np.isfinite(observed) and np.isfinite(null_median) else np.nan,
                "empirical_p": empirical_p,
                "support": bool(observed < null_median) if np.isfinite(observed) and np.isfinite(null_median) else "",
            }
        )
    return rows


def load_dataset(key: str) -> pd.DataFrame:
    info = DATASETS[key]
    path = info["path"]
    usecols = available_usecols(path, BASE_USECOLS + info["extra_usecols"])
    return pd.read_csv(path, usecols=usecols, encoding="utf-8-sig")


def summarize(rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return (
        df.groupby(["dataset", "tissue_site", "target_program", "null_type"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_observed_rho=("observed_rho", "median"),
            median_null=("null_median", "median"),
            median_delta=("delta", "median"),
            support_n=("support", lambda x: int(np.sum([v is True or str(v).lower() == "true" for v in x]))),
            support_fraction=("support", lambda x: float(np.mean([v is True or str(v).lower() == "true" for v in x]))),
            median_empirical_p=("empirical_p", "median"),
        )
        .reset_index()
    )


def write_report(summary: pd.DataFrame, n_perm: int, datasets: list[str]) -> None:
    lines = [
        "# Spatially contiguous random-core null",
        "",
        f"n_perm per sample: {n_perm}",
        f"datasets: {', '.join(datasets)}",
        "",
        "This analysis tests whether CAF-core gradients exceed random contiguous tissue regions of the same size.",
        "Negative observed-minus-null delta indicates stronger enrichment near the observed CAF-myeloid core than near a matched contiguous random core.",
        "",
        "## Summary",
        "",
    ]
    if summary.empty:
        lines.append("No results generated.")
    else:
        focus = summary[summary["target_program"].isin(["IFN/MHC", "immune-core", "tumor-aggressive", "SPP1/TAM", "TGF-beta/EMT"])]
        for _, row in focus.iterrows():
            lines.append(
                f"- {row['dataset']} {row['tissue_site']} {row['target_program']}: "
                f"median delta {row['median_delta']:.3f}, support {int(row['support_n'])}/{int(row['n_samples'])}."
            )
    lines.extend(
        [
            "",
            "## Claim boundary",
            "",
            "This is still an observational spatial null model. It controls for arbitrary contiguous tissue regions but does not establish causal CAF-to-immune or CAF-to-tumor signaling.",
            "",
        ]
    )
    (REPORTS / "stronger_null_contiguous_random_core_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=["mvp", "gse272362", "gse235315", "gse274557"], choices=sorted(DATASETS))
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--max-samples", type=int, default=0)
    args = parser.parse_args()

    ANALYSIS_OUT.mkdir(parents=True, exist_ok=True)
    SUPP_TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    all_rows: list[dict[str, object]] = []
    sample_counter = 0
    for key in args.datasets:
        info = DATASETS[key]
        print(f"Reading {key}: {info['path']}")
        df = load_dataset(key)
        for sample_id, sdf in df.groupby("sample_id", sort=True):
            if args.max_samples and sample_counter >= args.max_samples:
                break
            dataset_label = info["dataset_label"] or str(sdf["dataset_id"].iloc[0])
            specimen_col = info["specimen_col"]
            specimen = normalize_specimen(sdf[specimen_col].iloc[0]) if specimen_col and specimen_col in sdf.columns else "metadata_required"
            sample_counter += 1
            print(f"[{sample_counter}] {dataset_label} {sample_id} {specimen} n={len(sdf)}")
            rows = analyze_sample(sdf, dataset_label, specimen, args.n_perm)
            all_rows.extend(rows)
        if args.max_samples and sample_counter >= args.max_samples:
            break

    per_sample = pd.DataFrame(all_rows)
    summary = summarize(all_rows)
    per_sample_path = ANALYSIS_OUT / "stronger_null_contiguous_random_core_per_sample.csv"
    summary_path = ANALYSIS_OUT / "stronger_null_contiguous_random_core_summary.csv"
    supp_path = SUPP_TABLES / "Supplementary_Table_4_Stronger_Null_Sensitivity.csv"
    per_sample.to_csv(per_sample_path, index=False, encoding="utf-8")
    summary.to_csv(summary_path, index=False, encoding="utf-8")
    per_sample.to_csv(supp_path, index=False, encoding="utf-8")
    write_report(summary, args.n_perm, args.datasets)
    print(f"Wrote {per_sample_path}")
    print(f"Wrote {summary_path}")
    print(f"Wrote {supp_path}")


if __name__ == "__main__":
    main()
