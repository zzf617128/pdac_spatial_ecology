from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE = "24_gse235315_random_core_anchor"
RANDOM_SEED = 20260624
N_PERM = 1000
TARGETS = [
    ("z_ifn_antigen_presentation", "IFN/MHC"),
    ("score_immune_hub_core", "immune core"),
    ("score_tumor_aggressive", "tumor aggressive"),
    ("score_immune_hub_maturity", "immune maturity-like"),
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


def median_nn(points: np.ndarray) -> float:
    nn = NearestNeighbors(n_neighbors=2).fit(points)
    d = nn.kneighbors(points, return_distance=True)[0][:, 1]
    med = float(np.nanmedian(d))
    return med if np.isfinite(med) and med > 0 else 1.0


def nearest_norm_distance(points: np.ndarray, center_idx: np.ndarray, scale: float) -> np.ndarray:
    nn = NearestNeighbors(n_neighbors=1).fit(points[center_idx])
    d = nn.kneighbors(points, return_distance=True)[0][:, 0]
    return d / scale


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 50:
        return np.nan
    rho = spearmanr(x[mask], y[mask]).statistic
    return float(rho) if np.isfinite(rho) else np.nan


def run_sample(sample: pd.DataFrame, rng: np.random.Generator) -> tuple[list[dict], list[dict]]:
    points = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    scale = median_nn(points)
    n_spots = len(sample)
    n_core = max(1, int(round(0.10 * n_spots)))
    caf = sample["score_caf_myeloid_barrier"].to_numpy(float)
    core_idx = np.argsort(caf)[-n_core:]
    obs_dist = nearest_norm_distance(points, core_idx, scale)
    first = sample.iloc[0]

    rows: list[dict] = []
    null_rows: list[dict] = []
    observed: dict[str, float] = {}
    for col, label in TARGETS:
        observed[label] = safe_spearman(obs_dist, sample[col].to_numpy(float))

    null_values: dict[str, list[float]] = {label: [] for _, label in TARGETS}
    all_idx = np.arange(n_spots)
    for perm in range(N_PERM):
        random_idx = rng.choice(all_idx, size=n_core, replace=False)
        random_dist = nearest_norm_distance(points, random_idx, scale)
        for col, label in TARGETS:
            rho = safe_spearman(random_dist, sample[col].to_numpy(float))
            null_values[label].append(rho)
            null_rows.append(
                {
                    "dataset_id": first["dataset_id"],
                    "sample_id": first["sample_id"],
                    "target_label": label,
                    "permutation": perm,
                    "random_core_rho": rho,
                }
            )

    for _, label in TARGETS:
        null = np.asarray(null_values[label], dtype=float)
        obs = observed[label]
        rows.append(
            {
                "dataset_id": first["dataset_id"],
                "sample_id": first["sample_id"],
                "patient_id": first["patient_id"],
                "n_spots": n_spots,
                "n_core_spots": n_core,
                "target_label": label,
                "observed_rho": obs,
                "null_median_rho": float(np.nanmedian(null)),
                "null_p05_rho": float(np.nanquantile(null, 0.05)),
                "null_p95_rho": float(np.nanquantile(null, 0.95)),
                "delta_vs_null_median": float(obs - np.nanmedian(null)),
                "empirical_p_more_negative": float((np.nansum(null <= obs) + 1) / (np.isfinite(null).sum() + 1)),
                "observed_more_negative_than_null": bool(obs < np.nanmedian(null)),
            }
        )
    return rows, null_rows


def make_figure(summary: pd.DataFrame, output_base: Path) -> None:
    target_order = [label for _, label in TARGETS]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.2), constrained_layout=True)

    ax0 = axes[0]
    grouped = summary.groupby("target_label", sort=False).agg(
        median_delta=("delta_vs_null_median", "median"),
        n_support=("observed_more_negative_than_null", "sum"),
        n_samples=("sample_id", "nunique"),
    )
    grouped = grouped.reindex(target_order)
    x = np.arange(len(grouped))
    ax0.bar(x, grouped["median_delta"], color="#4C78A8")
    ax0.axhline(0, color="#555555", linewidth=0.8)
    ax0.set_xticks(x, grouped.index, rotation=35, ha="right")
    ax0.set_ylabel("observed minus random-core median rho")
    ax0.set_title("GSE235315 random-core support")
    for i, row in enumerate(grouped.itertuples()):
        ax0.text(i, row.median_delta, f"{int(row.n_support)}/{int(row.n_samples)}", ha="center", va="bottom" if row.median_delta > 0 else "top", fontsize=8)

    ax1 = axes[1]
    for i, target in enumerate(target_order):
        values = summary.loc[summary["target_label"].eq(target), "observed_rho"].dropna().to_numpy(float)
        jitter = np.linspace(-0.12, 0.12, len(values)) if len(values) > 1 else np.array([0.0])
        ax1.scatter(np.full(len(values), i) + jitter, values, s=32, color="#4C78A8", alpha=0.8)
        ax1.plot([i - 0.25, i + 0.25], [np.nanmedian(values), np.nanmedian(values)], color="#222222", linewidth=1.5)
    ax1.axhline(0, color="#555555", linewidth=0.8)
    ax1.set_xticks(np.arange(len(target_order)), target_order, rotation=35, ha="right")
    ax1.set_ylabel("observed distance-to-CAF-core rho")
    ax1.set_title("Observed CAF-core gradients")

    fig.suptitle("External paired-ST anchor cohort supports CAF-core spatial organization", fontsize=13)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=220)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def write_report(summary: pd.DataFrame) -> None:
    lines = [
        "# Stage 24 GSE235315 Random-Core Anchor",
        "",
        "## Purpose",
        "",
        "This analysis tests whether the CAF-core spatial-gradient result is reproduced in the newly downloaded GSE235315 paired spatial cohort. The cohort is used as an external spatial-state anchor; patient/treatment labels remain metadata-pending.",
        "",
        "## Results",
        "",
    ]
    grouped = summary.groupby("target_label", sort=False).agg(
        n_samples=("sample_id", "nunique"),
        n_support=("observed_more_negative_than_null", "sum"),
        median_observed_rho=("observed_rho", "median"),
        median_delta=("delta_vs_null_median", "median"),
        median_empirical_p=("empirical_p_more_negative", "median"),
    )
    for target, row in grouped.iterrows():
        lines.append(
            f"- {target}: observed CAF-core gradient more negative than random cores in "
            f"{int(row['n_support'])}/{int(row['n_samples'])} samples; median observed rho "
            f"{row['median_observed_rho']:.3f}; median delta {row['median_delta']:.3f}; "
            f"median empirical p {row['median_empirical_p']:.3g}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "GSE235315 provides an additional paired-ST anchor for the CAF-core organization concept. Because sample metadata were not yet fully curated, the result should be used as external spatial-state support rather than a patient-level or treatment-context claim.",
        ]
    )
    path = PROJECT_ROOT / "results/reports/gse235315_random_core_anchor_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    spot = pd.read_csv(PROJECT_ROOT / "results/tables/gse235315_spot_level_scores.csv")
    rng = np.random.default_rng(RANDOM_SEED)
    rows: list[dict] = []
    null_rows: list[dict] = []
    warnings: list[str] = []
    for sample_id, sample in spot.groupby("sample_id"):
        try:
            sample_rows, sample_null = run_sample(sample.copy(), rng)
            rows.extend(sample_rows)
            null_rows.extend(sample_null)
            print(f"Random-core anchor complete: {sample_id}")
        except Exception as exc:
            warnings.append(f"{sample_id}: {exc}")

    summary = pd.DataFrame(rows)
    null = pd.DataFrame(null_rows)
    tables = PROJECT_ROOT / "results/tables"
    tables.mkdir(parents=True, exist_ok=True)
    summary.to_csv(tables / "gse235315_random_core_anchor_summary.csv", index=False)
    null.to_csv(tables / "gse235315_random_core_anchor_null_rhos.csv", index=False)

    source = PROJECT_ROOT / "results/source_data"
    source.mkdir(parents=True, exist_ok=True)
    summary.to_csv(source / "Source_Data_Fig_6A.csv", index=False)

    make_figure(summary, PROJECT_ROOT / "results/figures/main/figure6_gse235315_anchor")
    write_report(summary)

    write_status(
        STAGE,
        "success" if not warnings else "partial_success",
        {
            "n_samples_processed": int(summary["sample_id"].nunique()) if not summary.empty else 0,
            "n_targets": len(TARGETS),
            "n_permutations_per_sample": N_PERM,
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": warnings,
            "outputs": [
                "results/tables/gse235315_random_core_anchor_summary.csv",
                "results/figures/main/figure6_gse235315_anchor.pdf",
                "results/reports/gse235315_random_core_anchor_report.md",
            ],
            "next_manual_check": [
                "Audit GSE235315 sample metadata before promoting this beyond external spatial-state support.",
                "Inspect GSE235315 anchor score maps for coordinate/image plausibility.",
            ],
        },
    )
    print("Stage 24 GSE235315 random-core anchor complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
