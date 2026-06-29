from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def bh_fdr(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    ranked = np.empty_like(p)
    n = len(p)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        idx = order[i]
        rank = i + 1
        q = min(prev, p[idx] * n / rank)
        ranked[idx] = q
        prev = q
    return ranked.tolist()


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
    gradients = pd.read_csv(PROJECT_ROOT / "results/tables/caf_myeloid_niche_gradient_stats.csv")
    rows: list[dict] = []
    for target, group in gradients.groupby("target"):
        values = group["rho_distance_to_caf_core"].dropna().to_numpy(float)
        if len(values) == 0:
            continue
        test = wilcoxon(values, zero_method="wilcox", alternative="two-sided")
        rows.append(
            {
                "analysis": "distance_to_caf_myeloid_core",
                "target": target,
                "analysis_level": "sample_section",
                "test_name": "Wilcoxon signed-rank vs zero",
                "n_samples": len(values),
                "median_rho": float(np.nanmedian(values)),
                "mean_rho": float(np.nanmean(values)),
                "n_negative": int((values < 0).sum()),
                "n_positive": int((values > 0).sum()),
                "effect_direction": "enriched_near_caf_core" if np.nanmedian(values) < 0 else "enriched_far_from_caf_core",
                "p_value": float(test.pvalue),
            }
        )
    q_values = bh_fdr([row["p_value"] for row in rows])
    for row, q in zip(rows, q_values):
        row["q_value"] = q
    out = pd.DataFrame(rows).sort_values("q_value")
    out.to_csv(PROJECT_ROOT / "results/tables/caf_myeloid_niche_gradient_model_stats.csv", index=False)
    write_status(
        "10_caf_myeloid_niche_stats",
        "success",
        {
            "n_datasets_processed": int(gradients["dataset_id"].nunique()),
            "n_samples_processed": int(gradients["sample_id"].nunique()),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": [
                "Sample sections are used as units; patient-level inference remains blocked until verified patient mapping."
            ],
            "next_manual_check": [
                "Inspect results/tables/caf_myeloid_niche_gradient_model_stats.csv.",
            ],
        },
    )
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

