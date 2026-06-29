from __future__ import annotations

import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def safe_spearman(x: pd.Series, y: pd.Series) -> float:
    mask = x.notna() & y.notna()
    if mask.sum() < 30:
        return math.nan
    if x[mask].nunique() < 2 or y[mask].nunique() < 2:
        return math.nan
    return float(spearmanr(x[mask], y[mask]).correlation)


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


def compute_spot_correlations() -> pd.DataFrame:
    usecols = [
        "dataset_id",
        "sample_id",
        "score_immune_hub_core",
        "score_immune_hub_maturity",
        "score_caf_myeloid_barrier",
        "score_tumor_aggressive",
        "score_ifn_antigen_presentation",
        "z_ifn_antigen_presentation",
    ]
    edge_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    if edge_path.exists():
        spot = pd.read_csv(edge_path, usecols=usecols + ["edge_or_background_risk"])
    else:
        spot = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv", usecols=usecols)
    rows: list[dict] = []
    for (dataset_id, sample_id), group in spot.groupby(["dataset_id", "sample_id"], sort=False):
        safe_group = group
        if "edge_or_background_risk" in group.columns:
            safe_group = group[~group["edge_or_background_risk"].fillna(False).astype(bool)]
        rows.append(
            {
                "dataset_id": dataset_id,
                "sample_id": sample_id,
                "n_spots": len(group),
                "n_safe_spots": len(safe_group),
                "rho_barrier_vs_immune_core": safe_spearman(
                    group["score_caf_myeloid_barrier"], group["score_immune_hub_core"]
                ),
                "rho_barrier_vs_ifn_mhc": safe_spearman(
                    group["score_caf_myeloid_barrier"], group["z_ifn_antigen_presentation"]
                ),
                "rho_barrier_vs_tumor_aggressive": safe_spearman(
                    group["score_caf_myeloid_barrier"], group["score_tumor_aggressive"]
                ),
                "rho_immune_core_vs_maturity": safe_spearman(
                    group["score_immune_hub_core"], group["score_immune_hub_maturity"]
                ),
                "rho_tumor_aggressive_vs_ifn_mhc": safe_spearman(
                    group["score_tumor_aggressive"], group["z_ifn_antigen_presentation"]
                ),
                "safe_rho_barrier_vs_immune_core": safe_spearman(
                    safe_group["score_caf_myeloid_barrier"], safe_group["score_immune_hub_core"]
                ),
                "safe_rho_barrier_vs_ifn_mhc": safe_spearman(
                    safe_group["score_caf_myeloid_barrier"], safe_group["z_ifn_antigen_presentation"]
                ),
                "safe_rho_barrier_vs_tumor_aggressive": safe_spearman(
                    safe_group["score_caf_myeloid_barrier"], safe_group["score_tumor_aggressive"]
                ),
            }
        )
    return pd.DataFrame(rows)


def summarize_metrics(sample: pd.DataFrame, corr: pd.DataFrame) -> dict:
    merged = sample.merge(corr, on=["dataset_id", "sample_id"], how="left")
    return {
        "n_samples": len(sample),
        "n_gse282302": int((sample["dataset_id"] == "GSE282302").sum()),
        "n_gse274103": int((sample["dataset_id"] == "GSE274103").sum()),
        "median_barrier_fraction_z_gt1": float(
            np.nanmedian(sample["fraction_caf_myeloid_barrier_z_gt1"])
        ),
        "median_immune_core_fraction_z_gt1": float(
            np.nanmedian(sample["fraction_immune_hub_core_z_gt1"])
        ),
        "median_tumor_aggressive_fraction_z_gt1": float(
            np.nanmedian(sample["fraction_tumor_aggressive_z_gt1"])
        ),
        "median_rho_barrier_vs_immune_core": float(
            np.nanmedian(merged["rho_barrier_vs_immune_core"])
        ),
        "median_rho_barrier_vs_ifn_mhc": float(np.nanmedian(merged["rho_barrier_vs_ifn_mhc"])),
        "median_rho_barrier_vs_tumor_aggressive": float(
            np.nanmedian(merged["rho_barrier_vs_tumor_aggressive"])
        ),
        "median_rho_immune_core_vs_maturity": float(
            np.nanmedian(merged["rho_immune_core_vs_maturity"])
        ),
        "median_safe_rho_barrier_vs_immune_core": float(
            np.nanmedian(merged["safe_rho_barrier_vs_immune_core"])
        )
        if "safe_rho_barrier_vs_immune_core" in merged
        else math.nan,
        "median_safe_rho_barrier_vs_ifn_mhc": float(
            np.nanmedian(merged["safe_rho_barrier_vs_ifn_mhc"])
        )
        if "safe_rho_barrier_vs_ifn_mhc" in merged
        else math.nan,
        "median_safe_rho_barrier_vs_tumor_aggressive": float(
            np.nanmedian(merged["safe_rho_barrier_vs_tumor_aggressive"])
        )
        if "safe_rho_barrier_vs_tumor_aggressive" in merged
        else math.nan,
        "n_detectable_immune_hub": int((sample["has_detectable_immune_hub"].astype(str) == "true").sum()),
        "n_mature_like_hub": int((sample["has_mature_like_hub"].astype(str) == "true").sum()),
    }


def build_story_table(metrics: dict) -> pd.DataFrame:
    rows: list[dict] = []

    barrier_support = 0
    if metrics["median_barrier_fraction_z_gt1"] > 0.03:
        barrier_support += 1
    if metrics["median_rho_barrier_vs_tumor_aggressive"] > 0.1:
        barrier_support += 1
    if metrics["median_rho_barrier_vs_immune_core"] < 0:
        barrier_support += 1
    rows.append(
        {
            "candidate_story": "CAF-myeloid inflammatory stromal niche with limited immune organization",
            "mvp_score_0_to_3": barrier_support,
            "evidence": (
                f"median barrier high-spot fraction={metrics['median_barrier_fraction_z_gt1']:.3f}; "
                f"median rho barrier~immune={metrics['median_rho_barrier_vs_immune_core']:.3f}; "
                f"median safe rho barrier~immune={metrics['median_safe_rho_barrier_vs_immune_core']:.3f}; "
                f"median safe rho barrier~tumor_aggressive={metrics['median_safe_rho_barrier_vs_tumor_aggressive']:.3f}"
            ),
            "caveat": "Current MVP shows positive barrier~immune/IFN correlation, so do not frame this as simple immune exclusion yet.",
            "next_action": "Deepen with H&E/spatial overlays, local gradients, and tests for mature immune organization.",
        }
    )

    immune_support = 0
    if metrics["n_detectable_immune_hub"] >= metrics["n_samples"] * 0.25:
        immune_support += 1
    if metrics["n_mature_like_hub"] < metrics["n_detectable_immune_hub"]:
        immune_support += 1
    if metrics["median_rho_immune_core_vs_maturity"] > 0.1:
        immune_support += 1
    rows.append(
        {
            "candidate_story": "Immune hub maturation arrest",
            "mvp_score_0_to_3": immune_support,
            "evidence": (
                f"detectable immune-hub-like signal in {metrics['n_detectable_immune_hub']}/{metrics['n_samples']} samples; "
                f"mature-like signal in {metrics['n_mature_like_hub']}/{metrics['n_samples']} samples; "
                f"median rho immune_core~maturity={metrics['median_rho_immune_core_vs_maturity']:.3f}"
            ),
            "caveat": "MVP detection uses signatures only; do not call mature TLS without B/T/FDC/GC morphology review.",
            "next_action": "Keep as secondary unless manual overlays show true lymphoid aggregates.",
        }
    )

    rows.append(
        {
            "candidate_story": "Primary-to-metastatic spatial ecology remodeling",
            "mvp_score_0_to_3": 1,
            "evidence": "GSE272362 is downloaded and has paired images/coordinates, but GEO supplementary lacks spot-level count matrices.",
            "caveat": "Needs Zenodo 11.6GB RDS, author processed object, or SRA processing before independent scoring.",
            "next_action": "Request/download Zenodo PDAC_Updated.rds only if this story becomes primary.",
        }
    )

    rows.append(
        {
            "candidate_story": "Post-neoadjuvant residual spatial ecology",
            "mvp_score_0_to_3": 2 if metrics["n_gse282302"] >= 100 else 1,
            "evidence": f"GSE282302 contributes {metrics['n_gse282302']} directly scored ST-H&E samples.",
            "caveat": "Treatment context and response metadata must be resolved before any therapy-response claim.",
            "next_action": "Curate GSE282302 sample/patient metadata and compare residual ecotype axes.",
        }
    )

    rows.append(
        {
            "candidate_story": "H&E-readable PDAC spatial ecotype",
            "mvp_score_0_to_3": 1,
            "evidence": "All directly scored samples have paired tissue images; model not trained yet.",
            "caveat": "No claim until patient/sample split H&E feature model is validated.",
            "next_action": "Train a small exploratory ResNet/timm feature model only after story axis is chosen.",
        }
    )

    rows.append(
        {
            "candidate_story": "Neural niche-associated immune exclusion",
            "mvp_score_0_to_3": 0,
            "evidence": "Neural niche has not been targeted in MVP; GSE202740 not yet added.",
            "caveat": "Only add if neural_schwann signal is strong in current scored samples or if project pivots.",
            "next_action": "Defer.",
        }
    )
    return pd.DataFrame(rows).sort_values("mvp_score_0_to_3", ascending=False)


def write_report(story: pd.DataFrame, metrics: dict) -> None:
    top = story.iloc[0]
    lines = [
        "# MVP Decision Report",
        "",
        f"Last updated UTC: {now_iso()}",
        "",
        "## Current Data Status",
        "",
        f"- Directly scored Visium samples: {metrics['n_samples']}",
        f"- GSE282302 scored samples: {metrics['n_gse282302']}",
        f"- GSE274103 scored samples: {metrics['n_gse274103']}",
        "- GSE272362 is downloaded but lacks spot-level count matrices in GEO supplementary files.",
        "",
        "## Preliminary Recommendation",
        "",
        f"Current top MVP direction: **{top['candidate_story']}**.",
        "",
        "This is a provisional direction. It should not be treated as a final claim until manual H&E/spatial overlay review and a stronger spatial-gradient/permutation analysis are complete.",
        "",
        "## Key Metrics",
        "",
    ]
    for key, value in metrics.items():
        lines.append(f"- {key}: {value}")
    lines.extend(["", "## Candidate Story Ranking", ""])
    table_cols = ["candidate_story", "mvp_score_0_to_3", "evidence", "caveat", "next_action"]
    lines.append("| candidate_story | mvp_score_0_to_3 | evidence | caveat | next_action |")
    lines.append("|---|---:|---|---|---|")
    for _, row in story[table_cols].iterrows():
        cells = []
        for col in table_cols:
            value = str(row[col]).replace("|", "/").replace("\n", " ")
            cells.append(value)
        lines.append(f"| {' | '.join(cells)} |")
    lines.extend(
        [
            "",
            "## Immediate Next Step",
            "",
            "1. Curate GSE282302 sample and patient metadata.",
        "2. Generate representative H&E + score overlays for high CAF-myeloid, high-immune, and high-tumor-aggressive samples.",
        "3. If the overlays are biologically plausible, deepen the CAF-myeloid inflammatory stromal niche / residual ecology story.",
            "4. Do not invest in GSE272362 primary-metastasis scoring unless the Zenodo RDS or equivalent count matrix is obtained.",
            "",
        ]
    )
    (PROJECT_ROOT / "docs/mvp_decision_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    sample_path = PROJECT_ROOT / "results/tables/mvp_sample_level_scores.csv"
    spot_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv"
    if not sample_path.exists() or not spot_path.exists():
        raise FileNotFoundError("Run scripts/03_mvp_score_visium.py before story screening.")

    sample = pd.read_csv(sample_path)
    corr = compute_spot_correlations()
    metrics = summarize_metrics(sample, corr)
    story = build_story_table(metrics)

    corr.to_csv(PROJECT_ROOT / "results/tables/mvp_sample_spot_correlations.csv", index=False)
    story.to_csv(PROJECT_ROOT / "results/tables/mvp_story_screening_table.csv", index=False)
    write_report(story, metrics)
    write_status(
        "04_mvp_story_screen",
        "success",
        {
            "n_datasets_processed": int(sample["dataset_id"].nunique()),
            "n_samples_processed": len(sample),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": [
                "GSE272362 requires additional spot-level count data before independent scoring."
            ],
            "next_manual_check": [
                "Inspect docs/mvp_decision_report.md.",
                "Review high-ranking sample overlays before choosing the final story.",
            ],
        },
    )
    print(f"Wrote story screening table for {len(sample)} samples")
    print(f"Top candidate: {story.iloc[0]['candidate_story']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
