from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def strip_score(name: str) -> str:
    return name.replace("score_", "").replace("_", " ")


def main() -> int:
    sample = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_sample_level_scores.csv")
    corr = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_sample_spot_correlations.csv")
    edge = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_edge_qc_sample_summary.csv")
    out_dir = PROJECT_ROOT / "results/figures/mvp/summary"
    out_dir.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 2, figsize=(13, 10), constrained_layout=True)

    ax = axes[0, 0]
    cols = [
        "fraction_caf_myeloid_barrier_z_gt1",
        "fraction_immune_hub_core_z_gt1",
        "fraction_tumor_aggressive_z_gt1",
        "fraction_ifn_antigen_presentation_z_gt1",
    ]
    data = [sample[col].dropna().to_numpy() for col in cols if col in sample]
    labels = [col.replace("fraction_", "").replace("_z_gt1", "").replace("_", "\n") for col in cols if col in sample]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(labels) + 1), labels)
    ax.set_ylabel("Fraction of spots with score > 1")
    ax.set_title("High-score spatial footprint")
    ax.tick_params(axis="x", labelsize=8)

    ax = axes[0, 1]
    corr_cols = [
        "rho_barrier_vs_immune_core",
        "safe_rho_barrier_vs_immune_core",
        "rho_barrier_vs_ifn_mhc",
        "safe_rho_barrier_vs_ifn_mhc",
        "rho_barrier_vs_tumor_aggressive",
        "safe_rho_barrier_vs_tumor_aggressive",
    ]
    data = [corr[col].dropna().to_numpy() for col in corr_cols if col in corr]
    labels = [col.replace("safe_", "safe\n").replace("rho_barrier_vs_", "").replace("_", "\n") for col in corr_cols if col in corr]
    ax.boxplot(data, showfliers=False)
    ax.set_xticks(range(1, len(labels) + 1), labels)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_ylabel("Within-sample Spearman rho")
    ax.set_title("CAF-myeloid co-localization persists after edge QC")
    ax.tick_params(axis="x", labelsize=8)

    ax = axes[1, 0]
    ax.scatter(
        edge["fraction_edge_or_background_risk"],
        edge["score_caf_myeloid_barrier_high_fraction_risk"],
        s=22,
        alpha=0.75,
        label="CAF-myeloid top spots",
    )
    ax.plot([0, 1], [0, 1], color="gray", lw=1, ls="--")
    ax.set_xlabel("All spots: edge/background risk fraction")
    ax.set_ylabel("Top CAF-myeloid spots: risk fraction")
    ax.set_title("High CAF-myeloid is not globally edge-dominated")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    story = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_story_screening_table.csv")
    story = story.sort_values("mvp_score_0_to_3")
    ax.barh(story["candidate_story"], story["mvp_score_0_to_3"], color="#4C78A8")
    ax.set_xlim(0, 3)
    ax.set_xlabel("MVP score")
    ax.set_title("Candidate story ranking")
    ax.tick_params(axis="y", labelsize=8)

    fig.suptitle("PDAC Spatial Ecology MVP Summary", fontsize=14)
    png = out_dir / "mvp_summary.png"
    pdf = out_dir / "mvp_summary.pdf"
    fig.savefig(png, dpi=180)
    fig.savefig(pdf)
    plt.close(fig)

    write_status(
        "08_mvp_summary_figures",
        "success",
        {
            "n_datasets_processed": int(sample["dataset_id"].nunique()),
            "n_samples_processed": len(sample),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": [],
            "next_manual_check": [
                "Inspect results/figures/mvp/summary/mvp_summary.png.",
                "Decide whether to deepen CAF-myeloid niche story before H&E modeling.",
            ],
        },
    )
    print(f"Wrote {png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
