from __future__ import annotations

import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data/raw/GSE235315"
STAGE = "23_gse235315_anchor"


def load_stage03():
    path = PROJECT_ROOT / "scripts/03_mvp_score_visium.py"
    spec = importlib.util.spec_from_file_location("stage03_mvp_score_visium", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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


def discover_samples() -> pd.DataFrame:
    tar_map: dict[str, str] = {}
    for tar_path in RAW_DIR.glob("GSM*_processed.tar.gz"):
        parts = tar_path.name.split("_")
        if len(parts) >= 2:
            tar_map[parts[1]] = parts[0]

    rows: list[dict] = []
    for sample_dir in sorted(RAW_DIR.rglob("SS*_processed")):
        ss_id = sample_dir.name.replace("_processed", "")
        gsm_id = tar_map.get(ss_id, "")
        rows.append(
            {
                "dataset_id": "GSE235315",
                "sample_id": f"{gsm_id}_{ss_id}" if gsm_id else ss_id,
                "gsm_id": gsm_id,
                "patient_id": ss_id,
                "specimen_type": "pdac_paired_st_metadata_pending",
                "treatment_context": "metadata_required",
                "expression_path": str(sample_dir / "filtered_feature_bc_matrix.h5"),
                "coordinates_path": str(sample_dir / "tissue_positions_list.csv"),
                "image_path": str(sample_dir / "tissue_hires_image.png"),
                "scalefactors_path": str(sample_dir / "scalefactors_json.json"),
                "sample_dir": str(sample_dir),
            }
        )
    return pd.DataFrame(rows)


def write_sample_figure(spot_scores: pd.DataFrame, output_path: Path) -> None:
    scores = [
        ("score_caf_myeloid_barrier", "CAF-myeloid"),
        ("score_tumor_aggressive", "tumor aggressive"),
        ("z_ifn_antigen_presentation", "IFN/MHC"),
        ("score_immune_hub_core", "immune core"),
    ]
    fig, axes = plt.subplots(1, len(scores), figsize=(4.2 * len(scores), 4.0), constrained_layout=True)
    for ax, (score, title) in zip(axes, scores):
        values = spot_scores[score].to_numpy(float)
        sc = ax.scatter(
            spot_scores["x_pixel"],
            -spot_scores["y_pixel"],
            c=values,
            s=5,
            cmap="viridis",
            linewidths=0,
            vmin=np.nanpercentile(values, 2),
            vmax=np.nanpercentile(values, 98),
        )
        ax.set_title(title)
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle(str(spot_scores["sample_id"].iloc[0]))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180)
    plt.close(fig)


def main() -> int:
    stage03 = load_stage03()
    signatures = stage03.parse_signatures(PROJECT_ROOT / "config/signatures.yaml")
    samples = discover_samples()
    warnings: list[str] = []
    all_spots: list[pd.DataFrame] = []
    coverage_rows: list[dict] = []
    sample_rows: list[dict] = []

    for _, row in samples.iterrows():
        expression_path = Path(row["expression_path"])
        coordinates_path = Path(row["coordinates_path"])
        if not expression_path.exists() or not coordinates_path.exists():
            warnings.append(f"Missing expression or coordinates for {row['sample_id']}")
            continue
        try:
            counts, genes, barcodes = stage03.read_10x_h5(expression_path)
            positions = stage03.read_positions(coordinates_path)
            spot_scores, cov = stage03.score_signatures(
                counts,
                genes,
                barcodes,
                positions,
                signatures,
                row["dataset_id"],
                row["sample_id"],
                row["patient_id"],
            )
            spot_scores.insert(3, "specimen_type", row["specimen_type"])
            all_spots.append(spot_scores)
            coverage_rows.extend(cov)
            summary = stage03.sample_summary(spot_scores, row)
            sample_rows.append(summary)
            safe = row["sample_id"].replace("/", "_").replace("\\", "_")
            write_sample_figure(
                spot_scores,
                PROJECT_ROOT / "results/figures/mvp/gse235315_anchor" / f"{safe}_scores.png",
            )
            print(f"Scored GSE235315 anchor sample: {row['sample_id']} ({len(spot_scores)} spots)")
        except Exception as exc:
            warnings.append(f"{row['sample_id']}: {exc}")

    out_dir = PROJECT_ROOT / "results/tables"
    out_dir.mkdir(parents=True, exist_ok=True)
    if all_spots:
        pd.concat(all_spots, ignore_index=True).to_csv(out_dir / "gse235315_spot_level_scores.csv", index=False)
    samples.to_csv(PROJECT_ROOT / "metadata/gse235315_sample_manifest.csv", index=False)
    pd.DataFrame(sample_rows).to_csv(out_dir / "gse235315_sample_level_scores.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(out_dir / "gse235315_signature_gene_coverage.csv", index=False)

    write_status(
        STAGE,
        "success" if not warnings else "partial_success",
        {
            "n_samples_discovered": int(len(samples)),
            "n_samples_processed": int(len(sample_rows)),
            "n_spots_scored": int(sum(len(df) for df in all_spots)),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": warnings,
            "outputs": [
                "metadata/gse235315_sample_manifest.csv",
                "results/tables/gse235315_spot_level_scores.csv",
                "results/tables/gse235315_sample_level_scores.csv",
                "results/tables/gse235315_signature_gene_coverage.csv",
            ],
            "next_manual_check": [
                "Audit GSE235315 sample metadata before using patient- or treatment-level labels.",
                "Use this cohort as external cell-state/spatial-state support unless metadata are fully curated.",
            ],
        },
    )
    print("Stage 23 GSE235315 anchor scoring complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
