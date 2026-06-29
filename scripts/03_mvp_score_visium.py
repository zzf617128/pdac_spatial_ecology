from __future__ import annotations

import argparse
import ast
import csv
import json
import math
import sys
from datetime import datetime, timezone
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse


PROJECT_ROOT = Path(__file__).resolve().parents[1]

COMPOSITES = {
    "immune_hub_core": ["b_cell", "t_cell", "dc_apc", "tls_chemokine"],
    "immune_hub_maturity": ["fdc_gc_like", "plasma_cell", "ifn_antigen_presentation"],
    "caf_myeloid_barrier": ["pan_caf", "myeloid", "spp1_tam", "tgfb_pathway", "icaf"],
    "tumor_aggressive": ["pdac_basal_like", "emt_invasion", "hypoxia", "proliferation"],
}

PLOT_SCORES = [
    "score_immune_hub_core",
    "score_immune_hub_maturity",
    "score_caf_myeloid_barrier",
    "score_tumor_aggressive",
    "z_ifn_antigen_presentation",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def parse_signatures(path: Path) -> dict[str, list[str]]:
    signatures: dict[str, list[str]] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if ": [" not in line:
            continue
        key, value = line.split(":", 1)
        genes = [gene.strip().upper() for gene in value.strip().strip("[]").split(",") if gene.strip()]
        signatures[key.strip()] = genes
    return signatures


def decode_array(values) -> list[str]:
    return [item.decode("utf-8") if isinstance(item, bytes) else str(item) for item in values]


def read_10x_h5(path: Path) -> tuple[sparse.csc_matrix, list[str], list[str]]:
    with h5py.File(path, "r") as handle:
        matrix = handle["matrix"]
        data = matrix["data"][()]
        indices = matrix["indices"][()]
        indptr = matrix["indptr"][()]
        shape = tuple(matrix["shape"][()])
        barcodes = decode_array(matrix["barcodes"][()])
        gene_names = decode_array(matrix["features"]["name"][()])
    x = sparse.csc_matrix((data, indices, indptr), shape=shape)
    return x, gene_names, barcodes


def zscore(values: np.ndarray) -> np.ndarray:
    values = values.astype(float)
    mean = np.nanmean(values)
    std = np.nanstd(values)
    if not np.isfinite(std) or std == 0:
        return np.zeros_like(values, dtype=float)
    return (values - mean) / std


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
    rename = {
        "pxl_col_in_fullres": "x_pixel",
        "pxl_row_in_fullres": "y_pixel",
    }
    positions = positions.rename(columns=rename)
    required = {"barcode", "in_tissue", "array_row", "array_col", "x_pixel", "y_pixel"}
    missing = required.difference(positions.columns)
    if missing:
        raise ValueError(f"Missing position columns in {path}: {sorted(missing)}")
    return positions


def score_signatures(
    counts: sparse.csc_matrix,
    gene_names: list[str],
    barcodes: list[str],
    positions: pd.DataFrame,
    signatures: dict[str, list[str]],
    dataset_id: str,
    sample_id: str,
    patient_id: str,
) -> tuple[pd.DataFrame, list[dict]]:
    gene_to_indices: dict[str, list[int]] = {}
    for idx, gene in enumerate(gene_names):
        gene_to_indices.setdefault(gene.upper(), []).append(idx)

    barcode_to_idx = {barcode: idx for idx, barcode in enumerate(barcodes)}
    in_tissue = positions[positions["in_tissue"].astype(int) == 1].copy()
    in_tissue["matrix_idx"] = in_tissue["barcode"].map(barcode_to_idx)
    in_tissue = in_tissue.dropna(subset=["matrix_idx"])
    spot_indices = in_tissue["matrix_idx"].astype(int).to_numpy()
    x = counts[:, spot_indices]

    n_counts = np.asarray(x.sum(axis=0)).ravel()
    n_genes = np.diff(x.tocsc().indptr)
    scale = np.divide(10000.0, n_counts, out=np.zeros_like(n_counts, dtype=float), where=n_counts > 0)

    spot_scores = pd.DataFrame(
        {
            "dataset_id": dataset_id,
            "sample_id": sample_id,
            "patient_id": patient_id,
            "barcode": in_tissue["barcode"].to_numpy(),
            "array_row": in_tissue["array_row"].to_numpy(),
            "array_col": in_tissue["array_col"].to_numpy(),
            "x_pixel": in_tissue["x_pixel"].to_numpy(),
            "y_pixel": in_tissue["y_pixel"].to_numpy(),
            "n_counts": n_counts,
            "n_genes": n_genes,
        }
    )

    coverage_rows: list[dict] = []
    for signature_name, genes in signatures.items():
        gene_indices: list[int] = []
        present_genes: list[str] = []
        missing_genes: list[str] = []
        for gene in genes:
            indices = gene_to_indices.get(gene, [])
            if indices:
                gene_indices.extend(indices)
                present_genes.append(gene)
            else:
                missing_genes.append(gene)
        if len(set(present_genes)) < 3:
            values = np.full(x.shape[1], np.nan)
            reliable = False
        else:
            sub = x[gene_indices, :].astype(float).tocsc()
            # Log-normalize only the signature genes to keep the MVP memory-light.
            sub = sub.multiply(scale)
            dense = sub.toarray()
            dense = np.log1p(dense)
            values = np.asarray(dense.mean(axis=0)).ravel()
            reliable = True
        spot_scores[f"score_{signature_name}"] = values
        spot_scores[f"z_{signature_name}"] = zscore(values) if reliable else np.nan
        coverage_rows.append(
            {
                "dataset_id": dataset_id,
                "sample_id": sample_id,
                "signature": signature_name,
                "n_genes_defined": len(genes),
                "n_genes_present": len(set(present_genes)),
                "n_genes_missing": len(set(missing_genes)),
                "present_genes": ";".join(sorted(set(present_genes))),
                "missing_genes": ";".join(sorted(set(missing_genes))),
                "reliable": str(reliable).lower(),
            }
        )

    for composite_name, components in COMPOSITES.items():
        component_cols = [f"z_{name}" for name in components if f"z_{name}" in spot_scores.columns]
        if component_cols:
            spot_scores[f"score_{composite_name}"] = spot_scores[component_cols].mean(axis=1)
        else:
            spot_scores[f"score_{composite_name}"] = np.nan
    return spot_scores, coverage_rows


def sample_summary(spot_scores: pd.DataFrame, manifest_row: pd.Series) -> dict:
    summary = {
        "dataset_id": manifest_row["dataset_id"],
        "sample_id": manifest_row["sample_id"],
        "patient_id": manifest_row["patient_id"],
        "specimen_type": manifest_row["specimen_type"],
        "treatment_context": manifest_row["treatment_context"],
        "n_spots_qc": len(spot_scores),
        "median_counts": float(np.nanmedian(spot_scores["n_counts"])),
        "median_genes": float(np.nanmedian(spot_scores["n_genes"])),
        "include_in_story_screen": "true",
        "exclude_reason": "",
    }
    summary_score_map = {
        "immune_hub_core": "score_immune_hub_core",
        "immune_hub_maturity": "score_immune_hub_maturity",
        "caf_myeloid_barrier": "score_caf_myeloid_barrier",
        "tumor_aggressive": "score_tumor_aggressive",
        "ifn_antigen_presentation": "z_ifn_antigen_presentation",
    }
    for prefix, score in summary_score_map.items():
        if score in spot_scores:
            summary[f"mean_{prefix}"] = float(np.nanmean(spot_scores[score]))
            summary[f"p90_{prefix}"] = float(np.nanpercentile(spot_scores[score], 90))
            summary[f"max_{prefix}"] = float(np.nanmax(spot_scores[score]))
            summary[f"fraction_{prefix}_z_gt1"] = float(np.nanmean(spot_scores[score] > 1.0))
        else:
            summary[f"mean_{prefix}"] = math.nan
            summary[f"p90_{prefix}"] = math.nan
            summary[f"max_{prefix}"] = math.nan
            summary[f"fraction_{prefix}_z_gt1"] = math.nan
    for raw_score in [
        "score_pan_caf",
        "score_myeloid",
        "score_spp1_tam",
        "score_tgfb_pathway",
        "score_icaf",
        "score_b_cell",
        "score_t_cell",
        "score_tls_chemokine",
        "score_dc_apc",
        "score_pdac_basal_like",
        "score_emt_invasion",
        "score_hypoxia",
    ]:
        if raw_score in spot_scores:
            prefix = raw_score.replace("score_", "")
            summary[f"mean_raw_{prefix}"] = float(np.nanmean(spot_scores[raw_score]))
            summary[f"p90_raw_{prefix}"] = float(np.nanpercentile(spot_scores[raw_score], 90))
    hub_core = spot_scores.get("score_immune_hub_core")
    hub_maturity = spot_scores.get("score_immune_hub_maturity")
    if hub_core is not None:
        summary["has_detectable_immune_hub"] = str(bool(np.nanmax(hub_core) > 1.0)).lower()
    else:
        summary["has_detectable_immune_hub"] = "false"
    if hub_core is not None and hub_maturity is not None:
        core_thr = np.nanpercentile(hub_core, 90)
        maturity_thr = np.nanpercentile(hub_maturity, 75)
        summary["has_mature_like_hub"] = str(bool(np.any((hub_core >= core_thr) & (hub_maturity >= maturity_thr)))).lower()
    else:
        summary["has_mature_like_hub"] = "false"
    summary["dominant_spatial_ecotype"] = "mvp_not_clustered_yet"
    return summary


def write_sample_figure(spot_scores: pd.DataFrame, output_dir: Path, sample_id: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    n = len(PLOT_SCORES)
    fig, axes = plt.subplots(1, n, figsize=(4 * n, 4), constrained_layout=True)
    if n == 1:
        axes = [axes]
    for ax, score in zip(axes, PLOT_SCORES):
        values = spot_scores[score] if score in spot_scores else np.zeros(len(spot_scores))
        sc = ax.scatter(
            spot_scores["x_pixel"],
            -spot_scores["y_pixel"],
            c=values,
            s=5,
            cmap="viridis",
            linewidths=0,
        )
        ax.set_title(score.replace("score_", "").replace("z_", "z_"))
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    fig.suptitle(sample_id)
    safe_id = sample_id.replace("/", "_").replace("\\", "_")
    fig.savefig(output_dir / f"{safe_id}_mvp_scores.png", dpi=180)
    plt.close(fig)


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


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run lightweight MVP signature scoring for Visium H5 samples.")
    parser.add_argument("--datasets", nargs="+", default=["GSE274103"])
    parser.add_argument("--max-samples", type=int, default=None)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    manifest = manifest[manifest["dataset_id"].isin(args.datasets)].copy()
    if args.max_samples is not None:
        manifest = manifest.head(args.max_samples)
    signatures = parse_signatures(PROJECT_ROOT / "config/signatures.yaml")

    output_dir = PROJECT_ROOT / "results/tables"
    figure_dir = PROJECT_ROOT / "results/figures/mvp"
    output_dir.mkdir(parents=True, exist_ok=True)

    all_spot_scores: list[pd.DataFrame] = []
    coverage_rows: list[dict] = []
    sample_rows: list[dict] = []
    qc_rows: list[dict] = []
    errors: list[str] = []

    for _, row in manifest.iterrows():
        expression_path = Path(row["expression_path"])
        coordinates_path = Path(row["coordinates_path"])
        if not expression_path.exists() or expression_path.suffix != ".h5":
            errors.append(f"Skipping non-H5 or missing expression file: {row['sample_id']}")
            continue
        if not coordinates_path.exists():
            errors.append(f"Skipping missing coordinates: {row['sample_id']}")
            continue
        try:
            counts, genes, barcodes = read_10x_h5(expression_path)
            positions = read_positions(coordinates_path)
            spot_scores, sample_coverage = score_signatures(
                counts,
                genes,
                barcodes,
                positions,
                signatures,
                row["dataset_id"],
                row["sample_id"],
                row["patient_id"],
            )
            all_spot_scores.append(spot_scores)
            coverage_rows.extend(sample_coverage)
            sample_rows.append(sample_summary(spot_scores, row))
            qc_rows.append(
                {
                    "dataset_id": row["dataset_id"],
                    "sample_id": row["sample_id"],
                    "spots_raw": len(barcodes),
                    "spots_in_tissue": len(spot_scores),
                    "spots_pass_qc": len(spot_scores),
                    "median_counts": float(np.nanmedian(spot_scores["n_counts"])),
                    "median_genes": float(np.nanmedian(spot_scores["n_genes"])),
                    "image_available": str(Path(row["image_path"]).exists()).lower(),
                    "include_after_qc": str(len(spot_scores) >= 100).lower(),
                    "qc_exclusion_reason": "" if len(spot_scores) >= 100 else "fewer_than_100_tissue_spots",
                }
            )
            write_sample_figure(spot_scores, figure_dir / row["dataset_id"], row["sample_id"])
        except Exception as exc:
            errors.append(f"{row['sample_id']}: {exc}")

    if all_spot_scores:
        pd.concat(all_spot_scores, ignore_index=True).to_csv(
            output_dir / "mvp_spot_level_scores.csv", index=False
        )
    pd.DataFrame(sample_rows).to_csv(output_dir / "mvp_sample_level_scores.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(output_dir / "signature_gene_coverage.csv", index=False)
    pd.DataFrame(qc_rows).to_csv(PROJECT_ROOT / "metadata/sample_qc_summary.csv", index=False)

    status = "success" if not errors else "partial_success"
    write_status(
        "03_mvp_score",
        status,
        {
            "n_datasets_processed": len(args.datasets),
            "n_samples_processed": len(sample_rows),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Inspect results/figures/mvp for spatial score plausibility.",
                "Review signature_gene_coverage.csv before interpreting missing signatures.",
                "Download GSE282302 if the GSE274103 pilot outputs look sane.",
            ],
        },
    )
    print(f"Scored {len(sample_rows)} samples")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
