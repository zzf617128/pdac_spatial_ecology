from __future__ import annotations

import importlib.util
import math
import tarfile
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_DATA = PROJECT / "data" / "external" / "GSE274557" / "pilot_selected"
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
REPORT_DIR = PROJECT / "results" / "reports"
TABLE_DIR.mkdir(parents=True, exist_ok=True)
FIG_DIR.mkdir(parents=True, exist_ok=True)
REPORT_DIR.mkdir(parents=True, exist_ok=True)

SIGNATURES = PROJECT / "config" / "signatures.yaml"
RNG = np.random.default_rng(20260627)


def load_stage03():
    path = PROJECT / "scripts" / "03_mvp_score_visium.py"
    spec = importlib.util.spec_from_file_location("stage03", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def find_positions(sample_dir: Path) -> Path:
    candidates = list(sample_dir.rglob("tissue_positions*.csv"))
    if not candidates:
        raise FileNotFoundError(f"No tissue_positions*.csv found under {sample_dir}")
    return candidates[0]


def median_nn_distance(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return 1.0
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    val = float(np.nanmedian(dists[:, 1]))
    return val if np.isfinite(val) and val > 0 else 1.0


def nearest_distance_to_core(xy: np.ndarray, core_mask: np.ndarray) -> np.ndarray:
    tree = cKDTree(xy[core_mask])
    dists, _ = tree.query(xy, k=1)
    return dists


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 4 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


def sample_from_filename(path: Path) -> str:
    return path.name.replace("_filtered_feature_bc_matrix.h5", "")


def build_composites(stage03, scores: pd.DataFrame) -> pd.DataFrame:
    for composite, components in stage03.COMPOSITES.items():
        z_cols = [f"z_{component}" for component in components if f"z_{component}" in scores.columns]
        if z_cols:
            scores[f"score_{composite}"] = scores[z_cols].mean(axis=1, skipna=True)
    return scores


def calc_sample_gradients(scores: pd.DataFrame, sample_meta: dict, n_random: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    xy = scores[["x_pixel", "y_pixel"]].to_numpy(float)
    core_score = scores["score_caf_myeloid_barrier"].to_numpy(float)
    n_core = max(5, int(math.ceil(0.10 * len(scores))))
    threshold = np.nanquantile(core_score, 1 - n_core / len(scores))
    core_mask = core_score >= threshold
    if core_mask.sum() < 3:
        core_mask[np.argsort(core_score)[-n_core:]] = True

    scale = median_nn_distance(xy)
    observed_dist = nearest_distance_to_core(xy, core_mask) / scale
    targets = {
        "IFN/MHC": "z_ifn_antigen_presentation",
        "immune_core": "score_immune_hub_core",
        "tumor_aggressive": "score_tumor_aggressive",
        "SPP1_TAM": "z_spp1_tam",
        "TGFb_EMT": "z_tgfb_pathway",
        "myCAF_matrix": "z_myCAF" if "z_myCAF" in scores.columns else "z_pan_caf",
    }

    rows = []
    core_rows = []
    for label, col in targets.items():
        if col not in scores.columns:
            continue
        vals = scores[col].to_numpy(float)
        obs_rho = safe_spearman(observed_dist, vals)
        random_rhos = []
        for _ in range(n_random):
            random_core = np.zeros(len(scores), dtype=bool)
            random_core[RNG.choice(len(scores), size=core_mask.sum(), replace=False)] = True
            random_dist = nearest_distance_to_core(xy, random_core) / scale
            random_rhos.append(safe_spearman(random_dist, vals))
        random_rhos = np.array(random_rhos, dtype=float)
        rows.append(
            {
                **sample_meta,
                "target_program": label,
                "n_spots": int(len(scores)),
                "n_caf_core_spots": int(core_mask.sum()),
                "observed_rho": obs_rho,
                "random_median_rho": float(np.nanmedian(random_rhos)),
                "delta_vs_random_median": obs_rho - float(np.nanmedian(random_rhos)),
                "observed_more_negative_than_random_median": bool(obs_rho < float(np.nanmedian(random_rhos))),
            }
        )
        core_rows.append(
            {
                **sample_meta,
                "target_program": label,
                "caf_core_enrichment": float(np.nanmean(vals[core_mask]) - np.nanmean(vals[~core_mask])),
                "caf_core_mean": float(np.nanmean(vals[core_mask])),
                "noncore_mean": float(np.nanmean(vals[~core_mask])),
            }
        )
    per_spot = scores.copy()
    per_spot["is_caf_core"] = core_mask
    per_spot["dist_to_caf_core_norm"] = observed_dist
    return pd.DataFrame(rows), pd.DataFrame(core_rows), per_spot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-prefix", default="gse274557_pilot")
    parser.add_argument("--n-random", type=int, default=300)
    args = parser.parse_args()
    input_dir = args.input_dir
    dataset_id = args.output_prefix.replace("gse", "GSE", 1)

    stage03 = load_stage03()
    signatures = stage03.parse_signatures(SIGNATURES)
    meta = pd.read_csv(TABLE_DIR / "gse274557_series_metadata.csv")
    meta_by_acc = meta.set_index("geo_accession").to_dict("index")

    gradient_all = []
    enrichment_all = []
    spot_all = []
    coverage_all = []

    h5_files = sorted(input_dir.glob("*_filtered_feature_bc_matrix.h5"))
    for h5_path in h5_files:
        sample_id = sample_from_filename(h5_path)
        accession = sample_id.split("_", 1)[0]
        spatial_dir = input_dir / sample_id
        positions_path = find_positions(spatial_dir)
        counts, genes, barcodes = stage03.read_10x_h5(h5_path)
        positions = stage03.read_positions(positions_path)
        row = meta_by_acc.get(accession, {})
        sample_meta = {
            "dataset_id": dataset_id,
            "geo_accession": accession,
            "sample_id": sample_id,
            "title": row.get("title", sample_id),
            "tissue": row.get("characteristic__tissue", "unknown"),
            "treatment": row.get("characteristic__treatment", "unknown"),
            "patient_id": row.get("title", sample_id).split("-")[1].split("_")[0] if "Pt-" in row.get("title", "") else row.get("title", sample_id),
        }
        scores, coverage = stage03.score_signatures(
            counts=counts,
            gene_names=genes,
            barcodes=barcodes,
            positions=positions,
            signatures=signatures,
            dataset_id=dataset_id,
            sample_id=sample_id,
            patient_id=sample_meta["patient_id"],
        )
        scores = build_composites(stage03, scores)
        for key, value in sample_meta.items():
            if key not in scores.columns:
                scores[key] = value
        gradients, enrichments, spot_scores = calc_sample_gradients(scores, sample_meta, args.n_random)
        gradient_all.append(gradients)
        enrichment_all.append(enrichments)
        spot_all.append(spot_scores)
        coverage_all.extend(coverage)

    gradients = pd.concat(gradient_all, ignore_index=True)
    enrichments = pd.concat(enrichment_all, ignore_index=True)
    spots = pd.concat(spot_all, ignore_index=True)
    coverage = pd.DataFrame(coverage_all)

    gradients.to_csv(TABLE_DIR / f"{args.output_prefix}_caf_core_gradients.csv", index=False)
    enrichments.to_csv(TABLE_DIR / f"{args.output_prefix}_caf_core_enrichment.csv", index=False)
    coverage.to_csv(TABLE_DIR / f"{args.output_prefix}_signature_coverage.csv", index=False)
    spots.to_csv(TABLE_DIR / f"{args.output_prefix}_spot_scores.csv", index=False)

    context = (
        gradients.groupby(["tissue", "target_program"])
        .agg(
            n_samples=("sample_id", "nunique"),
            median_observed_rho=("observed_rho", "median"),
            median_delta_vs_random=("delta_vs_random_median", "median"),
            n_more_negative=("observed_more_negative_than_random_median", "sum"),
        )
        .reset_index()
    )
    context.to_csv(TABLE_DIR / f"{args.output_prefix}_caf_core_context_summary.csv", index=False)

    plot_df = context[context["target_program"].isin(["IFN/MHC", "immune_core", "tumor_aggressive", "SPP1_TAM", "TGFb_EMT"])]
    tissue_order = ["Primary PDAC", "Liver metastasis", "Lung metastasis", "Peritoneal metastasis"]
    prog_order = ["IFN/MHC", "immune_core", "tumor_aggressive", "SPP1_TAM", "TGFb_EMT"]
    matrix = plot_df.pivot(index="target_program", columns="tissue", values="median_delta_vs_random").reindex(
        index=prog_order, columns=tissue_order
    )
    labels = plot_df.assign(label=lambda x: x["n_more_negative"].astype(int).astype(str) + "/" + x["n_samples"].astype(int).astype(str)).pivot(
        index="target_program", columns="tissue", values="label"
    ).reindex(index=prog_order, columns=tissue_order)

    mpl.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    fig, ax = plt.subplots(figsize=(8.8, 4.8))
    im = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax.set_xticks(range(len(tissue_order)))
    ax.set_xticklabels(tissue_order, rotation=25, ha="right")
    ax.set_yticks(range(len(prog_order)))
    ax.set_yticklabels(prog_order)
    title_prefix = "GSE274557 pilot" if "pilot" in args.output_prefix.lower() else "GSE274557 external validation"
    ax.set_title(f"{title_prefix}: CAF-core gradients vs random cores", loc="left", fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iat[i, j]
            lab = labels.iat[i, j]
            ax.text(j, i, "NA" if pd.isna(val) else f"{val:.2f}\n{lab}", ha="center", va="center", fontsize=8)
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label("median delta vs random rho")
    fig.tight_layout()
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(FIG_DIR / f"{args.output_prefix}_caf_core_external_validation.{ext}", dpi=300)
    plt.close(fig)

    lines = [
        "# GSE274557 CAF-Core External Validation",
        "",
        "Last updated: 2026-06-27",
        "",
        "## Scope",
        "",
        f"Downloaded and analyzed {gradients['sample_id'].nunique()} Nature 2025 metastatic PDAC Visium samples from GSE274557.",
        "",
        "## Outputs",
        "",
        f"- `results/tables/{args.output_prefix}_caf_core_gradients.csv`",
        f"- `results/tables/{args.output_prefix}_caf_core_enrichment.csv`",
        f"- `results/tables/{args.output_prefix}_caf_core_context_summary.csv`",
        f"- `results/tables/{args.output_prefix}_spot_scores.csv`",
        f"- `results/figures/submission/{args.output_prefix}_caf_core_external_validation.pdf`",
        "",
        "## Initial Interpretation",
        "",
        f"This run used {args.n_random} same-size random cores per sample. Negative delta values indicate that the observed CAF-core gradient is more target-program-centered than random cores.",
    ]
    for _, row in context.iterrows():
        lines.append(
            f"- {row.tissue} / {row.target_program}: median delta {row.median_delta_vs_random:.3f}, support {int(row.n_more_negative)}/{int(row.n_samples)} samples."
        )
    lines.append("")
    lines.append("## Next Step")
    lines.append("")
    if "pilot" in args.output_prefix.lower():
        lines.append("If pilot directions are biologically plausible, run the same analysis across the full non-PDX GSE274557 Visium set.")
    else:
        lines.append("Use this external cohort as Extended Data Figure 10. Frame the result as broad independent validation of CAF-core organization across primary and metastatic organ contexts, not as independent lymph-node validation or causal mechanism proof.")
    (REPORT_DIR / f"{args.output_prefix}_caf_core_external_validation.md").write_text("\n".join(lines), encoding="utf-8")

    print(f"Wrote {args.output_prefix} CAF-core outputs.")


if __name__ == "__main__":
    main()
