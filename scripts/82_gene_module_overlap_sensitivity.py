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
MODULE_TABLE = SUPP_TABLES / "Supplementary_Table_2_Gene_Modules.csv"
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

TARGETS = ["IFN/MHC", "immune-core", "tumor-aggressive", "SPP1/TAM", "TGF-beta/EMT", "myCAF/matrix"]
MODULE_NAME_MAP = {
    "IFN/MHC": "IFN/MHC antigen-presentation",
    "immune-core": "immune-core",
    "tumor-aggressive": "tumor-aggressive",
    "SPP1/TAM": "SPP1/TAM",
    "TGF-beta/EMT": "TGF-beta/EMT",
    "myCAF/matrix": "myCAF/matrix",
}


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


def load_modules() -> dict[str, set[str]]:
    rows = pd.read_csv(MODULE_TABLE, encoding="utf-8-sig")
    modules: dict[str, set[str]] = {}
    for module_name, group in rows.groupby("module_name"):
        modules[str(module_name)] = set(group["gene_symbol"].dropna().astype(str).str.upper())
    if "TGF-beta/EMT" not in modules:
        modules["TGF-beta/EMT"] = modules.get("TGF-beta", set()) | modules.get("EMT/invasion", set())
    return modules


def build_overlap_matrix(modules: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    names = sorted(modules)
    for a in names:
        genes_a = modules[a]
        for b in names:
            genes_b = modules[b]
            shared = genes_a & genes_b
            union = genes_a | genes_b
            rows.append(
                {
                    "module_a": a,
                    "module_b": b,
                    "module_a_genes_n": len(genes_a),
                    "module_b_genes_n": len(genes_b),
                    "shared_genes_n": len(shared),
                    "jaccard_index": len(shared) / len(union) if union else np.nan,
                    "percent_module_b_shared": len(shared) / len(genes_b) if genes_b else np.nan,
                    "shared_gene_symbols": ";".join(sorted(shared)),
                }
            )
    return pd.DataFrame(rows)


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 30 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
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


def target_values(df: pd.DataFrame, target: str) -> np.ndarray:
    if target == "IFN/MHC":
        return df["z_ifn_antigen_presentation"].to_numpy(float)
    if target == "immune-core":
        return df["score_immune_hub_core"].to_numpy(float)
    if target == "tumor-aggressive":
        return df["score_tumor_aggressive"].to_numpy(float)
    if target == "SPP1/TAM":
        return df["z_spp1_tam"].to_numpy(float)
    if target == "TGF-beta/EMT":
        return df[["z_tgfb_pathway", "z_emt_invasion"]].mean(axis=1, skipna=True).to_numpy(float)
    if target == "myCAF/matrix":
        return df[["z_mycaf", "z_pan_caf"]].mean(axis=1, skipna=True).to_numpy(float)
    raise KeyError(target)


def core_score_excluding_target_overlap(df: pd.DataFrame, target: str, shared_n: int) -> tuple[np.ndarray, str]:
    if shared_n == 0:
        return df["score_caf_myeloid_barrier"].to_numpy(float), "no_shared_genes_original_core"
    if target == "SPP1/TAM":
        return df[["z_mycaf", "z_pan_caf", "z_myeloid"]].mean(axis=1, skipna=True).to_numpy(float), "core_excluding_spp1_tam_component"
    if target == "myCAF/matrix":
        return df[["z_myeloid", "z_spp1_tam"]].mean(axis=1, skipna=True).to_numpy(float), "core_excluding_caf_matrix_component"
    return df["score_caf_myeloid_barrier"].to_numpy(float), "target_overlap_not_recomputable_from_existing_scores"


def analyze_sample(
    df: pd.DataFrame,
    dataset: str,
    tissue_site: str,
    n_perm: int,
    caf_genes: set[str],
    modules: dict[str, set[str]],
) -> list[dict[str, object]]:
    if "edge_or_background_risk" in df.columns:
        risk = df["edge_or_background_risk"].astype(str).str.lower().isin(["true", "1", "yes"])
        df = df[~risk].copy()
    df = df[np.isfinite(df["x_pixel"]) & np.isfinite(df["y_pixel"])].copy()
    if len(df) < 100:
        return []
    xy = df[["x_pixel", "y_pixel"]].to_numpy(float)
    scale = median_nn_scale(xy)
    n_core = max(10, int(math.ceil(0.10 * len(df))))
    random_indices = [RNG.choice(len(df), size=n_core, replace=False) for _ in range(n_perm)]
    sample_id = str(df["sample_id"].iloc[0])
    patient_id = str(df["patient_id"].iloc[0]) if "patient_id" in df.columns else ""
    rows: list[dict[str, object]] = []

    random_null: dict[str, np.ndarray] = {}
    for target in TARGETS:
        y = target_values(df, target)
        vals = []
        for ridx in random_indices:
            vals.append(safe_spearman(nearest_dist(xy, ridx, scale), y))
        random_null[target] = np.array(vals, dtype=float)

    for target in TARGETS:
        target_module = MODULE_NAME_MAP[target]
        target_genes = modules[target_module]
        shared = caf_genes & target_genes
        remaining = target_genes - caf_genes
        core_score, sensitivity_mode = core_score_excluding_target_overlap(df, target, len(shared))
        idx = np.argsort(core_score)[-n_core:]
        dist = nearest_dist(xy, idx, scale)
        y = target_values(df, target)
        observed = safe_spearman(dist, y)
        null = random_null[target]
        n_finite = int(np.isfinite(null).sum())
        median = float(np.nanmedian(null)) if n_finite else np.nan
        target_removed_status = "target_score_unchanged_no_overlap" if not shared else ("target_empty_after_overlap_removal" if not remaining else "target_overlap_removed_requires_raw_rescoring")
        rows.append(
            {
                "dataset": dataset,
                "sample_id": sample_id,
                "patient_id": patient_id,
                "tissue_site": tissue_site,
                "core_module": "CAF-myeloid",
                "target_module": target_module,
                "target_program": target,
                "shared_genes_n": len(shared),
                "target_genes_n": len(target_genes),
                "remaining_target_genes_n": len(remaining),
                "jaccard_index": len(shared) / len(caf_genes | target_genes) if (caf_genes | target_genes) else np.nan,
                "percent_target_shared_with_core": len(shared) / len(target_genes) if target_genes else np.nan,
                "shared_gene_symbols": ";".join(sorted(shared)),
                "sensitivity_mode": sensitivity_mode,
                "target_removed_status": target_removed_status,
                "n_spots": len(df),
                "n_core_spots": n_core,
                "n_perm": n_perm,
                "overlap_sensitive_observed_rho": observed,
                "overlap_sensitive_random_median_rho": median,
                "overlap_sensitive_delta": observed - median if np.isfinite(observed) and np.isfinite(median) else np.nan,
                "empirical_p": (1 + int(np.nansum(null <= observed))) / (1 + n_finite) if np.isfinite(observed) and n_finite else np.nan,
                "support": bool(observed < median) if np.isfinite(observed) and np.isfinite(median) else "",
            }
        )
    return rows


def load_dataset(key: str) -> pd.DataFrame:
    path, _, _, extras = DATASETS[key]
    cols = available_usecols(path, USECOLS + extras)
    return pd.read_csv(path, usecols=cols, encoding="utf-8-sig")


def summarize(per_sample: pd.DataFrame) -> pd.DataFrame:
    return (
        per_sample.groupby(["core_module", "target_module", "target_program", "sensitivity_mode", "target_removed_status"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            shared_genes_n=("shared_genes_n", "first"),
            target_genes_n=("target_genes_n", "first"),
            remaining_target_genes_n=("remaining_target_genes_n", "first"),
            jaccard_index=("jaccard_index", "first"),
            percent_target_shared_with_core=("percent_target_shared_with_core", "first"),
            median_overlap_sensitive_delta=("overlap_sensitive_delta", "median"),
            support_n=("support", lambda x: int(np.sum([v is True or str(v).lower() == "true" for v in x]))),
            support_fraction=("support", lambda x: float(np.mean([v is True or str(v).lower() == "true" for v in x]))),
            median_empirical_p=("empirical_p", "median"),
            shared_gene_symbols=("shared_gene_symbols", "first"),
        )
        .reset_index()
    )


def write_report(summary: pd.DataFrame, overlap_focus: pd.DataFrame, n_perm: int) -> None:
    lines = [
        "# Gene-module overlap sensitivity",
        "",
        f"n_perm per sample: {n_perm}",
        "",
        "This analysis quantifies marker overlap between CAF-myeloid core genes and target modules.",
        "For target modules without shared genes, the original target score is unchanged and overlap sensitivity is equivalent to the core-gradient analysis.",
        "For component modules that become empty after removing shared genes, the target-removed score is not estimable; instead, an anchor-side sensitivity was run by excluding the overlapping component from the CAF-myeloid core score where possible.",
        "",
        "## CAF-myeloid overlap summary",
        "",
    ]
    for _, row in overlap_focus.iterrows():
        lines.append(
            f"- CAF-myeloid vs {row['module_b']}: shared {int(row['shared_genes_n'])}/{int(row['module_b_genes_n'])} target genes "
            f"(Jaccard {row['jaccard_index']:.3f}); shared genes: {row['shared_gene_symbols'] or 'none'}."
        )
    lines.extend(["", "## Sensitivity summary", ""])
    for _, row in summary.iterrows():
        lines.append(
            f"- {row['target_program']}: mode={row['sensitivity_mode']}; status={row['target_removed_status']}; "
            f"median delta {row['median_overlap_sensitive_delta']:.3f}; support {int(row['support_n'])}/{int(row['n_samples'])}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The main IFN/MHC, immune-core and tumor-aggressive target modules do not share marker genes with the current CAF-myeloid core definition, so their CAF-core gradients are not explained by direct marker overlap.",
            "SPP1/TAM and myCAF/matrix are core-component interpretation modules. They should be described as components of the CAF-myeloid architecture rather than independent non-overlapping targets.",
            "",
        ]
    )
    (REPORTS / "gene_module_overlap_sensitivity_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--datasets", nargs="+", default=list(DATASETS), choices=sorted(DATASETS))
    parser.add_argument("--n-perm", type=int, default=1000)
    parser.add_argument("--max-samples", type=int, default=0)
    args = parser.parse_args()

    ANALYSIS_OUT.mkdir(parents=True, exist_ok=True)
    SUPP_TABLES.mkdir(parents=True, exist_ok=True)
    REPORTS.mkdir(parents=True, exist_ok=True)

    modules = load_modules()
    caf_genes = modules["CAF-myeloid"]
    overlap = build_overlap_matrix(modules)
    overlap.to_csv(ANALYSIS_OUT / "module_overlap_matrix_all_pairs.csv", index=False, encoding="utf-8")
    overlap_focus = overlap[(overlap["module_a"] == "CAF-myeloid") & (overlap["module_b"].isin(MODULE_NAME_MAP.values()))].copy()
    overlap_focus.to_csv(ANALYSIS_OUT / "caf_myeloid_target_module_overlap_matrix.csv", index=False, encoding="utf-8")

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
            rows.extend(analyze_sample(sdf, dataset, tissue_site, args.n_perm, caf_genes, modules))
        if args.max_samples and counter >= args.max_samples:
            break

    per_sample = pd.DataFrame(rows)
    summary = summarize(per_sample)
    per_sample.to_csv(ANALYSIS_OUT / "gene_module_overlap_sensitivity_per_sample.csv", index=False, encoding="utf-8")
    summary.to_csv(ANALYSIS_OUT / "gene_module_overlap_sensitivity_summary.csv", index=False, encoding="utf-8")
    summary.to_csv(SUPP_TABLES / "Supplementary_Table_5_Module_Overlap_Sensitivity.csv", index=False, encoding="utf-8")
    write_report(summary, overlap_focus, args.n_perm)
    print("Wrote gene-module overlap sensitivity outputs.", flush=True)


if __name__ == "__main__":
    main()
