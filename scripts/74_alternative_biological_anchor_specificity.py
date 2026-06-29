from __future__ import annotations

from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
REPORT_DIR = PROJECT / "results" / "reports"
SOURCE_DIR = PROJECT / "results" / "source_data"
for d in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260628)
N_RANDOM = 250
ANCHOR_FRACTION = 0.10


DATASETS = [
    ("GSE282302", "results/tables/mvp_spot_level_scores_with_edge_qc.csv", "post-NACT"),
    ("GSE272362", "results/tables/gse272362_rds_spot_level_scores.csv", None),
    ("GSE235315", "results/tables/gse235315_spot_level_scores.csv", "paired ST-H&E"),
    ("GSE274557", "results/tables/gse274557_full_spot_scores.csv", None),
]

ANCHORS = {
    "CAF-myeloid core": "score_caf_myeloid_barrier",
    "tumor-aggressive high": "score_tumor_aggressive",
    "tumor-epithelial high": "score_tumor_epithelial",
    "immune-core high": "score_immune_hub_core",
    "panCAF high": "score_pan_caf",
    "SPP1/TAM high": "score_spp1_tam",
}

TARGETS = {
    "IFN/MHC": "score_ifn_antigen_presentation",
    "immune core": "score_immune_hub_core",
    "tumor aggressive": "score_tumor_aggressive",
    "SPP1/TAM": "score_spp1_tam",
    "TGFb/EMT": "score_emt_invasion",
    "tumor epithelial": "score_tumor_epithelial",
}


def zscore_within_sample(df: pd.DataFrame, col: str) -> np.ndarray:
    values = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    mu = np.nanmean(values)
    sd = np.nanstd(values)
    if not np.isfinite(sd) or sd == 0:
        return np.full_like(values, np.nan)
    return (values - mu) / sd


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 20:
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


def nearest_dist_norm(xy: np.ndarray, mask: np.ndarray, scale: float) -> np.ndarray:
    tree = cKDTree(xy[mask])
    dists, _ = tree.query(xy, k=1)
    return dists / scale


def context_from_row(dataset_id: str, default_context: str | None, sample: pd.DataFrame) -> str:
    for col in ["specimen_type", "tissue", "treatment"]:
        if col in sample.columns:
            vals = sample[col].dropna().astype(str).unique()
            if len(vals):
                return vals[0]
    return default_context or dataset_id


def load_dataset(dataset_id: str, rel_path: str, default_context: str | None) -> pd.DataFrame:
    df = pd.read_csv(PROJECT / rel_path)
    if "edge_or_background_risk" in df.columns:
        df = df.loc[~df["edge_or_background_risk"].fillna(False).astype(bool)].copy()
    needed = {"sample_id", "x_pixel", "y_pixel", *ANCHORS.values(), *TARGETS.values()}
    missing = sorted(c for c in needed if c not in df.columns)
    if missing:
        raise ValueError(f"{dataset_id} missing columns: {missing}")
    df["analysis_dataset"] = dataset_id
    df["analysis_context"] = (
        df.groupby("sample_id", group_keys=False)
        .apply(lambda x: pd.Series([context_from_row(dataset_id, default_context, x)] * len(x), index=x.index), include_groups=False)
        .sort_index()
        .to_numpy()
    )
    return df


def analyze_sample(dataset_id: str, context: str, sample_id: str, sample: pd.DataFrame) -> tuple[list[dict], list[dict]]:
    sample = sample.dropna(subset=["x_pixel", "y_pixel"]).copy()
    if len(sample) < 80:
        return [], []
    xy = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    if not np.isfinite(xy).all():
        return [], []
    scale = median_nn_scale(xy)
    n_anchor = max(30, int(math.ceil(ANCHOR_FRACTION * len(sample))))
    n_anchor = min(n_anchor, len(sample) - 5)

    target_z = {label: zscore_within_sample(sample, col) for label, col in TARGETS.items()}
    rows: list[dict] = []
    curve_rows: list[dict] = []

    for anchor_label, anchor_col in ANCHORS.items():
        av = pd.to_numeric(sample[anchor_col], errors="coerce").to_numpy(float)
        if np.isfinite(av).sum() < n_anchor:
            continue
        cutoff = np.nanquantile(av, 1 - n_anchor / len(sample))
        mask = np.isfinite(av) & (av >= cutoff)
        if mask.sum() < 20:
            continue
        dist = nearest_dist_norm(xy, mask, scale)
        random_dists = []
        for _ in range(N_RANDOM):
            rand_mask = np.zeros(len(sample), dtype=bool)
            rand_mask[RNG.choice(len(sample), size=int(mask.sum()), replace=False)] = True
            random_dists.append(nearest_dist_norm(xy, rand_mask, scale))
        random_dists = np.asarray(random_dists)
        bins = pd.qcut(pd.Series(dist), q=5, labels=False, duplicates="drop")
        for target_label, y in target_z.items():
            observed = safe_spearman(dist, y)
            random_rhos = np.array([safe_spearman(rd, y) for rd in random_dists], dtype=float)
            rows.append(
                {
                    "dataset_id": dataset_id,
                    "context": context,
                    "sample_id": sample_id,
                    "anchor": anchor_label,
                    "target_program": target_label,
                    "n_spots": int(len(sample)),
                    "n_anchor_spots": int(mask.sum()),
                    "observed_rho": observed,
                    "random_median_rho": float(np.nanmedian(random_rhos)),
                    "delta_vs_random_median": observed - float(np.nanmedian(random_rhos)),
                    "observed_more_negative_than_random_median": bool(observed < float(np.nanmedian(random_rhos))),
                    "random_p05_rho": float(np.nanpercentile(random_rhos, 5)),
                    "random_p95_rho": float(np.nanpercentile(random_rhos, 95)),
                }
            )
            for b in sorted(pd.Series(bins).dropna().unique()):
                m = np.asarray(bins == b)
                curve_rows.append(
                    {
                        "dataset_id": dataset_id,
                        "context": context,
                        "sample_id": sample_id,
                        "anchor": anchor_label,
                        "target_program": target_label,
                        "distance_bin": int(b) + 1,
                        "mean_z": float(np.nanmean(y[m])),
                        "n_spots": int(m.sum()),
                    }
                )
    return rows, curve_rows


def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    all_rows: list[dict] = []
    all_curves: list[dict] = []
    for dataset_id, path, default_context in DATASETS:
        df = load_dataset(dataset_id, path, default_context)
        for sample_id, sample in df.groupby("sample_id", sort=False):
            context = str(sample["analysis_context"].iloc[0])
            rows, curves = analyze_sample(dataset_id, context, str(sample_id), sample)
            all_rows.extend(rows)
            all_curves.extend(curves)
            print(f"{dataset_id} {sample_id}: rows={len(rows)} curves={len(curves)}")
    detail = pd.DataFrame(all_rows)
    curves = pd.DataFrame(all_curves)
    summary = (
        detail.groupby(["anchor", "target_program"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_delta=("delta_vs_random_median", "median"),
            median_observed_rho=("observed_rho", "median"),
            support_n=("observed_more_negative_than_random_median", "sum"),
        )
    )
    summary["support_fraction"] = summary["support_n"] / summary["n_samples"]
    return detail, curves, summary


def make_figure(detail: pd.DataFrame, curves: pd.DataFrame, summary: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.8,
        }
    )
    anchor_order = list(ANCHORS.keys())
    target_order = ["IFN/MHC", "immune core", "tumor aggressive", "SPP1/TAM", "TGFb/EMT", "tumor epithelial"]
    mat = (
        summary.pivot(index="target_program", columns="anchor", values="median_delta")
        .reindex(index=target_order, columns=anchor_order)
    )
    support = (
        summary.pivot(index="target_program", columns="anchor", values="support_fraction")
        .reindex(index=target_order, columns=anchor_order)
    )

    fig = plt.figure(figsize=(13.8, 9.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.0, 1.0], height_ratios=[1.05, 1.0], wspace=0.48, hspace=0.55)
    fig.suptitle("Alternative biological anchors test CAF-myeloid core specificity", fontsize=16, fontweight="bold")

    ax_a = fig.add_subplot(gs[:, 0])
    im = ax_a.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax_a.set_title("A  Anchor-target specificity", loc="left", fontweight="bold")
    ax_a.set_xticks(range(len(anchor_order)))
    ax_a.set_xticklabels(anchor_order, rotation=45, ha="right", fontsize=8)
    ax_a.set_yticks(range(len(target_order)))
    ax_a.set_yticklabels(target_order)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            supp = support.iat[i, j]
            txt = "NA" if pd.isna(val) else f"{val:.2f}\n{supp:.0%}"
            ax_a.text(j, i, txt, ha="center", va="center", fontsize=7.5, color="white" if pd.notna(val) and abs(val) > 0.18 else "#111")
    cb = fig.colorbar(im, ax=ax_a, fraction=0.045, pad=0.03)
    cb.set_label("median observed-minus-random rho")

    ax_b = fig.add_subplot(gs[0, 1])
    focus = summary[summary["target_program"].isin(["IFN/MHC", "immune core", "tumor aggressive"])]
    caf = focus[focus["anchor"].eq("CAF-myeloid core")].set_index("target_program")
    others = focus[~focus["anchor"].eq("CAF-myeloid core")]
    best_other = others.sort_values("median_delta").groupby("target_program").first()
    xs = np.arange(len(caf))
    width = 0.36
    ax_b.bar(xs - width / 2, caf.loc[target_order[:3], "median_delta"], width=width, label="CAF-myeloid", color="#2D6A8E")
    ax_b.bar(xs + width / 2, best_other.loc[target_order[:3], "median_delta"], width=width, label="best alternative", color="#B57A3C")
    ax_b.axhline(0, color="#333", lw=0.8)
    ax_b.set_xticks(xs)
    ax_b.set_xticklabels(target_order[:3], rotation=25, ha="right", fontsize=8)
    ax_b.set_ylabel("median delta")
    ax_b.set_title("B  CAF core versus best alternative", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, fontsize=8)

    ax_c = fig.add_subplot(gs[0, 2])
    context = (
        detail[detail["anchor"].eq("CAF-myeloid core")]
        .groupby(["context", "target_program"], as_index=False)
        .agg(median_delta=("delta_vs_random_median", "median"))
    )
    keep_contexts = context.groupby("context")["median_delta"].size().sort_values(ascending=False).head(7).index
    cmat = (
        context[context["context"].isin(keep_contexts)]
        .pivot(index="context", columns="target_program", values="median_delta")
        .reindex(columns=target_order[:5])
    )
    im2 = ax_c.imshow(cmat.to_numpy(float), cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax_c.set_title("C  CAF-core specificity by context", loc="left", fontweight="bold")
    ax_c.set_xticks(range(cmat.shape[1]))
    ax_c.set_xticklabels(cmat.columns, rotation=35, ha="right", fontsize=8)
    ax_c.set_yticks(range(cmat.shape[0]))
    ax_c.set_yticklabels(cmat.index, fontsize=8)
    fig.colorbar(im2, ax=ax_c, fraction=0.046, pad=0.04)

    ax_d = fig.add_subplot(gs[1, 1])
    curve_focus = curves[
        curves["anchor"].isin(["CAF-myeloid core", "tumor-aggressive high", "immune-core high"])
        & curves["target_program"].eq("IFN/MHC")
    ]
    curve_summary = (
        curve_focus.groupby(["anchor", "distance_bin"], as_index=False)
        .agg(mean_z=("mean_z", "median"))
    )
    colors = {"CAF-myeloid core": "#2D6A8E", "tumor-aggressive high": "#B57A3C", "immune-core high": "#628F4E"}
    for anchor in ["CAF-myeloid core", "tumor-aggressive high", "immune-core high"]:
        sub = curve_summary[curve_summary["anchor"].eq(anchor)]
        ax_d.plot(sub["distance_bin"], sub["mean_z"], marker="o", lw=2, color=colors[anchor], label=anchor)
    ax_d.axhline(0, color="#333", lw=0.8)
    ax_d.set_xlabel("distance bin from anchor")
    ax_d.set_ylabel("median within-sample z")
    ax_d.set_title("D  IFN/MHC distance profiles", loc="left", fontweight="bold")
    ax_d.legend(frameon=False, fontsize=8)

    ax_e = fig.add_subplot(gs[1, 2])
    curve_focus = curves[
        curves["anchor"].isin(["CAF-myeloid core", "tumor-aggressive high", "SPP1/TAM high"])
        & curves["target_program"].eq("tumor aggressive")
    ]
    curve_summary = (
        curve_focus.groupby(["anchor", "distance_bin"], as_index=False)
        .agg(mean_z=("mean_z", "median"))
    )
    colors2 = {"CAF-myeloid core": "#2D6A8E", "tumor-aggressive high": "#B57A3C", "SPP1/TAM high": "#9C4E70"}
    for anchor in ["CAF-myeloid core", "tumor-aggressive high", "SPP1/TAM high"]:
        sub = curve_summary[curve_summary["anchor"].eq(anchor)]
        ax_e.plot(sub["distance_bin"], sub["mean_z"], marker="o", lw=2, color=colors2[anchor], label=anchor)
    ax_e.axhline(0, color="#333", lw=0.8)
    ax_e.set_xlabel("distance bin from anchor")
    ax_e.set_ylabel("median within-sample z")
    ax_e.set_title("E  Tumor-aggressive distance profiles", loc="left", fontweight="bold")
    ax_e.legend(frameon=False, fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / "extended_data_figure30_alternative_biological_anchor_specificity"
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(detail: pd.DataFrame, summary: pd.DataFrame) -> None:
    caf = summary[summary["anchor"].eq("CAF-myeloid core")].sort_values("median_delta")
    lines = [
        "# Alternative Biological Anchor Specificity",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Question",
        "",
        "This analysis tests whether CAF-myeloid cores are a specific spatial organizer rather than one of many high-scoring biological regions that would produce similar proximity gradients.",
        "",
        "## Method",
        "",
        f"For each section, the top {int(ANCHOR_FRACTION * 100)}% of spots for each anchor program were selected. Target-program proximity was measured as the Spearman correlation between normalized distance to the nearest anchor spot and within-sample target z-score. Each anchor was compared with {N_RANDOM} same-size random anchors per sample.",
        "",
        "## Key CAF-myeloid Core Results",
        "",
    ]
    for _, row in caf.iterrows():
        lines.append(
            f"- {row.target_program}: median observed-minus-random rho {row.median_delta:.3f}; support {int(row.support_n)}/{int(row.n_samples)} samples ({row.support_fraction:.1%})."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A negative observed-minus-random rho indicates that a target program is more centered on the anchor than expected for same-size random tissue regions. The analysis is designed as a biological-anchor specificity test, not a causal perturbation test.",
            "",
            "## Outputs",
            "",
            "- `results/tables/alternative_biological_anchor_specificity_per_sample.csv`",
            "- `results/tables/alternative_biological_anchor_specificity_summary.csv`",
            "- `results/tables/alternative_biological_anchor_distance_curves.csv`",
            "- `results/source_data/Source_Data_Extended_Data_Fig_30_alternative_anchor_specificity.csv`",
            "- `results/figures/submission/extended_data_figure30_alternative_biological_anchor_specificity.pdf`",
        ]
    )
    (REPORT_DIR / "alternative_biological_anchor_specificity_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    detail, curves, summary = run_analysis()
    detail.to_csv(TABLE_DIR / "alternative_biological_anchor_specificity_per_sample.csv", index=False)
    curves.to_csv(TABLE_DIR / "alternative_biological_anchor_distance_curves.csv", index=False)
    summary.to_csv(TABLE_DIR / "alternative_biological_anchor_specificity_summary.csv", index=False)
    summary.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_30_alternative_anchor_specificity.csv", index=False)
    curves.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_30_distance_curves.csv", index=False)
    make_figure(detail, curves, summary)
    write_report(detail, summary)
    print("Wrote alternative biological anchor specificity outputs.")


if __name__ == "__main__":
    main()
