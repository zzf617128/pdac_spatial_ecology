from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.7,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT = PROJECT_ROOT / "results/figures/submission/extended_data_figure5_external_anchor_robustness"
SOURCE_OUT = PROJECT_ROOT / "results/source_data/Source_Data_ED_Fig_5D_cross_context_support_fraction.csv"

COLORS = {
    "IFN/MHC": "#3B6FB6",
    "Immune core": "#2C7A51",
    "Tumor aggressive": "#B23A48",
    "Immune maturity-like": "#7B68A6",
}

TARGET_LABELS = {
    "ifn_mhc": "IFN/MHC",
    "immune_core": "Immune core",
    "tumor_aggressive": "Tumor aggressive",
    "immune_maturity": "Immune maturity-like",
    "immune_maturity_like": "Immune maturity-like",
}

CONTEXT_LABELS = {
    "GSE282302": "post-NACT\nGSE282302",
    "GSE274103": "treatment-naive\nGSE274103",
    "primary_tumor": "primary\nGSE272362",
    "liver_metastasis": "liver met\nGSE272362",
    "lymph_node_metastasis": "LN met\nGSE272362",
    "GSE235315": "external anchor\nGSE235315",
}


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.05, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E6E6E6", linewidth=0.55, zorder=0)


def build_support_table() -> pd.DataFrame:
    mvp = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_random_core_permutation_summary.csv")
    mvp = mvp.rename(
        columns={
            "dataset_id": "context",
            "target": "target_key",
            "median_delta_vs_null": "delta",
            "n_observed_more_negative_than_null": "n_support",
        }
    )
    gse = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_random_core_permutation_summary.csv")
    gse = gse.rename(
        columns={
            "specimen_type": "context",
            "target": "target_key",
            "median_delta_vs_null": "delta",
            "n_observed_more_negative_than_null": "n_support",
        }
    )
    anchor = pd.read_csv(PROJECT_ROOT / "results/tables/gse235315_random_core_anchor_summary.csv")
    anchor["context"] = "GSE235315"
    label_to_key = {
        "IFN/MHC": "ifn_mhc",
        "immune core": "immune_core",
        "tumor aggressive": "tumor_aggressive",
        "immune maturity-like": "immune_maturity",
    }
    anchor["target_key"] = anchor["target_label"].map(label_to_key)
    anchor = (
        anchor.groupby(["context", "target_key"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            delta=("delta_vs_null_median", "median"),
            n_support=("observed_more_negative_than_null", lambda s: int(s.astype(bool).sum())),
        )
    )
    keep_cols = ["context", "target_key", "n_samples", "delta", "n_support"]
    combined = pd.concat([mvp[keep_cols], gse[keep_cols], anchor[keep_cols]], ignore_index=True)
    combined["target"] = combined["target_key"].map(TARGET_LABELS)
    return combined[combined["target"].notna()].copy()


def plot_support_heatmap(ax: plt.Axes, support: pd.DataFrame) -> None:
    contexts = ["GSE282302", "GSE274103", "primary_tumor", "liver_metastasis", "lymph_node_metastasis", "GSE235315"]
    targets = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity-like"]
    mat = support.pivot_table(index="target", columns="context", values="delta", aggfunc="median").reindex(index=targets, columns=contexts)
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.42, vmax=0.42, aspect="auto")
    ax.set_xticks(np.arange(len(contexts)), [CONTEXT_LABELS[c] for c in contexts], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(targets)), targets)
    ax.set_title("Random-core support across cohorts", fontsize=9, fontweight="bold", loc="left")
    for i, target in enumerate(targets):
        for j, context in enumerate(contexts):
            val = mat.loc[target, context]
            if not np.isfinite(val):
                continue
            row = support[(support["target"].eq(target)) & (support["context"].eq(context))]
            label = f"{val:.2f}\n{int(row['n_support'].iloc[0])}/{int(row['n_samples'].iloc[0])}"
            ax.text(j, i, label, ha="center", va="center", fontsize=6.0)
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("observed rho - random-core median", fontsize=7)
    cb.ax.tick_params(labelsize=7)


def plot_threshold(ax: plt.Axes) -> None:
    df = pd.read_csv(PROJECT_ROOT / "results/tables/caf_core_threshold_sensitivity_summary.csv")
    df["target_label"] = df["target"].map(TARGET_LABELS)
    df = df[df["target_label"].isin(["IFN/MHC", "Immune core", "Tumor aggressive"])]
    df = df[df["specimen_type"].isin(["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis"])]
    df["threshold"] = df["core_label"].str.extract(r"(\d+)").astype(float)
    for target, group in df.groupby("target_label"):
        summary = group.groupby("threshold", as_index=False)["median_rho"].median().sort_values("threshold", ascending=False)
        ax.plot(summary["threshold"], summary["median_rho"], marker="o", ms=4, lw=1.8, color=COLORS[target], label=target)
    ax.axhline(0, color="#333333", lw=0.75)
    ax.set_xlim(16, 4)
    ax.set_xticks([15, 10, 5])
    ax.set_xlabel("CAF-core percentile", fontsize=8)
    ax.set_ylabel("median distance-to-core rho", fontsize=8)
    ax.set_title("Threshold sensitivity across contexts", fontsize=9, fontweight="bold", loc="left")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=7, loc="lower left")


def plot_anchor_sample_deltas(ax: plt.Axes) -> None:
    df = pd.read_csv(PROJECT_ROOT / "results/tables/gse235315_random_core_anchor_summary.csv")
    order = ["IFN/MHC", "immune core", "tumor aggressive", "immune maturity-like"]
    labels = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity-like"]
    positions = np.arange(len(order))
    for i, target in enumerate(order):
        sub = df[df["target_label"].eq(target)].copy()
        vals = sub["delta_vs_null_median"].astype(float).to_numpy()
        jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else []
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=22, color=COLORS[labels[i]], alpha=0.8, edgecolor="white", linewidth=0.4, zorder=3)
        ax.plot([i - 0.22, i + 0.22], [np.median(vals), np.median(vals)], color="#111111", lw=1.3, zorder=4)
    ax.axhline(0, color="#333333", lw=0.75)
    ax.set_xticks(positions, labels, rotation=28, ha="right")
    ax.set_ylabel("observed rho - random-core median", fontsize=8)
    ax.set_title("GSE235315 per-sample anchor", fontsize=9, fontweight="bold", loc="left")
    clean_axes(ax, axis="y")


def summarize_cross_context_support(support: pd.DataFrame) -> pd.DataFrame:
    order = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity-like"]
    summary = (
        support.groupby("target", as_index=False)
        .agg(
            n_support=("n_support", "sum"),
            n_samples=("n_samples", "sum"),
            median_delta=("delta", "median"),
        )
        .assign(support_fraction=lambda d: d["n_support"] / d["n_samples"])
    )
    summary["target"] = pd.Categorical(summary["target"], categories=order, ordered=True)
    return summary.sort_values("target").reset_index(drop=True)


def plot_support_fraction(ax: plt.Axes, support: pd.DataFrame) -> None:
    summary = summarize_cross_context_support(support)
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    summary.to_csv(SOURCE_OUT, index=False)

    y = np.arange(len(summary))
    colors = [COLORS[str(target)] for target in summary["target"]]
    ax.barh(y, summary["support_fraction"], color=colors, alpha=0.88, height=0.58, zorder=3)
    ax.set_yticks(y, summary["target"].astype(str))
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("samples beating random-core median", fontsize=8)
    ax.set_title("Cross-context support fraction", fontsize=9, fontweight="bold", loc="left")
    ax.invert_yaxis()
    clean_axes(ax, axis="x")
    ax.set_xticks([0, 0.5, 1.0], ["0", "0.5", "1.0"])
    for i, row in summary.iterrows():
        frac = float(row["support_fraction"])
        label = f"{int(row['n_support'])}/{int(row['n_samples'])}"
        ax.text(min(frac + 0.035, 0.92), i, label, va="center", ha="left", fontsize=7.2, color="#222222")


def main() -> int:
    support = build_support_table()
    fig = plt.figure(figsize=(13.2, 7.4))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.22, 1.0], width_ratios=[1.1, 1.0, 0.78], hspace=0.50, wspace=0.42)
    axA = fig.add_subplot(gs[0, :])
    panel_label(axA, "A")
    plot_support_heatmap(axA, support)
    axB = fig.add_subplot(gs[1, 0])
    panel_label(axB, "B")
    plot_threshold(axB)
    axC = fig.add_subplot(gs[1, 1])
    panel_label(axC, "C")
    plot_anchor_sample_deltas(axC)
    axD = fig.add_subplot(gs[1, 2])
    panel_label(axD, "D")
    plot_support_fraction(axD, support)
    fig.suptitle("External anchor and robustness of CAF-core spatial organization", fontsize=13, fontweight="bold", y=0.985)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = OUT.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
