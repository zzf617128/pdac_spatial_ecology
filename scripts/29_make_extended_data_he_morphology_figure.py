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
OUT = PROJECT_ROOT / "results/figures/submission/extended_data_figure4_he_morphology_bridge"

TARGET_ORDER = ["CAF-myeloid", "tumor aggressive", "IFN/MHC", "immune core"]
TARGET_LABELS = ["CAF-myeloid", "Tumor aggressive", "IFN/MHC", "Immune core"]
TARGET_COLORS = {
    "CAF-myeloid": "#8C2D04",
    "tumor aggressive": "#B23A48",
    "IFN/MHC": "#3B6FB6",
    "immune core": "#2C7A51",
}
FEATURE_ORDER = [
    "red OD",
    "stain density",
    "blue OD",
    "green OD",
    "brightness",
    "red-blue",
    "purple fraction",
    "saturation",
    "gray texture",
    "edge density",
]


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.13, 1.05, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E8E8E8", linewidth=0.55, zorder=0)


def load_cv_summary() -> pd.DataFrame:
    metrics = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_he_patch_grouped_cv_metrics.csv")
    observed = metrics[metrics["metric_scope"].eq("observed_all")].copy()
    null = (
        metrics[metrics["metric_scope"].eq("permuted_all")]
        .groupby(["target", "target_label"], as_index=False)
        .agg(
            null_median=("spearman_rho", "median"),
            null_q05=("spearman_rho", lambda s: s.quantile(0.05)),
            null_q95=("spearman_rho", lambda s: s.quantile(0.95)),
        )
    )
    return observed.merge(null, on=["target", "target_label"], how="left")


def plot_feature_heatmap(ax: plt.Axes) -> None:
    corr = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_he_patch_feature_correlation_summary.csv")
    corr = corr[corr["feature_label"].isin(FEATURE_ORDER)].copy()
    mat = (
        corr.pivot_table(index="feature_label", columns="target_label", values="median_rho", aggfunc="median")
        .reindex(index=FEATURE_ORDER, columns=TARGET_ORDER)
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.42, vmax=0.42, aspect="auto")
    ax.set_xticks(np.arange(len(TARGET_ORDER)), TARGET_LABELS, rotation=32, ha="right")
    ax.set_yticks(np.arange(len(FEATURE_ORDER)), FEATURE_ORDER)
    ax.set_title("H&E feature associations", fontsize=9, fontweight="bold", loc="left")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.3)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("median within-sample rho", fontsize=7)
    cb.ax.tick_params(labelsize=7)


def plot_cv(ax: plt.Axes) -> None:
    summary = load_cv_summary()
    summary["target_label"] = pd.Categorical(summary["target_label"], categories=TARGET_ORDER, ordered=True)
    summary = summary.sort_values("target_label")
    x = np.arange(len(summary))
    colors = [TARGET_COLORS[str(label)] for label in summary["target_label"]]
    ax.bar(x, summary["spearman_rho"], color=colors, width=0.62, alpha=0.9, zorder=3, label="observed")
    ax.errorbar(
        x,
        summary["null_median"],
        yerr=[
            summary["null_median"] - summary["null_q05"],
            summary["null_q95"] - summary["null_median"],
        ],
        fmt="o",
        color="#222222",
        markersize=4,
        capsize=3,
        lw=0.9,
        zorder=4,
        label="within-sample target shuffle",
    )
    for i, row in enumerate(summary.itertuples()):
        ax.text(i, float(row.spearman_rho) + 0.018, f"{float(row.spearman_rho):.2f}", ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color="#333333", lw=0.75)
    ax.set_xticks(x, TARGET_LABELS, rotation=28, ha="right")
    ax.set_ylabel("held-out grouped CV Spearman rho", fontsize=8)
    ax.set_ylim(-0.04, 0.48)
    ax.set_title("Grouped H&E feature model", fontsize=9, fontweight="bold", loc="left")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=7, loc="upper right")


def plot_scope(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_title("Patch-feature workflow", fontsize=9, fontweight="bold", loc="left")
    steps = [
        (0.10, 0.72, "H&E patch", "#E6C5CF"),
        (0.10, 0.43, "color +\ntexture", "#E8E8E8"),
        (0.10, 0.14, "grouped CV", "#DCE8F5"),
    ]
    for x0, y0, label, color in steps:
        ax.add_patch(
            plt.Rectangle((x0, y0), 0.58, 0.17, transform=ax.transAxes, facecolor=color, edgecolor="#555555", lw=0.8)
        )
        ax.text(x0 + 0.29, y0 + 0.085, label, ha="center", va="center", fontsize=8, fontweight="bold", transform=ax.transAxes)
    for y0, y1 in [(0.72, 0.60), (0.43, 0.31)]:
        ax.annotate(
            "",
            xy=(0.39, y1),
            xytext=(0.39, y0),
            xycoords=ax.transAxes,
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#666666"),
        )
    ax.add_patch(
        plt.Rectangle((0.72, 0.14), 0.23, 0.17, transform=ax.transAxes, facecolor="#F7F7F7", edgecolor="#888888", lw=0.8)
    )
    ax.text(0.835, 0.225, "target\nshuffle", ha="center", va="center", fontsize=7.3, transform=ax.transAxes)
    ax.annotate(
        "",
        xy=(0.72, 0.225),
        xytext=(0.68, 0.225),
        xycoords=ax.transAxes,
        arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#666666"),
    )
    ax.text(0.39, 0.02, "exploratory bridge", ha="center", fontsize=7.5, color="#333333", transform=ax.transAxes)


def main() -> int:
    fig = plt.figure(figsize=(12.4, 6.2))
    gs = GridSpec(1, 3, figure=fig, width_ratios=[1.35, 1.0, 0.78], wspace=0.38)
    axA = fig.add_subplot(gs[0, 0])
    panel_label(axA, "A")
    plot_feature_heatmap(axA)
    axB = fig.add_subplot(gs[0, 1])
    panel_label(axB, "B")
    plot_cv(axB)
    axC = fig.add_subplot(gs[0, 2])
    panel_label(axC, "C")
    plot_scope(axC)
    fig.suptitle("H&E patch morphology provides an exploratory bridge to the CAF-myeloid spatial axis", fontsize=12.5, fontweight="bold", y=0.98)
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
