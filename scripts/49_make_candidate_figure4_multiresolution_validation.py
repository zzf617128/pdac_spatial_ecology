from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
SOURCE_DIR = PROJECT / "results" / "source_data"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "figure4_submission_multiresolution_validation"

VISIUM_TISSUES = ["Primary PDAC", "Liver metastasis", "Lung metastasis", "Peritoneal metastasis"]
VISIUM_PROGRAMS = ["IFN/MHC", "immune_core", "tumor_aggressive", "SPP1_TAM", "TGFb_EMT"]
VISIUM_LABELS = {
    "IFN/MHC": "IFN/MHC",
    "immune_core": "immune core",
    "tumor_aggressive": "tumor aggressive",
    "SPP1_TAM": "SPP1/TAM",
    "TGFb_EMT": "TGF-beta/EMT",
}
VISIUM_COLORS = {
    "Primary PDAC": "#4C78A8",
    "Liver metastasis": "#59A14F",
    "Lung metastasis": "#E15759",
    "Peritoneal metastasis": "#B279A2",
}

XENIUM_TARGETS = ["SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT", "Tumor_epithelial", "SPP1_tumor_like"]
XENIUM_LABELS = {
    "SPP1_TAM": "SPP1/TAM",
    "IFN_APC": "IFN/APC",
    "T_NK": "T/NK",
    "TGFb_EMT": "TGF-beta/EMT",
    "Tumor_epithelial": "tumor epithelial",
    "SPP1_tumor_like": "SPP1 tumor-like",
}
ANCHOR_LABELS = {"CAF_APC": "CAF-APC", "CAF_SPP1TAM": "CAF-SPP1/TAM"}


def setup_style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def draw_cohort_scale(ax: plt.Axes, visium_gradients: pd.DataFrame, xenium_composition: pd.DataFrame) -> None:
    visium_counts = (
        visium_gradients[["sample_id", "tissue"]]
        .drop_duplicates()
        .groupby("tissue")
        .size()
        .reindex(VISIUM_TISSUES)
        .fillna(0)
        .astype(int)
    )
    xenium_counts = xenium_composition.groupby("treatment").size().reindex(["treatment-naive", "chemoradiotherapy-treated"]).fillna(0).astype(int)
    labels = ["Primary", "Liver", "Lung", "Peritoneal", "Xenium\nnaive", "Xenium\nCRT"]
    values = list(visium_counts.values) + list(xenium_counts.values)
    colors = [VISIUM_COLORS[t] for t in VISIUM_TISSUES] + ["#7A9E9F", "#9C755F"]
    x = np.arange(len(values))
    ax.bar(x, values, color=colors, edgecolor="#333333", linewidth=0.6)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=28, ha="right")
    ax.set_ylabel("sections")
    ax.set_title("A  External validation scale", loc="left", fontweight="bold", pad=8)
    ax.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax.set_axisbelow(True)
    for xi, val in zip(x, values):
        ax.text(xi, val + 0.45, str(int(val)), ha="center", va="bottom", fontsize=7.5)


def draw_visium_heatmap(ax: plt.Axes, context: pd.DataFrame) -> None:
    plot_df = context[context["target_program"].isin(VISIUM_PROGRAMS)].copy()
    matrix = plot_df.pivot(index="target_program", columns="tissue", values="median_delta_vs_random").reindex(
        index=VISIUM_PROGRAMS, columns=VISIUM_TISSUES
    )
    support = (
        plot_df.assign(
            label=lambda x: x["median_delta_vs_random"].map(lambda v: f"{v:.2f}")
            + "\n"
            + x["n_more_negative"].astype(int).astype(str)
            + "/"
            + x["n_samples"].astype(int).astype(str)
        )
        .pivot(index="target_program", columns="tissue", values="label")
        .reindex(index=VISIUM_PROGRAMS, columns=VISIUM_TISSUES)
    )
    im = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-0.50, vmax=0.20, aspect="auto")
    ax.set_title("B  GSE274557 Visium CAF-core validation", loc="left", fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(VISIUM_TISSUES)))
    ax.set_xticklabels(["Primary", "Liver", "Lung", "Peritoneal"], rotation=35, ha="right")
    ax.set_yticks(np.arange(len(VISIUM_PROGRAMS)))
    ax.set_yticklabels([VISIUM_LABELS[p] for p in VISIUM_PROGRAMS])
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iat[i, j]
            text_color = "white" if pd.notna(val) and val <= -0.32 else "#111111"
            ax.text(j, i, support.iat[i, j], ha="center", va="center", fontsize=7, color=text_color)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("median delta vs random")


def draw_xenium_summary(ax: plt.Axes, context: pd.DataFrame) -> None:
    context = context.copy()
    matrix = context.pivot(index="target_program", columns="anchor", values="median_delta_vs_random").reindex(
        index=XENIUM_TARGETS, columns=["CAF_APC", "CAF_SPP1TAM"]
    )
    support = (
        context.assign(
            label=lambda x: x["median_delta_vs_random"].map(lambda v: f"{v:.2f}")
            + "\n"
            + x["n_support"].astype(int).astype(str)
            + "/"
            + x["n_samples"].astype(int).astype(str)
        )
        .pivot(index="target_program", columns="anchor", values="label")
        .reindex(index=XENIUM_TARGETS, columns=["CAF_APC", "CAF_SPP1TAM"])
    )
    im = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax.set_title("C  GSE274673 Xenium fixed-anchor summary", loc="left", fontweight="bold", pad=8)
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels([ANCHOR_LABELS[a] for a in ["CAF_APC", "CAF_SPP1TAM"]], rotation=20, ha="right")
    ax.set_yticks(np.arange(len(XENIUM_TARGETS)))
    ax.set_yticklabels([XENIUM_LABELS[t] for t in XENIUM_TARGETS])
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iat[i, j]
            text_color = "white" if pd.notna(val) and abs(val) >= 0.34 else "#111111"
            ax.text(j, i, support.iat[i, j], ha="center", va="center", fontsize=7, color=text_color)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("median delta vs random")


def draw_xenium_sample_heatmap(ax: plt.Axes, gradients: pd.DataFrame) -> None:
    sample_order = ["GSM8454446", "GSM8454449", "GSM8454450", "GSM8454447", "GSM8454448", "GSM8454451"]
    sample_labels = ["Naive\nP1", "Naive\nP4", "Naive\nP5", "CRT\nP2", "CRT\nP3", "CRT\nP6"]
    target_order = ["SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT", "Tumor_epithelial", "SPP1_tumor_like"]
    sub = gradients[gradients["anchor"].eq("CAF_SPP1TAM")].copy()
    matrix = sub.pivot(index="target_program", columns="geo_accession", values="delta_vs_random_median").reindex(
        index=target_order, columns=sample_order
    )
    im = ax.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax.set_title("D  CAF-SPP1/TAM anchor across sections", loc="left", fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(sample_order)))
    ax.set_xticklabels(sample_labels)
    ax.set_yticks(np.arange(len(target_order)))
    ax.set_yticklabels([XENIUM_LABELS[t] for t in target_order])
    ax.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            val = matrix.iat[i, j]
            text_color = "white" if pd.notna(val) and abs(val) >= 0.34 else "#111111"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.8, color=text_color)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    cb.set_label("delta vs random")


def main() -> None:
    setup_style()
    visium_context = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_context_summary.csv")
    visium_gradients = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_gradients.csv")
    xenium_context = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_context_summary.csv")
    xenium_gradients = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_gradients.csv")
    xenium_composition = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")

    scale_source = pd.concat(
        [
            visium_gradients[["sample_id", "tissue"]].drop_duplicates().assign(dataset="GSE274557", context=lambda d: d["tissue"])[
                ["dataset", "context", "sample_id"]
            ],
            xenium_composition.assign(dataset="GSE274673", context=lambda d: d["treatment"], sample_id=lambda d: d["geo_accession"])[
                ["dataset", "context", "sample_id"]
            ],
        ],
        ignore_index=True,
    )
    scale_source.to_csv(SOURCE_DIR / "Source_Data_Fig_4A_multiresolution_scale.csv", index=False)
    visium_context.to_csv(SOURCE_DIR / "Source_Data_Fig_4B_GSE274557.csv", index=False)
    xenium_context.to_csv(SOURCE_DIR / "Source_Data_Fig_4C_GSE274673.csv", index=False)
    xenium_gradients[xenium_gradients["anchor"].eq("CAF_SPP1TAM")].to_csv(
        SOURCE_DIR / "Source_Data_Fig_4D_GSE274673.csv", index=False
    )

    fig = plt.figure(figsize=(13.2, 8.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[0.90, 1.45], height_ratios=[1.0, 1.0], wspace=0.45, hspace=0.55)
    fig.suptitle("Independent multi-resolution validation of CAF-domain organization", fontsize=15, fontweight="bold", y=0.965)

    ax_a = fig.add_subplot(gs[0, 0])
    draw_cohort_scale(ax_a, visium_gradients, xenium_composition)

    ax_b = fig.add_subplot(gs[0, 1])
    draw_visium_heatmap(ax_b, visium_context)

    ax_c = fig.add_subplot(gs[1, 0])
    draw_xenium_summary(ax_c, xenium_context)

    ax_d = fig.add_subplot(gs[1, 1])
    draw_xenium_sample_heatmap(ax_d, xenium_gradients)

    for ext in ["pdf", "png", "svg"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"Wrote {OUT}.pdf/.png/.svg")


if __name__ == "__main__":
    main()
