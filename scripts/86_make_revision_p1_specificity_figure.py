from __future__ import annotations

from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REVISION_ROOT = PROJECT_ROOT / "results" / "revision_2026_06_29"
OUT_BASE = REVISION_ROOT / "figures" / "Extended_Data_Figure_7_Specificity_Sensitivity"
SOURCE_OUT = REVISION_ROOT / "analysis_outputs" / "extended_data_figure7_specificity_sensitivity_source_data.csv"

TARGET_ORDER = ["IFN/MHC", "immune-core", "tumor-aggressive", "SPP1/TAM", "TGF-beta/EMT"]
ANCHOR_ORDER = ["CAF-myeloid combined", "CAF-only", "myeloid-only", "immune-high", "tumor-high"]
CONTEXT_ORDER = [
    ("GSE282302", "metadata_required", "GSE282302\npost-NACT"),
    ("GSE274103", "metadata_required", "GSE274103\ntreatment-naive"),
    ("GSE272362", "primary_tumor", "GSE272362\nprimary"),
    ("GSE272362", "liver_metastasis", "GSE272362\nliver met"),
    ("GSE272362", "lymph_node_metastasis", "GSE272362\nLN met"),
    ("GSE274557", "primary_tumor", "GSE274557\nprimary"),
    ("GSE274557", "liver_metastasis", "GSE274557\nliver met"),
    ("GSE274557", "peritoneal_metastasis", "GSE274557\nperitoneal"),
    ("GSE274557", "lung_metastasis", "GSE274557\nlung met"),
]
PANEL_LETTERS = dict(fontsize=13, fontweight="bold", ha="left", va="top")


def read_table(name: str) -> pd.DataFrame:
    return pd.read_csv(REVISION_ROOT / "analysis_outputs" / name)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, **PANEL_LETTERS)


def plot_contiguous_null(ax: plt.Axes, source_rows: list[pd.DataFrame]) -> None:
    df = read_table("stronger_null_contiguous_random_core_summary.csv")
    df = df[df["target_program"].isin(TARGET_ORDER)].copy()
    label_map = {(dataset, site): label for dataset, site, label in CONTEXT_ORDER}
    df["context_label"] = [label_map.get((r.dataset, r.tissue_site), None) for r in df.itertuples()]
    df = df[df["context_label"].notna()].copy()
    df["target_program"] = pd.Categorical(df["target_program"], categories=TARGET_ORDER, ordered=True)
    df["context_label"] = pd.Categorical(
        df["context_label"], categories=[item[2] for item in CONTEXT_ORDER], ordered=True
    )
    mat = df.pivot(index="context_label", columns="target_program", values="median_delta").reindex(
        index=[item[2] for item in CONTEXT_ORDER], columns=TARGET_ORDER
    )
    support = df.pivot(index="context_label", columns="target_program", values="support_n").reindex_like(mat)
    n_samples = df.pivot(index="context_label", columns="target_program", values="n_samples").reindex_like(mat)

    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.55, vmax=0.2, aspect="auto")
    ax.set_xticks(np.arange(len(TARGET_ORDER)), TARGET_ORDER, rotation=25, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(mat.index)), mat.index, fontsize=8)
    ax.set_title("Contiguous-null specificity", fontsize=10, pad=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if pd.isna(val):
                continue
            text = f"{support.iloc[i, j]:.0f}/{n_samples.iloc[i, j]:.0f}"
            ax.text(j, i, text, ha="center", va="center", fontsize=6.5, color="#111111")
    cbar = plt.colorbar(im, ax=ax, fraction=0.032, pad=0.015)
    cbar.set_label("observed rho - contiguous-null median", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    source = df.copy()
    source["panel"] = "A_contiguous_null"
    source_rows.append(source)


def plot_anchor_components(ax: plt.Axes, source_rows: list[pd.DataFrame]) -> None:
    df = read_table("caf_myeloid_component_anchor_comparison_summary.csv")
    df = df[df["anchor_type"].isin(ANCHOR_ORDER) & df["target_program"].isin(TARGET_ORDER + ["tumor epithelial"])].copy()
    target_order = TARGET_ORDER + ["tumor epithelial"]
    df["anchor_type"] = pd.Categorical(df["anchor_type"], categories=ANCHOR_ORDER, ordered=True)
    df["target_program"] = pd.Categorical(df["target_program"], categories=target_order, ordered=True)
    x_map = {target: i for i, target in enumerate(target_order)}
    y_map = {anchor: i for i, anchor in enumerate(ANCHOR_ORDER)}
    norm = Normalize(vmin=-0.55, vmax=0.15)
    sizes = 35 + 170 * df["support_fraction"].astype(float)
    sc = ax.scatter(
        df["target_program"].map(x_map),
        df["anchor_type"].map(y_map),
        c=df["median_delta"].astype(float),
        s=sizes,
        cmap="RdBu_r",
        norm=norm,
        edgecolor="#222222",
        linewidth=0.35,
    )
    ax.set_xticks(np.arange(len(target_order)), target_order, rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(ANCHOR_ORDER)), ANCHOR_ORDER, fontsize=8)
    ax.invert_yaxis()
    ax.set_title("Anchor-component dissection", fontsize=10, pad=8)
    ax.grid(axis="both", color="#E0E0E0", linewidth=0.45)
    cbar = plt.colorbar(sc, ax=ax, fraction=0.04, pad=0.015)
    cbar.set_label("delta vs random anchor", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    source = df.copy()
    source["panel"] = "B_anchor_components"
    source_rows.append(source)


def plot_overlap_sensitivity(ax: plt.Axes, source_rows: list[pd.DataFrame]) -> None:
    df = read_table("gene_module_overlap_sensitivity_summary.csv")
    df = df[df["target_program"].isin(TARGET_ORDER + ["myCAF/matrix"])].copy()
    order = ["IFN/MHC", "immune-core", "tumor-aggressive", "TGF-beta/EMT", "SPP1/TAM", "myCAF/matrix"]
    df["target_program"] = pd.Categorical(df["target_program"], categories=order, ordered=True)
    df = df.sort_values("target_program")
    y = np.arange(len(df))
    shared_pct = df["percent_target_shared_with_core"].astype(float)
    if shared_pct.max() <= 1.0:
        shared_pct = shared_pct * 100.0
    delta = df["median_overlap_sensitive_delta"].astype(float)
    colors = ["#4C78A8" if p == 0 else "#B04A37" for p in shared_pct]
    ax.barh(y, shared_pct, color=colors, alpha=0.85, height=0.62)
    ax.set_yticks(y, df["target_program"], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("% target genes shared with core", fontsize=8)
    ax.set_xlim(0, 105)
    ax.set_title("Marker-overlap sensitivity", fontsize=10, pad=8)
    ax.tick_params(axis="x", labelsize=7)
    ax2 = ax.twiny()
    ax2.patch.set_visible(False)
    ax2.scatter(delta, y, marker="D", s=34, color="#222222", zorder=3)
    ax2.axvline(0, color="#777777", linestyle="--", linewidth=0.8)
    ax2.set_xlim(-0.55, 0.05)
    ax2.set_xlabel("overlap-sensitive delta", fontsize=8)
    ax2.tick_params(axis="x", labelsize=7)
    for yi, row in enumerate(df.itertuples()):
        pct = float(row.percent_target_shared_with_core)
        if pct <= 1.0:
            pct *= 100.0
        text = "no shared genes" if pct == 0 else f"{int(row.shared_genes_n)}/{int(row.target_genes_n)} shared"
        ax.text(min(pct + 2, 84), yi, text, va="center", fontsize=7, color="#222222")
    source = df.copy()
    source["panel"] = "C_overlap_sensitivity"
    source_rows.append(source)


def plot_ln_leave_one_out(ax: plt.Axes, source_rows: list[pd.DataFrame]) -> None:
    ind = read_table("ln_metastasis_individual_sample_summary.csv")
    loo = read_table("ln_metastasis_leave_one_out_summary.csv")
    programs = [
        ("IFN/MHC", "delta_vs_null_median__IFN_MHC", "IFN_MHC_median_delta"),
        ("immune-core", "delta_vs_null_median__immune_core", "immune_core_median_delta"),
        ("tumor-aggressive", "delta_vs_null_median__tumor_aggressive", "tumor_aggressive_median_delta"),
    ]
    x = np.arange(len(programs))
    palette = {"IFN/MHC": "#4C78A8", "immune-core": "#59A14F", "tumor-aggressive": "#B04A37"}
    for i, (label, ind_col, loo_col) in enumerate(programs):
        vals = ind[ind_col].astype(float).to_numpy()
        jitter = np.linspace(-0.14, 0.14, len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, s=32, color=palette[label], alpha=0.85, edgecolor="white", linewidth=0.4)
        all_val = float(loo.loc[loo["analysis"].eq("all_LN_samples"), loo_col].iloc[0])
        loo_vals = loo.loc[loo["analysis"].eq("leave_one_out"), loo_col].astype(float).to_numpy()
        ax.plot([i - 0.20, i + 0.20], [all_val, all_val], color="#111111", linewidth=2.0)
        ax.vlines(i, loo_vals.min(), loo_vals.max(), color="#111111", linewidth=1.2)
    ax.axhline(0, color="#777777", linestyle="--", linewidth=0.8)
    ax.set_xticks(x, [item[0] for item in programs], rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("delta vs contiguous null", fontsize=8)
    ax.set_title("LN subset leave-one-out", fontsize=10, pad=8)
    ax.tick_params(axis="y", labelsize=7)
    ax.text(
        0.02,
        0.04,
        "points: individual LN samples\nblack bar: all-LN median\nvertical line: leave-one-out range",
        transform=ax.transAxes,
        fontsize=7,
        color="#333333",
        va="bottom",
    )
    source_rows.append(ind.assign(panel="D_LN_individual"))
    source_rows.append(loo.assign(panel="D_LN_leave_one_out"))


def plot_nmf_rank(ax: plt.Axes, source_rows: list[pd.DataFrame]) -> None:
    summary = read_table("nmf_rank_stability_summary.csv")
    reference = read_table("nmf_rank_nndsvda_reference.csv")
    ax.plot(summary["rank"], summary["explained_fraction_mean"], marker="o", color="#246A73", linewidth=2.0, label="randomized mean")
    ax.plot(reference["rank"], reference["explained_fraction"], marker="^", color="#6D597A", linewidth=1.6, label="NNDSVDa reference")
    ax.axvline(4, color="#222222", linestyle="--", linewidth=1.0)
    ax.set_xlabel("NMF rank", fontsize=8)
    ax.set_ylabel("explained fraction", fontsize=8)
    ax.set_ylim(0.955, 0.996)
    ax.tick_params(labelsize=7)
    ax.set_title("NMF rank stability", fontsize=10, pad=8)
    ax2 = ax.twinx()
    ax2.plot(summary["rank"], summary["ari_mean"], marker="s", color="#B04A37", linewidth=1.5, label="ARI")
    ax2.plot(summary["rank"], 1 - summary["consensus_pac_0.1_0.9"], marker="d", color="#59A14F", linewidth=1.5, label="1 - PAC")
    ax2.set_ylabel("assignment stability", fontsize=8)
    ax2.set_ylim(0.94, 1.005)
    ax2.tick_params(labelsize=7)
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines + lines2, labels + labels2, frameon=False, fontsize=7, loc="lower right")
    source_rows.append(summary.assign(panel="E_NMF_rank_randomized"))
    source_rows.append(reference.assign(panel="E_NMF_rank_reference"))


def make_figure() -> None:
    source_rows: list[pd.DataFrame] = []
    fig = plt.figure(figsize=(15.6, 11.2), constrained_layout=False)
    gs = fig.add_gridspec(
        3,
        6,
        height_ratios=[1.28, 1.0, 1.0],
        width_ratios=[1, 1, 1, 1, 1, 1],
        hspace=1.02,
        wspace=0.90,
    )
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, :3])
    ax_c = fig.add_subplot(gs[1, 3:])
    ax_d = fig.add_subplot(gs[2, :3])
    ax_e = fig.add_subplot(gs[2, 3:])

    plot_contiguous_null(ax_a, source_rows)
    plot_anchor_components(ax_b, source_rows)
    plot_overlap_sensitivity(ax_c, source_rows)
    plot_ln_leave_one_out(ax_d, source_rows)
    plot_nmf_rank(ax_e, source_rows)

    for label, ax in zip("ABCDE", [ax_a, ax_b, ax_c, ax_d, ax_e]):
        panel_label(ax, label)

    fig.suptitle(
        "Revision specificity suite for CAF-myeloid spatial architecture",
        fontsize=15,
        fontweight="bold",
        y=0.985,
    )
    fig.text(
        0.5,
        0.012,
        "P1 sensitivity analyses test five reviewer-risk alternatives: arbitrary contiguous tissue regions, single-component anchors, marker overlap, single LN-sample dominance and arbitrary NMF rank. In panels A-B, more negative deltas indicate stronger enrichment near the observed anchor than the matched null.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color="#333333",
    )
    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=320, bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    source = pd.concat(source_rows, ignore_index=True, sort=False)
    source.to_csv(SOURCE_OUT, index=False, encoding="utf-8")
    print(OUT_BASE.with_suffix(".pdf"))
    print(SOURCE_OUT)


if __name__ == "__main__":
    make_figure()
