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
        "axes.linewidth": 0.65,
        "xtick.major.width": 0.55,
        "ytick.major.width": 0.55,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [SOURCE_DIR, FIG_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "extended_data_figure25_spatial_robustness_module_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Extended_Data_Fig_25_spatial_robustness_module.csv"
REPORT_OUT = REPORT_DIR / "extended_data_figure25_spatial_robustness_module_notes.md"

TARGET_LABELS = {
    "ifn_mhc": "IFN/MHC",
    "immune_core": "Immune core",
    "immune_maturity": "Immune maturity",
    "immune_maturity_like": "Immune maturity",
    "tumor_aggressive": "Tumor aggressive",
    "SPP1_TAM": "SPP1/TAM",
    "TGFb_EMT": "TGFb/EMT",
    "DC_APC": "DC/APC",
    "T_NK": "T/NK",
    "B_plasma": "B/plasma",
    "IFN_APC": "IFN/APC",
}

PROGRAM_COLORS = {
    "IFN/MHC": "#3B6FB6",
    "Immune core": "#2C7A51",
    "Immune maturity": "#7B68A6",
    "Tumor aggressive": "#B23A48",
    "SPP1/TAM": "#8A5A44",
    "TGFb/EMT": "#C77C2D",
    "IFN/APC": "#4C78A8",
    "DC/APC": "#72B7B2",
    "T/NK": "#2C7A51",
    "B/plasma": "#7B68A6",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def context_label(value: str) -> str:
    labels = {
        "GSE282302": "post-NACT\nGSE282302",
        "GSE274103": "treat-naive\nGSE274103",
        "primary_tumor": "primary\nGSE272362",
        "liver_metastasis": "liver met\nGSE272362",
        "lymph_node_metastasis": "LN met\nGSE272362",
        "GSE235315": "external\nGSE235315",
        "post-NACT": "post-NACT",
        "treatment-naive": "treat-naive",
        "primary": "primary",
        "liver met": "liver met",
        "LN met": "LN met",
    }
    return labels.get(value, str(value).replace("_", " "))


def load_random_support() -> pd.DataFrame:
    mvp = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Random_Core_MVP.csv")
    mvp = mvp.assign(context=mvp["dataset_id"])
    mvp["target_label"] = mvp["target"].map(TARGET_LABELS)
    mvp = mvp.rename(columns={"median_delta_vs_null": "delta", "n_observed_more_negative_than_null": "n_support"})

    gse = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Random_Core_GSE272362.csv")
    gse = gse.assign(context=gse["specimen_type"])
    gse["target_label"] = gse["target"].map(TARGET_LABELS)
    gse = gse.rename(columns={"median_delta_vs_null": "delta", "n_observed_more_negative_than_null": "n_support"})

    anchor = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_6A.csv")
    anchor = anchor.assign(context="GSE235315")
    anchor["target_label"] = anchor["target_label"].replace({"immune core": "Immune core", "tumor aggressive": "Tumor aggressive"})
    anchor = (
        anchor.groupby(["context", "target_label"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            delta=("delta_vs_null_median", "median"),
            n_support=("observed_more_negative_than_null", lambda s: int(s.astype(bool).sum())),
        )
    )

    cols = ["context", "target_label", "n_samples", "delta", "n_support"]
    out = pd.concat([mvp[cols], gse[cols], anchor[cols]], ignore_index=True)
    return out[out["target_label"].isin(["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"])].copy()


def plot_random_heatmap(ax: plt.Axes, data: pd.DataFrame) -> None:
    contexts = ["GSE282302", "GSE274103", "primary_tumor", "liver_metastasis", "lymph_node_metastasis", "GSE235315"]
    targets = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"]
    mat = data.pivot_table(index="target_label", columns="context", values="delta", aggfunc="median").reindex(targets, columns=contexts)
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax.set_xticks(np.arange(len(contexts)))
    ax.set_xticklabels([context_label(c) for c in contexts], rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(targets)))
    ax.set_yticklabels(targets, fontsize=7)
    ax.set_title("Random-core specificity across cohorts", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "A")
    for i, target in enumerate(targets):
        for j, context in enumerate(contexts):
            val = mat.loc[target, context]
            if not np.isfinite(val):
                continue
            row = data[(data["target_label"].eq(target)) & (data["context"].eq(context))]
            label = f"{val:.2f}\n{int(row['n_support'].iloc[0])}/{int(row['n_samples'].iloc[0])}"
            ax.text(j, i, label, ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > 0.25 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("median delta vs random", fontsize=6)
    cbar.ax.tick_params(labelsize=6)


def plot_threshold(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Threshold_Sensitivity.csv")
    df["target_label"] = df["target"].map(TARGET_LABELS)
    df = df[df["target_label"].isin(["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"])].copy()
    order = ["top_15", "top_10", "top_5"]
    x = np.arange(len(order))
    summary = df.groupby(["target_label", "core_label"], as_index=False).agg(median_rho=("median_rho", "median"))
    for target in ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"]:
        sub = summary[summary["target_label"].eq(target)].set_index("core_label").reindex(order)
        ax.plot(x, sub["median_rho"], marker="o", lw=1.5, ms=3.5, color=PROGRAM_COLORS[target], label=target)
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels(["top 15%", "top 10%", "top 5%"], fontsize=7)
    ax.set_ylabel("median rho to CAF core", fontsize=7)
    ax.set_title("Core-threshold sensitivity", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "B")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=6, loc="lower left", ncol=1)
    return summary


def plot_distance(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv")
    df = df[df["target_label"].isin(["IFN/MHC", "immune core", "tumor aggressive", "immune maturity-like"])].copy()
    df["target_label"] = df["target_label"].replace(
        {"immune core": "Immune core", "tumor aggressive": "Tumor aggressive", "immune maturity-like": "Immune maturity"}
    )
    df = df[df["panel_context"].eq("Discovery/support cohorts")]
    order = ["Core", "Near", "Mid", "Far"]
    x = np.arange(len(order))
    for target in ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"]:
        sub = df[df["target_label"].eq(target)].sort_values("bin_index")
        ax.plot(x, sub["median"], marker="o", lw=1.5, ms=3.5, color=PROGRAM_COLORS[target], label=target)
        ax.fill_between(x, sub["q25"], sub["q75"], color=PROGRAM_COLORS[target], alpha=0.12, linewidth=0)
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels(order, fontsize=7)
    ax.set_ylabel("module score", fontsize=7)
    ax.set_title("Distance-gradient dynamics", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "C")
    clean_axes(ax)
    return df


def plot_core_far_heatmap(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_12D_core_to_far.csv")
    df["target_label"] = df["target_label"].replace(
        {"immune core": "Immune core", "tumor aggressive": "Tumor aggressive", "immune maturity-like": "Immune maturity"}
    )
    df = df[df["target_label"].isin(["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"])].copy()
    contexts = ["post-NACT", "treatment-naive", "primary", "liver met", "LN met"]
    targets = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity"]
    mat = df.pivot_table(index="target_label", columns="context", values="core_minus_far_median", aggfunc="median").reindex(targets, columns=contexts)
    im = ax.imshow(mat.to_numpy(float), cmap="YlGnBu", vmin=0, vmax=max(0.8, np.nanmax(mat.to_numpy(float))), aspect="auto")
    ax.set_xticks(np.arange(len(contexts)))
    ax.set_xticklabels([context_label(c) for c in contexts], rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(targets)))
    ax.set_yticklabels(targets, fontsize=7)
    ax.set_title("Core-to-far effect by context", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "D")
    for i, target in enumerate(targets):
        for j, context in enumerate(contexts):
            val = mat.loc[target, context]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="#FFFFFF" if val > 0.45 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("core - far", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_gse274557(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_4B_GSE274557.csv")
    df = df[df["target_program"].isin(["IFN/MHC", "SPP1_TAM", "TGFb_EMT", "DC_APC", "T_NK", "B_plasma"])].copy()
    df["target_label"] = df["target_program"].map(TARGET_LABELS).fillna(df["target_program"])
    summary = df.groupby("target_label", as_index=False).agg(
        median_delta=("median_delta_vs_random", "median"),
        n_support=("n_more_negative", "sum"),
        n_samples=("n_samples", "sum"),
    )
    summary = summary.sort_values("median_delta")
    y = np.arange(len(summary))
    ax.axvline(0, color="#333333", lw=0.65)
    colors = [PROGRAM_COLORS.get(x, "#8A8F98") for x in summary["target_label"]]
    ax.hlines(y, 0, summary["median_delta"], color=colors, lw=1.8)
    ax.scatter(summary["median_delta"], y, s=34, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(summary["target_label"], fontsize=6.7)
    ax.set_xlim(min(-0.42, float(summary["median_delta"].min()) - 0.05), 0.04)
    ax.set_xlabel("median delta vs random", fontsize=7)
    ax.set_title("GSE274557 external Visium", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax)
    for i, row in summary.reset_index(drop=True).iterrows():
        ax.text(row["median_delta"] + 0.035, i, f"{int(row['n_support'])}/{int(row['n_samples'])}", ha="left", va="center", fontsize=6.2)
    return summary


def plot_gse274673(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_4C_GSE274673.csv")
    df = df[df["target_program"].isin(["IFN_APC", "SPP1_TAM", "TGFb_EMT", "DC_APC", "T_NK", "B_plasma", "SPP1_tumor_like"])].copy()
    df["target_label"] = df["target_program"].map(TARGET_LABELS).fillna(df["target_program"].str.replace("_", " "))
    pivot = df.pivot_table(index="target_label", columns="anchor", values="median_delta_vs_random", aggfunc="median")
    order = ["IFN/APC", "SPP1/TAM", "TGFb/EMT", "DC/APC", "T/NK", "B/plasma", "SPP1 tumor like"]
    pivot = pivot.reindex([x for x in order if x in pivot.index])
    im = ax.imshow(pivot.to_numpy(float), cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(np.arange(len(pivot.columns)))
    ax.set_xticklabels([str(c).replace("_", "\n") for c in pivot.columns], fontsize=7)
    ax.set_yticks(np.arange(len(pivot.index)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_title("GSE274673 Xenium cell-resolution", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "F")
    for i, idx in enumerate(pivot.index):
        for j, col in enumerate(pivot.columns):
            val = pivot.loc[idx, col]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="#FFFFFF" if abs(val) > 0.28 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("delta vs random", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_support_fraction(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_ED_Fig_5D_cross_context_support_fraction.csv")
    df["target"] = df["target"].replace({"immune core": "Immune core", "tumor aggressive": "Tumor aggressive"})
    short = {
        "Tumor aggressive": "Tumor aggr.",
        "IFN/MHC": "IFN/MHC",
        "Immune core": "Immune core",
        "Immune maturity-like": "Immune mat.",
    }
    df = df.sort_values("support_fraction")
    y = np.arange(len(df))
    colors = [PROGRAM_COLORS.get(x, "#8A8F98") for x in df["target"]]
    ax.barh(y, df["support_fraction"], color=colors, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels([short.get(x, x) for x in df["target"]], fontsize=6.6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("support fraction", fontsize=7)
    ax.set_title("Cross-context support", loc="left", fontsize=9.5, fontweight="bold")
    panel_label(ax, "G", x=-0.24)
    clean_axes(ax)
    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(row["support_fraction"] + 0.025, i, f"{int(row['n_support'])}/{int(row['n_samples'])}", va="center", fontsize=6.5)
    return df


def write_source_data(random_support: pd.DataFrame, threshold: pd.DataFrame, distance: pd.DataFrame, core_far: pd.DataFrame, g274557: pd.DataFrame, g274673: pd.DataFrame, support_fraction: pd.DataFrame) -> None:
    rows: list[dict[str, object]] = []
    for _, row in random_support.iterrows():
        rows.append({"panel": "A", "source": "Source_Data_Extended_Random_Core_* plus Source_Data_Fig_6A", "item": f"{row['context']}|{row['target_label']}", "metric": "median_delta_vs_random", "value": row["delta"]})
    for _, row in threshold.iterrows():
        rows.append({"panel": "B", "source": "Source_Data_Extended_Threshold_Sensitivity.csv", "item": f"{row['target_label']}|{row['core_label']}", "metric": "median_rho_across_contexts", "value": row["median_rho"]})
    for _, row in distance.iterrows():
        rows.append({"panel": "C", "source": "Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv", "item": f"{row['target_label']}|{row['bin_label']}", "metric": "median_module_score", "value": row["median"]})
    for _, row in core_far.iterrows():
        rows.append({"panel": "D", "source": "Source_Data_Extended_Data_Fig_12D_core_to_far.csv", "item": f"{row['context']}|{row['target_label']}", "metric": "core_minus_far_median", "value": row["core_minus_far_median"]})
    for _, row in g274557.iterrows():
        rows.append({"panel": "E", "source": "Source_Data_Fig_4B_GSE274557.csv", "item": row["target_label"], "metric": "median_delta_vs_random", "value": row["median_delta"]})
    for _, row in g274673.iterrows():
        rows.append({"panel": "F", "source": "Source_Data_Fig_4C_GSE274673.csv", "item": f"{row['anchor']}|{row['target_label']}", "metric": "median_delta_vs_random", "value": row["median_delta_vs_random"]})
    for _, row in support_fraction.iterrows():
        rows.append({"panel": "G", "source": "Source_Data_ED_Fig_5D_cross_context_support_fraction.csv", "item": row["target"], "metric": "support_fraction", "value": row["support_fraction"]})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Extended Data Figure 25 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "NC-style robustness module for the CAF-myeloid spatial niche model. The figure consolidates random-core specificity, threshold stability, distance-gradient behavior, core-to-far effects and independent Visium/Xenium validation summaries.\n\n"
        "## Panel contract\n\n"
        "- A: cross-cohort random-core support heatmap; each cell shows median delta and supported samples.\n"
        "- B: robustness of the core-distance association under top 15%, top 10% and top 5% CAF-core definitions.\n"
        "- C: discovery/support distance-gradient dynamics from core to far regions.\n"
        "- D: context-level core-to-far effect sizes.\n"
        "- E: independent GSE274557 Visium support summarized across tissue contexts.\n"
        "- F: GSE274673 Xenium cell-resolution support by CAF-domain anchor.\n"
        "- G: cross-context support fraction for the main programs.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf')}`\n"
        f"- `{OUT.with_suffix('.svg')}`\n"
        f"- `{OUT.with_suffix('.png')}`\n"
        f"- `{SOURCE_OUT}`\n",
        encoding="utf-8",
    )


def main() -> None:
    random_support = load_random_support()
    fig = plt.figure(figsize=(14.8, 13.2), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[1.15, 1.0, 1.0, 0.92], hspace=0.95, wspace=0.9)
    fig.suptitle("Spatial robustness module for the CAF-myeloid niche model", fontsize=15, fontweight="bold", y=0.985)

    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 4:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:6])
    ax_e = fig.add_subplot(gs[2, 0:2])
    ax_f = fig.add_subplot(gs[2, 2:5])
    ax_g = fig.add_subplot(gs[2, 5:6])
    ax_note = fig.add_subplot(gs[3, :])

    plot_random_heatmap(ax_a, random_support)
    threshold = plot_threshold(ax_b)
    distance = plot_distance(ax_c)
    core_far = plot_core_far_heatmap(ax_d)
    g274557 = plot_gse274557(ax_e)
    g274673 = plot_gse274673(ax_f)
    support_fraction = plot_support_fraction(ax_g)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.82,
        "Interpretation: the CAF-myeloid spatial niche signal is not a single-threshold or single-cohort artifact. "
        "Random-core controls, threshold sensitivity, distance gradients and independent Visium/Xenium summaries support the organizing-core model, "
        "while remaining observational and program-defined.",
        fontsize=8.2,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.35,
        "Boundary: this module supports spatial robustness and external consistency; it does not establish perturbational causality, clinical prediction or histology-annotated compartment identity.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_source_data(random_support, threshold, distance, core_far, g274557, g274673, support_fraction)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
