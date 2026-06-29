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
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from PIL import Image


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
SOURCE_DIR = PROJECT / "results" / "source_data"

DISTANCE_OUT = FIG_DIR / "extended_data_figure12_distance_to_caf_core_dynamics"
XENIUM_OUT = FIG_DIR / "extended_data_figure13_xenium_cell_domain_maps"
ATLAS_OUT = FIG_DIR / "extended_data_figure14_spatial_atlas_overview"

TARGET_ORDER = ["ifn_mhc", "immune_core", "tumor_aggressive", "immune_maturity"]
TARGET_LABELS = {
    "ifn_mhc": "IFN/MHC",
    "immune_core": "immune core",
    "tumor_aggressive": "tumor aggressive",
    "immune_maturity": "immune maturity-like",
}
TARGET_COLORS = {
    "ifn_mhc": "#3B6FB6",
    "immune_core": "#2C7A51",
    "tumor_aggressive": "#B23A48",
    "immune_maturity": "#7B68A6",
}
BIN_ORDER = ["core_0_1.5", "near_1.5_3", "mid_3_6", "far_gt6"]
BIN_LABELS = ["Core", "Near", "Mid", "Far"]


def clean_axes(ax: plt.Axes, grid_axis: str = "y") -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=8)
    ax.grid(axis=grid_axis, color="#E6E6E6", linewidth=0.6, zorder=0)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.10, 1.07, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def save_figure(fig: plt.Figure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = out.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def prepare_distance_data() -> pd.DataFrame:
    mvp = pd.read_csv(TABLE_DIR / "caf_myeloid_niche_distance_bins.csv")
    mvp["context"] = mvp["dataset_id"].map(
        {"GSE282302": "post-NACT", "GSE274103": "treatment-naive"}
    )
    mvp["source_group"] = "Discovery/support"
    gse = pd.read_csv(TABLE_DIR / "gse272362_rds_caf_myeloid_distance_bins.csv")
    gse["context"] = gse["specimen_type"].map(
        {
            "primary_tumor": "primary",
            "liver_metastasis": "liver met",
            "lymph_node_metastasis": "LN met",
            "normal_pancreas": "normal pancreas",
        }
    )
    gse["source_group"] = "GSE272362"
    cols = [
        "dataset_id",
        "sample_id",
        "context",
        "source_group",
        "target",
        "distance_bin",
        "n_spots_bin",
        "median_score",
    ]
    combined = pd.concat([mvp[cols], gse[cols]], ignore_index=True)
    combined = combined[combined["target"].isin(TARGET_ORDER)]
    combined = combined[combined["distance_bin"].isin(BIN_ORDER)]
    combined["target_label"] = combined["target"].map(TARGET_LABELS)
    combined["bin_index"] = combined["distance_bin"].map({b: i for i, b in enumerate(BIN_ORDER)})
    combined["bin_label"] = combined["distance_bin"].map(dict(zip(BIN_ORDER, BIN_LABELS)))
    return combined


def summarize_distance_lines(df: pd.DataFrame, contexts: list[str]) -> pd.DataFrame:
    sub = df[df["context"].isin(contexts)].copy()
    return (
        sub.groupby(["target", "target_label", "distance_bin", "bin_index", "bin_label"], as_index=False)
        .agg(
            median=("median_score", "median"),
            q25=("median_score", lambda s: np.nanquantile(s, 0.25)),
            q75=("median_score", lambda s: np.nanquantile(s, 0.75)),
            n_samples=("sample_id", "nunique"),
        )
        .sort_values(["target", "bin_index"])
    )


def plot_distance_lines(ax: plt.Axes, summary: pd.DataFrame, title: str) -> None:
    for target in TARGET_ORDER:
        g = summary[summary["target"].eq(target)].sort_values("bin_index")
        if g.empty:
            continue
        x = g["bin_index"].to_numpy(float)
        y = g["median"].to_numpy(float)
        ax.plot(
            x,
            y,
            marker="o",
            ms=4,
            lw=1.8,
            color=TARGET_COLORS[target],
            label=TARGET_LABELS[target],
            zorder=3,
        )
        ax.fill_between(
            x,
            g["q25"].to_numpy(float),
            g["q75"].to_numpy(float),
            color=TARGET_COLORS[target],
            alpha=0.12,
            linewidth=0,
            zorder=2,
        )
    ax.axhline(0, color="#333333", lw=0.75)
    ax.set_xticks(np.arange(len(BIN_LABELS)), BIN_LABELS)
    ax.set_ylabel("median program score", fontsize=8)
    ax.set_title(title, fontsize=9.5, fontweight="bold", loc="left")
    clean_axes(ax, grid_axis="y")


def plot_core_to_far_heatmap(ax: plt.Axes, df: pd.DataFrame) -> pd.DataFrame:
    contexts = ["post-NACT", "treatment-naive", "primary", "liver met", "LN met"]
    rows = []
    for context in contexts:
        sub = df[df["context"].eq(context)]
        for target in TARGET_ORDER:
            g = sub[sub["target"].eq(target)]
            core = g[g["distance_bin"].eq("core_0_1.5")]["median_score"].median()
            far = g[g["distance_bin"].eq("far_gt6")]["median_score"].median()
            rows.append(
                {
                    "context": context,
                    "target": target,
                    "target_label": TARGET_LABELS[target],
                    "core_minus_far_median": core - far,
                    "n_samples": g["sample_id"].nunique(),
                }
            )
    out = pd.DataFrame(rows)
    mat = out.pivot(index="target_label", columns="context", values="core_minus_far_median").reindex(
        index=[TARGET_LABELS[t] for t in TARGET_ORDER], columns=contexts
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.25, vmax=0.75, aspect="auto")
    ax.set_xticks(np.arange(len(contexts)), contexts, rotation=30, ha="right")
    ax.set_yticks(np.arange(len(mat.index)), mat.index)
    ax.set_title("Core-to-far enrichment", fontsize=9.5, fontweight="bold", loc="left")
    ax.tick_params(length=0, labelsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.set_label("core minus far median", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    return out


def make_distance_figure() -> None:
    df = prepare_distance_data()
    summary_all = []
    panel_defs = {
        "A": ("Discovery/support cohorts", ["post-NACT", "treatment-naive"]),
        "B": ("GSE272362 primary + liver", ["primary", "liver met"]),
        "C": ("GSE272362 lymph node", ["LN met"]),
    }
    fig = plt.figure(figsize=(12.6, 7.2))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.0, 1.03], hspace=0.55, wspace=0.42)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax, (label, (title, contexts)) in zip(axes, panel_defs.items()):
        panel_label(ax, label)
        summary = summarize_distance_lines(df, contexts)
        summary["panel"] = label
        summary["panel_context"] = title
        summary_all.append(summary)
        plot_distance_lines(ax, summary, title)
    axes[0].legend(frameon=False, fontsize=7, loc="lower left", ncol=1)
    ax_d = fig.add_subplot(gs[1, :])
    panel_label(ax_d, "D")
    heat = plot_core_to_far_heatmap(ax_d, df)
    fig.suptitle("Distance-to-CAF-core dynamics of spatial programs", fontsize=13.5, fontweight="bold", y=0.99)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    pd.concat(summary_all, ignore_index=True).to_csv(
        SOURCE_DIR / "Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv", index=False
    )
    heat.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_12D_core_to_far.csv", index=False)
    save_figure(fig, DISTANCE_OUT)


def scale_values(values: pd.Series) -> np.ndarray:
    arr = values.astype(float).to_numpy()
    lo, hi = np.nanpercentile(arr, [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        lo, hi = np.nanmin(arr), np.nanmax(arr)
    if lo == hi:
        return np.zeros_like(arr)
    return np.clip((arr - lo) / (hi - lo), 0, 1)


def plot_cell_map(ax: plt.Axes, df: pd.DataFrame, column: str, title: str, cmap: str = "viridis") -> None:
    plot_df = df.copy()
    if len(plot_df) > 80000:
        plot_df = plot_df.sample(n=80000, random_state=17)
    colors = scale_values(plot_df[column])
    ax.scatter(
        plot_df["x_centroid"],
        plot_df["y_centroid"],
        c=colors,
        s=0.28,
        cmap=cmap,
        linewidths=0,
        alpha=0.82,
        rasterized=True,
    )
    if column == "anchor_CAF_SPP1TAM":
        cutoff = df[column].quantile(0.90)
        anchors = df[df[column].ge(cutoff)]
        if len(anchors) > 18000:
            anchors = anchors.sample(n=18000, random_state=31)
        ax.scatter(
            anchors["x_centroid"],
            anchors["y_centroid"],
            s=0.7,
            facecolors="none",
            edgecolors="#111111",
            linewidths=0.12,
            alpha=0.45,
            rasterized=True,
        )
    ax.set_title(title, fontsize=8.8, fontweight="bold", loc="left")
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def make_xenium_figure() -> None:
    cells = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_cell_scores.csv")
    comp = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")
    selected = (
        comp.sort_values("n_cells", ascending=False)
        .groupby("treatment", as_index=False)
        .head(1)
        .sort_values("treatment", ascending=False)
    )
    selected_ids = selected["geo_accession"].tolist()
    source = cells[cells["geo_accession"].isin(selected_ids)].copy()
    source.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_13_selected_xenium_cell_scores.csv", index=False)
    fig, axes = plt.subplots(2, 4, figsize=(13.0, 7.2), constrained_layout=False)
    columns = [
        ("anchor_CAF_SPP1TAM", "CAF-SPP1/TAM anchor"),
        ("score_SPP1_TAM", "SPP1/TAM program"),
        ("score_IFN_APC", "IFN/APC program"),
        ("score_Tumor_epithelial", "tumor epithelial"),
    ]
    row_labels = []
    for row, sample_id in enumerate(selected_ids):
        sub = source[source["geo_accession"].eq(sample_id)].copy()
        title = sub["title"].iloc[0].replace("Patient ", "P")
        treatment = sub["treatment"].iloc[0]
        row_labels.append(f"{title} | {treatment}")
        for col, (column, label) in enumerate(columns):
            plot_cell_map(axes[row, col], sub, column, label)
            if row == 0:
                panel_label(axes[row, col], chr(ord("A") + col))
        axes[row, 0].text(
            -0.08,
            0.5,
            row_labels[-1],
            transform=axes[row, 0].transAxes,
            rotation=90,
            ha="center",
            va="center",
            fontsize=8.5,
            fontweight="bold",
        )
    fig.suptitle("Cell-resolution CAF-domain maps in GSE274673 Xenium sections", fontsize=13.5, fontweight="bold", y=0.985)
    fig.text(0.50, 0.03, "Cell colors are per-panel scores clipped at the 2nd-98th percentiles; black outlines mark top 10% CAF-SPP1/TAM anchor cells.", ha="center", fontsize=8)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.91, bottom=0.08, wspace=0.05, hspace=0.08)
    save_figure(fig, XENIUM_OUT)


def build_atlas_counts() -> pd.DataFrame:
    rows = []
    summary = pd.read_csv(TABLE_DIR / "submission_cohort_summary.csv")
    for _, row in summary.iterrows():
        rows.append(
            {
                "cohort": row["cohort"],
                "context": row["specimen_group"],
                "technology": "Visium/ST-H&E",
                "n_sections": int(row["n_samples"]),
                "n_spots_or_cells": int(row["n_spots"]),
                "role": row["analysis_role"],
            }
        )
    g274557 = pd.read_csv(TABLE_DIR / "gse274557_full_spot_scores.csv", low_memory=False)
    for tissue, sub in g274557.groupby("tissue"):
        rows.append(
            {
                "cohort": "GSE274557",
                "context": tissue,
                "technology": "Visium",
                "n_sections": sub["sample_id"].nunique(),
                "n_spots_or_cells": len(sub),
                "role": "Independent broad external metastatic validation.",
            }
        )
    comp = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")
    for treatment, sub in comp.groupby("treatment"):
        rows.append(
            {
                "cohort": "GSE274673",
                "context": treatment,
                "technology": "Xenium",
                "n_sections": len(sub),
                "n_spots_or_cells": int(sub["n_cells"].sum()),
                "role": "Cell-resolution CAF-domain validation.",
            }
        )
    return pd.DataFrame(rows)


def plot_atlas_counts(ax: plt.Axes, atlas: pd.DataFrame) -> None:
    order = [
        "GSE282302",
        "GSE274103",
        "GSE272362",
        "GSE235315",
        "GSE274557",
        "GSE274673",
    ]
    totals = atlas.groupby("cohort", as_index=False).agg(n_sections=("n_sections", "sum"), n_units=("n_spots_or_cells", "sum"))
    totals["cohort"] = pd.Categorical(totals["cohort"], categories=order, ordered=True)
    totals = totals.sort_values("cohort").reset_index(drop=True)
    colors = ["#4C78A8", "#72B7B2", "#59A14F", "#B279A2", "#E15759", "#9C755F"]
    y = np.arange(len(totals))
    ax.barh(y, totals["n_sections"], color=colors[: len(totals)], alpha=0.90)
    ax.set_yticks(y, totals["cohort"].astype(str))
    ax.invert_yaxis()
    ax.set_xlabel("sections / samples", fontsize=8)
    ax.set_title("Cohort scale", fontsize=9.5, fontweight="bold", loc="left")
    clean_axes(ax, grid_axis="x")
    for i, row in totals.iterrows():
        ax.text(row["n_sections"] + 1.0, i, str(int(row["n_sections"])), va="center", fontsize=7.5)


def plot_evidence_matrix(ax: plt.Axes) -> None:
    rows = ["CAF-core discovery", "metastatic contrast", "external Visium", "Xenium cells", "H&E bridge"]
    cols = ["GSE282302", "GSE274103", "GSE272362", "GSE274557", "GSE274673", "GSE235315"]
    val = np.array(
        [
            [2, 1, 1, 0, 0, 1],
            [0, 0, 2, 0, 0, 0],
            [0, 0, 0, 2, 0, 0],
            [0, 0, 0, 0, 2, 0],
            [2, 1, 0, 0, 0, 0],
        ],
        dtype=float,
    )
    cmap = ListedColormap(["#F5F5F5", "#9CC9E2", "#2F6C9F"])
    ax.imshow(val, cmap=cmap, vmin=0, vmax=2, aspect="auto")
    ax.set_xticks(np.arange(len(cols)), cols, rotation=35, ha="right")
    ax.set_yticks(np.arange(len(rows)), rows)
    ax.tick_params(length=0, labelsize=7.5)
    ax.set_title("Evidence role by cohort", fontsize=9.5, fontweight="bold", loc="left")
    for i in range(val.shape[0]):
        for j in range(val.shape[1]):
            if val[i, j] > 0:
                ax.text(j, i, "main" if val[i, j] == 2 else "support", ha="center", va="center", fontsize=6.8)


def crop_to_square(img: Image.Image) -> Image.Image:
    w, h = img.size
    side = min(w, h)
    left = (w - side) // 2
    top = (h - side) // 2
    return img.crop((left, top, left + side, top + side))


def plot_thumbnail_strip(ax: plt.Axes) -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    manifest = pd.read_csv(TABLE_DIR / "gse272362_rds_overlay_manifest.csv")
    picks = [
        ("primary", "primary_tumor"),
        ("liver", "liver_metastasis"),
        ("LN", "lymph_node_metastasis"),
    ]
    paths = []
    for label, specimen in picks:
        sub = manifest[manifest["specimen_type"].eq(specimen)]
        if len(sub):
            paths.append((label, Path(sub.iloc[0]["image_path"])))
    mvp = pd.read_csv(TABLE_DIR / "mvp_overlay_manifest.csv")
    if len(mvp):
        paths.insert(0, ("post-NACT", Path(mvp.iloc[0]["image_path"])))
    x = 0.02
    for label, path in paths:
        if not path.exists():
            continue
        img = Image.open(path).convert("RGB")
        img = crop_to_square(img).resize((170, 170), Image.Resampling.LANCZOS)
        ax.imshow(img, extent=(x, x + 0.22, 0.12, 0.86), aspect="auto", zorder=2)
        ax.add_patch(plt.Rectangle((x, 0.12), 0.22, 0.74, fill=False, edgecolor="#333333", linewidth=0.5, zorder=3))
        ax.text(x + 0.11, 0.91, label, ha="center", va="bottom", fontsize=7.5, fontweight="bold")
        x += 0.24
    ax.set_title("Representative H&E sections", fontsize=9.5, fontweight="bold", loc="left")


def plot_ecotype_composition(ax: plt.Axes) -> None:
    df = pd.read_csv(TABLE_DIR / "spatial_ecotype_context_counts.csv")
    contexts = ["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
    labels = ["post-NACT", "treatment-naive", "primary", "liver met", "LN met"]
    ecotypes = ["NMF1", "NMF2", "NMF3", "NMF4"]
    colors = {"NMF1": "#B23A48", "NMF2": "#59A14F", "NMF3": "#3B6FB6", "NMF4": "#9C755F"}
    mat = (
        df.pivot_table(index="cohort_context", columns="dominant_nmf_ecotype", values="n_samples", aggfunc="sum")
        .reindex(index=contexts, columns=ecotypes)
        .fillna(0)
    )
    frac = mat.div(mat.sum(axis=1), axis=0).fillna(0)
    bottom = np.zeros(len(frac))
    x = np.arange(len(frac))
    for ecotype in ecotypes:
        ax.bar(x, frac[ecotype], bottom=bottom, color=colors[ecotype], label=ecotype, width=0.72)
        bottom += frac[ecotype].to_numpy()
    ax.set_xticks(x, labels, rotation=28, ha="right")
    ax.set_ylim(0, 1)
    ax.set_ylabel("fraction of samples", fontsize=8)
    ax.set_title("CAF-core ecotype composition", fontsize=9.5, fontweight="bold", loc="left")
    clean_axes(ax, grid_axis="y")
    ax.legend(frameon=False, fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.18))


def make_atlas_figure() -> None:
    atlas = build_atlas_counts()
    atlas.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_14A_atlas_counts.csv", index=False)
    fig = plt.figure(figsize=(13.2, 8.2))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05], width_ratios=[0.9, 1.35], hspace=0.56, wspace=0.40)
    ax_a = fig.add_subplot(gs[0, 0])
    panel_label(ax_a, "A")
    plot_atlas_counts(ax_a, atlas)
    ax_b = fig.add_subplot(gs[0, 1])
    panel_label(ax_b, "B")
    plot_evidence_matrix(ax_b)
    ax_c = fig.add_subplot(gs[1, 0])
    panel_label(ax_c, "C")
    plot_ecotype_composition(ax_c)
    ax_d = fig.add_subplot(gs[1, 1])
    panel_label(ax_d, "D")
    plot_thumbnail_strip(ax_d)
    fig.suptitle("PDAC spatial ecology atlas and figure grammar summary", fontsize=13.5, fontweight="bold", y=0.985)
    save_figure(fig, ATLAS_OUT)


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    make_distance_figure()
    make_xenium_figure()
    make_atlas_figure()
    print(f"Wrote {DISTANCE_OUT.with_suffix('.pdf')}")
    print(f"Wrote {XENIUM_OUT.with_suffix('.pdf')}")
    print(f"Wrote {ATLAS_OUT.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
