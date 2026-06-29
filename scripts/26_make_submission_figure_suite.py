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
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.8,
        "ytick.major.size": 2.8,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from matplotlib import patheffects as pe
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "results/figures/submission"
OUT_DIR.mkdir(parents=True, exist_ok=True)

COLORS = {
    "caf": "#8C2D04",
    "ifn": "#3B6FB6",
    "immune": "#2C7A51",
    "tumor": "#B23A48",
    "maturity": "#7B68A6",
    "primary": "#1B9E77",
    "liver": "#D95F02",
    "ln": "#7570B3",
    "normal": "#666666",
    "post": "#4C78A8",
    "naive": "#72B7B2",
    "anchor": "#E45756",
}

TARGET_LABELS = {
    "ifn_mhc": "IFN/MHC",
    "immune_core": "Immune core",
    "tumor_aggressive": "Tumor aggressive",
    "immune_maturity": "Immune maturity-like",
}
TARGET_COLORS = {
    "IFN/MHC": COLORS["ifn"],
    "Immune core": COLORS["immune"],
    "Tumor aggressive": COLORS["tumor"],
    "Immune maturity-like": COLORS["maturity"],
}
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "external_paired_st_anchor": "GSE235315",
}
CONTEXT_COLORS = {
    "post_neoadjuvant_sections": COLORS["post"],
    "treatment_naive_primary": COLORS["naive"],
    "primary_tumor": COLORS["primary"],
    "liver_metastasis": COLORS["liver"],
    "lymph_node_metastasis": COLORS["ln"],
    "normal_pancreas": COLORS["normal"],
    "external_paired_st_anchor": COLORS["anchor"],
}


def clean_axes(ax: plt.Axes) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis="x", color="#E6E6E6", linewidth=0.5, zorder=0)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.04,
        label,
        transform=ax.transAxes,
        fontsize=11,
        fontweight="bold",
        va="bottom",
        ha="left",
    )


def save_figure(fig: plt.Figure, name: str) -> None:
    for ext in ["pdf", "svg", "png"]:
        path = OUT_DIR / f"{name}.{ext}"
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def draw_workflow(ax: plt.Axes) -> None:
    ax.set_axis_off()
    panel_label(ax, "A")
    nodes = [
        ("Public PDAC\nST-H&E cohorts", "GSE282302 108\nGSE274103 5\nGSE272362 30\nGSE235315 7", "#F0F4F8"),
        ("CAF-myeloid\ncore definition", "top 10% within\nsample; 15/5%\nsensitivity", "#FFF4E6"),
        ("Spatial specificity", "1,000 random\nsame-size cores\nper sample", "#EAF6F0"),
        ("Biology readout", "primary/liver\nvalidation; LN\nimmune decoupling", "#F8EEF5"),
        ("Candidate axes", "SPP1-TAM/matrix\nTGF-beta/EMT\nH&E bridge", "#EEF2FF"),
    ]
    xs = np.linspace(0.05, 0.95, len(nodes))
    y = 0.58
    for i, (title, body, color) in enumerate(nodes):
        ax.annotate(
            "",
            xy=(xs[i] - 0.08, y),
            xytext=(xs[i - 1] + 0.08, y),
            arrowprops=dict(arrowstyle="->", color="#555555", lw=1.0),
            annotation_clip=False,
        ) if i > 0 else None
        ax.text(
            xs[i],
            y + 0.12,
            title,
            ha="center",
            va="center",
            fontsize=8.2,
            fontweight="bold",
            bbox=dict(boxstyle="round,pad=0.35,rounding_size=0.08", fc=color, ec="#666666", lw=0.6),
        )
        ax.text(xs[i], y - 0.08, body, ha="center", va="top", fontsize=7.0, color="#333333")
    ax.text(
        0.5,
        0.08,
        "Design principle: spatial core claims require same-sample random-core controls and explicit claim boundaries.",
        ha="center",
        va="center",
        fontsize=7.3,
        color="#333333",
    )


def plot_random_core_bar(ax: plt.Axes, df: pd.DataFrame, title: str, group_order: list[str]) -> None:
    sub = df[df["panel"].eq("B")].copy()
    sub["target"] = sub["metric"].str.replace(" delta_vs_null", "", regex=False).map(TARGET_LABELS)
    sub = sub[sub["group"].isin(group_order)]
    sub["group"] = pd.Categorical(sub["group"], categories=group_order, ordered=True)
    targets = ["IFN/MHC", "Immune core", "Tumor aggressive", "Immune maturity-like"]
    y_positions = []
    y_labels = []
    colors = []
    values = []
    supports = []
    y = 0
    for group in group_order:
        g = sub[sub["group"].eq(group)]
        for target in targets:
            row = g[g["target"].eq(target)]
            if row.empty:
                continue
            y_positions.append(y)
            y_labels.append(f"{group}\n{target}" if target == targets[0] else f"  {target}")
            values.append(float(row["value"].iloc[0]))
            supports.append(row["support"].iloc[0])
            colors.append(TARGET_COLORS[target])
            y += 1
        y += 0.6
    ax.barh(y_positions, values, color=colors, height=0.68, zorder=3)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_xlim(-0.39, 0.08)
    ax.set_yticks(y_positions, y_labels)
    ax.invert_yaxis()
    ax.set_xlabel("observed rho - random-core median", fontsize=8)
    ax.set_title(title, fontsize=9, fontweight="bold", loc="left")
    clean_axes(ax)
    for yp, sup in zip(y_positions, supports):
        ax.text(0.012, yp, str(sup), va="center", ha="left", fontsize=6.8)


def plot_threshold_lines(ax: plt.Axes, df: pd.DataFrame) -> None:
    sub = df[df["panel"].eq("C")].copy()
    sub["target"] = sub["metric"].str.replace(" median_rho", "", regex=False).map(TARGET_LABELS)
    # threshold is stored in core_label in current source table.
    if "core_label" in sub.columns:
        sub["threshold"] = sub["core_label"].str.extract(r"top\\s*(\\d+)%", expand=False)
    else:
        sub["threshold"] = np.nan
    if sub["threshold"].isna().all():
        # Fallback: rows are ordered 15/10/5 within each group/target.
        sub["threshold"] = sub.groupby(["group", "target"]).cumcount().map({0: "15", 1: "10", 2: "5"})
    sub["threshold"] = sub["threshold"].astype(float)
    sub = sub[sub["target"].isin(["IFN/MHC", "Immune core", "Tumor aggressive"])]
    for target, g in sub.groupby("target"):
        summary = g.groupby("threshold", as_index=False)["value"].median().sort_values("threshold", ascending=False)
        ax.plot(summary["threshold"], summary["value"].astype(float), marker="o", lw=1.8, ms=4, color=TARGET_COLORS[target], label=target)
    ax.axhline(0, color="#333333", lw=0.7)
    ax.set_xticks([15, 10, 5])
    ax.set_xlim(16, 4)
    ax.set_xlabel("CAF-core percentile", fontsize=8)
    ax.set_ylabel("median distance-to-core rho", fontsize=8)
    ax.set_title("Threshold-stable gradients", fontsize=9, fontweight="bold", loc="left")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=7, loc="lower left")


def add_image_panel(ax: plt.Axes, path: Path, title: str) -> None:
    ax.set_axis_off()
    img = Image.open(path).convert("RGB")
    arr = np.asarray(img)
    mask = np.any(arr < 248, axis=2)
    if mask.any():
        ys, xs = np.where(mask)
        pad = 8
        left = max(0, xs.min() - pad)
        right = min(img.width, xs.max() + pad)
        top = max(0, ys.min() - pad)
        bottom = min(img.height, ys.max() + pad)
        img = img.crop((left, top, right, bottom))
    ax.imshow(img)
    ax.set_title(title, fontsize=8.5, fontweight="bold", loc="left")


def read_hires_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0))


def make_submission_spatial_examples() -> None:
    spots = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv")
    manifest = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_overlay_manifest.csv")
    rows = [
        ("IU_PDA_T1", "Primary tumor"),
        ("IU_PDA_HM10", "Liver metastasis"),
        ("IU_PDA_LNM7", "Lymph node metastasis"),
    ]
    score_cols = [
        ("score_caf_myeloid_barrier", "CAF-myeloid", COLORS["caf"]),
        ("z_ifn_antigen_presentation", "IFN/MHC", COLORS["ifn"]),
        ("score_tumor_aggressive", "Tumor aggressive", COLORS["tumor"]),
        ("score_immune_hub_core", "Immune core", COLORS["immune"]),
    ]
    fig, axes = plt.subplots(len(rows), 1 + len(score_cols), figsize=(13.8, 8.9), constrained_layout=False)
    fig.subplots_adjust(left=0.055, right=0.995, top=0.90, bottom=0.075, wspace=0.035, hspace=0.12)
    for r, (sample_id, row_label) in enumerate(rows):
        meta = manifest[manifest["sample_id"].eq(sample_id)].iloc[0]
        sample = spots[spots["sample_id"].eq(sample_id)].copy()
        image = Image.open(meta["image_path"]).convert("RGB")
        scale_path = Path(str(meta["image_path"]).replace("_tissue_hires_image.png", "_scalefactors_json.json"))
        scale = read_hires_scale(scale_path)
        x = sample["x_pixel"].to_numpy(float) * scale
        y = sample["y_pixel"].to_numpy(float) * scale
        caf = sample["score_caf_myeloid_barrier"].to_numpy(float)
        core = caf >= np.nanpercentile(caf, 90)
        pad = 0.08
        xmin, xmax = np.nanquantile(x, [pad, 1 - pad])
        ymin, ymax = np.nanquantile(y, [pad, 1 - pad])
        dx = xmax - xmin
        dy = ymax - ymin
        xmin = max(0, xmin - 0.15 * dx)
        xmax = min(image.width, xmax + 0.15 * dx)
        ymin = max(0, ymin - 0.15 * dy)
        ymax = min(image.height, ymax + 0.15 * dy)

        ax = axes[r, 0]
        ax.imshow(image)
        ax.scatter(x, y, s=4, c="#222222", alpha=0.14, linewidths=0)
        ax.scatter(x[core], y[core], s=16, facecolors="none", edgecolors="#D62728", linewidths=0.55, alpha=0.95)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("H&E + CAF core" if r == 0 else "", fontsize=8.3, fontweight="bold")
        ax.set_ylabel(row_label, fontsize=8.4, fontweight="bold")

        for c, (col, title, color) in enumerate(score_cols, start=1):
            ax = axes[r, c]
            values = sample[col].to_numpy(float)
            lo, hi = np.nanpercentile(values, [2, 98])
            ax.imshow(image)
            ax.scatter(x, y, c=values, s=9, cmap="viridis", vmin=lo, vmax=hi, alpha=0.88, linewidths=0)
            ax.scatter(x[core], y[core], s=17, facecolors="none", edgecolors="#D62728", linewidths=0.45, alpha=0.85)
            ax.set_xlim(xmin, xmax)
            ax.set_ylim(ymax, ymin)
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(title if r == 0 else "", fontsize=8.3, fontweight="bold")
    fig.suptitle("Representative GSE272362 spatial programs around CAF-myeloid cores", fontsize=13, fontweight="bold")
    fig.text(
        0.5,
        0.028,
        "Colored spots are scaled within each sample and program by the 2nd-98th percentiles; red open circles mark top 10% CAF-myeloid core spots.",
        ha="center",
        va="bottom",
        fontsize=7.5,
        color="#333333",
    )
    save_figure(fig, "figure2_supplement_submission_spatial_examples")


def read_sample_spots(path: Path, sample_id: str, usecols: list[str]) -> pd.DataFrame:
    chunks: list[pd.DataFrame] = []
    for chunk in pd.read_csv(path, usecols=usecols, chunksize=250_000):
        hit = chunk[chunk["sample_id"].eq(sample_id)]
        if not hit.empty:
            chunks.append(hit)
    if not chunks:
        raise ValueError(f"No spots found for {sample_id} in {path}")
    return pd.concat(chunks, ignore_index=True)


def make_submission_post_nact_spatial_example() -> None:
    sample_id = "GSM8641105_C3_D8_ROI3"
    score_cols = [
        ("score_caf_myeloid_barrier", "CAF-myeloid", COLORS["caf"]),
        ("z_ifn_antigen_presentation", "IFN/MHC", COLORS["ifn"]),
        ("score_tumor_aggressive", "Tumor aggressive", COLORS["tumor"]),
        ("score_immune_hub_core", "Immune core", COLORS["immune"]),
    ]
    usecols = ["dataset_id", "sample_id", "x_pixel", "y_pixel"] + [col for col, _, _ in score_cols]
    spot_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv"
    if not spot_path.exists():
        spot_path = PROJECT_ROOT / "results/tables/mvp_spot_level_scores.csv"
    sample = read_sample_spots(spot_path, sample_id, usecols)
    manifest = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_overlay_manifest.csv")
    meta = manifest[manifest["sample_id"].eq(sample_id)].iloc[0]
    image = Image.open(meta["image_path"]).convert("RGB")
    scale_path = Path(str(meta["image_path"]).replace("_tissue_hires_image.png", "_scalefactors_json.json"))
    scale = read_hires_scale(scale_path)
    x = sample["x_pixel"].to_numpy(float) * scale
    y = sample["y_pixel"].to_numpy(float) * scale
    caf = sample["score_caf_myeloid_barrier"].to_numpy(float)
    core = caf >= np.nanpercentile(caf, 90)
    xmin, xmax = np.nanquantile(x, [0.04, 0.96])
    ymin, ymax = np.nanquantile(y, [0.04, 0.96])
    dx = xmax - xmin
    dy = ymax - ymin
    xmin = max(0, xmin - 0.12 * dx)
    xmax = min(image.width, xmax + 0.12 * dx)
    ymin = max(0, ymin - 0.12 * dy)
    ymax = min(image.height, ymax + 0.12 * dy)

    fig, axes = plt.subplots(1, 1 + len(score_cols), figsize=(13.8, 3.2), constrained_layout=False)
    fig.subplots_adjust(left=0.02, right=0.995, top=0.82, bottom=0.13, wspace=0.035)
    ax = axes[0]
    ax.imshow(image)
    ax.scatter(x, y, s=3.2, c="#222222", alpha=0.13, linewidths=0)
    ax.scatter(x[core], y[core], s=15, facecolors="none", edgecolors="#D62728", linewidths=0.55, alpha=0.92)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("H&E + CAF core", fontsize=8.6, fontweight="bold")
    for ax, (col, title, color) in zip(axes[1:], score_cols):
        values = sample[col].to_numpy(float)
        lo, hi = np.nanpercentile(values, [2, 98])
        ax.imshow(image)
        ax.scatter(x, y, c=values, s=8, cmap="viridis", vmin=lo, vmax=hi, alpha=0.88, linewidths=0)
        ax.scatter(x[core], y[core], s=15, facecolors="none", edgecolors="#D62728", linewidths=0.45, alpha=0.82)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymax, ymin)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(title, fontsize=8.6, fontweight="bold")
    fig.suptitle("Representative post-NACT PDAC section (GSE282302)", fontsize=12, fontweight="bold")
    fig.text(
        0.5,
        0.035,
        "Colored spots are scaled within this sample and program by the 2nd-98th percentiles; red open circles mark top 10% CAF-myeloid core spots.",
        ha="center",
        va="bottom",
        fontsize=7.2,
        color="#333333",
    )
    save_figure(fig, "figure1_supplement_submission_post_nact_spatial_example")


def make_submission_figure1() -> None:
    fig1 = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_1.csv")
    gse235315 = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_6A.csv")
    gse235315 = gse235315.rename(columns={"target": "metric"})
    if "metric" in gse235315.columns:
        gse235315["metric"] = gse235315["metric"].str.replace("_", " ") + " delta_vs_null"
    extra_rows = []
    target_map = {
        "ifn_mhc": "ifn_mhc delta_vs_null",
        "immune_core": "immune_core delta_vs_null",
        "tumor_aggressive": "tumor_aggressive delta_vs_null",
        "immune_maturity": "immune_maturity delta_vs_null",
    }
    raw6 = pd.read_csv(PROJECT_ROOT / "results/tables/gse235315_random_core_anchor_summary.csv")
    label_to_metric = {
        "IFN/MHC": "ifn_mhc delta_vs_null",
        "immune core": "immune_core delta_vs_null",
        "tumor aggressive": "tumor_aggressive delta_vs_null",
        "immune maturity-like": "immune_maturity delta_vs_null",
    }
    raw6_summary = (
        raw6.groupby("target_label", dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            value=("delta_vs_null_median", "median"),
            n_support=("observed_more_negative_than_null", lambda s: int(s.astype(str).str.lower().eq("true").sum())),
        )
        .reset_index()
    )
    for _, row in raw6_summary.iterrows():
        extra_rows.append(
            {
                "figure": "Figure 1",
                "panel": "B",
                "source_file": "gse235315_random_core_anchor_summary.csv",
                "metric": label_to_metric.get(row["target_label"], f"{row['target_label']} delta_vs_null"),
                "group": "GSE235315",
                "n_samples": row["n_samples"],
                "value": row["value"],
                "support": f"{row['n_support']}/{row['n_samples']}",
                "core_label": "",
            }
        )
    fig1_aug = pd.concat([fig1, pd.DataFrame(extra_rows)], ignore_index=True)

    fig = plt.figure(figsize=(13.2, 8.2))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[0.74, 1.08, 1.08], width_ratios=[1.15, 1.05, 1.45], hspace=0.58, wspace=0.46)
    axA = fig.add_subplot(gs[0, :])
    draw_workflow(axA)

    axB = fig.add_subplot(gs[1:, 0])
    panel_label(axB, "B")
    plot_random_core_bar(axB, fig1_aug, "CAF-core specificity across cohorts", ["GSE282302", "GSE274103", "GSE235315"])

    axC = fig.add_subplot(gs[1, 1])
    panel_label(axC, "C")
    plot_threshold_lines(axC, fig1)

    axD = fig.add_subplot(gs[2, 1:])
    panel_label(axD, "D")
    add_image_panel(
        axD,
        PROJECT_ROOT / "results/figures/submission/figure1_supplement_submission_post_nact_spatial_example.png",
        "Representative post-NACT section",
    )

    axE = fig.add_subplot(gs[1, 2])
    axE.set_axis_off()
    panel_label(axE, "E")
    axE.text(0.02, 0.88, "Evidence ladder", fontsize=10, fontweight="bold", transform=axE.transAxes)
    ladder = [
        ("1", "Same-size\nrandom cores", "1,000/sample", COLORS["caf"]),
        ("2", "Core-threshold\nsensitivity", "15% | 10% | 5%", COLORS["ifn"]),
        ("3", "External\npaired-ST anchor", "GSE235315", COLORS["immune"]),
    ]
    y_positions = [0.68, 0.43, 0.18]
    for (num, label, detail, color), y in zip(ladder, y_positions):
        axE.add_patch(
            plt.Circle((0.12, y + 0.03), 0.055, transform=axE.transAxes, facecolor=color, edgecolor="#222222", lw=0.6)
        )
        axE.text(0.12, y + 0.03, num, color="white", fontsize=8.5, fontweight="bold", ha="center", va="center", transform=axE.transAxes)
        axE.add_patch(
            plt.Rectangle((0.22, y - 0.04), 0.70, 0.14, transform=axE.transAxes, facecolor="#F7F7F7", edgecolor="#BBBBBB", lw=0.6)
        )
        axE.text(0.25, y + 0.045, label, fontsize=7.5, fontweight="bold", va="center", transform=axE.transAxes)
        axE.text(0.63, y + 0.045, detail, fontsize=7.2, color="#333333", va="center", transform=axE.transAxes)
    for y0, y1 in zip(y_positions[:-1], y_positions[1:]):
        axE.annotate(
            "",
            xy=(0.12, y1 + 0.095),
            xytext=(0.12, y0 - 0.025),
            xycoords=axE.transAxes,
            arrowprops=dict(arrowstyle="-|>", lw=0.8, color="#666666"),
        )
    fig.suptitle("CAF-myeloid cores define reproducible spatial organizing regions in PDAC", fontsize=13, fontweight="bold", y=0.985)
    save_figure(fig, "figure1_submission_spatial_specificity")


def make_submission_figure2() -> None:
    fig2 = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_2.csv")
    fig3a = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_3A.csv")
    dec = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_5C.csv")

    fig = plt.figure(figsize=(13.2, 6.4))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[0.85, 1.2], width_ratios=[0.95, 1.15, 1.35], hspace=0.55, wspace=0.42)

    axA = fig.add_subplot(gs[0, 0])
    panel_label(axA, "A")
    counts = fig2[(fig2["panel"].eq("A")) & (fig2["metric"].eq("site sample and spot counts"))].copy()
    order = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis", "normal_pancreas"]
    counts["group"] = pd.Categorical(counts["group"], categories=order, ordered=True)
    counts = counts.sort_values("group")
    axA.bar(np.arange(len(counts)), counts["n_samples"].astype(int), color=[CONTEXT_COLORS.get(g, "#999999") for g in counts["group"]])
    axA.set_xticks(np.arange(len(counts)), [CONTEXT_LABELS.get(g, g) for g in counts["group"]], rotation=35, ha="right")
    axA.set_ylabel("specimens", fontsize=8)
    axA.set_title("GSE272362 validation atlas", fontsize=9, fontweight="bold", loc="left")
    clean_axes(axA)
    axA.set_ylim(0, max(counts["n_samples"].astype(int)) + 3.0)
    for i, row in enumerate(counts.itertuples()):
        axA.text(i, int(row.n_samples) + 0.45, f"{int(row.value)/1000:.1f}k\nspots", ha="center", va="bottom", fontsize=6.4)

    axB = fig.add_subplot(gs[0, 1:])
    panel_label(axB, "B")
    sub = fig2[fig2["panel"].eq("B")].copy()
    sub["target"] = sub["metric"].str.replace(" delta_vs_null", "", regex=False).map(TARGET_LABELS)
    sub = sub[sub["target"].isin(["IFN/MHC", "Immune core", "Tumor aggressive"])]
    x_order = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
    offsets = {"IFN/MHC": -0.22, "Immune core": 0.0, "Tumor aggressive": 0.22}
    for target, group in sub.groupby("target"):
        xs = [x_order.index(g) + offsets[target] for g in group["group"] if g in x_order]
        vals = [float(v) for g, v in zip(group["group"], group["value"]) if g in x_order]
        bars = axB.bar(xs, vals, width=0.20, color=TARGET_COLORS[target], label=target, zorder=3)
        for bar, sup in zip(bars, [s for g, s in zip(group["group"], group["support"]) if g in x_order]):
            val = bar.get_height()
            axB.text(
                bar.get_x() + bar.get_width() / 2,
                val - 0.025 if val < 0 else val + 0.025,
                sup,
                ha="center",
                va="top" if val < 0 else "bottom",
                fontsize=6.5,
            )
    axB.axhline(0, color="#333333", lw=0.8)
    axB.set_xticks(range(len(x_order)), [CONTEXT_LABELS[g] for g in x_order])
    axB.set_ylabel("observed rho - random-core median", fontsize=8)
    axB.set_title("Primary/liver validation; lymph-node immune divergence", fontsize=9, fontweight="bold", loc="left")
    clean_axes(axB)
    axB.legend(frameon=False, fontsize=7, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.28))

    axC = fig.add_subplot(gs[1, :2])
    panel_label(axC, "C")
    programs = ["myCAF", "Myeloid", "SPP1/TREM2 TAM", "TGF-beta", "EMT/invasion", "Tumor-aggressive", "IFN/MHC", "Immune core", "T cell", "B cell", "DC/APC"]
    sites = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
    pivot = fig3a.pivot_table(index="program_label", columns="specimen_type", values="median_rho", aggfunc="median").reindex(index=programs, columns=sites)
    im = axC.imshow(pivot.to_numpy(float), cmap="RdBu_r", vmin=-0.55, vmax=0.55, aspect="auto")
    axC.set_xticks(np.arange(len(sites)), [CONTEXT_LABELS[s] for s in sites], rotation=20, ha="right")
    axC.set_yticks(np.arange(len(programs)), programs)
    axC.set_title("CAF-core subprogram decomposition", fontsize=9, fontweight="bold", loc="left")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            if np.isfinite(val):
                axC.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6.2)
    cb = fig.colorbar(im, ax=axC, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("rho to CAF-core distance", fontsize=7)

    axD = fig.add_subplot(gs[1, 2])
    panel_label(axD, "D")
    order_dec = ["treatment_naive_primary", "primary_tumor", "liver_metastasis", "post_neoadjuvant_sections", "lymph_node_metastasis"]
    dec = dec[dec["cohort_context"].isin(order_dec)].copy()
    dec["cohort_context"] = pd.Categorical(dec["cohort_context"], categories=order_dec, ordered=True)
    dec = dec.sort_values("cohort_context")
    vals = dec["median_immune_decoupling_index"].astype(float)
    axD.barh(np.arange(len(dec)), vals, color=[CONTEXT_COLORS[str(c)] for c in dec["cohort_context"]])
    axD.axvline(0, color="#333333", lw=0.8)
    axD.set_yticks(np.arange(len(dec)), [CONTEXT_LABELS[str(c)] for c in dec["cohort_context"]])
    axD.set_xlabel("immune-decoupling index", fontsize=8)
    axD.set_title("LN metastases decouple immune programs", fontsize=9, fontweight="bold", loc="left")
    clean_axes(axD)

    fig.suptitle("Metastatic site remodels CAF-core immune coupling", fontsize=13, fontweight="bold", y=0.985)
    save_figure(fig, "figure2_submission_metastatic_decoupling")


def make_submission_figure3() -> None:
    nmf = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_5A.csv")
    counts = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_5B.csv")
    axes = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_7A_B.csv")
    corr = pd.read_csv(PROJECT_ROOT / "results/source_data/Source_Data_Fig_7C.csv")

    fig = plt.figure(figsize=(13.4, 9.5))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.05], width_ratios=[1.1, 1.0], hspace=0.42, wspace=0.42)

    axA = fig.add_subplot(gs[0, 0])
    panel_label(axA, "A")
    programs = ["myCAF", "iCAF", "apCAF", "myeloid", "SPP1/TREM2 TAM", "TGF-beta", "EMT/invasion", "hypoxia", "basal-like", "tumor aggressive", "IFN/MHC", "immune core", "T cell", "B cell", "DC/APC", "plasma cell"]
    mat = nmf.set_index("nmf_ecotype")[programs]
    mat = mat.div(mat.max(axis=1), axis=0)
    im = axA.imshow(mat.to_numpy(float), cmap="viridis", aspect="auto", vmin=0, vmax=1)
    axA.set_xticks(np.arange(len(programs)), programs, rotation=45, ha="right", fontsize=7)
    axA.set_yticks(np.arange(len(mat)), mat.index)
    axA.set_title("CAF-core ecotype loadings", fontsize=9, fontweight="bold", loc="left")
    cb = fig.colorbar(im, ax=axA, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=7)

    axB = fig.add_subplot(gs[0, 1])
    panel_label(axB, "B")
    context_order = ["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
    ecotypes = sorted(counts["dominant_nmf_ecotype"].unique())
    pivot = counts.pivot_table(index="cohort_context", columns="dominant_nmf_ecotype", values="n_samples", aggfunc="sum").fillna(0).reindex(context_order)
    frac = pivot.div(pivot.sum(axis=1), axis=0).fillna(0)
    bottom = np.zeros(len(frac))
    pal = ["#4C78A8", "#F58518", "#54A24B", "#B279A2"]
    for i, ec in enumerate(ecotypes):
        if ec not in frac.columns:
            continue
        axB.bar(np.arange(len(frac)), frac[ec], bottom=bottom, color=pal[i % len(pal)], label=ec, width=0.72)
        bottom += frac[ec].to_numpy(float)
    axB.set_xticks(np.arange(len(frac)), [CONTEXT_LABELS.get(c, c) for c in frac.index], rotation=30, ha="right")
    axB.set_ylim(0, 1)
    axB.set_ylabel("fraction of samples", fontsize=8)
    axB.set_title("Dominant CAF-core ecotypes by context", fontsize=9, fontweight="bold", loc="left")
    clean_axes(axB)
    axB.legend(frameon=False, fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.25))

    axC = fig.add_subplot(gs[1, 0])
    panel_label(axC, "C")
    context_order2 = ["post_neoadjuvant_sections", "treatment_naive_primary", "primary_tumor", "liver_metastasis", "lymph_node_metastasis", "external_paired_st_anchor"]
    axis_order = ["SPP1-TAM/matrix", "TGF-beta/EMT invasive", "IFN/APC antigen", "B/plasma lymphoid", "T cell/checkpoint", "Basal-hypoxic tumor"]
    core = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_core_enrichment", aggfunc="median").reindex(index=axis_order, columns=context_order2)
    interface = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_interface_enrichment", aggfunc="median").reindex(index=axis_order, columns=context_order2)
    combined = np.concatenate([core.to_numpy(float), interface.to_numpy(float)], axis=1)
    col_labels = [CONTEXT_LABELS[c] for c in context_order2] + [CONTEXT_LABELS[c] for c in context_order2]
    im2 = axC.imshow(combined, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    axC.axvline(len(context_order2) - 0.5, color="#333333", lw=1.0)
    axC.set_xticks(np.arange(len(col_labels)), col_labels, rotation=45, ha="right", fontsize=6.6)
    axC.set_yticks(np.arange(len(axis_order)), axis_order)
    axC.text(0.22, 1.04, "CAF-core enrichment", transform=axC.transAxes, ha="center", fontsize=8, fontweight="bold")
    axC.text(0.75, 1.04, "interface enrichment", transform=axC.transAxes, ha="center", fontsize=8, fontweight="bold")
    axC.set_title("Candidate-axis enrichment", fontsize=9, fontweight="bold", loc="left", pad=24)
    for i in range(combined.shape[0]):
        for j in range(combined.shape[1]):
            val = combined[i, j]
            if np.isfinite(val):
                axC.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8)
    cb2 = fig.colorbar(im2, ax=axC, fraction=0.035, pad=0.02)
    cb2.ax.tick_params(labelsize=7)

    axD = fig.add_subplot(gs[1, 1])
    panel_label(axD, "D")
    c = corr[corr["metric"].eq("core_enrichment")].set_index("axis_label").reindex(axis_order)
    vals = c["rho_with_immune_decoupling_index"].astype(float)
    y = np.arange(len(c))
    axD.barh(y, vals, color=["#4C78A8" if v > 0 else "#B279A2" for v in vals])
    axD.axvline(0, color="#333333", lw=0.8)
    axD.set_yticks(y, c.index)
    axD.set_xlabel("rho with immune-decoupling index", fontsize=8)
    axD.set_title("Immune decoupling separates invasive and immune axes", fontsize=9, fontweight="bold", loc="left")
    clean_axes(axD)
    axD.set_xlim(-0.75, 0.55)
    for yy, val in zip(y, vals):
        if abs(val) > 0.18:
            axD.text(val - 0.02 if val > 0 else val + 0.02, yy, f"{val:.2f}", va="center", ha="right" if val > 0 else "left", fontsize=7, color="white")
        else:
            axD.text(val + (0.02 if val >= 0 else -0.02), yy, f"{val:.2f}", va="center", ha="left" if val >= 0 else "right", fontsize=7)

    fig.suptitle("CAF-core ecotypes resolve candidate invasive-interface axes", fontsize=13, fontweight="bold", y=0.985)
    save_figure(fig, "figure3_submission_ecotypes_mechanism_axes")


def main() -> int:
    make_submission_post_nact_spatial_example()
    make_submission_figure1()
    make_submission_figure2()
    make_submission_figure3()
    make_submission_spatial_examples()
    print(f"Wrote submission figure suite to {OUT_DIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
