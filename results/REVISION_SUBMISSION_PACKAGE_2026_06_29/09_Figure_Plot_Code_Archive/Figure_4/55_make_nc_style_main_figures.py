from __future__ import annotations

import importlib.util
import json
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
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
REVISION_ANALYSIS = ROOT / "results" / "revision_2026_06_29" / "analysis_outputs"
FIG_DIR.mkdir(parents=True, exist_ok=True)
SOURCE_DIR.mkdir(parents=True, exist_ok=True)

TARGET_LABELS = {
    "ifn_mhc": "IFN/MHC",
    "immune_core": "Immune core",
    "tumor_aggressive": "Tumor aggressive",
    "immune_maturity": "Immune maturity-like",
}
TARGET_COLORS = {
    "IFN/MHC": "#3B6FB6",
    "Immune core": "#2C7A51",
    "Tumor aggressive": "#B23A48",
    "Immune maturity-like": "#7B68A6",
}
CONTEXT_LABELS = {
    "GSE282302": "post-NACT\nGSE282302",
    "GSE274103": "treatment-naive\nGSE274103",
    "GSE235315": "external\nGSE235315",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
}
CONTEXT_COLORS = {
    "primary_tumor": "#1B9E77",
    "liver_metastasis": "#D95F02",
    "lymph_node_metastasis": "#7570B3",
    "normal_pancreas": "#666666",
    "post_neoadjuvant_sections": "#4C78A8",
    "treatment_naive_primary": "#72B7B2",
    "GSE282302": "#4C78A8",
    "GSE274103": "#72B7B2",
    "GSE235315": "#E45756",
}


def load_module(script_name: str, name: str):
    path = ROOT / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


FIG26 = load_module("26_make_submission_figure_suite.py", "fig26")
FIG49 = load_module("49_make_candidate_figure4_multiresolution_validation.py", "fig49")
FIG50 = load_module("50_make_reference_style_supplement_figures.py", "fig50")


def save(fig: plt.Figure, base: str) -> None:
    for ext in ["pdf", "svg", "png"]:
        path = FIG_DIR / f"{base}.{ext}"
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom", ha="left")


def clean_axes(ax: plt.Axes, grid_axis: str | None = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    if grid_axis:
        ax.grid(axis=grid_axis, color="#E6E6E6", linewidth=0.45, zorder=0)
        ax.set_axisbelow(True)


def read_hires_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0))


def spatial_bounds(x: np.ndarray, y: np.ndarray, image: Image.Image, q: float = 0.05) -> tuple[float, float, float, float]:
    xmin, xmax = np.nanquantile(x, [q, 1 - q])
    ymin, ymax = np.nanquantile(y, [q, 1 - q])
    dx, dy = xmax - xmin, ymax - ymin
    return (
        max(0, xmin - 0.12 * dx),
        min(image.width, xmax + 0.12 * dx),
        max(0, ymin - 0.12 * dy),
        min(image.height, ymax + 0.12 * dy),
    )


def plot_spatial_program(
    ax: plt.Axes,
    spots: pd.DataFrame,
    image: Image.Image,
    scale: float,
    column: str | None,
    title: str,
    core_col: str = "score_caf_myeloid_barrier",
    s: float = 5.5,
) -> None:
    x = spots["x_pixel"].to_numpy(float) * scale
    y = spots["y_pixel"].to_numpy(float) * scale
    core = spots[core_col].to_numpy(float) >= np.nanpercentile(spots[core_col].to_numpy(float), 90)
    xmin, xmax, ymin, ymax = spatial_bounds(x, y, image)
    ax.imshow(image)
    if column is None:
        ax.scatter(x, y, s=s * 0.55, c="#222222", alpha=0.12, linewidths=0)
    else:
        vals = spots[column].to_numpy(float)
        lo, hi = np.nanpercentile(vals, [2, 98])
        ax.scatter(x, y, c=vals, s=s, cmap="viridis", vmin=lo, vmax=hi, alpha=0.88, linewidths=0)
    ax.scatter(x[core], y[core], s=s * 2.0, facecolors="none", edgecolors="#D73027", linewidths=0.42, alpha=0.85)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, fontsize=7.2, fontweight="bold", pad=2)


def plot_cohort_scale(ax: plt.Axes) -> None:
    labels = ["GSE282302", "GSE274103", "GSE272362", "GSE235315", "GSE274557", "GSE274673"]
    values = [108, 5, 30, 7, 55, 6]
    colors = ["#4C78A8", "#72B7B2", "#7570B3", "#E45756", "#59A14F", "#9C755F"]
    ax.bar(np.arange(len(values)), values, color=colors, edgecolor="#333333", linewidth=0.4)
    ax.set_xticks(np.arange(len(values)), ["GSE282302", "GSE274103", "GSE272362", "GSE235315", "GSE274557", "GSE274673"], rotation=35, ha="right")
    ax.set_ylabel("sections / samples", fontsize=7.5)
    ax.set_title("Cohort scale and evidence roles", loc="left", fontsize=9.2, fontweight="bold")
    for i, v in enumerate(values):
        ax.text(i, v + 2.0, str(v), ha="center", va="bottom", fontsize=6.8)
    clean_axes(ax, "y")


def plot_random_specificity(ax: plt.Axes) -> None:
    fig1 = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_1.csv")
    raw6 = pd.read_csv(TABLE_DIR / "gse235315_random_core_anchor_summary.csv")
    label_to_metric = {
        "IFN/MHC": "ifn_mhc delta_vs_null",
        "immune core": "immune_core delta_vs_null",
        "tumor aggressive": "tumor_aggressive delta_vs_null",
        "immune maturity-like": "immune_maturity delta_vs_null",
    }
    rows = []
    for _, row in raw6.groupby("target_label", dropna=False).agg(
        n_samples=("sample_id", "nunique"),
        value=("delta_vs_null_median", "median"),
        n_support=("observed_more_negative_than_null", lambda s: int(s.astype(str).str.lower().eq("true").sum())),
    ).reset_index().iterrows():
        rows.append(
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
    df = pd.concat([fig1, pd.DataFrame(rows)], ignore_index=True)
    FIG26.plot_random_core_bar(ax, df, "Same-size random-core specificity", ["GSE282302", "GSE274103", "GSE235315"])


def plot_null_diagnostics(ax: plt.Axes) -> None:
    nulls = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_19_random_core_null_diagnostics.csv")
    order = [
        ("GSE282302", "metadata_required", "ifn_mhc", "post-NACT\nIFN/MHC"),
        ("GSE272362", "liver_metastasis", "immune_core", "liver met\nimmune core"),
        ("GSE272362", "lymph_node_metastasis", "tumor_aggressive", "LN met\ntumor aggressive"),
    ]
    y = np.arange(len(order))
    obs, lo, hi, med, labels = [], [], [], [], []
    for dataset, specimen_type, target, label in order:
        row = nulls[(nulls["dataset_id"].eq(dataset)) & (nulls["specimen_type"].eq(specimen_type)) & (nulls["target"].eq(target))].iloc[0]
        obs.append(float(row["observed_rho"]))
        med.append(float(row["null_median_rho"]))
        lo.append(float(row["null_p05_rho"]))
        hi.append(float(row["null_p95_rho"]))
        labels.append(label)
    ax.hlines(y, lo, hi, color="#AAAAAA", lw=5, alpha=0.8, label="random 5-95%")
    ax.scatter(med, y, s=24, color="#555555", zorder=3, label="random median")
    ax.scatter(obs, y, s=34, color="#D73027", zorder=4, label="observed")
    ax.axvline(0, color="#333333", lw=0.7)
    ax.set_yticks(y, labels)
    ax.invert_yaxis()
    ax.set_xlabel("distance-to-core rho", fontsize=7.5)
    ax.set_title("Representative random-core nulls", loc="left", fontsize=9.2, fontweight="bold")
    clean_axes(ax, "x")
    ax.legend(frameon=False, fontsize=6.3, loc="lower right")


def plot_contiguous_null_summary(ax: plt.Axes) -> None:
    path = REVISION_ANALYSIS / "stronger_null_contiguous_random_core_summary.csv"
    if not path.exists():
        ax.text(0.5, 0.5, "contiguous-null\nsummary missing", ha="center", va="center", fontsize=8)
        ax.axis("off")
        return
    data = pd.read_csv(path)
    contexts = [
        ("GSE282302", "post_neoadjuvant_sections", "post-NACT"),
        ("GSE272362", "primary_tumor", "primary"),
        ("GSE272362", "liver_metastasis", "liver met"),
        ("GSE272362", "lymph_node_metastasis", "LN met"),
    ]
    targets = ["IFN/MHC", "immune-core", "tumor-aggressive"]
    rows = []
    for dataset, tissue, label in contexts:
        for target in targets:
            sub = data[
                data["dataset"].eq(dataset)
                & data["tissue_site"].eq(tissue)
                & data["target_program"].eq(target)
            ]
            if sub.empty:
                continue
            row = sub.iloc[0]
            rows.append(
                {
                    "context": label,
                    "target": target,
                    "delta": float(row["median_delta"]),
                    "support": float(row["support_fraction"]),
                    "support_label": f"{int(row['support_n'])}/{int(row['n_samples'])}",
                }
            )
    plot_df = pd.DataFrame(rows)
    x_labels = [c[2] for c in contexts]
    y_labels = targets
    matrix = np.full((len(y_labels), len(x_labels)), np.nan)
    labels = [["" for _ in x_labels] for _ in y_labels]
    support = np.full_like(matrix, np.nan, dtype=float)
    for _, row in plot_df.iterrows():
        i = y_labels.index(row["target"])
        j = x_labels.index(row["context"])
        matrix[i, j] = row["delta"]
        labels[i][j] = row["support_label"]
        support[i, j] = row["support"]
    im = ax.imshow(matrix, cmap="RdBu_r", vmin=-0.20, vmax=0.20, aspect="auto")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if np.isfinite(matrix[i, j]):
                color = "white" if abs(matrix[i, j]) > 0.12 else "#222222"
                ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=6.9, color=color)
                ax.scatter(j, i + 0.31, s=16 + 85 * support[i, j], facecolors="none", edgecolors=color, linewidths=0.7)
    ax.set_xticks(np.arange(len(x_labels)), x_labels, rotation=25, ha="right")
    ax.set_yticks(np.arange(len(y_labels)), ["IFN/MHC", "immune-core", "tumor-aggressive"])
    ax.tick_params(labelsize=7)
    ax.set_title("Spatially contiguous random-core specificity", loc="left", fontsize=9.2, fontweight="bold")
    ax.text(
        0.0,
        -0.28,
        "labels: support n/N; more negative = stronger CAF-core centering",
        transform=ax.transAxes,
        fontsize=6.3,
        color="#444444",
        ha="left",
        va="top",
    )
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    cbar.ax.tick_params(labelsize=6.5)
    cbar.set_label("observed - contiguous-null rho", fontsize=6.8)


def plot_distance_bins(ax: plt.Axes) -> None:
    bins = pd.read_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_12A_C_distance_dynamics.csv")
    sub = bins[bins["panel_context"].eq("Discovery/support cohorts")].copy()
    order = ["Core", "Near", "Mid", "Far"]
    for target in ["IFN/MHC", "immune core", "Tumor aggressive"]:
        g = sub[sub["target_label"].eq(target)].sort_values("bin_index")
        color = TARGET_COLORS["Immune core" if target == "immune core" else target]
        ax.plot(np.arange(len(g)), g["median"].astype(float), marker="o", lw=1.6, ms=3.8, color=color, label=target)
        ax.fill_between(np.arange(len(g)), g["q25"].astype(float), g["q75"].astype(float), color=color, alpha=0.12, linewidth=0)
    ax.axhline(0, color="#333333", lw=0.7)
    ax.set_xticks(np.arange(len(order)), order)
    ax.set_ylabel("within-sample z score", fontsize=7.5)
    ax.set_title("Program decay from CAF cores", loc="left", fontsize=9.2, fontweight="bold")
    clean_axes(ax, "y")
    ax.legend(frameon=False, fontsize=6.5, loc="upper right")


def make_figure1_nc() -> None:
    fig = plt.figure(figsize=(15.6, 10.6))
    gs = GridSpec(4, 4, figure=fig, height_ratios=[0.95, 1.0, 1.15, 1.15], width_ratios=[1.1, 1.15, 1.15, 1.15], hspace=0.55, wspace=0.42)
    fig.suptitle("CAF-myeloid cores mark reproducible spatial architectures in PDAC", fontsize=15.2, fontweight="bold", y=0.985)

    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0:2, 1]),
        fig.add_subplot(gs[0:2, 2]),
        fig.add_subplot(gs[0:2, 3]),
        fig.add_subplot(gs[1, 0]),
    ]
    letters = list("ABCDE")
    for ax, letter in zip(axes, letters):
        panel_label(ax, letter)
    plot_cohort_scale(axes[0])
    plot_random_specificity(axes[1])
    FIG26.plot_threshold_lines(axes[2], pd.read_csv(SOURCE_DIR / "Source_Data_Fig_1.csv"))
    plot_contiguous_null_summary(axes[3])
    plot_distance_bins(axes[4])

    sub = gs[2:, :].subgridspec(2, 5, hspace=0.10, wspace=0.04)
    sample_id = "GSM8641105_C3_D8_ROI3"
    usecols = ["dataset_id", "sample_id", "x_pixel", "y_pixel", "score_caf_myeloid_barrier", "z_ifn_antigen_presentation", "score_tumor_aggressive", "score_immune_hub_core"]
    spot_path = TABLE_DIR / "mvp_spot_level_scores_with_edge_qc.csv"
    if not spot_path.exists():
        spot_path = TABLE_DIR / "mvp_spot_level_scores.csv"
    sample = FIG26.read_sample_spots(spot_path, sample_id, usecols)
    manifest = pd.read_csv(TABLE_DIR / "mvp_overlay_manifest.csv")
    meta = manifest[manifest["sample_id"].eq(sample_id)].iloc[0]
    image = Image.open(meta["image_path"]).convert("RGB")
    scale = read_hires_scale(Path(str(meta["image_path"]).replace("_tissue_hires_image.png", "_scalefactors_json.json")))
    map_specs = [
        (None, "H&E + CAF core"),
        ("score_caf_myeloid_barrier", "CAF-myeloid"),
        ("z_ifn_antigen_presentation", "IFN/MHC"),
        ("score_tumor_aggressive", "tumor aggressive"),
        ("score_immune_hub_core", "immune core"),
    ]
    for i, (col, title) in enumerate(map_specs):
        ax = fig.add_subplot(sub[0:, i])
        panel_label(ax, chr(ord("F") + i), x=-0.04, y=1.02)
        plot_spatial_program(ax, sample, image, scale, col, title, s=5.0)
    save(fig, "figure1_submission_spatial_specificity_nc_style")


def make_figure2_nc() -> None:
    fig2 = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_2.csv")
    fig3a = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_3A.csv")
    dec = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_5C.csv")
    fig = plt.figure(figsize=(15.8, 13.2))
    gs = GridSpec(4, 4, figure=fig, height_ratios=[0.86, 0.98, 1.12, 1.12], width_ratios=[0.95, 1.18, 1.18, 1.18], hspace=0.68, wspace=0.40)
    fig.suptitle("Metastatic site remodels CAF-core immune coupling", fontsize=15.2, fontweight="bold", y=0.986)

    axA = fig.add_subplot(gs[0, 0])
    panel_label(axA, "A")
    counts = fig2[(fig2["panel"].eq("A")) & (fig2["metric"].eq("site sample and spot counts"))].copy()
    order = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis", "normal_pancreas"]
    counts["group"] = pd.Categorical(counts["group"], categories=order, ordered=True)
    counts = counts.sort_values("group")
    axA.bar(np.arange(len(counts)), counts["n_samples"].astype(int), color=[CONTEXT_COLORS[str(g)] for g in counts["group"]], edgecolor="#333333", linewidth=0.45)
    axA.set_xticks(np.arange(len(counts)), [CONTEXT_LABELS[str(g)] for g in counts["group"]], rotation=35, ha="right")
    axA.set_ylabel("specimens", fontsize=7.5)
    axA.set_title("GSE272362 validation atlas", fontsize=9.2, fontweight="bold", loc="left")
    for i, row in enumerate(counts.itertuples()):
        axA.text(i, int(row.n_samples) + 0.45, f"{int(row.value)/1000:.1f}k\nspots", ha="center", va="bottom", fontsize=6.5)
    clean_axes(axA, "y")

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
            axB.text(bar.get_x() + bar.get_width() / 2, val - 0.025 if val < 0 else val + 0.025, sup, ha="center", va="top" if val < 0 else "bottom", fontsize=6.5)
    axB.axhline(0, color="#333333", lw=0.8)
    axB.set_xticks(range(len(x_order)), [CONTEXT_LABELS[g] for g in x_order])
    axB.set_ylabel("observed rho - random-core median", fontsize=7.5)
    axB.set_title("Primary/liver validation; lymph-node immune divergence", fontsize=9.2, fontweight="bold", loc="left")
    clean_axes(axB, "y")
    axB.legend(frameon=False, fontsize=7, ncol=3, loc="lower center", bbox_to_anchor=(0.5, -0.28))

    axC = fig.add_subplot(gs[1, 0:2])
    panel_label(axC, "C")
    programs = ["myCAF", "Myeloid", "SPP1/TREM2 TAM", "TGF-beta", "EMT/invasion", "Tumor-aggressive", "IFN/MHC", "Immune core", "T cell", "B cell", "DC/APC"]
    sites = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
    pivot = fig3a.pivot_table(index="program_label", columns="specimen_type", values="median_rho", aggfunc="median").reindex(index=programs, columns=sites)
    im = axC.imshow(pivot.to_numpy(float), cmap="RdBu_r", vmin=-0.55, vmax=0.55, aspect="auto")
    axC.set_xticks(np.arange(len(sites)), [CONTEXT_LABELS[s] for s in sites], rotation=20, ha="right")
    axC.set_yticks(np.arange(len(programs)), programs, fontsize=7)
    axC.set_title("CAF-core subprogram decomposition", fontsize=9.2, fontweight="bold", loc="left")
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            val = pivot.iloc[i, j]
            if np.isfinite(val):
                axC.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8)
    cb = fig.colorbar(im, ax=axC, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=6.5)
    cb.set_label("rho to CAF-core distance", fontsize=6.5)

    axD = fig.add_subplot(gs[1, 2:])
    panel_label(axD, "D")
    order_dec = ["treatment_naive_primary", "primary_tumor", "liver_metastasis", "post_neoadjuvant_sections", "lymph_node_metastasis"]
    d = dec[dec["cohort_context"].isin(order_dec)].copy()
    d["cohort_context"] = pd.Categorical(d["cohort_context"], categories=order_dec, ordered=True)
    d = d.sort_values("cohort_context")
    axD.barh(np.arange(len(d)), d["median_immune_decoupling_index"].astype(float), color=[CONTEXT_COLORS[str(c)] for c in d["cohort_context"]])
    axD.axvline(0, color="#333333", lw=0.8)
    axD.set_yticks(np.arange(len(d)), [CONTEXT_LABELS[str(c)] for c in d["cohort_context"]])
    axD.set_xlabel("immune-decoupling index", fontsize=7.5)
    axD.set_title("LN subset suggests immune decoupling", fontsize=9.2, fontweight="bold", loc="left")
    clean_axes(axD, "x")

    map_grid = gs[2:, :].subgridspec(3, 4, hspace=0.28, wspace=0.06)
    spots = pd.read_csv(TABLE_DIR / "gse272362_rds_spot_level_scores.csv")
    manifest = pd.read_csv(TABLE_DIR / "gse272362_rds_overlay_manifest.csv")
    sample_rows = [("IU_PDA_T1", "primary"), ("IU_PDA_HM10", "liver met"), ("IU_PDA_LNM7", "LN met")]
    map_specs = [
        (None, "H&E + CAF core"),
        ("z_ifn_antigen_presentation", "IFN/MHC"),
        ("score_immune_hub_core", "immune core"),
        ("score_tumor_aggressive", "tumor aggressive"),
    ]
    letter = ord("E")
    for r, (sample_id, row_label) in enumerate(sample_rows):
        meta = manifest[manifest["sample_id"].eq(sample_id)].iloc[0]
        sample = spots[spots["sample_id"].eq(sample_id)].copy()
        image = Image.open(meta["image_path"]).convert("RGB")
        scale = read_hires_scale(Path(str(meta["image_path"]).replace("_tissue_hires_image.png", "_scalefactors_json.json")))
        for c, (column, title) in enumerate(map_specs):
            ax = fig.add_subplot(map_grid[r, c])
            panel_label(ax, chr(letter), x=-0.04, y=1.02)
            letter += 1
            plot_spatial_program(ax, sample, image, scale, column, f"{row_label} {title}", s=4.8)
    save(fig, "figure2_submission_metastatic_decoupling_nc_style")


def make_figure4_nc() -> None:
    visium_context = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_context_summary.csv")
    visium_gradients = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_gradients.csv")
    xenium_context = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_context_summary.csv")
    xenium_gradients = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_gradients.csv")
    xenium_composition = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")
    cells = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_cell_scores.csv")

    fig = plt.figure(figsize=(17.0, 13.0))
    gs = GridSpec(4, 4, figure=fig, height_ratios=[0.88, 1.05, 1.10, 1.10], width_ratios=[0.85, 1.25, 1.22, 1.22], hspace=0.70, wspace=0.72)
    fig.suptitle("Independent multi-resolution validation of CAF-domain organization", fontsize=15.2, fontweight="bold", y=0.986)

    axA = fig.add_subplot(gs[0, 0])
    panel_label(axA, "A", x=-0.16)
    FIG49.draw_cohort_scale(axA, visium_gradients, xenium_composition)
    axA.set_title("External validation scale", loc="left", fontsize=9.2, fontweight="bold")
    axA.set_xticklabels(["Primary", "Liver", "Lung", "Perit.", "Xenium\nnaive", "Xenium\nCRT"], rotation=35, ha="right", fontsize=7)

    axB = fig.add_subplot(gs[0, 1:])
    panel_label(axB, "B", x=-0.07)
    FIG49.draw_visium_heatmap(axB, visium_context)
    axB.set_title("GSE274557 Visium CAF-core validation", loc="left", fontsize=9.2, fontweight="bold")

    axC = fig.add_subplot(gs[1, 0:2])
    panel_label(axC, "C")
    FIG49.draw_xenium_summary(axC, xenium_context)
    axC.set_title("GSE274673 Xenium fixed-anchor summary", loc="left", fontsize=9.2, fontweight="bold")
    axC.text(
        0.01,
        -0.22,
        "more negative delta = stronger centering around CAF-domain anchor",
        transform=axC.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="#444444",
    )

    axD = fig.add_subplot(gs[1, 2:])
    panel_label(axD, "D")
    FIG49.draw_xenium_sample_heatmap(axD, xenium_gradients)
    axD.set_title("CAF-SPP1/TAM anchor across sections", loc="left", fontsize=9.2, fontweight="bold")
    axD.text(
        0.01,
        -0.18,
        "more negative = stronger CAF-domain centering",
        transform=axD.transAxes,
        ha="left",
        va="top",
        fontsize=7,
        color="#444444",
    )

    selected = (
        xenium_composition.sort_values("n_cells", ascending=False)
        .groupby("treatment", as_index=False)
        .head(1)
        .sort_values("treatment", ascending=False)
    )
    selected_ids = selected["geo_accession"].tolist()
    source = cells[cells["geo_accession"].isin(selected_ids)].copy()
    source.to_csv(SOURCE_DIR / "Source_Data_Fig_4_NC_style_selected_xenium_cell_scores.csv", index=False)
    map_grid = gs[2:, :].subgridspec(2, 4, hspace=0.20, wspace=0.10)
    columns = [
        ("anchor_CAF_SPP1TAM", "CAF-SPP1/TAM anchor"),
        ("score_SPP1_TAM", "SPP1/TAM"),
        ("score_IFN_APC", "IFN/APC"),
        ("score_Tumor_epithelial", "tumor epithelial"),
    ]
    letter = ord("E")
    for r, sample_id in enumerate(selected_ids):
        sub = source[source["geo_accession"].eq(sample_id)].copy()
        treatment = sub["treatment"].iloc[0].replace("chemoradiotherapy-treated", "CRT")
        title_prefix = sub["title"].iloc[0].replace("Patient ", "P")
        for c, (column, title) in enumerate(columns):
            ax = fig.add_subplot(map_grid[r, c])
            panel_label(ax, chr(letter), x=-0.09, y=1.08)
            letter += 1
            FIG50.plot_cell_map(ax, sub, column, title if r == 0 else "")
            ax.set_title("", loc="left")
            ax.set_title(title if r == 0 else "", fontsize=8.0, fontweight="bold", loc="center", pad=2)
            if c == 0:
                ax.text(
                    -0.08,
                    0.5,
                    f"{title_prefix}\n{treatment}",
                    transform=ax.transAxes,
                    rotation=90,
                    ha="center",
                    va="center",
                    fontsize=7.2,
                    fontweight="bold",
                )
    save(fig, "figure4_submission_multiresolution_validation_nc_style")


def main() -> int:
    make_figure1_nc()
    make_figure2_nc()
    make_figure4_nc()
    panel_index = pd.DataFrame(
        [
            ("Figure 1", "A-E", "cohort scale, random-core specificity, threshold sensitivity, contiguous-null specificity and distance gradients"),
            ("Figure 1", "F-J", "representative post-NACT H&E and program maps"),
            ("Figure 2", "A-D", "GSE272362 validation, CAF-core subprogram decomposition and immune decoupling"),
            ("Figure 2", "E-P", "primary/liver/LN representative H&E and program maps"),
            ("Figure 4", "A-D", "external Visium and Xenium cohort-level validation"),
            ("Figure 4", "E-L", "representative Xenium cell-level CAF-domain maps"),
        ],
        columns=["figure", "panels", "content"],
    )
    panel_index.to_csv(SOURCE_DIR / "Source_Data_Main_Figures_NC_style_panel_index.csv", index=False)
    print("Wrote NC-style main Figures 1, 2 and 4")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
