from __future__ import annotations

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
        "axes.linewidth": 0.7,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from PIL import Image, ImageOps
from scipy.spatial import cKDTree


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
SOURCE_DIR = PROJECT / "results" / "source_data"
MANIFEST = PROJECT / "metadata" / "dataset_manifest_curated.csv"

ED15_OUT = FIG_DIR / "extended_data_figure15_local_spatial_program_maps"
ED16_OUT = FIG_DIR / "extended_data_figure16_interface_compartment_maps"
ED17_OUT = FIG_DIR / "extended_data_figure17_he_patch_examples"
ED18_OUT = FIG_DIR / "extended_data_figure18_xenium_program_neighborhoods"
ED19_OUT = FIG_DIR / "extended_data_figure19_random_core_null_diagnostics"

PROGRAMS = {
    "CAF-myeloid": "score_caf_myeloid_barrier",
    "IFN/MHC": "z_ifn_antigen_presentation",
    "immune core": "score_immune_hub_core",
    "tumor aggressive": "score_tumor_aggressive",
    "SPP1/TAM": "z_spp1_tam",
    "TGF-beta/EMT": "z_tgfb_pathway",
}
DOMINANT_PROGRAMS = ["IFN/MHC", "immune core", "tumor aggressive", "SPP1/TAM"]
PROGRAM_COLORS = {
    "IFN/MHC": "#3B6FB6",
    "immune core": "#2C7A51",
    "tumor aggressive": "#B23A48",
    "SPP1/TAM": "#9C755F",
    "CAF core": "#D73027",
    "tumor-high": "#F0A202",
    "interface": "#7B68A6",
    "other": "#D8D8D8",
}
XENIUM_PROGRAMS = {
    "CAF/matrix": "score_CAF_matrix",
    "SPP1/TAM": "score_SPP1_TAM",
    "IFN/APC": "score_IFN_APC",
    "T/NK": "score_T_NK",
    "TGF-beta/EMT": "score_TGFb_EMT",
    "tumor epithelial": "score_Tumor_epithelial",
}


def read_scale(path: Path) -> float:
    if not path.exists():
        return 1.0
    return float(json.loads(path.read_text(encoding="utf-8")).get("tissue_hires_scalef", 1.0))


def save_figure(fig: plt.Figure, out: Path) -> None:
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        if ext == "png":
            fig.savefig(out.with_suffix(f".{ext}"), dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(out.with_suffix(f".{ext}"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.08, 1.05, label, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def image_paths_for_sample(dataset_id: str, sample_id: str, geo_base: str | None = None) -> tuple[Path | None, Path | None]:
    manifest = pd.read_csv(MANIFEST)
    if dataset_id == "GSE272362" and geo_base:
        row = manifest[(manifest["dataset_id"].eq("GSE272362")) & (manifest["sample_id"].astype(str).eq(geo_base))]
    else:
        row = manifest[(manifest["dataset_id"].eq(dataset_id)) & (manifest["sample_id"].astype(str).eq(sample_id))]
    if row.empty:
        return None, None
    image = Path(str(row.iloc[0]["image_path"]))
    scale = Path(str(row.iloc[0]["scalefactors_path"]))
    return image, scale


def selected_spatial_samples() -> pd.DataFrame:
    mvp_manifest = pd.read_csv(TABLE_DIR / "mvp_overlay_manifest.csv")
    gse_manifest = pd.read_csv(TABLE_DIR / "gse272362_rds_overlay_manifest.csv")
    rows = []
    if len(mvp_manifest):
        row = mvp_manifest.iloc[0]
        rows.append({"dataset_id": "GSE282302", "sample_id": row["sample_id"], "label": "post-NACT", "geo_base": ""})
    for label, specimen in [("primary", "primary_tumor"), ("liver met", "liver_metastasis"), ("LN met", "lymph_node_metastasis")]:
        sub = gse_manifest[gse_manifest["specimen_type"].eq(specimen)]
        if len(sub):
            row = sub.iloc[0]
            rows.append({"dataset_id": "GSE272362", "sample_id": row["sample_id"], "label": label, "geo_base": row["geo_base"]})
    return pd.DataFrame(rows)


def load_spots(dataset_id: str, sample_id: str) -> pd.DataFrame:
    if dataset_id == "GSE272362":
        df = pd.read_csv(TABLE_DIR / "gse272362_rds_spot_level_scores.csv")
    else:
        df = pd.read_csv(TABLE_DIR / "mvp_spot_level_scores_with_edge_qc.csv")
    return df[df["sample_id"].astype(str).eq(str(sample_id))].copy()


def classify_spot_compartments(spots: pd.DataFrame) -> pd.DataFrame:
    spots = spots.copy()
    x = spots[["x_pixel", "y_pixel"]].to_numpy(float)
    caf = spots["score_caf_myeloid_barrier"].astype(float)
    tumor = spots["score_tumor_aggressive"].astype(float)
    spots["is_caf_core"] = caf >= caf.quantile(0.90)
    spots["is_tumor_high"] = tumor >= tumor.quantile(0.80)
    core_xy = x[spots["is_caf_core"].to_numpy()]
    tumor_xy = x[spots["is_tumor_high"].to_numpy()]
    if len(core_xy) and len(tumor_xy):
        med_step = np.nanmedian(cKDTree(x).query(x, k=2)[0][:, 1])
        radius = max(med_step * 2.5, 1.0)
        d_core = cKDTree(core_xy).query(x, k=1)[0]
        d_tumor = cKDTree(tumor_xy).query(x, k=1)[0]
        near_core = d_core <= radius
        near_tumor = d_tumor <= radius
    else:
        near_core = np.zeros(len(spots), dtype=bool)
        near_tumor = np.zeros(len(spots), dtype=bool)
    comp = np.full(len(spots), "other", dtype=object)
    comp[near_tumor] = "tumor-high"
    comp[near_core] = "CAF core"
    comp[near_core & near_tumor] = "interface"
    spots["spatial_compartment"] = comp
    zcols = {label: col for label, col in PROGRAMS.items() if col in spots.columns and label in DOMINANT_PROGRAMS}
    zmat = np.column_stack([spots[col].astype(float).to_numpy() for col in zcols.values()])
    zmat = (zmat - np.nanmean(zmat, axis=0)) / np.nanstd(zmat, axis=0)
    labels = list(zcols.keys())
    spots["dominant_local_program"] = [labels[i] for i in np.nanargmax(zmat, axis=1)]
    return spots


def plot_he(ax: plt.Axes, image: Image.Image, spots: pd.DataFrame, scale: float, title: str) -> None:
    ax.imshow(image)
    core = spots["is_caf_core"].to_numpy(bool)
    x = spots["x_pixel"].to_numpy(float) * scale
    y = spots["y_pixel"].to_numpy(float) * scale
    ax.scatter(x, y, s=1.2, c="#333333", alpha=0.10, linewidths=0)
    ax.scatter(x[core], y[core], s=8, facecolors="none", edgecolors="#D73027", linewidths=0.35, alpha=0.75)
    ax.set_title(title, fontsize=8.5, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])


def plot_categorical_spots(ax: plt.Axes, image: Image.Image, spots: pd.DataFrame, scale: float, column: str, title: str, colors: dict[str, str]) -> None:
    ax.imshow(ImageOps.grayscale(image).convert("RGB"), alpha=0.38)
    x = spots["x_pixel"].to_numpy(float) * scale
    y = spots["y_pixel"].to_numpy(float) * scale
    for label, color in colors.items():
        sub = spots[spots[column].eq(label)]
        if len(sub):
            ax.scatter(
                sub["x_pixel"].to_numpy(float) * scale,
                sub["y_pixel"].to_numpy(float) * scale,
                s=4.0,
                c=color,
                alpha=0.78,
                linewidths=0,
                label=label,
            )
    ax.set_title(title, fontsize=8.5, fontweight="bold")
    ax.set_xticks([])
    ax.set_yticks([])


def make_ed15_local_program_maps() -> None:
    rows = selected_spatial_samples()
    source_rows = []
    fig, axes = plt.subplots(len(rows), 3, figsize=(10.8, 3.1 * len(rows)))
    for r, row in rows.iterrows():
        spots = classify_spot_compartments(load_spots(row["dataset_id"], row["sample_id"]))
        image_path, scale_path = image_paths_for_sample(row["dataset_id"], row["sample_id"], row.get("geo_base", ""))
        if image_path is None or not image_path.exists():
            continue
        scale = read_scale(scale_path or Path(""))
        image = Image.open(image_path).convert("RGB")
        plot_he(axes[r, 0], image, spots, scale, f"{row['label']} | H&E + CAF core")
        plot_categorical_spots(
            axes[r, 1],
            image,
            spots,
            scale,
            "dominant_local_program",
            "dominant local program",
            {k: PROGRAM_COLORS[k] for k in DOMINANT_PROGRAMS},
        )
        plot_categorical_spots(
            axes[r, 2],
            image,
            spots,
            scale,
            "spatial_compartment",
            "CAF/tumor compartment",
            {k: PROGRAM_COLORS[k] for k in ["CAF core", "tumor-high", "interface", "other"]},
        )
        if r == 0:
            panel_label(axes[r, 0], "A")
            panel_label(axes[r, 1], "B")
            panel_label(axes[r, 2], "C")
        source_rows.append(
            spots[
                [
                    "dataset_id",
                    "sample_id",
                    "barcode",
                    "x_pixel",
                    "y_pixel",
                    "dominant_local_program",
                    "spatial_compartment",
                    "is_caf_core",
                    "is_tumor_high",
                ]
            ]
        )
    axes[0, 1].legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.24), ncol=4)
    axes[0, 2].legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.24), ncol=4)
    fig.suptitle("Local spatial program and CAF/tumor compartment maps", fontsize=13.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    if source_rows:
        pd.concat(source_rows, ignore_index=True).to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_15_local_program_maps.csv", index=False)
    save_figure(fig, ED15_OUT)


def make_ed16_interface_maps() -> None:
    rows = selected_spatial_samples().tail(3).reset_index(drop=True)
    fig, axes = plt.subplots(len(rows), 4, figsize=(13.2, 3.2 * len(rows)))
    source_rows = []
    for r, row in rows.iterrows():
        spots = classify_spot_compartments(load_spots(row["dataset_id"], row["sample_id"]))
        image_path, scale_path = image_paths_for_sample(row["dataset_id"], row["sample_id"], row.get("geo_base", ""))
        if image_path is None or not image_path.exists():
            continue
        scale = read_scale(scale_path or Path(""))
        image = Image.open(image_path).convert("RGB")
        plot_categorical_spots(
            axes[r, 0],
            image,
            spots,
            scale,
            "spatial_compartment",
            f"{row['label']} compartment",
            {k: PROGRAM_COLORS[k] for k in ["CAF core", "tumor-high", "interface", "other"]},
        )
        for c, (col, title) in enumerate(
            [
                ("z_spp1_tam", "SPP1/TAM"),
                ("z_tgfb_pathway", "TGF-beta"),
                ("score_tumor_aggressive", "tumor aggressive"),
            ],
            start=1,
        ):
            ax = axes[r, c]
            ax.imshow(ImageOps.grayscale(image).convert("RGB"), alpha=0.35)
            vals = spots[col].astype(float).to_numpy()
            lo, hi = np.nanpercentile(vals, [2, 98])
            sc = ax.scatter(
                spots["x_pixel"].to_numpy(float) * scale,
                spots["y_pixel"].to_numpy(float) * scale,
                c=vals,
                cmap="magma",
                vmin=lo,
                vmax=hi,
                s=4,
                linewidths=0,
                alpha=0.78,
            )
            interface = spots["spatial_compartment"].eq("interface").to_numpy()
            ax.scatter(
                spots.loc[interface, "x_pixel"].to_numpy(float) * scale,
                spots.loc[interface, "y_pixel"].to_numpy(float) * scale,
                facecolors="none",
                edgecolors="#55FFFF",
                s=12,
                linewidths=0.35,
                alpha=0.75,
            )
            ax.set_title(title, fontsize=8.5, fontweight="bold")
            ax.set_xticks([])
            ax.set_yticks([])
            if r == len(rows) - 1:
                fig.colorbar(sc, ax=ax, fraction=0.035, pad=0.02)
        if r == 0:
            for c, label in enumerate(["A", "B", "C", "D"]):
                panel_label(axes[r, c], label)
        source_rows.append(
            spots[["dataset_id", "sample_id", "barcode", "x_pixel", "y_pixel", "spatial_compartment", "z_spp1_tam", "z_tgfb_pathway", "score_tumor_aggressive"]]
        )
    axes[0, 0].legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, 1.24), ncol=4)
    fig.suptitle("Tumor-stroma interface compartments and candidate programs", fontsize=13.5, fontweight="bold", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    if source_rows:
        pd.concat(source_rows, ignore_index=True).to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_16_interface_compartments.csv", index=False)
    save_figure(fig, ED16_OUT)


def load_image_for_patch(row: pd.Series, manifest: pd.DataFrame) -> tuple[Image.Image | None, float]:
    match = manifest[(manifest["dataset_id"].eq(row["dataset_id"])) & (manifest["sample_id"].eq(row["sample_id"]))]
    if match.empty:
        return None, 1.0
    image_path = Path(str(match.iloc[0]["image_path"]))
    scale_path = Path(str(match.iloc[0]["scalefactors_path"]))
    if not image_path.exists():
        return None, 1.0
    return Image.open(image_path).convert("RGB"), read_scale(scale_path)


def crop_patch(image: Image.Image, x: float, y: float, radius: int = 48) -> Image.Image:
    left = int(max(0, x - radius))
    upper = int(max(0, y - radius))
    right = int(min(image.width, x + radius))
    lower = int(min(image.height, y + radius))
    patch = image.crop((left, upper, right, lower))
    return ImageOps.pad(patch, (96, 96), color="white")


def make_ed17_he_patch_examples() -> None:
    patch = pd.read_csv(TABLE_DIR / "mvp_he_patch_morphology_features.csv")
    patch = patch[patch["analysis_eligible"].astype(bool)].copy()
    manifest = pd.read_csv(MANIFEST)
    categories = [
        ("CAF-myeloid high", "score_caf_myeloid_barrier", False),
        ("tumor-aggressive high", "score_tumor_aggressive", False),
        ("IFN/MHC high", "z_ifn_antigen_presentation", False),
        ("CAF-myeloid low", "score_caf_myeloid_barrier", True),
    ]
    fig, axes = plt.subplots(len(categories), 8, figsize=(12.5, 6.8))
    source_rows = []
    used = set()
    for r, (label, col, ascending) in enumerate(categories):
        ranked = patch.sort_values(col, ascending=ascending)
        picked = []
        for _, row in ranked.iterrows():
            key = (row["dataset_id"], row["sample_id"], row["barcode"])
            if key in used:
                continue
            image, scale = load_image_for_patch(row, manifest)
            if image is None:
                continue
            picked.append((row, image, scale))
            used.add(key)
            if len(picked) == 8:
                break
        for c, (row, image, scale) in enumerate(picked):
            axes[r, c].imshow(crop_patch(image, float(row["x_pixel"]) * scale, float(row["y_pixel"]) * scale))
            axes[r, c].set_xticks([])
            axes[r, c].set_yticks([])
            axes[r, c].set_title(f"{row[col]:.2f}", fontsize=6.5)
            source_rows.append(row)
        axes[r, 0].set_ylabel(label, fontsize=8.5, fontweight="bold")
    fig.suptitle("Representative H&E patches from program-high and program-low regions", fontsize=13.5, fontweight="bold", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if source_rows:
        pd.DataFrame(source_rows).to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_17_he_patch_examples.csv", index=False)
    save_figure(fig, ED17_OUT)


def assign_xenium_program_labels(cells: pd.DataFrame) -> pd.DataFrame:
    cells = cells.copy()
    label = np.full(len(cells), "other", dtype=object)
    scores = np.column_stack([cells[col].astype(float).to_numpy() for col in XENIUM_PROGRAMS.values()])
    names = list(XENIUM_PROGRAMS.keys())
    for sample_id, idx in cells.groupby("geo_accession").groups.items():
        idx = np.array(list(idx))
        sample_scores = scores[idx]
        max_idx = np.nanargmax(sample_scores, axis=1)
        max_val = sample_scores[np.arange(len(idx)), max_idx]
        cut = np.nanquantile(sample_scores, 0.80)
        sample_label = np.array(["other"] * len(idx), dtype=object)
        sample_label[max_val >= cut] = [names[i] for i in max_idx[max_val >= cut]]
        label[idx] = sample_label
    cells["program_defined_label"] = label
    return cells


def neighborhood_composition_for_sample(cells: pd.DataFrame, anchor_col: str) -> pd.DataFrame:
    xy = cells[["x_centroid", "y_centroid"]].to_numpy(float)
    tree = cKDTree(xy)
    nn = tree.query(xy[:: max(1, len(xy) // 50000)], k=2)[0][:, 1]
    radius = float(max(np.nanmedian(nn) * 8, 35.0))
    labels = list(XENIUM_PROGRAMS.keys()) + ["other"]
    rows = []
    rng = np.random.default_rng(20260628)
    for anchor_type in [anchor_col, "random"]:
        if anchor_type == "random":
            anchor_idx = rng.choice(np.arange(len(cells)), size=max(1000, int(len(cells) * 0.10)), replace=False)
        else:
            cutoff = cells[anchor_col].quantile(0.90)
            anchor_idx = np.flatnonzero(cells[anchor_col].to_numpy(float) >= cutoff)
            if len(anchor_idx) > 6000:
                anchor_idx = rng.choice(anchor_idx, size=6000, replace=False)
        neighbor_lists = tree.query_ball_point(xy[anchor_idx], r=radius)
        neighbor_idx = np.unique(np.concatenate([np.asarray(v, dtype=int) for v in neighbor_lists if len(v)]))
        counts = cells.iloc[neighbor_idx]["program_defined_label"].value_counts().reindex(labels, fill_value=0)
        total = counts.sum()
        for label, count in counts.items():
            rows.append(
                {
                    "geo_accession": cells["geo_accession"].iloc[0],
                    "title": cells["title"].iloc[0],
                    "treatment": cells["treatment"].iloc[0],
                    "anchor_type": "CAF-SPP1/TAM" if anchor_type != "random" else "random cells",
                    "radius": radius,
                    "program_defined_label": label,
                    "n_neighbor_cells": int(count),
                    "fraction": float(count / total) if total else np.nan,
                }
            )
    return pd.DataFrame(rows)


def make_ed18_xenium_neighborhoods() -> None:
    cells = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_cell_scores.csv")
    cells = assign_xenium_program_labels(cells)
    rows = []
    for _, sub in cells.groupby("geo_accession"):
        rows.append(neighborhood_composition_for_sample(sub.reset_index(drop=True), "anchor_CAF_SPP1TAM"))
    comp = pd.concat(rows, ignore_index=True)
    comp.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_18_xenium_program_neighborhoods.csv", index=False)
    labels = list(XENIUM_PROGRAMS.keys()) + ["other"]
    colors = {
        "CAF/matrix": "#8C6D31",
        "SPP1/TAM": "#9C755F",
        "IFN/APC": "#3B6FB6",
        "T/NK": "#59A14F",
        "TGF-beta/EMT": "#B279A2",
        "tumor epithelial": "#E15759",
        "other": "#D8D8D8",
    }
    fig = plt.figure(figsize=(13.2, 7.4))
    gs = GridSpec(2, 2, figure=fig, height_ratios=[1.0, 1.0], width_ratios=[1.35, 1.0], hspace=0.48, wspace=0.38)
    ax_a = fig.add_subplot(gs[0, :])
    panel_label(ax_a, "A")
    sample_order = comp[["geo_accession", "title", "treatment"]].drop_duplicates().sort_values(["treatment", "geo_accession"])
    xlabels = []
    x = 0
    ticks = []
    for _, sample in sample_order.iterrows():
        for anchor in ["CAF-SPP1/TAM", "random cells"]:
            sub = comp[(comp["geo_accession"].eq(sample["geo_accession"])) & (comp["anchor_type"].eq(anchor))]
            bottom = 0
            for label in labels:
                frac = float(sub[sub["program_defined_label"].eq(label)]["fraction"].iloc[0])
                ax_a.bar(x, frac, bottom=bottom, color=colors[label], width=0.82)
                bottom += frac
            ticks.append(x)
            xlabels.append(("CAF" if anchor.startswith("CAF") else "rand") + "\n" + sample["geo_accession"].replace("GSM845", ""))
            x += 1
        x += 0.35
    ax_a.set_xticks(ticks, xlabels)
    ax_a.set_ylabel("neighbor fraction", fontsize=8)
    ax_a.set_ylim(0, 1)
    ax_a.set_title("Program-defined neighborhood composition around CAF-SPP1/TAM anchors", fontsize=9.5, fontweight="bold", loc="left")
    ax_a.legend(labels, frameon=False, fontsize=7, ncol=7, loc="upper center", bbox_to_anchor=(0.5, 1.20))
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_b = fig.add_subplot(gs[1, 0])
    panel_label(ax_b, "B")
    pivot = comp.pivot_table(index=["geo_accession", "anchor_type"], columns="program_defined_label", values="fraction").reindex(columns=labels)
    enrich_rows = []
    for sample in sample_order["geo_accession"]:
        caf = pivot.loc[(sample, "CAF-SPP1/TAM")]
        rand = pivot.loc[(sample, "random cells")]
        enrich_rows.append((caf - rand).rename(sample))
    enrich = pd.DataFrame(enrich_rows)
    im = ax_b.imshow(enrich.to_numpy(float), cmap="RdBu_r", vmin=-0.25, vmax=0.25, aspect="auto")
    ax_b.set_yticks(np.arange(len(enrich.index)), [s.replace("GSM845", "") for s in enrich.index])
    ax_b.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
    ax_b.set_title("CAF-anchor minus random-neighbor fraction", fontsize=9.5, fontweight="bold", loc="left")
    for i in range(enrich.shape[0]):
        for j in range(enrich.shape[1]):
            ax_b.text(j, i, f"{enrich.iat[i, j]:.2f}", ha="center", va="center", fontsize=6.5)
    plt.colorbar(im, ax=ax_b, fraction=0.035, pad=0.02)
    ax_c = fig.add_subplot(gs[1, 1])
    panel_label(ax_c, "C")
    treatment_summary = (
        comp[comp["anchor_type"].eq("CAF-SPP1/TAM")]
        .groupby(["treatment", "program_defined_label"], as_index=False)["fraction"]
        .median()
    )
    for treatment, sub in treatment_summary.groupby("treatment"):
        vals = sub.set_index("program_defined_label").reindex(labels)["fraction"].fillna(0)
        ax_c.plot(labels, vals, marker="o", lw=1.8, label=treatment)
    ax_c.set_xticks(np.arange(len(labels)), labels, rotation=35, ha="right")
    ax_c.set_ylabel("median neighbor fraction", fontsize=8)
    ax_c.set_title("Treatment-context neighborhood summary", fontsize=9.5, fontweight="bold", loc="left")
    ax_c.legend(frameon=False, fontsize=7)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.grid(axis="y", color="#E6E6E6", linewidth=0.6)
    fig.suptitle("Program-defined neighborhoods around Xenium CAF-SPP1/TAM domains", fontsize=13.5, fontweight="bold", y=0.985)
    save_figure(fig, ED18_OUT)


def make_ed19_null_diagnostics() -> None:
    mvp_stats = pd.read_csv(TABLE_DIR / "mvp_random_core_permutation_sample_stats.csv")
    mvp_null = pd.read_csv(TABLE_DIR / "mvp_random_core_permutation_null_rhos.csv")
    gse_stats = pd.read_csv(TABLE_DIR / "gse272362_rds_random_core_permutation_sample_stats.csv")
    gse_null = pd.read_csv(TABLE_DIR / "gse272362_rds_random_core_permutation_null_rhos.csv")
    picks = [
        ("GSE282302", mvp_stats, mvp_null, "GSE282302 strong IFN/MHC"),
        ("GSE272362 primary", gse_stats[gse_stats["specimen_type"].eq("primary_tumor")], gse_null, "primary tumor"),
        ("GSE272362 liver", gse_stats[gse_stats["specimen_type"].eq("liver_metastasis")], gse_null, "liver metastasis"),
        ("GSE272362 LN", gse_stats[gse_stats["specimen_type"].eq("lymph_node_metastasis")], gse_null, "lymph node metastasis"),
    ]
    targets = ["ifn_mhc", "immune_core", "tumor_aggressive"]
    fig, axes = plt.subplots(len(picks), len(targets), figsize=(12.8, 9.0), sharex=False, sharey=False)
    source_rows = []
    for r, (_, stats, nulls, title) in enumerate(picks):
        for c, target in enumerate(targets):
            ax = axes[r, c]
            sub = stats[stats["target"].eq(target)].sort_values("delta_vs_null_median").head(1)
            if sub.empty:
                continue
            sample = sub.iloc[0]["sample_id"]
            observed = float(sub.iloc[0]["observed_rho"])
            n = nulls[(nulls["sample_id"].eq(sample)) & (nulls["target"].eq(target))]["null_rho"].astype(float)
            ax.hist(n, bins=35, color="#B8C7D9", edgecolor="white")
            ax.axvline(observed, color="#B23A48", lw=2.0, label="observed")
            ax.axvline(n.median(), color="#333333", lw=1.2, ls="--", label="null median")
            ax.set_title(f"{title}\n{target} | {sample}", fontsize=7.5)
            ax.tick_params(labelsize=7)
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            if r == 0 and c == 0:
                panel_label(ax, "A")
                ax.legend(frameon=False, fontsize=6)
            source_rows.append(sub.assign(null_median=float(n.median()), null_n=len(n)))
    fig.suptitle("Representative same-size random-core null distributions", fontsize=13.5, fontweight="bold", y=0.985)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    if source_rows:
        pd.concat(source_rows, ignore_index=True).to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_19_random_core_null_diagnostics.csv", index=False)
    save_figure(fig, ED19_OUT)


def main() -> int:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)
    make_ed15_local_program_maps()
    make_ed16_interface_maps()
    make_ed17_he_patch_examples()
    make_ed18_xenium_neighborhoods()
    make_ed19_null_diagnostics()
    for out in [ED15_OUT, ED16_OUT, ED17_OUT, ED18_OUT, ED19_OUT]:
        print(f"Wrote {out.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
