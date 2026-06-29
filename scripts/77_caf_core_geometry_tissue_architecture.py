from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import ConvexHull, cKDTree
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
REPORT_DIR = PROJECT / "results" / "reports"
SOURCE_DIR = PROJECT / "results" / "source_data"
for d in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

CAF_CORE_FRACTION = 0.10
TUMOR_HIGH_FRACTION = 0.20
NEIGHBOR_RADIUS_SCALE = 1.55
INTERFACE_RADIUS_SCALE = 2.0

DATASETS = [
    ("GSE282302", "results/tables/mvp_spot_level_scores_with_edge_qc.csv", "post-NACT"),
    ("GSE272362", "results/tables/gse272362_rds_spot_level_scores.csv", None),
    ("GSE235315", "results/tables/gse235315_spot_level_scores.csv", "paired ST-H&E"),
    ("GSE274557", "results/tables/gse274557_full_spot_scores.csv", None),
]

CONTEXT_ORDER = [
    "post-NACT",
    "treatment-naive primary",
    "primary tumor",
    "liver metastasis",
    "lymph-node metastasis",
    "other metastasis",
    "paired ST-H&E",
]
CONTEXT_COLORS = {
    "post-NACT": "#8C6D31",
    "treatment-naive primary": "#B55A30",
    "primary tumor": "#4E79A7",
    "liver metastasis": "#59A14F",
    "lymph-node metastasis": "#9C755F",
    "other metastasis": "#6F4E7C",
    "paired ST-H&E": "#7F7F7F",
}


def context_label(dataset_id: str, sample: pd.DataFrame, default_context: str | None) -> str:
    if dataset_id == "GSE282302":
        if sample["dataset_id"].astype(str).eq("GSE274103").any():
            return "treatment-naive primary"
        return "post-NACT"
    if dataset_id == "GSE235315":
        return "paired ST-H&E"
    if dataset_id == "GSE272362":
        value = str(sample["specimen_type"].dropna().iloc[0]) if "specimen_type" in sample and sample["specimen_type"].notna().any() else ""
        return {
            "primary_tumor": "primary tumor",
            "liver_metastasis": "liver metastasis",
            "lymph_node_metastasis": "lymph-node metastasis",
            "normal_pancreas": "normal pancreas",
        }.get(value, value or "GSE272362")
    if dataset_id == "GSE274557":
        tissue = str(sample["tissue"].dropna().iloc[0]) if "tissue" in sample and sample["tissue"].notna().any() else ""
        if "Primary" in tissue:
            return "primary tumor"
        if "Liver" in tissue:
            return "liver metastasis"
        if "Lung" in tissue or "Peritoneal" in tissue:
            return "other metastasis"
        return tissue or "GSE274557"
    return default_context or dataset_id


def load_spots() -> pd.DataFrame:
    frames = []
    keep_base = {
        "dataset_id",
        "sample_id",
        "patient_id",
        "barcode",
        "x_pixel",
        "y_pixel",
        "score_caf_myeloid_barrier",
        "score_tumor_epithelial",
        "score_tumor_aggressive",
        "score_spp1_tam",
        "score_ifn_antigen_presentation",
        "score_immune_hub_core",
        "score_mycaf",
        "score_tgfb_pathway",
        "score_emt_invasion",
    }
    for dataset_id, rel_path, default_context in DATASETS:
        df = pd.read_csv(PROJECT / rel_path)
        if "edge_or_background_risk" in df.columns:
            df = df.loc[~df["edge_or_background_risk"].fillna(False).astype(bool)].copy()
        keep = [c for c in df.columns if c in keep_base or c in {"specimen_type", "tissue", "treatment", "title", "geo_accession"}]
        df = df[keep].copy()
        df["analysis_dataset"] = dataset_id
        df["default_context"] = default_context or dataset_id
        frames.append(df)
    spots = pd.concat(frames, ignore_index=True)
    spots = spots[np.isfinite(spots["x_pixel"]) & np.isfinite(spots["y_pixel"])].copy()
    return spots


def median_nn_scale(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return 1.0
    tree = cKDTree(xy)
    d, _ = tree.query(xy, k=2)
    val = float(np.nanmedian(d[:, 1]))
    return val if np.isfinite(val) and val > 0 else 1.0


def nearest_distance(xy: np.ndarray, mask: np.ndarray) -> np.ndarray:
    tree = cKDTree(xy[mask])
    d, _ = tree.query(xy, k=1)
    return d


def hull_compactness(points: np.ndarray) -> tuple[float, float, float]:
    if len(points) < 4:
        return np.nan, np.nan, np.nan
    try:
        hull = ConvexHull(points)
        area = float(hull.volume)
        perimeter = float(hull.area)
        compact = float((4 * np.pi * area) / (perimeter * perimeter)) if perimeter > 0 else np.nan
        return area, perimeter, compact
    except Exception:
        return np.nan, np.nan, np.nan


def component_labels(xy: np.ndarray, core_mask: np.ndarray, radius: float) -> tuple[np.ndarray, int, np.ndarray]:
    core_idx = np.where(core_mask)[0]
    labels_full = np.full(len(xy), -1, dtype=int)
    if len(core_idx) == 0:
        return labels_full, 0, np.array([], dtype=int)
    if len(core_idx) == 1:
        labels_full[core_idx[0]] = 0
        return labels_full, 1, np.array([1], dtype=int)
    core_xy = xy[core_idx]
    pairs = list(cKDTree(core_xy).query_pairs(radius))
    if pairs:
        row = [p[0] for p in pairs] + [p[1] for p in pairs]
        col = [p[1] for p in pairs] + [p[0] for p in pairs]
        data = np.ones(len(row), dtype=int)
        graph = coo_matrix((data, (row, col)), shape=(len(core_idx), len(core_idx)))
    else:
        graph = coo_matrix((len(core_idx), len(core_idx)))
    n_comp, labels = connected_components(graph, directed=False, return_labels=True)
    labels_full[core_idx] = labels
    sizes = np.bincount(labels, minlength=n_comp)
    return labels_full, int(n_comp), sizes


def analyze_sample(dataset_id: str, sample_id: str, sample: pd.DataFrame) -> tuple[dict | None, pd.DataFrame]:
    sample = sample.dropna(subset=["x_pixel", "y_pixel", "score_caf_myeloid_barrier", "score_tumor_epithelial"]).copy()
    if len(sample) < 120:
        return None, pd.DataFrame()
    xy = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    scale = median_nn_scale(xy)
    caf = pd.to_numeric(sample["score_caf_myeloid_barrier"], errors="coerce").to_numpy(float)
    tumor = pd.to_numeric(sample["score_tumor_epithelial"], errors="coerce").to_numpy(float)
    core_mask = caf >= np.nanquantile(caf, 1 - CAF_CORE_FRACTION)
    tumor_mask = tumor >= np.nanquantile(tumor, 1 - TUMOR_HIGH_FRACTION)
    if core_mask.sum() < 20 or tumor_mask.sum() < 20:
        return None, pd.DataFrame()

    neighbor_radius = NEIGHBOR_RADIUS_SCALE * scale
    interface_radius = INTERFACE_RADIUS_SCALE * scale
    labels, n_components, sizes = component_labels(xy, core_mask, neighbor_radius)
    largest_label = int(np.argmax(sizes)) if len(sizes) else -1
    largest_mask = labels == largest_label
    largest_area, largest_perim, largest_compact = hull_compactness(xy[largest_mask])

    all_tree = cKDTree(xy)
    core_idx = np.where(core_mask)[0]
    boundary_contacts = 0
    tumor_boundary_contacts = 0
    for idx in core_idx:
        neighbors = all_tree.query_ball_point(xy[idx], r=neighbor_radius)
        neighbors = [n for n in neighbors if n != idx]
        boundary_contacts += int(np.sum(~core_mask[neighbors]))
        tumor_boundary_contacts += int(np.sum(tumor_mask[neighbors]))

    d_tumor = nearest_distance(xy, tumor_mask) / scale
    d_core = nearest_distance(xy, core_mask) / scale
    interface_core = core_mask & (d_tumor <= INTERFACE_RADIUS_SCALE)
    tumor_near_core = tumor_mask & (d_core <= INTERFACE_RADIUS_SCALE)

    context = context_label(dataset_id, sample, str(sample["default_context"].iloc[0]))
    core_area, core_perim, core_compact = hull_compactness(xy[core_mask])
    row = {
        "dataset_id": dataset_id,
        "context": context,
        "sample_id": sample_id,
        "n_spots": int(len(sample)),
        "n_core_spots": int(core_mask.sum()),
        "n_tumor_high_spots": int(tumor_mask.sum()),
        "core_fraction": float(core_mask.mean()),
        "median_nn_scale": float(scale),
        "n_core_components": int(n_components),
        "fragmentation_per_100_core_spots": float(n_components / core_mask.sum() * 100),
        "largest_component_fraction": float(sizes.max() / core_mask.sum()) if len(sizes) else np.nan,
        "core_convex_hull_area_scaled": float(core_area / (scale * scale)) if np.isfinite(core_area) else np.nan,
        "core_convex_hull_compactness": float(core_compact),
        "largest_component_area_scaled": float(largest_area / (scale * scale)) if np.isfinite(largest_area) else np.nan,
        "largest_component_compactness": float(largest_compact),
        "boundary_contacts_per_core_spot": float(boundary_contacts / core_mask.sum()),
        "tumor_boundary_contacts_per_core_spot": float(tumor_boundary_contacts / core_mask.sum()),
        "interface_core_fraction": float(interface_core.sum() / core_mask.sum()),
        "tumor_near_core_fraction": float(tumor_near_core.sum() / tumor_mask.sum()),
        "median_core_to_tumor_distance": float(np.nanmedian(d_tumor[core_mask])),
        "median_tumor_to_core_distance": float(np.nanmedian(d_core[tumor_mask])),
    }

    map_df = sample[["sample_id", "x_pixel", "y_pixel"]].copy()
    map_df["dataset_id"] = dataset_id
    map_df["context"] = context
    map_df["is_caf_core"] = core_mask
    map_df["is_tumor_high"] = tumor_mask
    map_df["is_core_tumor_interface"] = interface_core
    map_df["core_component"] = labels
    map_df["dist_to_tumor_high_norm"] = d_tumor
    map_df["dist_to_caf_core_norm"] = d_core
    return row, map_df


def select_map_samples(map_df: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for context in ["primary tumor", "liver metastasis", "lymph-node metastasis"]:
        sub = map_df[map_df["context"].eq(context)]
        if sub.empty:
            continue
        counts = sub.groupby("sample_id").size().sort_values(ascending=False)
        chosen = None
        for sample_id in counts.index:
            s = sub[sub["sample_id"].eq(sample_id)]
            if s["is_caf_core"].sum() > 30 and s["is_tumor_high"].sum() > 50:
                chosen = s.copy()
                break
        if chosen is not None:
            if len(chosen) > 5000:
                chosen = chosen.sample(5000, random_state=20260628)
            frames.append(chosen)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize(metrics: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    context = (
        metrics.groupby("context", as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_components=("n_core_components", "median"),
            median_fragmentation=("fragmentation_per_100_core_spots", "median"),
            median_largest_component_fraction=("largest_component_fraction", "median"),
            median_compactness=("largest_component_compactness", "median"),
            median_interface_core_fraction=("interface_core_fraction", "median"),
            median_core_to_tumor_distance=("median_core_to_tumor_distance", "median"),
            median_tumor_boundary_contacts=("tumor_boundary_contacts_per_core_spot", "median"),
        )
    )

    bio = pd.read_csv(TABLE_DIR / "mechanism_candidate_axis_sample_summary.csv")
    transition = pd.read_csv(TABLE_DIR / "core_to_interface_transition_sample_summary.csv")
    transition_pivot = transition.pivot_table(index="sample_id", columns="program", values="rho_with_core_to_tumor_axis", aggfunc="median")
    transition_pivot = transition_pivot.add_prefix("transition_rho__").reset_index()
    merged = metrics.merge(
        bio[["sample_id", "stromal_tumor_core_coupling", "immune_core_coupling", "immune_decoupling_index"]],
        on="sample_id",
        how="left",
    ).merge(transition_pivot, on="sample_id", how="left")

    geom_cols = [
        "fragmentation_per_100_core_spots",
        "largest_component_fraction",
        "largest_component_compactness",
        "interface_core_fraction",
        "tumor_boundary_contacts_per_core_spot",
        "median_core_to_tumor_distance",
    ]
    bio_cols = [
        "immune_decoupling_index",
        "stromal_tumor_core_coupling",
        "immune_core_coupling",
        "transition_rho__SPP1/TAM",
        "transition_rho__TGFb/EMT",
        "transition_rho__tumor aggressive",
    ]
    rows = []
    for g in geom_cols:
        for b in bio_cols:
            tmp = merged[[g, b]].replace([np.inf, -np.inf], np.nan).dropna()
            rho = np.nan
            p = np.nan
            if len(tmp) >= 10 and tmp[g].std() > 0 and tmp[b].std() > 0:
                rho, p = spearmanr(tmp[g], tmp[b])
            rows.append({"geometry_metric": g, "biological_readout": b, "n_samples": len(tmp), "spearman_rho": rho, "p_value": p})
    corr = pd.DataFrame(rows)
    return context, corr


def draw_architecture_schematic(ax: plt.Axes) -> None:
    ax.set_title("A  CAF-core tissue architecture metrics", loc="left", fontweight="bold")
    ax.axis("off")
    rng = np.random.default_rng(1)
    core1 = rng.normal([0.30, 0.58], [0.06, 0.08], size=(70, 2))
    core2 = rng.normal([0.56, 0.44], [0.04, 0.05], size=(30, 2))
    tumor = rng.normal([0.74, 0.56], [0.07, 0.10], size=(85, 2))
    ax.scatter(tumor[:, 0], tumor[:, 1], s=16, color="#6F5C99", alpha=0.65, label="tumor-high")
    ax.scatter(core1[:, 0], core1[:, 1], s=18, color="#2D6A8E", alpha=0.9, label="CAF-core")
    ax.scatter(core2[:, 0], core2[:, 1], s=18, color="#2D6A8E", alpha=0.9)
    ax.plot([0.61, 0.66], [0.48, 0.52], color="#B57A3C", lw=4, solid_capstyle="round")
    ax.text(0.12, 0.82, "components\nfragmentation\ncompactness", fontsize=9, ha="left")
    ax.text(0.61, 0.28, "interface\ncontacts", fontsize=9, ha="center", color="#7A4A1C")
    ax.text(0.72, 0.78, "core-tumor\ndistance", fontsize=9, ha="center")
    ax.annotate("", xy=(0.68, 0.66), xytext=(0.39, 0.62), arrowprops=dict(arrowstyle="<->", lw=1.2, color="#333"))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def plot_sample_map(ax: plt.Axes, sample: pd.DataFrame, title: str) -> None:
    ax.scatter(sample["x_pixel"], sample["y_pixel"], s=4, color="#D8D8D8", alpha=0.65, linewidth=0)
    tumor = sample[sample["is_tumor_high"]]
    core = sample[sample["is_caf_core"]]
    interface = sample[sample["is_core_tumor_interface"]]
    ax.scatter(tumor["x_pixel"], tumor["y_pixel"], s=5, color="#6F5C99", alpha=0.7, linewidth=0)
    ax.scatter(core["x_pixel"], core["y_pixel"], s=6, color="#2D6A8E", alpha=0.9, linewidth=0)
    ax.scatter(interface["x_pixel"], interface["y_pixel"], s=8, color="#B57A3C", alpha=0.95, linewidth=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")


def make_figure(metrics: pd.DataFrame, context: pd.DataFrame, corr: pd.DataFrame, map_df: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.8,
        }
    )
    fig = plt.figure(figsize=(16.6, 11.0))
    gs = fig.add_gridspec(3, 4, height_ratios=[0.92, 1.0, 1.0], width_ratios=[1.0, 1.0, 1.0, 1.0], hspace=0.58, wspace=0.48)
    fig.suptitle("CAF-core geometry reveals tissue architecture behind spatial niche organization", fontsize=17, fontweight="bold")

    ax_a = fig.add_subplot(gs[0, 0])
    draw_architecture_schematic(ax_a)

    for i, ctx in enumerate(["primary tumor", "liver metastasis", "lymph-node metastasis"]):
        ax = fig.add_subplot(gs[0, i + 1])
        sub = map_df[map_df["context"].eq(ctx)]
        if sub.empty:
            ax.axis("off")
        else:
            plot_sample_map(ax, sub, ("B  " if i == 0 else "") + ctx)

    metric_plots = [
        ("largest_component_fraction", "C  largest component fraction", "fraction"),
        ("fragmentation_per_100_core_spots", "D  fragmentation", "components / 100 core spots"),
        ("interface_core_fraction", "E  core spots near tumor-high", "fraction"),
        ("median_core_to_tumor_distance", "F  core-to-tumor distance", "spot-neighbor units"),
    ]
    for idx, (metric, title, ylabel) in enumerate(metric_plots):
        ax = fig.add_subplot(gs[1, idx])
        contexts = [c for c in CONTEXT_ORDER if c in set(metrics["context"])]
        vals = [metrics.loc[metrics["context"].eq(c), metric].dropna().to_numpy(float) for c in contexts]
        parts = ax.violinplot(vals, positions=np.arange(len(contexts)), showmeans=False, showmedians=True, widths=0.75)
        for body, ctx in zip(parts["bodies"], contexts):
            body.set_facecolor(CONTEXT_COLORS.get(ctx, "#777"))
            body.set_alpha(0.32)
            body.set_edgecolor("none")
        for key in ["cmedians", "cbars", "cmins", "cmaxes"]:
            if key in parts:
                parts[key].set_color("#333333")
                parts[key].set_linewidth(0.8)
        for x, ctx in enumerate(contexts):
            y = vals[x]
            if len(y):
                jitter = np.random.default_rng(42 + idx + x).normal(x, 0.045, len(y))
                ax.scatter(jitter, y, s=9, color=CONTEXT_COLORS.get(ctx, "#777"), alpha=0.55, linewidth=0)
        ax.set_xticks(np.arange(len(contexts)))
        ax.set_xticklabels(contexts, rotation=35, ha="right", fontsize=8)
        ax.set_title(title, loc="left", fontweight="bold")
        ax.set_ylabel(ylabel)
        ax.grid(axis="y", color="#E5E5E5", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)

    ax_g = fig.add_subplot(gs[2, 0:2])
    merged = metrics.merge(
        pd.read_csv(TABLE_DIR / "mechanism_candidate_axis_sample_summary.csv")[["sample_id", "immune_decoupling_index"]],
        on="sample_id",
        how="left",
    )
    for ctx, sub in merged.dropna(subset=["interface_core_fraction", "immune_decoupling_index"]).groupby("context"):
        if ctx not in CONTEXT_ORDER:
            continue
        ax_g.scatter(
            sub["interface_core_fraction"],
            sub["immune_decoupling_index"],
            s=22,
            color=CONTEXT_COLORS.get(ctx, "#777"),
            alpha=0.70,
            edgecolor="white",
            linewidth=0.25,
            label=ctx,
        )
    tmp = merged[["interface_core_fraction", "immune_decoupling_index"]].dropna()
    if len(tmp) >= 10:
        rho, p = spearmanr(tmp["interface_core_fraction"], tmp["immune_decoupling_index"])
        coeff = np.polyfit(tmp["interface_core_fraction"], tmp["immune_decoupling_index"], deg=1)
        xx = np.linspace(tmp["interface_core_fraction"].min(), tmp["interface_core_fraction"].max(), 100)
        ax_g.plot(xx, coeff[0] * xx + coeff[1], color="#222222", lw=1)
        ax_g.text(0.03, 0.97, f"rho = {rho:.2f}\np = {p:.2g}", transform=ax_g.transAxes, ha="left", va="top", fontsize=9, bbox=dict(fc="white", ec="#DDDDDD", boxstyle="round,pad=0.25", lw=0.6))
    ax_g.set_title("G  Interface geometry versus immune decoupling", loc="left", fontweight="bold")
    ax_g.set_xlabel("CAF-core fraction within 2 spot-neighbor units of tumor-high")
    ax_g.set_ylabel("immune-decoupling index")
    ax_g.legend(frameon=False, fontsize=7, ncol=2)
    ax_g.grid(color="#E5E5E5", lw=0.6)
    ax_g.spines[["top", "right"]].set_visible(False)

    ax_h = fig.add_subplot(gs[2, 2:4])
    geom_order = [
        "fragmentation_per_100_core_spots",
        "largest_component_fraction",
        "largest_component_compactness",
        "interface_core_fraction",
        "median_core_to_tumor_distance",
    ]
    bio_order = [
        "immune_decoupling_index",
        "stromal_tumor_core_coupling",
        "immune_core_coupling",
        "transition_rho__SPP1/TAM",
        "transition_rho__tumor aggressive",
    ]
    mat = corr.pivot(index="geometry_metric", columns="biological_readout", values="spearman_rho").reindex(index=geom_order, columns=bio_order)
    im = ax_h.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.55, vmax=0.55, aspect="auto")
    ax_h.set_yticks(np.arange(len(geom_order)), [g.replace("_", " ").replace("per 100 core spots", "/100 core") for g in geom_order], fontsize=8)
    ax_h.set_xticks(np.arange(len(bio_order)), [b.replace("transition_rho__", "rho ").replace("_", " ") for b in bio_order], rotation=35, ha="right", fontsize=8)
    ax_h.set_title("H  Geometry-biological readout associations", loc="left", fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.notna(val):
                ax_h.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="white" if abs(val) > 0.38 else "#111")
    cb = fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.03)
    cb.set_label("Spearman rho", fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / "extended_data_figure33_caf_core_geometry_tissue_architecture"
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(metrics: pd.DataFrame, context: pd.DataFrame, corr: pd.DataFrame) -> None:
    lines = [
        "# CAF-Core Geometry And Tissue Architecture",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Question",
        "",
        "This analysis asks whether CAF-myeloid cores have interpretable tissue architecture rather than behaving only as high-scoring spots.",
        "",
        "## Method",
        "",
        f"CAF-core spots were defined as the top {int(CAF_CORE_FRACTION * 100)}% CAF-myeloid-barrier score and tumor-high spots as the top {int(TUMOR_HIGH_FRACTION * 100)}% tumor-epithelial score. CAF-core connected components were built using a radius of {NEIGHBOR_RADIUS_SCALE} median spot-neighbor distances. Interface contacts were defined as CAF-core spots within {INTERFACE_RADIUS_SCALE} median spot-neighbor distances of tumor-high spots.",
        "",
        "## Context-Level Geometry",
        "",
    ]
    for _, row in context.sort_values("context").iterrows():
        lines.append(
            f"- {row.context}: n={int(row.n_samples)}, median components {row.median_components:.1f}, "
            f"largest-component fraction {row.median_largest_component_fraction:.2f}, "
            f"interface-core fraction {row.median_interface_core_fraction:.2f}, "
            f"core-to-tumor distance {row.median_core_to_tumor_distance:.2f}."
        )
    lines.extend(["", "## Strongest Geometry Associations", ""])
    top = corr.dropna(subset=["spearman_rho"]).assign(abs_rho=lambda x: x["spearman_rho"].abs()).sort_values("abs_rho", ascending=False).head(10)
    for _, row in top.iterrows():
        lines.append(
            f"- {row.geometry_metric} vs {row.biological_readout}: rho {row.spearman_rho:.2f}, p={row.p_value:.2g}, n={int(row.n_samples)}."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This module adds a tissue-architecture layer to the CAF-myeloid niche model. It supports interpretation of CAF cores as spatial structures with fragmentation, compactness and tumor-interface geometry. It remains expression-defined and should not be described as pathologist-annotated histologic compartments.",
            "",
            "## Outputs",
            "",
            "- `results/tables/caf_core_geometry_metrics_per_sample.csv`",
            "- `results/tables/caf_core_geometry_context_summary.csv`",
            "- `results/tables/caf_core_geometry_biological_correlations.csv`",
            "- `results/source_data/Source_Data_Extended_Data_Fig_33_caf_core_geometry.csv`",
            "- `results/figures/submission/extended_data_figure33_caf_core_geometry_tissue_architecture.pdf`",
        ]
    )
    (REPORT_DIR / "caf_core_geometry_tissue_architecture_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    spots = load_spots()
    rows = []
    map_frames = []
    for (dataset_id, sample_id), sample in spots.groupby(["analysis_dataset", "sample_id"], sort=False):
        row, sample_map = analyze_sample(str(dataset_id), str(sample_id), sample)
        if row is not None:
            rows.append(row)
            map_frames.append(sample_map)
        print(f"{dataset_id} {sample_id}: {'ok' if row is not None else 'skipped'}")
    metrics = pd.DataFrame(rows)
    maps = pd.concat(map_frames, ignore_index=True) if map_frames else pd.DataFrame()
    selected_maps = select_map_samples(maps)
    context, corr = summarize(metrics)

    metrics.to_csv(TABLE_DIR / "caf_core_geometry_metrics_per_sample.csv", index=False)
    context.to_csv(TABLE_DIR / "caf_core_geometry_context_summary.csv", index=False)
    corr.to_csv(TABLE_DIR / "caf_core_geometry_biological_correlations.csv", index=False)
    selected_maps.to_csv(TABLE_DIR / "caf_core_geometry_selected_spot_maps.csv", index=False)

    metrics.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_33_caf_core_geometry.csv", index=False)
    context.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_33_context_summary.csv", index=False)
    corr.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_33_biological_correlations.csv", index=False)
    selected_maps.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_33_selected_spot_maps.csv", index=False)

    make_figure(metrics, context, corr, selected_maps)
    write_report(metrics, context, corr)
    print("Wrote CAF-core geometry tissue-architecture outputs.")


if __name__ == "__main__":
    main()
