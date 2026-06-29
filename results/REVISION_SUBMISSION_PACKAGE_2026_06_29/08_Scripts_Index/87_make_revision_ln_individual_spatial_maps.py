from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.6,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from PIL import Image


PROJECT_ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = PROJECT_ROOT / "data" / "raw" / "GSE272362"
REVISION_ROOT = PROJECT_ROOT / "results" / "revision_2026_06_29"
OUT_BASE = REVISION_ROOT / "figures" / "Extended_Data_Figure_8_LN_Individual_Spatial_Maps"
SOURCE_OUT = REVISION_ROOT / "analysis_outputs" / "extended_data_figure8_ln_individual_spatial_maps_source_data.csv"

PROGRAMS = [
    ("z_ifn_antigen_presentation", "IFN/MHC", "#4C78A8"),
    ("score_immune_hub_core", "immune core", "#59A14F"),
    ("score_tumor_aggressive", "tumor aggressive", "#B04A37"),
]

METRIC_COLUMNS = [
    "delta_vs_null_median__IFN_MHC",
    "delta_vs_null_median__immune_core",
    "delta_vs_null_median__tumor_aggressive",
    "immune_decoupling_index",
    "dominant_nmf_ecotype_label",
]


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def read_json(path: Path) -> dict:
    with open_text(path) as handle:
        return json.load(handle)


def read_image(path: Path) -> Image.Image:
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as handle:
            return Image.open(handle).convert("RGB")
    return Image.open(path).convert("RGB")


def bare_barcode(value: str) -> str:
    return str(value).split("_")[-1]


def base_from_position_path(path: Path) -> str:
    name = path.name
    if name.endswith("_tissue_positions.csv.gz"):
        return name[: -len("_tissue_positions.csv.gz")]
    return name[: -len("_tissue_positions.csv")]


def paired_path(base: str, stem: str) -> Path:
    plain = RAW_DIR / f"{base}_{stem}"
    gz = RAW_DIR / f"{base}_{stem}.gz"
    return plain if plain.exists() else gz


def read_positions(path: Path) -> pd.DataFrame:
    with open_text(path) as handle:
        df = pd.read_csv(handle)
    if "barcode" not in df.columns and df.shape[1] == 6:
        with open_text(path) as handle:
            df = pd.read_csv(
                handle,
                header=None,
                names=[
                    "barcode",
                    "in_tissue",
                    "array_row",
                    "array_col",
                    "pxl_row_in_fullres",
                    "pxl_col_in_fullres",
                ],
            )
    return df[["barcode", "in_tissue", "pxl_row_in_fullres", "pxl_col_in_fullres"]].copy()


def build_position_cache() -> dict[str, dict]:
    cache: dict[str, dict] = {}
    paths = sorted(RAW_DIR.glob("*_tissue_positions.csv")) + sorted(RAW_DIR.glob("*_tissue_positions.csv.gz"))
    seen: set[str] = set()
    for pos_path in paths:
        base = base_from_position_path(pos_path)
        if base in seen:
            continue
        seen.add(base)
        positions = read_positions(pos_path)
        positions["bare_barcode"] = positions["barcode"].map(bare_barcode)
        cache[base] = {
            "positions_path": pos_path,
            "positions": positions,
            "image_path": paired_path(base, "tissue_hires_image.png"),
            "scalefactors_path": paired_path(base, "scalefactors_json.json"),
        }
    return cache


def read_hires_scale(path: Path) -> float:
    payload = read_json(path) if path.exists() else {}
    return float(payload.get("tissue_hires_scalef", 1.0))


def match_sample_to_geo(sample_spots: pd.DataFrame, position_cache: dict[str, dict]) -> dict:
    query = sample_spots[["barcode", "x_pixel", "y_pixel"]].copy()
    query["bare_barcode"] = query["barcode"].map(bare_barcode)
    if len(query) > 600:
        query = query.sample(n=600, random_state=20260629)

    best: dict | None = None
    for base, payload in position_cache.items():
        merged = query.merge(payload["positions"], on="bare_barcode", how="inner")
        if merged.empty:
            continue
        exact = (
            (merged["x_pixel"].astype(float) - merged["pxl_col_in_fullres"].astype(float)).abs() <= 1
        ) & ((merged["y_pixel"].astype(float) - merged["pxl_row_in_fullres"].astype(float)).abs() <= 1)
        candidate = {
            **payload,
            "geo_base": base,
            "n_barcode_matches": int(len(merged)),
            "n_exact_coordinate_matches": int(exact.sum()),
            "match_fraction": float(exact.sum() / max(1, len(query))),
        }
        if best is None or candidate["n_exact_coordinate_matches"] > best["n_exact_coordinate_matches"]:
            best = candidate

    if best is None or best["n_exact_coordinate_matches"] < 20:
        sample_id = sample_spots["sample_id"].iloc[0]
        raise ValueError(f"Could not confidently match {sample_id} to a GEO H&E image.")
    return best


def make_tinted_cmap(hex_color: str) -> LinearSegmentedColormap:
    return LinearSegmentedColormap.from_list("program_tint", ["#F3F3F3", hex_color], N=256)


def robust_norm(values: pd.Series) -> Normalize:
    lo, hi = np.nanpercentile(values.to_numpy(float), [2, 98])
    if not np.isfinite(lo) or not np.isfinite(hi) or np.isclose(lo, hi):
        lo, hi = float(np.nanmin(values)), float(np.nanmax(values))
    if np.isclose(lo, hi):
        hi = lo + 1.0
    return Normalize(vmin=lo, vmax=hi)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.025, 1.04, label, transform=ax.transAxes, fontsize=11, fontweight="bold", ha="left", va="bottom")


def show_background(ax: plt.Axes, image: Image.Image, x: np.ndarray, y: np.ndarray, pad: float = 500.0) -> None:
    ax.imshow(image, alpha=0.92)
    ax.set_xlim(float(np.nanmin(x) - pad), float(np.nanmax(x) + pad))
    ax.set_ylim(float(np.nanmax(y) + pad), float(np.nanmin(y) - pad))
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def short_label(text: str) -> str:
    text = str(text)
    return text.replace(" / ", "/").replace("basal-like/", "basal/").replace("plasma cell", "plasma")


def make_figure() -> None:
    spots = pd.read_csv(PROJECT_ROOT / "results" / "tables" / "gse272362_rds_spot_level_scores.csv")
    ln_stats = pd.read_csv(REVISION_ROOT / "analysis_outputs" / "ln_metastasis_individual_sample_summary.csv")
    ln_stats = ln_stats.sort_values("immune_decoupling_index", ascending=False).reset_index(drop=True)
    sample_order = ln_stats["sample_id"].tolist()
    spots = spots[spots["sample_id"].isin(sample_order)].copy()

    position_cache = build_position_cache()
    matches = {sample_id: match_sample_to_geo(spots[spots["sample_id"].eq(sample_id)], position_cache) for sample_id in sample_order}
    images = {sample_id: read_image(Path(matches[sample_id]["image_path"])) for sample_id in sample_order}
    scales = {sample_id: read_hires_scale(Path(matches[sample_id]["scalefactors_path"])) for sample_id in sample_order}

    source_rows: list[pd.DataFrame] = []
    fig = plt.figure(figsize=(14.8, 11.4), constrained_layout=False)
    gs = fig.add_gridspec(
        nrows=len(sample_order) + 1,
        ncols=5,
        height_ratios=[0.18] + [1.0] * len(sample_order),
        width_ratios=[1.06, 1.0, 1.0, 1.0, 0.88],
        hspace=0.02,
        wspace=0.05,
    )

    headers = ["H&E + CAF core", "IFN/MHC", "immune core", "tumor aggressive", "per-sample metrics"]
    for col, header in enumerate(headers):
        ax = fig.add_subplot(gs[0, col])
        ax.axis("off")
        ax.text(0.5, 0.25, header, ha="center", va="center", fontsize=10.5, fontweight="bold")

    panel_counter = 0
    for row_idx, sample_id in enumerate(sample_order, start=1):
        sample = spots[spots["sample_id"].eq(sample_id)].copy()
        sample["caf_core_top10"] = sample["score_caf_myeloid_barrier"] >= np.nanpercentile(
            sample["score_caf_myeloid_barrier"].to_numpy(float), 90
        )
        x = sample["x_pixel"].to_numpy(float) * scales[sample_id]
        y = sample["y_pixel"].to_numpy(float) * scales[sample_id]
        core = sample["caf_core_top10"].to_numpy(bool)

        ax = fig.add_subplot(gs[row_idx, 0])
        show_background(ax, images[sample_id], x, y)
        ax.scatter(x, y, s=2.0, c="#181818", alpha=0.16, linewidths=0)
        ax.scatter(x[core], y[core], s=12, facecolors="none", edgecolors="#D64535", linewidths=0.55, alpha=0.9)
        ax.text(
            0.015,
            0.985,
            sample_id,
            transform=ax.transAxes,
            ha="left",
            va="top",
            fontsize=8.6,
            fontweight="bold",
            bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.72, "pad": 1.2},
        )
        panel_counter += 1
        add_panel_label(ax, chr(ord("a") + panel_counter - 1))

        for col_idx, (score_col, _, color) in enumerate(PROGRAMS, start=1):
            ax = fig.add_subplot(gs[row_idx, col_idx])
            show_background(ax, images[sample_id], x, y)
            values = sample[score_col].astype(float)
            sc = ax.scatter(
                x,
                y,
                c=values,
                s=3.2,
                cmap=make_tinted_cmap(color),
                norm=robust_norm(values),
                alpha=0.82,
                linewidths=0,
            )
            ax.scatter(x[core], y[core], s=10, facecolors="none", edgecolors="#111111", linewidths=0.35, alpha=0.42)
            if row_idx == len(sample_order):
                cbar = fig.colorbar(sc, ax=ax, orientation="horizontal", fraction=0.07, pad=0.025)
                cbar.ax.tick_params(labelsize=6, length=2)
            panel_counter += 1
            add_panel_label(ax, chr(ord("a") + panel_counter - 1))

        metric_ax = fig.add_subplot(gs[row_idx, 4])
        metric_ax.axis("off")
        row = ln_stats[ln_stats["sample_id"].eq(sample_id)].iloc[0]
        lines = [
            f"IFN/MHC delta  {row['delta_vs_null_median__IFN_MHC']:+.2f}",
            f"immune delta   {row['delta_vs_null_median__immune_core']:+.2f}",
            f"tumor delta    {row['delta_vs_null_median__tumor_aggressive']:+.2f}",
            f"decoupling     {row['immune_decoupling_index']:.2f}",
            short_label(row["dominant_nmf_ecotype_label"]),
            f"spots n={len(sample):,}",
            f"image match={matches[sample_id]['n_exact_coordinate_matches']}",
        ]
        metric_ax.text(
            0.02,
            0.5,
            "\n".join(lines),
            ha="left",
            va="center",
            fontsize=8.2,
            linespacing=1.55,
        )

        source = sample[
            [
                "dataset_id",
                "sample_id",
                "patient_id",
                "specimen_type",
                "barcode",
                "x_pixel",
                "y_pixel",
                "score_caf_myeloid_barrier",
                "z_ifn_antigen_presentation",
                "score_immune_hub_core",
                "score_tumor_aggressive",
                "caf_core_top10",
            ]
        ].copy()
        for metric in METRIC_COLUMNS:
            source[metric] = row[metric]
        source["geo_base"] = matches[sample_id]["geo_base"]
        source["image_match_fraction"] = matches[sample_id]["match_fraction"]
        source_rows.append(source)

    fig.suptitle(
        "Individual lymph-node metastasis maps show retained tumor coupling with heterogeneous immune organization",
        fontsize=14,
        fontweight="bold",
        y=0.991,
    )
    fig.text(
        0.5,
        0.008,
        "Rows show all five GSE272362 lymph-node metastasis samples, ordered by immune-decoupling index. Red rings mark top 10% CAF-myeloid spots; black rings in program maps mark the same CAF core. Deltas are observed distance-to-core Spearman rho minus same-size random-core median, where negative values indicate stronger CAF-core centering.",
        ha="center",
        va="bottom",
        fontsize=8.4,
        color="#333333",
    )

    OUT_BASE.parent.mkdir(parents=True, exist_ok=True)
    SOURCE_OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_BASE.with_suffix(".png"), dpi=320, bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT_BASE.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)

    pd.concat(source_rows, ignore_index=True).to_csv(SOURCE_OUT, index=False, encoding="utf-8")
    print(OUT_BASE.with_suffix(".pdf"))
    print(SOURCE_OUT)


if __name__ == "__main__":
    make_figure()
