from __future__ import annotations

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


PROJECT = Path(__file__).resolve().parents[1]
DATA_DIR = PROJECT / "data" / "external" / "GSE277782"
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "internal"
REPORT_DIR = PROJECT / "results" / "reports"
SOURCE_DIR = PROJECT / "results" / "source_data"
for path in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR]:
    path.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260627)
N_RANDOM = 1000


def read_slide_metadata(slide: int) -> pd.DataFrame:
    path = DATA_DIR / f"GSE277782_Run5961_PDAC_Slide{slide}_metadata_file.csv.gz"
    df = pd.read_csv(path)
    # The integrated metadata uses Run5961.PDAC.Slide{slide}_{cell_ID}_{fov}.
    df["global_cell_id"] = (
        "Run5961.PDAC.Slide"
        + str(slide)
        + "_"
        + df["cell_ID"].astype(int).astype(str)
        + "_"
        + df["fov"].astype(int).astype(str)
    )
    return df[
        [
            "global_cell_id",
            "CenterX_global_px",
            "CenterY_global_px",
            "Mean.PanCK",
            "Mean.CD68",
            "Mean.CD45",
            "Area",
        ]
    ]


def median_nn_scale(xy: np.ndarray) -> float:
    if len(xy) < 3:
        return 1.0
    tree = cKDTree(xy)
    dists, _ = tree.query(xy, k=2)
    scale = float(np.nanmedian(dists[:, 1]))
    return scale if np.isfinite(scale) and scale > 0 else 1.0


def nearest_dist_to_anchor(target_xy: np.ndarray, anchor_xy: np.ndarray) -> np.ndarray:
    tree = cKDTree(anchor_xy)
    dists, _ = tree.query(target_xy, k=1)
    return dists


def tissue_group(value: str) -> str:
    value = str(value)
    if value == "Pri":
        return "Primary"
    if value.startswith("LiM"):
        return "Liver metastasis"
    return value


def sample_sort_key(sample: str) -> tuple[int, str]:
    match = re.search(r"Pt-(\d+)", sample)
    pt = int(match.group(1)) if match else 999
    return pt, sample


def main() -> None:
    meta = pd.read_csv(DATA_DIR / "GSE277782_Meta.data_CoxMx.csv.gz").rename(columns={"Unnamed: 0": "global_cell_id"})
    coords = pd.concat([read_slide_metadata(1), read_slide_metadata(2)], ignore_index=True)
    cells = meta.merge(coords, on="global_cell_id", how="inner")
    cells["Tissue_group"] = cells["Tissue"].map(tissue_group)

    composition = (
        cells.groupby(["Sample_name", "Pt", "Tissue", "Tissue_group", "Annotation_main"], as_index=False)
        .agg(n_cells=("global_cell_id", "count"))
    )
    composition.to_csv(TABLE_DIR / "gse277782_cosmx_annotation_composition.csv", index=False)

    rows = []
    target_annotations = ["Tumor cell", "Macrophage", "Plasma cell", "Pancreatic islet", "Hepatocyte"]
    for sample_name, sdf in cells.groupby("Sample_name", sort=False):
        sdf = sdf.reset_index(drop=True)
        xy = sdf[["CenterX_global_px", "CenterY_global_px"]].to_numpy(float)
        scale = median_nn_scale(xy)
        caf_idx = np.flatnonzero(sdf["Annotation_main"].eq("CAF").to_numpy())
        if len(caf_idx) < 10:
            continue
        caf_xy = xy[caf_idx]
        for target in target_annotations:
            target_idx = np.flatnonzero(sdf["Annotation_main"].eq(target).to_numpy())
            if len(target_idx) < 20:
                continue
            target_xy = xy[target_idx]
            observed = float(np.nanmedian(nearest_dist_to_anchor(target_xy, caf_xy) / scale))
            random_medians = []
            for _ in range(N_RANDOM):
                random_idx = RNG.choice(len(sdf), size=len(caf_idx), replace=False)
                random_medians.append(
                    float(np.nanmedian(nearest_dist_to_anchor(target_xy, xy[random_idx]) / scale))
                )
            random_medians = np.array(random_medians)
            random_median = float(np.nanmedian(random_medians))
            empirical_p_closer = float((np.sum(random_medians <= observed) + 1) / (len(random_medians) + 1))
            rows.append(
                {
                    "dataset_id": "GSE277782_CosMx",
                    "sample_name": sample_name,
                    "patient_id": sdf["Pt"].iloc[0],
                    "tissue": sdf["Tissue"].iloc[0],
                    "tissue_group": sdf["Tissue_group"].iloc[0],
                    "target_annotation": target,
                    "n_cells": int(len(sdf)),
                    "n_caf": int(len(caf_idx)),
                    "n_target": int(len(target_idx)),
                    "median_target_to_caf_nn": observed,
                    "random_median_target_to_anchor_nn": random_median,
                    "delta_vs_random": observed - random_median,
                    "target_closer_to_caf_than_random": bool(observed < random_median),
                    "empirical_p_closer": empirical_p_closer,
                }
            )

    proximity = pd.DataFrame(rows)
    proximity.to_csv(TABLE_DIR / "gse277782_cosmx_caf_neighborhood_proximity.csv", index=False)

    context = (
        proximity.groupby(["tissue_group", "target_annotation"], as_index=False)
        .agg(
            n_samples=("sample_name", "nunique"),
            median_delta_vs_random=("delta_vs_random", "median"),
            n_closer=("target_closer_to_caf_than_random", "sum"),
            median_empirical_p=("empirical_p_closer", "median"),
        )
    )
    context.to_csv(TABLE_DIR / "gse277782_cosmx_caf_neighborhood_context_summary.csv", index=False)

    SOURCE_DIR.joinpath("Source_Data_Extended_Data_Fig_11A.csv").write_text(
        composition.to_csv(index=False), encoding="utf-8"
    )
    SOURCE_DIR.joinpath("Source_Data_Extended_Data_Fig_11B_C.csv").write_text(
        context.to_csv(index=False), encoding="utf-8"
    )

    make_figure(composition, context)
    write_report(composition, proximity, context)
    print("Wrote GSE277782 CosMx cell-neighborhood validation outputs.")


def make_figure(composition: pd.DataFrame, context: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    fig = plt.figure(figsize=(14.8, 7.6), constrained_layout=False)
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.35, 1.05], height_ratios=[1, 1], wspace=0.75, hspace=0.62)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[:, 2])
    fig.suptitle("Candidate Extended Data Fig. 11 | GSE277782 CosMx CAF-neighborhood stress test", x=0.5, y=0.98, fontsize=16, fontweight="bold")

    comp = composition.copy()
    sample_order = sorted(comp["Sample_name"].unique(), key=sample_sort_key)
    ann_order = ["CAF", "Macrophage", "Tumor cell", "Plasma cell", "Pancreatic islet", "Hepatocyte", "Undetermined"]
    palette = {
        "CAF": "#5AA05A",
        "Macrophage": "#4C78A8",
        "Tumor cell": "#D55E5E",
        "Plasma cell": "#B07AA1",
        "Pancreatic islet": "#E6A141",
        "Hepatocyte": "#7F7F7F",
        "Undetermined": "#C7C7C7",
    }
    bottom = np.zeros(len(sample_order))
    for ann in ann_order:
        vals = (
            comp[comp["Annotation_main"].eq(ann)]
            .set_index("Sample_name")["n_cells"]
            .reindex(sample_order)
            .fillna(0)
            .to_numpy()
        )
        ax_a.bar(range(len(sample_order)), vals, bottom=bottom, color=palette.get(ann, "#999999"), label=ann, width=0.8)
        bottom += vals
    ax_a.set_title("A  Cell annotations by sample", loc="left", fontweight="bold")
    ax_a.set_ylabel("cells")
    ax_a.set_xticks(range(len(sample_order)))
    ax_a.set_xticklabels([s.replace("Slide_", "S").replace("_", "\n") for s in sample_order], rotation=0, fontsize=8)
    ax_a.legend(frameon=False, fontsize=8, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.16))

    tissue_order = ["Primary", "Liver metastasis"]
    target_order = ["Macrophage", "Tumor cell", "Plasma cell", "Pancreatic islet", "Hepatocyte"]
    matrix = (
        context.pivot(index="target_annotation", columns="tissue_group", values="median_delta_vs_random")
        .reindex(index=target_order, columns=tissue_order)
    )
    labels = (
        context.assign(label=lambda x: x["median_delta_vs_random"].map(lambda v: f"{v:.2f}") + "\n" + x["n_closer"].astype(int).astype(str) + "/" + x["n_samples"].astype(int).astype(str))
        .pivot(index="target_annotation", columns="tissue_group", values="label")
        .reindex(index=target_order, columns=tissue_order)
    )
    im = ax_b.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-1.0, vmax=1.0, aspect="auto")
    ax_b.set_title("B  Target-cell distance to CAF vs random", loc="left", fontweight="bold")
    ax_b.set_xticks(range(len(tissue_order)))
    ax_b.set_xticklabels(tissue_order)
    ax_b.set_yticks(range(len(target_order)))
    ax_b.set_yticklabels(target_order)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax_b.text(j, i, labels.iat[i, j] if pd.notna(labels.iat[i, j]) else "NA", ha="center", va="center", fontsize=8)
    cb = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.06)
    cb.set_label("delta vs random")

    macrophage = context[context["target_annotation"].eq("Macrophage")].set_index("tissue_group").reindex(tissue_order)
    ax_c.bar(range(len(tissue_order)), macrophage["median_delta_vs_random"], color=["#4C78A8", "#72B7B2"], width=0.6)
    ax_c.axhline(0, color="#333333", lw=0.8)
    ax_c.set_title("C  Macrophage-to-CAF proximity", loc="left", fontweight="bold")
    ax_c.set_xticks(range(len(tissue_order)))
    ax_c.set_xticklabels(tissue_order)
    ax_c.set_ylabel("median delta vs random")
    for i, (_, row) in enumerate(macrophage.iterrows()):
        if pd.notna(row["median_delta_vs_random"]):
            ax_c.text(i, row["median_delta_vs_random"], f"{row['n_closer']:.0f}/{row['n_samples']:.0f}", ha="center", va="bottom" if row["median_delta_vs_random"] >= 0 else "top", fontsize=9)

    ax_d.axis("off")
    ax_d.set_title("D  Submission-use boundary", loc="left", fontweight="bold")
    ax_d.text(
        0,
        0.93,
        "Cell-level interpretation:\n\n"
        "GSE277782 provides CosMx cell annotations and\n"
        "global cell coordinates for paired primary and\n"
        "liver-metastatic PDAC samples.\n\n"
        "The analysis asks whether annotated macrophage,\n"
        "tumor and other cell classes are closer to CAFs\n"
        "than to same-size random cell anchors within\n"
        "each sample.\n\n"
        "In the current annotation-level test, most target\n"
        "classes are not closer to CAFs than to random\n"
        "anchors. This should be treated as a stress test\n"
        "or reviewer-response reserve, not promoted as\n"
        "positive validation without a better-defined\n"
        "neighborhood model.",
        va="top",
        ha="left",
        fontsize=10,
        linespacing=1.2,
    )

    out = FIG_DIR / "extended_data_figure11_gse277782_cosmx_cell_neighborhoods"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(composition: pd.DataFrame, proximity: pd.DataFrame, context: pd.DataFrame) -> None:
    lines = [
        "# GSE277782 CosMx Cell-Neighborhood Validation",
        "",
        "Last updated: 2026-06-27",
        "",
        "## Scope",
        "",
        f"Analyzed {composition['Sample_name'].nunique()} CosMx samples and {composition['n_cells'].sum():,} annotated cells from GSE277782.",
        "",
        "## Outputs",
        "",
        "- `results/tables/gse277782_cosmx_annotation_composition.csv`",
        "- `results/tables/gse277782_cosmx_caf_neighborhood_proximity.csv`",
        "- `results/tables/gse277782_cosmx_caf_neighborhood_context_summary.csv`",
        "- `results/source_data/Source_Data_Extended_Data_Fig_11A.csv`",
        "- `results/source_data/Source_Data_Extended_Data_Fig_11B_C.csv`",
        "- `results/figures/internal/extended_data_figure11_gse277782_cosmx_cell_neighborhoods.pdf`",
        "",
        "## Initial Interpretation",
        "",
        "Negative delta values indicate that the target cell class is closer to annotated CAFs than to same-size random cell anchors within the same sample. Positive values indicate that annotated CAFs are not closer than random anchors under this stringent annotation-level nearest-neighbor test.",
    ]
    for _, row in context.sort_values(["tissue_group", "target_annotation"]).iterrows():
        lines.append(
            f"- {row.tissue_group} / {row.target_annotation}: median delta {row.median_delta_vs_random:.3f}, support {int(row.n_closer)}/{int(row.n_samples)} samples, median empirical p {row.median_empirical_p:.3f}."
        )
    lines.extend(
        [
            "",
            "## Claim Boundary",
            "",
            "This analysis should not be promoted as positive validation in its current form. It is useful as a stress test showing that simple annotation-level nearest-neighbor proximity in GSE277782 does not recapitulate the Visium CAF-core organization. A stronger CosMx analysis would need expression-defined neighborhoods, region masks, or manuscript-specific spatial domains before inclusion.",
        ]
    )
    (REPORT_DIR / "gse277782_cosmx_cell_neighborhood_validation.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
