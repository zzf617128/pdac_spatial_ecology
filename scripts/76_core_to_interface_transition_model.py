from __future__ import annotations

from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
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
N_BINS = 6

DATASETS = [
    ("GSE282302", "results/tables/mvp_spot_level_scores_with_edge_qc.csv", "post_neoadjuvant"),
    ("GSE272362", "results/tables/gse272362_rds_spot_level_scores.csv", None),
    ("GSE235315", "results/tables/gse235315_spot_level_scores.csv", "paired_ST_HE"),
    ("GSE274557", "results/tables/gse274557_full_spot_scores.csv", None),
]

PROGRAMS = {
    "CAF/myCAF": ["score_mycaf", "score_pan_caf"],
    "SPP1/TAM": ["score_spp1_tam", "score_myeloid"],
    "TGFb/EMT": ["score_tgfb_pathway", "score_emt_invasion"],
    "IFN/APC": ["score_ifn_antigen_presentation", "score_dc_apc"],
    "T/NK": ["score_t_cell", "score_cd8_effector"],
    "tumor aggressive": ["score_tumor_aggressive"],
    "tumor epithelial": ["score_tumor_epithelial"],
}

PROGRAM_ORDER = ["CAF/myCAF", "SPP1/TAM", "TGFb/EMT", "IFN/APC", "T/NK", "tumor aggressive", "tumor epithelial"]

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
PROGRAM_COLORS = {
    "CAF/myCAF": "#2D6A8E",
    "SPP1/TAM": "#9C4E70",
    "TGFb/EMT": "#B57A3C",
    "IFN/APC": "#628F4E",
    "T/NK": "#4F7EC0",
    "tumor aggressive": "#C44E52",
    "tumor epithelial": "#6F5C99",
}


def zscore(values: pd.Series | np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    mu = np.nanmean(arr)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return np.zeros(len(arr), dtype=float)
    return (arr - mu) / sd


def safe_spearman(x: np.ndarray, y: np.ndarray) -> float:
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 30 or np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(spearmanr(x[mask], y[mask]).statistic)


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
    needed = {
        "dataset_id",
        "sample_id",
        "patient_id",
        "barcode",
        "x_pixel",
        "y_pixel",
        "score_caf_myeloid_barrier",
        "score_tumor_epithelial",
        *[c for cols in PROGRAMS.values() for c in cols],
    }
    for dataset_id, rel_path, default_context in DATASETS:
        path = PROJECT / rel_path
        df = pd.read_csv(path)
        if "edge_or_background_risk" in df.columns:
            df = df.loc[~df["edge_or_background_risk"].fillna(False).astype(bool)].copy()
        keep = [c for c in df.columns if c in needed or c in {"specimen_type", "tissue", "treatment", "title", "geo_accession"}]
        df = df[keep].copy()
        df["analysis_dataset"] = dataset_id
        df["default_context"] = default_context or dataset_id
        frames.append(df)
    spots = pd.concat(frames, ignore_index=True)
    spots = spots[np.isfinite(spots["x_pixel"]) & np.isfinite(spots["y_pixel"])].copy()
    return spots


def add_program_scores(sample: pd.DataFrame) -> pd.DataFrame:
    out = sample.copy()
    for program, cols in PROGRAMS.items():
        z_cols = []
        for col in cols:
            if col in out.columns:
                z_col = f"z__{col}"
                out[z_col] = zscore(pd.to_numeric(out[col], errors="coerce"))
                z_cols.append(z_col)
        out[f"program__{program}"] = out[z_cols].mean(axis=1) if z_cols else np.nan
    return out


def analyze_sample(dataset_id: str, sample_id: str, sample: pd.DataFrame) -> tuple[list[dict], list[dict], pd.DataFrame]:
    sample = sample.dropna(subset=["x_pixel", "y_pixel", "score_caf_myeloid_barrier", "score_tumor_epithelial"]).copy()
    if len(sample) < 120:
        return [], [], pd.DataFrame()
    sample = add_program_scores(sample)
    xy = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    scale = median_nn_scale(xy)

    caf_values = pd.to_numeric(sample["score_caf_myeloid_barrier"], errors="coerce").to_numpy(float)
    tumor_values = pd.to_numeric(sample["score_tumor_epithelial"], errors="coerce").to_numpy(float)
    caf_mask = caf_values >= np.nanquantile(caf_values, 1 - CAF_CORE_FRACTION)
    tumor_mask = tumor_values >= np.nanquantile(tumor_values, 1 - TUMOR_HIGH_FRACTION)
    if caf_mask.sum() < 20 or tumor_mask.sum() < 20:
        return [], [], pd.DataFrame()

    d_caf = nearest_distance(xy, caf_mask) / scale
    d_tumor = nearest_distance(xy, tumor_mask) / scale
    denom = d_caf + d_tumor
    axis = np.divide(d_caf, denom, out=np.full_like(d_caf, 0.5), where=denom > 0)
    axis = np.clip(axis, 0, 1)
    bin_id = np.minimum(np.floor(np.nan_to_num(axis, nan=0.5) * N_BINS).astype(int), N_BINS - 1) + 1

    context = context_label(dataset_id, sample, str(sample["default_context"].iloc[0]))
    meta = {
        "dataset_id": dataset_id,
        "context": context,
        "sample_id": sample_id,
        "n_spots": int(len(sample)),
        "n_caf_core": int(caf_mask.sum()),
        "n_tumor_high": int(tumor_mask.sum()),
        "median_nn_scale": float(scale),
    }

    bin_rows = []
    sample_rows = []
    for program in PROGRAM_ORDER:
        values = sample[f"program__{program}"].to_numpy(float)
        if not np.isfinite(values).any():
            continue
        for b in range(1, N_BINS + 1):
            mask = bin_id == b
            if mask.sum() == 0:
                continue
            bin_rows.append(
                {
                    **meta,
                    "program": program,
                    "axis_bin": b,
                    "axis_midpoint": (b - 0.5) / N_BINS,
                    "mean_z": float(np.nanmean(values[mask])),
                    "median_z": float(np.nanmedian(values[mask])),
                    "n_spots_bin": int(mask.sum()),
                }
            )
        bin_means = pd.DataFrame([r for r in bin_rows if r["program"] == program and r["sample_id"] == sample_id])
        if not bin_means.empty:
            peak_row = bin_means.loc[bin_means["mean_z"].idxmax()]
            early = bin_means.loc[bin_means["axis_bin"].isin([1, 2]), "mean_z"].mean()
            mid = bin_means.loc[bin_means["axis_bin"].isin([3, 4]), "mean_z"].mean()
            late = bin_means.loc[bin_means["axis_bin"].isin([5, 6]), "mean_z"].mean()
            sample_rows.append(
                {
                    **meta,
                    "program": program,
                    "rho_with_core_to_tumor_axis": safe_spearman(axis, values),
                    "peak_bin": int(peak_row["axis_bin"]),
                    "peak_axis_midpoint": float(peak_row["axis_midpoint"]),
                    "core_zone_mean_z": float(early),
                    "interface_zone_mean_z": float(mid),
                    "tumor_zone_mean_z": float(late),
                    "interface_minus_core": float(mid - early),
                    "tumor_minus_core": float(late - early),
                }
            )

    selected = sample[["sample_id", "x_pixel", "y_pixel"]].copy()
    selected["dataset_id"] = dataset_id
    selected["context"] = context
    selected["core_to_tumor_axis"] = axis
    selected["is_caf_core"] = caf_mask
    selected["is_tumor_high"] = tumor_mask
    for program in ["SPP1/TAM", "IFN/APC", "TGFb/EMT", "tumor aggressive"]:
        selected[f"program__{program}"] = sample[f"program__{program}"].to_numpy(float)
    return bin_rows, sample_rows, selected


def choose_map_samples(selected: pd.DataFrame) -> pd.DataFrame:
    keep_contexts = ["primary tumor", "liver metastasis", "lymph-node metastasis"]
    frames = []
    for context in keep_contexts:
        sub = selected[selected["context"].eq(context)]
        if sub.empty:
            continue
        counts = sub.groupby("sample_id").size().sort_values(ascending=False)
        for sample_id in counts.index:
            sample = sub[sub["sample_id"].eq(sample_id)].copy()
            if sample["is_caf_core"].sum() > 20 and sample["is_tumor_high"].sum() > 20:
                if len(sample) > 4500:
                    sample = sample.sample(4500, random_state=20260628)
                frames.append(sample)
                break
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def summarize(bin_df: pd.DataFrame, sample_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    context_summary = (
        bin_df.groupby(["context", "program", "axis_bin", "axis_midpoint"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_mean_z=("mean_z", "median"),
            q25_mean_z=("mean_z", lambda x: float(np.nanpercentile(x, 25))),
            q75_mean_z=("mean_z", lambda x: float(np.nanpercentile(x, 75))),
        )
    )
    sample_summary = (
        sample_df.groupby(["context", "program"], as_index=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_peak_axis=("peak_axis_midpoint", "median"),
            median_rho=("rho_with_core_to_tumor_axis", "median"),
            median_interface_minus_core=("interface_minus_core", "median"),
            median_tumor_minus_core=("tumor_minus_core", "median"),
        )
    )
    wide = sample_summary.pivot_table(index="program", columns="context", values="median_rho", aggfunc="median").reset_index()
    return context_summary, sample_summary, wide


def draw_schematic(ax: plt.Axes) -> None:
    ax.set_title("A  Core-to-interface transition model", loc="left", fontweight="bold")
    ax.axis("off")
    caf = Circle((0.18, 0.52), 0.13, color="#2D6A8E", alpha=0.9)
    interface = Circle((0.50, 0.52), 0.12, color="#B57A3C", alpha=0.75)
    tumor = Circle((0.82, 0.52), 0.14, color="#6F5C99", alpha=0.9)
    for patch in [caf, interface, tumor]:
        ax.add_patch(patch)
    ax.text(0.18, 0.52, "CAF\ncore", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    ax.text(0.50, 0.52, "interface\nzone", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    ax.text(0.82, 0.52, "tumor\nhigh", ha="center", va="center", color="white", fontsize=10, fontweight="bold")
    ax.add_patch(FancyArrowPatch((0.30, 0.52), (0.70, 0.52), arrowstyle="->", mutation_scale=18, lw=1.6, color="#333333"))
    ax.text(0.50, 0.75, "core-to-tumor coordinate", ha="center", fontsize=10)
    ax.text(0.18, 0.25, "0", ha="center", fontsize=9)
    ax.text(0.50, 0.25, "0.5", ha="center", fontsize=9)
    ax.text(0.82, 0.25, "1", ha="center", fontsize=9)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)


def make_figure(context_summary: pd.DataFrame, sample_summary: pd.DataFrame, selected: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.8,
        }
    )
    fig = plt.figure(figsize=(16.8, 11.2))
    gs = fig.add_gridspec(
        3,
        4,
        height_ratios=[0.84, 1.08, 1.08],
        width_ratios=[1.0, 1.0, 1.0, 0.88],
        wspace=0.46,
        hspace=0.62,
    )
    fig.suptitle("Core-to-interface transition trajectories deepen the CAF-domain model", fontsize=17, fontweight="bold")

    ax_a = fig.add_subplot(gs[0, 0])
    draw_schematic(ax_a)

    map_contexts = ["primary tumor", "liver metastasis", "lymph-node metastasis"]
    for i, context in enumerate(map_contexts):
        ax = fig.add_subplot(gs[0, i + 1])
        sub = selected[selected["context"].eq(context)]
        if sub.empty:
            ax.axis("off")
            continue
        im = ax.scatter(sub["x_pixel"], sub["y_pixel"], c=sub["core_to_tumor_axis"], s=5, cmap="viridis", vmin=0, vmax=1, linewidth=0, alpha=0.9)
        core = sub[sub["is_caf_core"]]
        tumor = sub[sub["is_tumor_high"]]
        ax.scatter(core["x_pixel"], core["y_pixel"], s=4, color="#2D6A8E", alpha=0.55, linewidth=0)
        ax.scatter(tumor["x_pixel"], tumor["y_pixel"], s=4, color="#6F5C99", alpha=0.45, linewidth=0)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(("B" if i == 0 else "") + f"  {context}", loc="left", fontweight="bold" if i == 0 else "normal")
        if i == 2:
            cb = fig.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
            cb.set_label("core-to-tumor axis", fontsize=8)

    plot_programs = ["SPP1/TAM", "TGFb/EMT", "IFN/APC", "T/NK", "tumor aggressive", "tumor epithelial"]
    plot_contexts = ["primary tumor", "liver metastasis", "lymph-node metastasis"]
    for idx, program in enumerate(plot_programs):
        row = 1 + idx // 3
        col = idx % 3
        ax = fig.add_subplot(gs[row, col])
        for context in plot_contexts:
            sub = context_summary[context_summary["context"].eq(context) & context_summary["program"].eq(program)].sort_values("axis_bin")
            if sub.empty:
                continue
            x = sub["axis_midpoint"].to_numpy(float)
            y = sub["median_mean_z"].to_numpy(float)
            lo = sub["q25_mean_z"].to_numpy(float)
            hi = sub["q75_mean_z"].to_numpy(float)
            ax.plot(x, y, marker="o", lw=2, color=CONTEXT_COLORS[context], label=context)
            ax.fill_between(x, lo, hi, color=CONTEXT_COLORS[context], alpha=0.13, linewidth=0)
        ax.axhline(0, color="#333333", lw=0.8)
        ax.axvspan(0.34, 0.66, color="#F0E5D5", alpha=0.45, zorder=0)
        ax.set_title(("C  " if idx == 0 else "") + program, loc="left", fontweight="bold" if idx == 0 else "normal")
        ax.set_xlabel("CAF-core to tumor-high coordinate")
        ax.set_ylabel("within-sample z")
        ax.set_xlim(0, 1)
        ax.grid(color="#E5E5E5", lw=0.6)
        ax.spines[["top", "right"]].set_visible(False)
        if idx == 0:
            ax.legend(frameon=False, fontsize=8)

    ax_g = fig.add_subplot(gs[1:, 3])
    peak_contexts = ["primary tumor", "liver metastasis", "lymph-node metastasis"]
    peak_programs = ["SPP1/TAM", "TGFb/EMT", "IFN/APC", "T/NK", "tumor aggressive", "tumor epithelial"]
    y = np.arange(len(peak_programs))
    offsets = {"primary tumor": -0.18, "liver metastasis": 0.0, "lymph-node metastasis": 0.18}
    for context in peak_contexts:
        sub = sample_summary[sample_summary["context"].eq(context)].set_index("program")
        vals = sub.reindex(peak_programs)["median_peak_axis"].to_numpy(float)
        ax_g.scatter(vals, y + offsets[context], s=42, color=CONTEXT_COLORS[context], label=context, zorder=3)
        for i, val in enumerate(vals):
            if np.isfinite(val):
                ax_g.plot([0.5, val], [i + offsets[context], i + offsets[context]], color="#BBBBBB", lw=0.7, zorder=0)
    ax_g.axvline(0.5, color="#333333", lw=0.8)
    ax_g.set_yticks(y, peak_programs, fontsize=8)
    ax_g.set_xlabel("median peak coordinate")
    ax_g.set_xlim(0, 1)
    ax_g.invert_yaxis()
    ax_g.set_title("D  Program peak along axis", loc="left", fontweight="bold")
    ax_g.legend(frameon=False, fontsize=7, loc="lower right")
    ax_g.grid(axis="x", color="#E5E5E5", lw=0.6)
    ax_g.spines[["top", "right"]].set_visible(False)

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / "extended_data_figure32_core_to_interface_transition_model"
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(context_summary: pd.DataFrame, sample_summary: pd.DataFrame) -> None:
    focus = sample_summary[sample_summary["context"].isin(["primary tumor", "liver metastasis", "lymph-node metastasis"])]
    lines = [
        "# Core-To-Interface Transition Model",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Question",
        "",
        "This analysis asks whether candidate programs follow reproducible trajectories along a pseudo-spatial axis from CAF-myeloid cores to tumor-high regions, rather than appearing only as static proximity enrichments.",
        "",
        "## Method",
        "",
        f"Within each section, CAF-core spots were defined as the top {int(CAF_CORE_FRACTION * 100)}% of CAF-myeloid-barrier score and tumor-high spots as the top {int(TUMOR_HIGH_FRACTION * 100)}% of tumor-epithelial score. Each spot received a coordinate d(CAF)/(d(CAF)+d(tumor)), where 0 indicates CAF-core proximity and 1 indicates tumor-high proximity. Program scores were within-sample z-scored and summarized across {N_BINS} axis bins.",
        "",
        "## Key Trajectory Summaries",
        "",
    ]
    for program in PROGRAM_ORDER:
        sub = focus[focus["program"].eq(program)].sort_values("context")
        if sub.empty:
            continue
        parts = []
        for _, row in sub.iterrows():
            parts.append(
                f"{row.context}: rho {row.median_rho:.2f}, peak {row.median_peak_axis:.2f}, interface-core {row.median_interface_minus_core:.2f}"
            )
        lines.append(f"- {program}: " + "; ".join(parts) + ".")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The transition model should be described as a spatial trajectory or gradient analysis, not a temporal trajectory and not a causal process. It is useful because it separates CAF-core-proximal, interface-like and tumor-high-proximal behavior across tissue contexts.",
            "",
            "## Outputs",
            "",
            "- `results/tables/core_to_interface_transition_bin_summary.csv`",
            "- `results/tables/core_to_interface_transition_sample_summary.csv`",
            "- `results/tables/core_to_interface_transition_context_summary.csv`",
            "- `results/source_data/Source_Data_Extended_Data_Fig_32_core_to_interface_transition.csv`",
            "- `results/figures/submission/extended_data_figure32_core_to_interface_transition_model.pdf`",
        ]
    )
    (REPORT_DIR / "core_to_interface_transition_model_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    spots = load_spots()
    all_bin_rows: list[dict] = []
    all_sample_rows: list[dict] = []
    selected_frames = []
    for (dataset_id, sample_id), sample in spots.groupby(["analysis_dataset", "sample_id"], sort=False):
        bin_rows, sample_rows, selected = analyze_sample(str(dataset_id), str(sample_id), sample)
        all_bin_rows.extend(bin_rows)
        all_sample_rows.extend(sample_rows)
        if not selected.empty:
            selected_frames.append(selected)
        print(f"{dataset_id} {sample_id}: bin_rows={len(bin_rows)} sample_rows={len(sample_rows)}")
    bin_df = pd.DataFrame(all_bin_rows)
    sample_df = pd.DataFrame(all_sample_rows)
    selected_all = pd.concat(selected_frames, ignore_index=True) if selected_frames else pd.DataFrame()
    selected_maps = choose_map_samples(selected_all)
    context_summary, sample_summary, wide = summarize(bin_df, sample_df)

    bin_df.to_csv(TABLE_DIR / "core_to_interface_transition_bin_summary.csv", index=False)
    sample_df.to_csv(TABLE_DIR / "core_to_interface_transition_sample_summary.csv", index=False)
    context_summary.to_csv(TABLE_DIR / "core_to_interface_transition_context_summary.csv", index=False)
    wide.to_csv(TABLE_DIR / "core_to_interface_transition_rho_matrix.csv", index=False)
    selected_maps.to_csv(TABLE_DIR / "core_to_interface_transition_selected_spot_maps.csv", index=False)

    context_summary.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_32_core_to_interface_transition.csv", index=False)
    sample_summary.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_32_sample_summary.csv", index=False)
    selected_maps.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_32_selected_spot_maps.csv", index=False)

    make_figure(context_summary, sample_summary, selected_maps)
    write_report(context_summary, sample_summary)
    print("Wrote core-to-interface transition model outputs.")


if __name__ == "__main__":
    main()
