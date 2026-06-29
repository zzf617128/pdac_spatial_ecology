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

OUT = FIG_DIR / "extended_data_figure26_cell_state_reference_xenium_module_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Extended_Data_Fig_26_cell_state_reference_xenium_module.csv"
REPORT_OUT = REPORT_DIR / "extended_data_figure26_cell_state_reference_xenium_module_notes.md"

CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
]
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treat-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
}
MARKER_STATES = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK", "B/plasma cell", "epithelial/tumor"]
REF_STATES = ["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "B_plasma", "epithelial_tumor"]
STATE_LABELS = {
    "myCAF/matrix": "myCAF/matrix",
    "SPP1/TAM": "SPP1/TAM",
    "DC/APC": "DC/APC",
    "T/NK": "T/NK",
    "B/plasma cell": "B/plasma",
    "epithelial/tumor": "epithelial",
    "myCAF_matrix": "myCAF/matrix",
    "SPP1_TAM": "SPP1/TAM",
    "DC_APC": "DC/APC",
    "T_NK": "T/NK",
    "B_plasma": "B/plasma",
    "epithelial_tumor": "epithelial",
    "IFN_APC": "IFN/APC",
    "TGFb_EMT": "TGFb/EMT",
    "SPP1_tumor_like": "SPP1 tumor",
    "Tumor_epithelial": "tumor epi.",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E8E8E8", linewidth=0.55, zorder=0)


def heatmap(ax: plt.Axes, mat: pd.DataFrame, title: str, panel: str, cmap: str, vmin: float, vmax: float, cbar_label: str) -> None:
    im = ax.imshow(mat.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels([CONTEXT_LABELS.get(c, c) for c in mat.columns], rotation=35, ha="right", fontsize=6.7)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels([STATE_LABELS.get(x, x) for x in mat.index], fontsize=7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, panel)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if not np.isfinite(val):
                continue
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > (0.6 * max(abs(vmin), abs(vmax))) else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
    cbar.set_label(cbar_label, fontsize=6)
    cbar.ax.tick_params(labelsize=6)


def plot_marker_core(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(MARKER_STATES)].copy()
    mat = (
        df.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
        .reindex(index=MARKER_STATES, columns=CONTEXT_ORDER)
    )
    heatmap(ax, mat, "Marker-state CAF-core enrichment", "A", "YlGnBu", -0.2, 0.85, "core - noncore")
    return df


def plot_marker_spearman(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(MARKER_STATES)].copy()
    mat = (
        df.pivot_table(index="cell_state", columns="cohort_context", values="median_spot_spearman_with_caf_myeloid", aggfunc="median")
        .reindex(index=MARKER_STATES, columns=CONTEXT_ORDER)
    )
    heatmap(ax, mat, "Marker-state correlation with CAF-myeloid score", "B", "RdBu_r", -0.35, 0.75, "median Spearman rho")
    return df


def plot_decoupling_correlations(ax: plt.Axes) -> pd.DataFrame:
    marker = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_correlations.csv")
    marker = marker[marker["target"].eq("immune_decoupling_index") & marker["cell_state"].isin(MARKER_STATES)].copy()
    marker["source"] = "marker"
    marker["state_label"] = marker["cell_state"].map(STATE_LABELS)

    ref = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_correlations.csv")
    ref = ref[ref["target"].eq("immune_decoupling_index") & ref["cell_state"].isin(REF_STATES)].copy()
    ref["source"] = "full GSE202051 projection"
    ref["state_label"] = ref["cell_state"].map(STATE_LABELS)

    keep = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK", "B/plasma", "epithelial"]
    combined = pd.concat(
        [
            marker[["source", "state_label", "spearman_rho", "p_value", "n_samples"]],
            ref[["source", "state_label", "spearman_rho", "p_value", "n_samples"]],
        ],
        ignore_index=True,
    )
    combined = combined[combined["state_label"].isin(keep)].copy()
    combined["state_label"] = pd.Categorical(combined["state_label"], categories=keep, ordered=True)
    combined = combined.sort_values(["state_label", "source"])

    y_base = np.arange(len(keep))
    offsets = {"marker": -0.16, "full GSE202051 projection": 0.16}
    colors = {"marker": "#4C78A8", "full GSE202051 projection": "#B23A48"}
    ax.axvline(0, color="#333333", lw=0.65)
    for source in ["marker", "full GSE202051 projection"]:
        sub = combined[combined["source"].eq(source)]
        y = np.array([keep.index(x) for x in sub["state_label"].astype(str)]) + offsets[source]
        ax.scatter(sub["spearman_rho"], y, s=38, color=colors[source], edgecolor="white", linewidth=0.5, label=source, zorder=3)
        ax.hlines(y, 0, sub["spearman_rho"], color=colors[source], lw=1.4)
    ax.set_yticks(y_base)
    ax.set_yticklabels(keep, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(-0.82, 0.82)
    ax.set_xlabel("rho with immune-decoupling index", fontsize=7)
    ax.set_title("Cell-state links to immune decoupling", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "C")
    clean_axes(ax, axis="x")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    return combined


def plot_full_reference(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(REF_STATES)].copy()
    mat = (
        df.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
        .reindex(index=REF_STATES, columns=CONTEXT_ORDER)
    )
    heatmap(ax, mat, "Full GSE202051 projection enrichment", "D", "RdBu_r", -0.08, 0.08, "projection core - noncore")
    return df


def plot_reference_stability(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_reference_projection_small_vs_full_comparison.csv")
    df = df[df["cell_state"].isin(REF_STATES)].copy()
    df["state_label"] = df["cell_state"].map(STATE_LABELS)
    order = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK", "B/plasma", "epithelial"]
    df = df[df["state_label"].isin(order)].copy()
    df["state_label"] = pd.Categorical(df["state_label"], categories=order, ordered=True)
    df = df.sort_values("state_label")
    y = np.arange(len(df))
    colors = np.where(df["same_direction"], "#2C7A51", "#8A8F98")
    ax.barh(y, df["spearman_rho"], color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(df["state_label"].astype(str), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("small vs full per-sample rho", fontsize=7)
    ax.set_title("Projection stability across reference size", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax)
    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(row["spearman_rho"] + 0.025, i, f"{row['same_direction_fraction']:.2f}", va="center", fontsize=6.5)
    return df


def plot_xenium_anchor(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_context_summary.csv")
    targets = ["IFN_APC", "SPP1_TAM", "TGFb_EMT", "T_NK", "SPP1_tumor_like", "Tumor_epithelial"]
    df = df[df["target_program"].isin(targets)].copy()
    df["target_label"] = df["target_program"].map(STATE_LABELS).fillna(df["target_program"])
    mat = df.pivot_table(index="target_label", columns="anchor", values="median_delta_vs_random", aggfunc="median")
    order = ["IFN/APC", "SPP1/TAM", "TGFb/EMT", "T/NK", "SPP1 tumor", "tumor epi."]
    mat = mat.reindex([x for x in order if x in mat.index])
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.5, vmax=0.5, aspect="auto")
    ax.set_xticks(np.arange(len(mat.columns)))
    ax.set_xticklabels([str(c).replace("_", "\n") for c in mat.columns], fontsize=7)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_title("Xenium CAF-domain target centering", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "F")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=6, color="#FFFFFF" if abs(val) > 0.28 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label("delta vs random", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_xenium_scale_coverage(ax: plt.Axes) -> tuple[pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")
    cov = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_signature_coverage.csv")
    comp = comp.sort_values("n_cells", ascending=True).copy()
    y = np.arange(len(comp))
    colors = np.where(comp["treatment"].str.contains("chemo", case=False, na=False), "#7B68A6", "#4C78A8")
    ax.barh(y, comp["n_cells"] / 1000, color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels([x.replace("Patient ", "P").replace(" PDAC_", "-") for x in comp["title"]], fontsize=6.6)
    ax.set_xlabel("Xenium cells (thousand)", fontsize=7)
    ax.set_title("Xenium scale and panel coverage", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "G")
    clean_axes(ax)

    ax2 = ax.inset_axes([0.57, 0.16, 0.39, 0.58])
    sig_order = ["CAF_matrix", "SPP1_TAM", "IFN_APC", "T_NK", "TGFb_EMT", "Tumor_epithelial"]
    cov = cov[cov["signature"].isin(sig_order)].copy()
    cov["coverage"] = cov["n_present"] / cov["n_genes"]
    cov_summary = cov.groupby("signature", as_index=False).agg(coverage=("coverage", "median"))
    cov_summary["signature"] = pd.Categorical(cov_summary["signature"], categories=sig_order, ordered=True)
    cov_summary = cov_summary.sort_values("signature")
    yy = np.arange(len(cov_summary))
    ax2.barh(yy, cov_summary["coverage"], color="#8A8F98", height=0.55)
    ax2.set_xlim(0, 1)
    ax2.set_yticks(yy)
    ax2.set_yticklabels([STATE_LABELS.get(x, x.replace("_", "/")) for x in cov_summary["signature"].astype(str)], fontsize=5.8)
    ax2.set_xticks([0, 0.5, 1.0])
    ax2.tick_params(axis="x", labelsize=5.5)
    ax2.set_title("gene coverage", fontsize=5.8, loc="left", pad=1.5)
    ax2.spines[["top", "right"]].set_visible(False)
    return comp, cov_summary


def write_source_data(
    marker_core: pd.DataFrame,
    marker_spearman: pd.DataFrame,
    decoupling: pd.DataFrame,
    ref_context: pd.DataFrame,
    stability: pd.DataFrame,
    xenium_anchor: pd.DataFrame,
    xenium_comp: pd.DataFrame,
    xenium_cov: pd.DataFrame,
) -> None:
    rows: list[dict[str, object]] = []
    for _, row in marker_core.iterrows():
        rows.append({"panel": "A", "source": "gap2_cell_state_marker_attribution_context_summary.csv", "item": f"{row['cohort_context']}|{row['cell_state']}", "metric": "median_core_enrichment", "value": row["median_core_enrichment"]})
    for _, row in marker_spearman.iterrows():
        rows.append({"panel": "B", "source": "gap2_cell_state_marker_attribution_context_summary.csv", "item": f"{row['cohort_context']}|{row['cell_state']}", "metric": "median_spot_spearman_with_caf_myeloid", "value": row["median_spot_spearman_with_caf_myeloid"]})
    for _, row in decoupling.iterrows():
        rows.append({"panel": "C", "source": row["source"], "item": row["state_label"], "metric": "spearman_rho_with_immune_decoupling", "value": row["spearman_rho"]})
    for _, row in ref_context.iterrows():
        rows.append({"panel": "D", "source": "gap2_full_reference_projection_deconvolution_context_summary.csv", "item": f"{row['cohort_context']}|{row['cell_state']}", "metric": "median_core_enrichment", "value": row["median_core_enrichment"]})
    for _, row in stability.iterrows():
        rows.append({"panel": "E", "source": "gap2_reference_projection_small_vs_full_comparison.csv", "item": row["state_label"], "metric": "small_vs_full_spearman_rho", "value": row["spearman_rho"]})
    for _, row in xenium_anchor.iterrows():
        rows.append({"panel": "F", "source": "gse274673_xenium_fixed_anchor_context_summary.csv", "item": f"{row['anchor']}|{row['target_program']}", "metric": "median_delta_vs_random", "value": row["median_delta_vs_random"]})
    for _, row in xenium_comp.iterrows():
        rows.append({"panel": "G", "source": "gse274673_xenium_fixed_anchor_sample_composition.csv", "item": row["geo_accession"], "metric": "n_cells", "value": row["n_cells"]})
    for _, row in xenium_cov.iterrows():
        rows.append({"panel": "G", "source": "gse274673_xenium_fixed_anchor_signature_coverage.csv", "item": row["signature"], "metric": "median_gene_coverage", "value": row["coverage"]})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Extended Data Figure 26 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "NC-style cell-state triangulation module. The figure consolidates marker-level spatial attribution, GSE202051 reference-projection support and GSE274673 Xenium cell-resolution support.\n\n"
        "## Panel contract\n\n"
        "- A-B: marker-state CAF-core enrichment and spot-level coupling to CAF-myeloid score.\n"
        "- C: marker and full-reference links between CAF-core cell-state enrichment and immune decoupling.\n"
        "- D: full GSE202051 reference-projection CAF-core enrichment across contexts.\n"
        "- E: small-reference versus full-reference projection stability.\n"
        "- F: Xenium CAF-domain target-program centering around CAF-APC and CAF-SPP1/TAM anchors.\n"
        "- G: Xenium sample scale and signature gene coverage.\n\n"
        "## Boundary\n\n"
        "This module supports cell-state interpretation but is not formal image segmentation, immunostaining, final validated deconvolution or causal cell-cell interaction inference.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf')}`\n"
        f"- `{OUT.with_suffix('.svg')}`\n"
        f"- `{OUT.with_suffix('.png')}`\n"
        f"- `{SOURCE_OUT}`\n",
        encoding="utf-8",
    )


def main() -> None:
    fig = plt.figure(figsize=(14.8, 13.6), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[1.05, 1.0, 1.0, 0.72], hspace=0.92, wspace=0.88)
    fig.suptitle("Cell-state, reference-projection and Xenium support module", fontsize=15, fontweight="bold", y=0.986)

    ax_a = fig.add_subplot(gs[0, 0:3])
    ax_b = fig.add_subplot(gs[0, 3:6])
    ax_c = fig.add_subplot(gs[1, 0:2])
    ax_d = fig.add_subplot(gs[1, 2:5])
    ax_e = fig.add_subplot(gs[1, 5:6])
    ax_f = fig.add_subplot(gs[2, 0:3])
    ax_g = fig.add_subplot(gs[2, 3:6])
    ax_note = fig.add_subplot(gs[3, :])

    marker_core = plot_marker_core(ax_a)
    marker_spearman = plot_marker_spearman(ax_b)
    decoupling = plot_decoupling_correlations(ax_c)
    ref_context = plot_full_reference(ax_d)
    stability = plot_reference_stability(ax_e)
    xenium_anchor = plot_xenium_anchor(ax_f)
    xenium_comp, xenium_cov = plot_xenium_scale_coverage(ax_g)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.78,
        "Interpretation: CAF-myeloid cores are supported by recognizable myCAF/matrix and SPP1/TAM marker states, "
        "with attenuated immune-state coupling in immune-decoupled contexts. Full GSE202051 reference projection preserves the same direction, "
        "and GSE274673 Xenium data show cell-resolution CAF-domain centering of immune/myeloid programs.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.32,
        "Boundary: these analyses support cell-state interpretation of expression-defined spatial programs; they are not final deconvolution, immunostaining, image segmentation or causal interaction evidence.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)

    write_source_data(marker_core, marker_spearman, decoupling, ref_context, stability, xenium_anchor, xenium_comp, xenium_cov)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
