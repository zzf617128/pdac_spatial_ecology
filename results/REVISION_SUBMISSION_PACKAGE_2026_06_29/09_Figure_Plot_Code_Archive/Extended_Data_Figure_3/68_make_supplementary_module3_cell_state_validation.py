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
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
REPORT_DIR = ROOT / "results" / "reports"

OUT = FIG_DIR / "supplementary_module3_cell_state_multiresolution_validation"
SOURCE_OUT = SOURCE_DIR / "Source_Data_supplementary_module3_cell_state_multiresolution_validation.csv"
REPORT_OUT = REPORT_DIR / "supplementary_module3_cell_state_multiresolution_validation_notes.md"

for d in (SOURCE_DIR, FIG_DIR, REPORT_DIR):
    d.mkdir(parents=True, exist_ok=True)

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
NNLS_STATES = ["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "B_plasma", "epithelial_tumor"]
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
    "CAF_matrix": "CAF/matrix",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E8E8E8", linewidth=0.55, zorder=0)


def heatmap(
    ax: plt.Axes,
    mat: pd.DataFrame,
    title: str,
    panel: str,
    cmap: str,
    vmin: float,
    vmax: float,
    cbar_label: str,
    context_axis: bool = True,
) -> None:
    im = ax.imshow(mat.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(mat.columns)))
    if context_axis:
        xlabels = [CONTEXT_LABELS.get(c, c) for c in mat.columns]
    else:
        xlabels = [str(c).replace("_", "\n") for c in mat.columns]
    ax.set_xticklabels(xlabels, rotation=35 if len(xlabels) > 2 else 0, ha="right" if len(xlabels) > 2 else "center", fontsize=6.7)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels([STATE_LABELS.get(x, x) for x in mat.index], fontsize=7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, panel)
    threshold = 0.58 * max(abs(vmin), abs(vmax))
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > threshold else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.036, pad=0.02)
    cbar.set_label(cbar_label, fontsize=6)
    cbar.ax.tick_params(labelsize=6)


def marker_enrichment(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(MARKER_STATES)].copy()
    mat = df.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
    mat = mat.reindex(index=MARKER_STATES, columns=CONTEXT_ORDER)
    heatmap(ax, mat, "Marker-state CAF-core enrichment", "A", "YlGnBu", -0.2, 0.85, "core - noncore")
    return df


def reference_enrichment(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(REF_STATES)].copy()
    mat = df.pivot_table(index="cell_state", columns="cohort_context", values="median_core_enrichment", aggfunc="median")
    mat = mat.reindex(index=REF_STATES, columns=CONTEXT_ORDER)
    heatmap(ax, mat, "Full GSE202051 projection", "B", "RdBu_r", -0.08, 0.08, "projection core - noncore")
    return df


def nnls_enrichment(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "strict_nnls_reference_deconvolution_context_summary.csv")
    df = df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(NNLS_STATES)].copy()
    mat = df.pivot_table(index="cell_state", columns="cohort_context", values="median_nnls_core_enrichment", aggfunc="median")
    mat = mat.reindex(index=NNLS_STATES, columns=CONTEXT_ORDER)
    heatmap(ax, mat, "Strict NNLS CAF-core enrichment", "C", "RdBu_r", -0.08, 0.08, "NNLS core - noncore")
    return df


def decoupling_correlations(ax: plt.Axes) -> pd.DataFrame:
    marker = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_correlations.csv")
    marker = marker[marker["target"].eq("immune_decoupling_index") & marker["cell_state"].isin(MARKER_STATES)].copy()
    marker["layer"] = "marker"
    marker["state_label"] = marker["cell_state"].map(STATE_LABELS)

    ref = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_correlations.csv")
    ref = ref[ref["target"].eq("immune_decoupling_index") & ref["cell_state"].isin(REF_STATES)].copy()
    ref["layer"] = "projection"
    ref["state_label"] = ref["cell_state"].map(STATE_LABELS)

    nnls = pd.read_csv(TABLE_DIR / "strict_nnls_reference_deconvolution_correlations.csv")
    nnls = nnls[nnls["target"].eq("immune_decoupling_index") & nnls["cell_state"].isin(NNLS_STATES)].copy()
    nnls["layer"] = "NNLS"
    nnls["state_label"] = nnls["cell_state"].map(STATE_LABELS)

    combined = pd.concat(
        [
            marker[["layer", "state_label", "spearman_rho", "p_value", "n_samples"]],
            ref[["layer", "state_label", "spearman_rho", "p_value", "n_samples"]],
            nnls[["layer", "state_label", "spearman_rho", "p_value", "n_samples"]],
        ],
        ignore_index=True,
    )
    order = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK", "B/plasma", "epithelial"]
    combined = combined[combined["state_label"].isin(order)].copy()
    y_base = np.arange(len(order))
    offsets = {"marker": -0.22, "projection": 0.0, "NNLS": 0.22}
    colors = {"marker": "#4C78A8", "projection": "#B23A48", "NNLS": "#2C7A51"}
    ax.axvline(0, color="#333333", lw=0.65)
    for layer in ["marker", "projection", "NNLS"]:
        sub = combined[combined["layer"].eq(layer)]
        y = np.array([order.index(x) for x in sub["state_label"]]) + offsets[layer]
        ax.scatter(sub["spearman_rho"], y, s=34, color=colors[layer], edgecolor="white", linewidth=0.45, label=layer, zorder=3)
        ax.hlines(y, 0, sub["spearman_rho"], color=colors[layer], lw=1.2, alpha=0.85)
    ax.set_yticks(y_base)
    ax.set_yticklabels(order, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(-0.85, 0.85)
    ax.set_xlabel("rho with immune-decoupling index", fontsize=7)
    ax.set_title("Cell-state links to immune decoupling", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "D")
    clean_axes(ax, axis="x")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    return combined


def nnls_projection_agreement(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "strict_nnls_vs_projection_comparison.csv")
    df = df[df["cell_state"].isin(NNLS_STATES)].copy()
    df["state_label"] = df["cell_state"].map(STATE_LABELS)
    order = ["myCAF/matrix", "SPP1/TAM", "DC/APC", "T/NK", "B/plasma", "epithelial"]
    df["state_label"] = pd.Categorical(df["state_label"], categories=order, ordered=True)
    df = df.sort_values("state_label")
    y = np.arange(len(df))
    ax.barh(y, df["nnls_vs_projection_spearman_rho"], color="#2C7A51", height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(df["state_label"].astype(str), fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.05)
    ax.set_xlabel("NNLS vs projection rho", fontsize=7)
    ax.set_title("Strict NNLS agrees with projection", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax)
    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(row["nnls_vs_projection_spearman_rho"] + 0.025, i, f"{row['same_direction_fraction']:.2f}", va="center", fontsize=6.5)
    return df


def xenium_anchor(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_context_summary.csv")
    targets = ["IFN_APC", "SPP1_TAM", "TGFb_EMT", "T_NK", "SPP1_tumor_like", "Tumor_epithelial"]
    df = df[df["target_program"].isin(targets)].copy()
    df["target_label"] = df["target_program"].map(STATE_LABELS).fillna(df["target_program"])
    mat = df.pivot_table(index="target_label", columns="anchor", values="median_delta_vs_random", aggfunc="median")
    order = ["IFN/APC", "SPP1/TAM", "TGFb/EMT", "T/NK", "SPP1 tumor", "tumor epi."]
    mat = mat.reindex([x for x in order if x in mat.index])
    heatmap(ax, mat, "Xenium CAF-domain target centering", "F", "RdBu_r", -0.5, 0.5, "delta vs random", context_axis=False)
    return df


def xenium_scale_coverage(ax: plt.Axes) -> tuple[pd.DataFrame, pd.DataFrame]:
    comp = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_sample_composition.csv")
    cov = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_signature_coverage.csv")
    comp = comp.sort_values("n_cells", ascending=True).copy()
    y = np.arange(len(comp))
    colors = np.where(comp["treatment"].str.contains("chemo", case=False, na=False), "#7B68A6", "#4C78A8")
    ax.barh(y, comp["n_cells"] / 1000, color=colors, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels([x.replace("Patient ", "P").replace(" PDAC_", "-") for x in comp["title"]], fontsize=6.5)
    ax.set_xlabel("cells (thousand)", fontsize=7)
    ax.set_title("Xenium scale and gene coverage", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "G")
    clean_axes(ax)

    ax2 = ax.inset_axes([0.55, 0.14, 0.40, 0.62])
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
    ax2.set_yticklabels([STATE_LABELS.get(x, x.replace("_", "/")) for x in cov_summary["signature"].astype(str)], fontsize=5.7)
    ax2.set_xticks([0, 0.5, 1.0])
    ax2.tick_params(axis="x", labelsize=5.5)
    ax2.set_title("median gene coverage", fontsize=5.8, loc="left", pad=1.5)
    ax2.spines[["top", "right"]].set_visible(False)
    return comp, cov_summary


def claim_boundary(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "H", x=-0.04, y=1.02)
    ax.text(0.02, 0.98, "Interpretation boundary", fontsize=10, fontweight="bold", va="top")
    rows = [
        ("myCAF/matrix + SPP1/TAM", "supported across marker/projection/NNLS", "#2C7A51"),
        ("Immune-state attenuation", "tracks immune decoupling", "#4C78A8"),
        ("Xenium CAF-domain layer", "supports immune/myeloid centering", "#7B68A6"),
        ("Single-cell abundance", "not claimed", "#8A8F98"),
        ("Causal interaction", "not claimed", "#8A8F98"),
    ]
    for i, (claim, scope, color) in enumerate(rows):
        y = 0.82 - i * 0.145
        ax.add_patch(Rectangle((0.02, y - 0.055), 0.96, 0.105, facecolor="#F6F7F9", edgecolor="#DDDDDD", linewidth=0.5))
        ax.add_patch(Rectangle((0.02, y - 0.055), 0.025, 0.105, facecolor=color, edgecolor=color, linewidth=0))
        ax.text(0.06, y + 0.018, claim, va="center", fontsize=6.8)
        ax.text(0.06, y - 0.024, scope, va="center", fontsize=6.2, color=color, fontweight="bold")
    return pd.DataFrame(rows, columns=["claim", "scope", "color"])


def write_source_data(*dfs: pd.DataFrame) -> None:
    panel_names = ["A_marker", "B_projection", "C_nnls", "D_decoupling", "E_agreement", "F_xenium", "G_xenium_scale", "G_xenium_coverage", "H_boundary"]
    rows: list[dict[str, object]] = []
    for panel, df in zip(panel_names, dfs):
        for _, row in df.iterrows():
            record = {"panel": panel}
            for col in df.columns:
                if col != "color":
                    record[col] = row[col]
            rows.append(record)
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Supplementary Module 3 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Title\n\nCell-state interpretation and multi-resolution validation\n\n"
        "## Role\n\n"
        "This module-style supplementary figure consolidates marker-level attribution, GSE202051 full-reference projection, strict NNLS sensitivity and GSE274673 Xenium cell-resolution support.\n\n"
        "## Core claim\n\n"
        "CAF-myeloid cores are supported by recognizable myCAF/matrix and SPP1/TAM-rich states, while immune-decoupled contexts show attenuated immune-state coupling; Xenium data support cell-resolution CAF-domain centering of immune/myeloid programs.\n\n"
        "## Boundary\n\n"
        "This module supports cell-state interpretation of expression-defined spatial programs. It is not immunostaining, image segmentation, causal cell-cell interaction evidence or single-cell-resolved spatial ground truth.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf').relative_to(ROOT)}`\n"
        f"- `{OUT.with_suffix('.svg').relative_to(ROOT)}`\n"
        f"- `{OUT.with_suffix('.png').relative_to(ROOT)}`\n"
        f"- `{SOURCE_OUT.relative_to(ROOT)}`\n",
        encoding="utf-8",
    )


def main() -> None:
    fig = plt.figure(figsize=(15.2, 15.0), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[1.05, 1.0, 1.0, 0.76], hspace=0.88, wspace=0.88)
    fig.suptitle("Cell-state interpretation and multi-resolution validation", fontsize=15.5, fontweight="bold", y=0.986)

    ax_a = fig.add_subplot(gs[0, 0:3])
    ax_b = fig.add_subplot(gs[0, 3:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:5])
    ax_e = fig.add_subplot(gs[1, 5:6])
    ax_f = fig.add_subplot(gs[2, 0:3])
    ax_g = fig.add_subplot(gs[2, 3:5])
    ax_h = fig.add_subplot(gs[2, 5:6])
    ax_note = fig.add_subplot(gs[3, :])

    a = marker_enrichment(ax_a)
    b = reference_enrichment(ax_b)
    c = nnls_enrichment(ax_c)
    d = decoupling_correlations(ax_d)
    e = nnls_projection_agreement(ax_e)
    f = xenium_anchor(ax_f)
    g1, g2 = xenium_scale_coverage(ax_g)
    h = claim_boundary(ax_h)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.76,
        "Interpretation: marker-state enrichment, full GSE202051 projection, strict NNLS sensitivity and Xenium cell-resolution anchors converge on a cellular interpretation of the CAF-myeloid niche: "
        "myCAF/matrix and SPP1/TAM-rich states align with CAF cores, while immune-state coupling is attenuated in immune-decoupled contexts.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.30,
        "Boundary: these are expression-derived and reference-dependent support layers. They should not be described as immunostaining, image segmentation, direct abundance measurement or causal cell-cell interaction evidence.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    write_source_data(a, b, c, d, e, f, g1, g2, h)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
