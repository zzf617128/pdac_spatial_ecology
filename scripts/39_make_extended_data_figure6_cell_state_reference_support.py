from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT_STEM = FIG_DIR / "extended_data_figure6_cell_state_reference_support"

CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
]
CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
}
MARKER_STATES = [
    "myCAF/matrix",
    "SPP1/TAM",
    "DC/APC",
    "T/NK cell",
    "B/plasma cell",
    "epithelial/tumor",
]
REFERENCE_STATES = [
    "myCAF_matrix",
    "SPP1_TAM",
    "DC_APC",
    "T_NK",
    "B_plasma",
    "epithelial_tumor",
]
DISPLAY_LABELS = {
    "myCAF/matrix": "myCAF/matrix",
    "SPP1/TAM": "SPP1/TAM",
    "DC/APC": "DC/APC",
    "T/NK cell": "T/NK",
    "B/plasma cell": "B/plasma",
    "epithelial/tumor": "epithelial/tumor",
    "myCAF_matrix": "myCAF/matrix",
    "SPP1_TAM": "SPP1/TAM",
    "DC_APP": "DC/APC",
    "DC_APC": "DC/APC",
    "T_NK": "T/NK",
    "B_plasma": "B/plasma",
    "epithelial_tumor": "epithelial/tumor",
}


def pivot_summary(df, states):
    values = (
        df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(states)]
        .pivot(index="cell_state", columns="cohort_context", values="median_core_enrichment")
        .reindex(index=states, columns=CONTEXT_ORDER)
    )
    positives = (
        df[df["cohort_context"].isin(CONTEXT_ORDER) & df["cell_state"].isin(states)]
        .assign(label=lambda x: x["n_core_positive"].astype(int).astype(str) + "/" + x["n_samples"].astype(int).astype(str))
        .pivot(index="cell_state", columns="cohort_context", values="label")
        .reindex(index=states, columns=CONTEXT_ORDER)
    )
    return values, positives


def draw_heatmap(ax, data, labels, title, cmap, vmin, vmax, fmt):
    im = ax.imshow(data.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", pad=8)
    ax.set_xticks(np.arange(len(CONTEXT_ORDER)))
    ax.set_xticklabels([CONTEXT_LABELS[c] for c in CONTEXT_ORDER], rotation=35, ha="right", fontsize=8)
    ax.set_yticks(np.arange(len(data.index)))
    ax.set_yticklabels([DISPLAY_LABELS.get(s, s) for s in data.index], fontsize=8)
    ax.tick_params(length=0)
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            val = data.iat[i, j]
            label = labels.iat[i, j]
            if pd.isna(val):
                text = "NA"
            else:
                text = fmt.format(val) + "\n" + str(label)
            ax.text(j, i, text, ha="center", va="center", fontsize=6.6, color="#1a1a1a")
    for spine in ax.spines.values():
        spine.set_linewidth(0.8)
        spine.set_color("#333333")
    return im


def main():
    marker = pd.read_csv(TABLE_DIR / "gap2_cell_state_marker_attribution_context_summary.csv")
    full_ref = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_context_summary.csv")
    corr = pd.read_csv(TABLE_DIR / "gap2_full_reference_projection_deconvolution_correlations.csv")
    comp = pd.read_csv(TABLE_DIR / "gap2_reference_projection_small_vs_full_comparison.csv")

    marker_data, marker_labels = pivot_summary(marker, MARKER_STATES)
    ref_data, ref_labels = pivot_summary(full_ref, REFERENCE_STATES)

    corr_order = ["T_NK", "DC_APC", "B_plasma", "myCAF_matrix", "SPP1_TAM", "epithelial_tumor"]
    corr_plot = (
        corr[(corr["metric"] == "core_enrichment") & (corr["target"] == "immune_decoupling_index")]
        .set_index("cell_state")
        .reindex(corr_order)
    )
    comp_plot = comp.set_index("cell_state").reindex(corr_order)

    mpl.rcParams.update(
        {
            "font.family": "Arial",
            "font.size": 8,
            "axes.titlesize": 10,
            "axes.labelsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )

    fig = plt.figure(figsize=(13.8, 8.2), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.28, 1.28, 0.92],
        height_ratios=[1.0, 0.78],
        left=0.06,
        right=0.985,
        top=0.90,
        bottom=0.10,
        wspace=0.44,
        hspace=0.50,
    )
    fig.suptitle(
        "Extended Data Fig. 6 | Cell-state support for CAF-myeloid core biology",
        fontsize=16,
        fontweight="bold",
        y=0.965,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, :])

    im_a = draw_heatmap(
        ax_a,
        marker_data,
        marker_labels,
        "A  Marker-level CAF-core enrichment",
        "RdBu_r",
        -0.85,
        0.85,
        "{:.2f}",
    )
    cbar_a = fig.colorbar(im_a, ax=ax_a, fraction=0.046, pad=0.025)
    cbar_a.set_label("median enrichment")

    im_b = draw_heatmap(
        ax_b,
        ref_data,
        ref_labels,
        "B  Full-reference GSE202051 projection",
        "RdBu_r",
        -0.06,
        0.06,
        "{:.2f}",
    )
    cbar_b = fig.colorbar(im_b, ax=ax_b, fraction=0.046, pad=0.025)
    cbar_b.set_label("median projection enrichment")

    colors = ["#4C78A8" if x < 0 else "#B279A2" for x in corr_plot["spearman_rho"]]
    y = np.arange(len(corr_plot))
    ax_c.barh(y, corr_plot["spearman_rho"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax_c.axvline(0, color="#333333", linewidth=0.8)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels([DISPLAY_LABELS.get(s, s) for s in corr_plot.index])
    ax_c.invert_yaxis()
    ax_c.set_xlabel("rho with immune-decoupling index")
    ax_c.set_title("C  Association with immune decoupling", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_c.set_xlim(-0.75, 0.65)
    ax_c.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax_c.set_axisbelow(True)
    for yi, row in enumerate(corr_plot.itertuples()):
        x = row.spearman_rho
        if x < 0:
            xpos = x + 0.05
            ha = "left"
            color = "white"
        else:
            xpos = x + 0.035
            ha = "left"
            color = "#1a1a1a"
        ax_c.text(
            xpos,
            yi,
            f"{x:.2f}",
            va="center",
            ha=ha,
            fontsize=8,
            color=color,
        )

    x = np.arange(len(comp_plot))
    width = 0.38
    ax_d.bar(
        x - width / 2,
        comp_plot["spearman_rho"],
        width=width,
        color="#59A14F",
        edgecolor="#333333",
        linewidth=0.5,
        label="per-sample small-vs-full rho",
    )
    ax_d.bar(
        x + width / 2,
        comp_plot["same_direction_fraction"],
        width=width,
        color="#F28E2B",
        edgecolor="#333333",
        linewidth=0.5,
        label="same-direction fraction",
    )
    ax_d.set_xticks(x)
    ax_d.set_xticklabels([DISPLAY_LABELS.get(s, s) for s in comp_plot.index], rotation=0)
    ax_d.set_ylim(0, 1.05)
    ax_d.set_ylabel("agreement")
    ax_d.set_title("D  Small-reference versus full-reference projection stability", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_d.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_d.set_axisbelow(True)
    ax_d.legend(frameon=False, loc="upper left", ncol=2)
    for xi, row in zip(x, comp_plot.itertuples()):
        ax_d.text(xi - width / 2, row.spearman_rho + 0.025, f"{row.spearman_rho:.2f}", ha="center", va="bottom", fontsize=7)
        ax_d.text(
            xi + width / 2,
            row.same_direction_fraction + 0.025,
            f"{row.same_direction_fraction:.2f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )

    note = (
        "Marker attribution and GSE202051 projection support cellular interpretation of CAF-myeloid cores; "
        "projection coefficients are not treated as definitive spatial cell abundance."
    )
    fig.text(0.06, 0.035, note, ha="left", va="bottom", fontsize=8, color="#333333")

    for ext in ["pdf", "png", "svg"]:
        fig.savefig(OUT_STEM.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)
    print(f"Wrote {OUT_STEM}.pdf/.png/.svg")


if __name__ == "__main__":
    main()
