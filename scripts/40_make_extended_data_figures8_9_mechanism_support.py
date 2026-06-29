from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
FIG_DIR.mkdir(parents=True, exist_ok=True)

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
CONTEXT_COLORS = {
    "post_neoadjuvant_sections": "#9C7A3E",
    "treatment_naive_primary": "#C7794F",
    "primary_tumor": "#4C78A8",
    "liver_metastasis": "#59A14F",
    "lymph_node_metastasis": "#9C755F",
}


def setup_style():
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


def annotated_heatmap(ax, values, labels, title, cmap="RdBu_r", vmin=-1, vmax=1, cbar_label=None):
    im = ax.imshow(values.to_numpy(float), cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold", pad=8)
    ax.set_xticks(np.arange(values.shape[1]))
    ax.set_xticklabels([CONTEXT_LABELS.get(c, c) for c in values.columns], rotation=35, ha="right")
    ax.set_yticks(np.arange(values.shape[0]))
    ax.set_yticklabels(values.index)
    ax.tick_params(length=0)
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            val = values.iat[i, j]
            lab = labels.iat[i, j]
            ax.text(j, i, f"{val:.2f}\n{lab}", ha="center", va="center", fontsize=6.5, color="#1a1a1a")
    for spine in ax.spines.values():
        spine.set_color("#333333")
        spine.set_linewidth(0.8)
    cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.025)
    if cbar_label:
        cbar.set_label(cbar_label)
    return im


def format_p(p):
    if p < 1e-4:
        return f"{p:.1e}"
    return f"{p:.3f}"


def make_figure8():
    context = pd.read_csv(TABLE_DIR / "gap1_cxcl9_spp1_polarity_context_summary.csv")
    sample = pd.read_csv(TABLE_DIR / "gap1_cxcl9_spp1_polarity_per_sample.csv")
    corr = pd.read_csv(TABLE_DIR / "gap1_cxcl9_spp1_polarity_decoupling_correlations.csv")

    features = [
        "SPP1",
        "APOE",
        "TREM2",
        "CD68",
        "SPP1 TAM program",
        "CXCL9",
        "CXCL10",
        "CXCL9 IFN program",
        "SPP1-high/CXCL9-low polarity",
    ]
    display = {
        "SPP1 TAM program": "SPP1 TAM\nprogram",
        "CXCL9 IFN program": "CXCL9 IFN\nprogram",
        "SPP1-high/CXCL9-low polarity": "SPP1-high /\nCXCL9-low",
    }
    h = context[context["feature"].isin(features) & context["cohort_context"].isin(CONTEXT_ORDER)].copy()
    h["feature_display"] = h["feature"].map(lambda x: display.get(x, x))
    h["label"] = h["n_core_positive"].astype(int).astype(str) + "/" + h["n_samples"].astype(int).astype(str)
    order_display = [display.get(x, x) for x in features]
    values = h.pivot(index="feature_display", columns="cohort_context", values="median_core_enrichment").reindex(
        index=order_display, columns=CONTEXT_ORDER
    )
    labels = h.pivot(index="feature_display", columns="cohort_context", values="label").reindex(
        index=order_display, columns=CONTEXT_ORDER
    )

    plot_context = (
        context[
            (context["feature"] == "SPP1-high/CXCL9-low polarity")
            & context["cohort_context"].isin(CONTEXT_ORDER)
        ]
        .set_index("cohort_context")
        .reindex(CONTEXT_ORDER)
    )
    corr_features = ["CXCL9 IFN program", "CXCL9", "SPP1 TAM program", "SPP1", "SPP1-high/CXCL9-low polarity"]
    corr_plot = (
        corr[(corr["analysis_group"] == "all_non_normal") & corr["feature"].isin(corr_features)]
        .set_index("feature")
        .reindex(corr_features)
    )

    fig = plt.figure(figsize=(13.2, 7.6), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.45, 1.0, 1.0],
        height_ratios=[1.0, 0.82],
        left=0.07,
        right=0.98,
        top=0.90,
        bottom=0.11,
        wspace=0.42,
        hspace=0.48,
    )
    fig.suptitle(
        "Extended Data Fig. 8 | SPP1/TAM enrichment and CXCL9/IFN attenuation around CAF-myeloid cores",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[0, 2])
    ax_d = fig.add_subplot(gs[1, 1:])

    annotated_heatmap(ax_a, values, labels, "A  CAF-core marker/program enrichment", vmin=-0.85, vmax=0.85, cbar_label="median enrichment")

    x = np.arange(len(CONTEXT_ORDER))
    bars = plot_context["median_core_enrichment"].to_numpy(float)
    ax_b.bar(
        x,
        bars,
        color=[CONTEXT_COLORS[c] for c in CONTEXT_ORDER],
        edgecolor="#333333",
        linewidth=0.5,
    )
    for i, context_id in enumerate(CONTEXT_ORDER):
        vals = sample[sample["cohort_context"] == context_id]["core_enrichment__SPP1-high/CXCL9-low polarity"].dropna()
        jitter = np.linspace(-0.13, 0.13, len(vals)) if len(vals) > 1 else np.array([0.0])
        ax_b.scatter(
            np.full(len(vals), i) + jitter,
            vals,
            s=10,
            color="#2f2f2f",
            alpha=0.35,
            linewidths=0,
            zorder=3,
        )
    ax_b.axhline(0, color="#333333", linewidth=0.8)
    ax_b.set_xticks(x)
    ax_b.set_xticklabels([CONTEXT_LABELS[c] for c in CONTEXT_ORDER], rotation=35, ha="right")
    ax_b.set_ylabel("CAF-core polarity enrichment")
    ax_b.set_title("B  SPP1-high / CXCL9-low polarity", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_b.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_b.set_axisbelow(True)

    corr_labels = ["CXCL9 IFN\nprogram", "CXCL9", "SPP1 TAM\nprogram", "SPP1", "SPP1-high /\nCXCL9-low"]
    y = np.arange(len(corr_plot))
    colors = ["#4C78A8" if v < 0 else "#B279A2" for v in corr_plot["spearman_rho_vs_immune_decoupling"]]
    ax_c.barh(y, corr_plot["spearman_rho_vs_immune_decoupling"], color=colors, edgecolor="#333333", linewidth=0.5)
    ax_c.axvline(0, color="#333333", linewidth=0.8)
    ax_c.set_yticks(y)
    ax_c.set_yticklabels(corr_labels)
    ax_c.invert_yaxis()
    ax_c.set_xlim(-0.50, 0.20)
    ax_c.set_xlabel("rho with immune decoupling")
    ax_c.set_title("C  Decoupling correlation", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_c.grid(axis="x", color="#dddddd", linewidth=0.6)
    ax_c.set_axisbelow(True)
    for yi, row in enumerate(corr_plot.itertuples()):
        rho = row.spearman_rho_vs_immune_decoupling
        label = f"{rho:.2f}\np={format_p(row.p_value)}"
        if rho < 0:
            xpos = rho + 0.025
            ha = "left"
            color = "white"
        else:
            xpos = rho + 0.015
            ha = "left"
            color = "#1a1a1a"
        ax_c.text(
            xpos,
            yi,
            label,
            va="center",
            ha=ha,
            fontsize=7,
            color=color,
        )

    xvals = sample["core_enrichment__SPP1-high/CXCL9-low polarity"].astype(float)
    yvals = sample["immune_decoupling_index"].astype(float)
    for context_id in CONTEXT_ORDER:
        sub = sample[sample["cohort_context"].eq(context_id)].copy()
        ax_d.scatter(
            sub["core_enrichment__SPP1-high/CXCL9-low polarity"].astype(float),
            sub["immune_decoupling_index"].astype(float),
            s=18,
            color=CONTEXT_COLORS[context_id],
            edgecolor="#333333",
            linewidth=0.25,
            alpha=0.62,
            label=CONTEXT_LABELS[context_id],
        )
    mask = np.isfinite(xvals) & np.isfinite(yvals)
    if mask.sum() > 2:
        coef = np.polyfit(xvals[mask], yvals[mask], 1)
        xs = np.linspace(float(xvals[mask].min()), float(xvals[mask].max()), 100)
        ax_d.plot(xs, coef[0] * xs + coef[1], color="#222222", lw=1.1)
    rho_row = corr_plot.loc["SPP1-high/CXCL9-low polarity"]
    ax_d.text(
        0.03,
        0.92,
        f"rho={rho_row.spearman_rho_vs_immune_decoupling:.2f}, p={format_p(rho_row.p_value)}",
        transform=ax_d.transAxes,
        fontsize=8,
        ha="left",
        va="top",
    )
    ax_d.axhline(0, color="#999999", lw=0.7, ls="--")
    ax_d.axvline(0, color="#999999", lw=0.7, ls="--")
    ax_d.set_xlabel("CAF-core SPP1-high / CXCL9-low enrichment")
    ax_d.set_ylabel("immune decoupling index")
    ax_d.set_title("D  Polarity versus immune decoupling", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_d.grid(color="#dddddd", linewidth=0.5)
    ax_d.set_axisbelow(True)
    ax_d.legend(frameon=False, ncol=3, fontsize=7, loc="lower right")

    out = FIG_DIR / "extended_data_figure8_cxcl9_spp1_polarity"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)
    print(f"Wrote {out}.pdf/.png/.svg")


def make_figure9():
    context = pd.read_csv(TABLE_DIR / "gap3_focused_lr_interface_context_summary.csv")
    corr = pd.read_csv(TABLE_DIR / "gap3_focused_lr_interface_correlations.csv")

    axes = ["SPP1-CD44/integrin", "TGF-beta/TGFBR", "matrix-integrin", "IL6-OSM/LIF-JAKSTAT"]
    metrics = ["ligand_core_enrichment", "receptor_interface_enrichment", "response_interface_enrichment"]
    metric_labels = {
        "ligand_core_enrichment": "ligand core",
        "receptor_interface_enrichment": "receptor interface",
        "response_interface_enrichment": "response interface",
    }
    h = context[
        context["axis"].isin(axes) & context["metric"].isin(metrics) & context["cohort_context"].isin(CONTEXT_ORDER)
    ].copy()
    h["row"] = h["axis"] + "\n" + h["metric"].map(metric_labels)
    h["label"] = h["n_positive"].astype(int).astype(str) + "/" + h["n_samples"].astype(int).astype(str)
    row_order = [axis + "\n" + metric_labels[m] for axis in axes for m in metrics]
    values = h.pivot(index="row", columns="cohort_context", values="median_value").reindex(index=row_order, columns=CONTEXT_ORDER)
    labels = h.pivot(index="row", columns="cohort_context", values="label").reindex(index=row_order, columns=CONTEXT_ORDER)

    scores = context[
        (context["axis"].isin(axes))
        & (context["metric"] == "directional_core_to_interface_score")
        & context["cohort_context"].isin(CONTEXT_ORDER)
    ].copy()
    score_pivot = scores.pivot(index="axis", columns="cohort_context", values="median_value").reindex(index=axes, columns=CONTEXT_ORDER)

    corr_plot = (
        corr[(corr["metric"] == "directional_core_to_interface_score") & corr["axis"].isin(axes)]
        .pivot(index="axis", columns="target", values="spearman_rho")
        .reindex(index=axes)
    )
    targets = ["immune_decoupling_index", "stromal_tumor_core_coupling", "immune_core_coupling"]
    target_labels = {
        "immune_decoupling_index": "immune\ndecoupling",
        "stromal_tumor_core_coupling": "stromal-tumor\ncoupling",
        "immune_core_coupling": "immune-core\ncoupling",
    }

    fig = plt.figure(figsize=(14.0, 8.0), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        3,
        width_ratios=[1.40, 1.05, 1.05],
        height_ratios=[1.0, 0.80],
        left=0.18,
        right=0.985,
        top=0.90,
        bottom=0.12,
        wspace=0.48,
        hspace=0.52,
    )
    fig.suptitle(
        "Extended Data Fig. 9 | Focused CAF-core/interface candidate communication axes",
        fontsize=15,
        fontweight="bold",
        y=0.965,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1:])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[1, 2])

    annotated_heatmap(ax_a, values, labels, "A  Directional component enrichment", vmin=-0.85, vmax=0.85, cbar_label="median enrichment")

    x = np.arange(len(axes))
    width = 0.15
    for k, context_id in enumerate(CONTEXT_ORDER):
        ax_b.bar(
            x + (k - 2) * width,
            score_pivot[context_id],
            width=width,
            color=CONTEXT_COLORS[context_id],
            edgecolor="#333333",
            linewidth=0.4,
            label=CONTEXT_LABELS[context_id],
        )
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(axes, rotation=18, ha="right")
    ax_b.set_ylabel("ligand core + receptor/response interface")
    ax_b.set_title("B  Combined directional nomination score", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_b.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_b.set_axisbelow(True)
    ax_b.legend(frameon=False, ncol=3, fontsize=7, loc="upper left")

    corr_heat = corr_plot[targets]
    im = ax_c.imshow(corr_heat.to_numpy(float), cmap="RdBu_r", vmin=-0.50, vmax=0.50, aspect="auto")
    ax_c.set_xticks(np.arange(len(targets)))
    ax_c.set_xticklabels([target_labels[t] for t in targets], rotation=0)
    ax_c.set_yticks(np.arange(len(axes)))
    ax_c.set_yticklabels(axes)
    ax_c.set_title("C  Correlation of nomination score", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_c.tick_params(length=0)
    for i in range(corr_heat.shape[0]):
        for j in range(corr_heat.shape[1]):
            ax_c.text(j, i, f"{corr_heat.iat[i, j]:.2f}", ha="center", va="center", fontsize=7)
    plt.colorbar(im, ax=ax_c, fraction=0.046, pad=0.025, label="Spearman rho")

    ax_d.axis("off")
    ax_d.set_title("D  Interface-score schematic", loc="left", fontsize=10, fontweight="bold", pad=8)
    boxes = [
        (0.08, 0.66, "CAF core", "ligand", COLORS["caf"] if "COLORS" in globals() else "#8C2D04"),
        (0.56, 0.66, "Tumor-stroma\ninterface", "receptor", "#B23A48"),
        (0.56, 0.25, "Tumor-stroma\ninterface", "response", "#4C78A8"),
    ]
    for x0, y0, title, subtitle, color in boxes:
        ax_d.add_patch(
            plt.Rectangle((x0, y0), 0.34, 0.20, transform=ax_d.transAxes, facecolor="#F7F7F7", edgecolor=color, lw=1.2)
        )
        ax_d.text(x0 + 0.17, y0 + 0.125, title, ha="center", va="center", fontsize=7.5, fontweight="bold", transform=ax_d.transAxes)
        ax_d.text(x0 + 0.17, y0 + 0.055, subtitle, ha="center", va="center", fontsize=7, color="#333333", transform=ax_d.transAxes)
    ax_d.annotate(
        "",
        xy=(0.56, 0.76),
        xytext=(0.42, 0.76),
        xycoords=ax_d.transAxes,
        arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#666666"),
    )
    ax_d.annotate(
        "",
        xy=(0.70, 0.45),
        xytext=(0.70, 0.66),
        xycoords=ax_d.transAxes,
        arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#666666"),
    )
    ax_d.text(0.44, 0.81, "candidate\naxis", ha="center", va="bottom", fontsize=7, color="#333333", transform=ax_d.transAxes)
    ax_d.text(0.50, 0.08, "score = ligand-core + receptor-interface + response-interface", fontsize=7.2, ha="center", transform=ax_d.transAxes)

    out = FIG_DIR / "extended_data_figure9_focused_interface_axes"
    for ext in ["pdf", "png", "svg"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)
    print(f"Wrote {out}.pdf/.png/.svg")


def main():
    setup_style()
    make_figure8()
    make_figure9()


if __name__ == "__main__":
    main()
