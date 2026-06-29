from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"

OUT = FIG_DIR / "extended_data_figure20_ecotype_context_flow"
SOURCE = SOURCE_DIR / "Source_Data_Extended_Data_Fig_20_ecotype_context_flow.csv"

CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal pancreas",
}

ECOTYPE_COLORS = {
    "NMF1": "#c84e4e",
    "NMF2": "#4c9f70",
    "NMF3": "#4c78a8",
    "NMF4": "#a8795b",
}


def clean_context(value: str) -> str:
    return CONTEXT_LABELS.get(value, value.replace("_", " "))


def ribbon(ax, x0, x1, y0a, y0b, y1a, y1b, color, alpha=0.72):
    dx = x1 - x0
    verts = [
        (x0, y0a),
        (x0 + 0.45 * dx, y0a),
        (x1 - 0.45 * dx, y1a),
        (x1, y1a),
        (x1, y1b),
        (x1 - 0.45 * dx, y1b),
        (x0 + 0.45 * dx, y0b),
        (x0, y0b),
        (x0, y0a),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.LINETO,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CURVE4,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(PathPatch(MplPath(verts, codes), facecolor=color, edgecolor="none", alpha=alpha))


def stacked_positions(keys, sizes, gap=0.015, y0=0.0, y1=1.0):
    total = sum(sizes)
    height = y1 - y0
    usable = height - gap * (len(keys) - 1)
    scale = usable / total if total else 1.0
    y = y0
    positions = {}
    for key, size in zip(keys, sizes):
        h = size * scale
        positions[key] = (y, y + h)
        y += h + gap
    return positions


def main() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    SOURCE_DIR.mkdir(parents=True, exist_ok=True)

    counts = pd.read_csv(TABLE_DIR / "spatial_ecotype_context_counts.csv")
    samples = pd.read_csv(TABLE_DIR / "spatial_ecotype_sample_summary.csv")

    contexts = [
        c
        for c in [
            "post_neoadjuvant_sections",
            "treatment_naive_primary",
            "primary_tumor",
            "liver_metastasis",
            "lymph_node_metastasis",
            "normal_pancreas",
        ]
        if c in set(counts["cohort_context"])
    ]
    ecotypes = ["NMF1", "NMF2", "NMF3", "NMF4"]

    counts = counts[counts["cohort_context"].isin(contexts)].copy()
    counts["context_label"] = counts["cohort_context"].map(clean_context)
    counts["ecotype_label_short"] = counts["dominant_nmf_ecotype"].map(
        {
            "NMF1": "basal/tumor",
            "NMF2": "lymphoid",
            "NMF3": "IFN/APC",
            "NMF4": "EMT/myCAF",
        }
    )
    counts.to_csv(SOURCE, index=False)

    context_sizes = [int(counts.loc[counts["cohort_context"] == c, "n_samples"].sum()) for c in contexts]
    ecotype_sizes = [
        int(counts.loc[counts["dominant_nmf_ecotype"] == e, "n_samples"].sum()) for e in ecotypes
    ]
    left_pos = stacked_positions(contexts, context_sizes, gap=0.018, y0=0.02, y1=0.91)
    right_pos = stacked_positions(ecotypes, ecotype_sizes, gap=0.035, y0=0.02, y1=0.91)

    fig = plt.figure(figsize=(12.8, 8.6), constrained_layout=False)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 0.95], width_ratios=[1.15, 1.0], hspace=0.34, wspace=0.28)
    ax_flow = fig.add_subplot(gs[0, :])
    ax_dec = fig.add_subplot(gs[1, 0])
    ax_frac = fig.add_subplot(gs[1, 1])

    x_left, x_right = 0.12, 0.88
    left_cursor = {k: left_pos[k][0] for k in contexts}
    right_cursor = {k: right_pos[k][0] for k in ecotypes}
    total = counts["n_samples"].sum()
    usable_left = 0.89 - 0.018 * (len(contexts) - 1)
    flow_scale = usable_left / total

    for context in contexts:
        sub = counts[counts["cohort_context"] == context].set_index("dominant_nmf_ecotype")
        for eco in ecotypes:
            if eco not in sub.index:
                continue
            n = float(sub.loc[eco, "n_samples"])
            if n <= 0:
                continue
            h = n * flow_scale
            y0a, y0b = left_cursor[context], left_cursor[context] + h
            y1a, y1b = right_cursor[eco], right_cursor[eco] + h
            ribbon(ax_flow, x_left + 0.03, x_right - 0.03, y0a, y0b, y1a, y1b, ECOTYPE_COLORS[eco])
            left_cursor[context] += h
            right_cursor[eco] += h

    for context, (y0, y1) in left_pos.items():
        ax_flow.add_patch(Rectangle((x_left - 0.03, y0), 0.035, y1 - y0, color="#303030"))
        ax_flow.text(
            x_left - 0.045,
            (y0 + y1) / 2,
            f"{clean_context(context)} (n={context_sizes[contexts.index(context)]})",
            ha="right",
            va="center",
            fontsize=8.0,
        )
    for eco, (y0, y1) in right_pos.items():
        label = counts.loc[counts["dominant_nmf_ecotype"] == eco, "dominant_nmf_ecotype_label"].dropna()
        label_text = label.iloc[0] if len(label) else eco
        ax_flow.add_patch(Rectangle((x_right - 0.005, y0), 0.035, y1 - y0, color=ECOTYPE_COLORS[eco]))
        ax_flow.text(x_right + 0.045, (y0 + y1) / 2, f"{eco}: {label_text}\n(n={ecotype_sizes[ecotypes.index(eco)]})",
                     ha="left", va="center", fontsize=8.5)
    ax_flow.text(0.02, 0.98, "A", transform=ax_flow.transAxes, fontsize=14, fontweight="bold")
    ax_flow.text(
        0.06,
        0.985,
        "Context-to-CAF-core ecotype flow",
        transform=ax_flow.transAxes,
        fontsize=11,
        fontweight="bold",
        va="top",
    )
    ax_flow.set_xlim(0, 1)
    ax_flow.set_ylim(0, 1)
    ax_flow.axis("off")

    sample_sub = samples[samples["dominant_nmf_ecotype"].isin(ecotypes)].copy()
    data = [
        sample_sub.loc[sample_sub["dominant_nmf_ecotype"] == eco, "immune_decoupling_index"].dropna().values
        for eco in ecotypes
    ]
    parts = ax_dec.violinplot(data, positions=np.arange(len(ecotypes)), showmedians=False, widths=0.78)
    for idx, body in enumerate(parts["bodies"]):
        body.set_facecolor(ECOTYPE_COLORS[ecotypes[idx]])
        body.set_edgecolor("none")
        body.set_alpha(0.34)
    rng = np.random.default_rng(7)
    for idx, vals in enumerate(data):
        x = idx + rng.uniform(-0.12, 0.12, size=len(vals))
        ax_dec.scatter(x, vals, s=18, color=ECOTYPE_COLORS[ecotypes[idx]], alpha=0.72, edgecolor="white", linewidth=0.25)
        if len(vals):
            ax_dec.plot([idx - 0.25, idx + 0.25], [np.median(vals), np.median(vals)], color="black", lw=1.2)
    ax_dec.axhline(0, color="#777777", lw=0.8, ls=":")
    ax_dec.set_xticks(np.arange(len(ecotypes)))
    ax_dec.set_xticklabels(["basal/\ntumor", "lymphoid", "IFN/APC", "EMT/\nmyCAF"], fontsize=8.5)
    ax_dec.set_ylabel("immune-decoupling index")
    ax_dec.set_title("B  Immune-decoupling by dominant ecotype", loc="left", fontsize=10, fontweight="bold")
    ax_dec.spines[["top", "right"]].set_visible(False)

    frac = (
        counts.pivot_table(index="context_label", columns="ecotype_label_short", values="n_samples", aggfunc="sum", fill_value=0)
        .reindex([clean_context(c) for c in contexts])
    )
    frac = frac.div(frac.sum(axis=1), axis=0)
    order = ["basal/tumor", "lymphoid", "IFN/APC", "EMT/myCAF"]
    bottom = np.zeros(len(frac))
    for label in order:
        vals = frac[label].values if label in frac.columns else np.zeros(len(frac))
        eco = {"basal/tumor": "NMF1", "lymphoid": "NMF2", "IFN/APC": "NMF3", "EMT/myCAF": "NMF4"}[label]
        ax_frac.barh(np.arange(len(frac)), vals, left=bottom, color=ECOTYPE_COLORS[eco], label=label, height=0.72)
        bottom += vals
    ax_frac.set_yticks(np.arange(len(frac)))
    ax_frac.set_yticklabels(frac.index, fontsize=8.5)
    ax_frac.set_xlim(0, 1)
    ax_frac.set_xlabel("fraction of samples")
    ax_frac.set_title("C  Ecotype composition by context", loc="left", fontsize=10, fontweight="bold")
    ax_frac.legend(frameon=False, fontsize=8, loc="lower center", bbox_to_anchor=(0.5, -0.33), ncol=2)
    ax_frac.spines[["top", "right"]].set_visible(False)

    fig.suptitle("CAF-core ecotype architecture across PDAC contexts", fontsize=13.5, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.015,
        "Ecotypes are NMF-derived expression programs in CAF-core-enriched regions; flows summarize sample-level dominant ecotypes and do not imply temporal transitions.",
        ha="center",
        fontsize=8.5,
    )
    fig.savefig(OUT.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(OUT.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)
    print(f"Wrote {OUT.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
