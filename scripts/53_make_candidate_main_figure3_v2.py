from __future__ import annotations

import importlib.util
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
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.8,
        "ytick.major.size": 2.8,
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
OUT = FIG_DIR / "figure3_candidate_v2_ecotype_interface_story"


def load_ecotype_helpers():
    path = ROOT / "scripts" / "52_make_ecotype_flow_figure.py"
    spec = importlib.util.spec_from_file_location("ecotype_flow", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.06, 1.04, label, transform=ax.transAxes, fontsize=13, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def plot_flow(ax: plt.Axes, helper) -> None:
    counts = pd.read_csv(TABLE_DIR / "spatial_ecotype_context_counts.csv")
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
    context_sizes = [int(counts.loc[counts["cohort_context"].eq(c), "n_samples"].sum()) for c in contexts]
    ecotype_sizes = [int(counts.loc[counts["dominant_nmf_ecotype"].eq(e), "n_samples"].sum()) for e in ecotypes]
    left_pos = helper.stacked_positions(contexts, context_sizes, gap=0.018, y0=0.04, y1=0.88)
    right_pos = helper.stacked_positions(ecotypes, ecotype_sizes, gap=0.035, y0=0.04, y1=0.88)
    x_left, x_right = 0.14, 0.83
    left_cursor = {k: left_pos[k][0] for k in contexts}
    right_cursor = {k: right_pos[k][0] for k in ecotypes}
    flow_scale = (0.84 - 0.018 * (len(contexts) - 1)) / counts["n_samples"].sum()

    for context in contexts:
        sub = counts[counts["cohort_context"].eq(context)].set_index("dominant_nmf_ecotype")
        for eco in ecotypes:
            if eco not in sub.index:
                continue
            n = float(sub.loc[eco, "n_samples"])
            if n <= 0:
                continue
            h = n * flow_scale
            helper.ribbon(
                ax,
                x_left + 0.03,
                x_right - 0.03,
                left_cursor[context],
                left_cursor[context] + h,
                right_cursor[eco],
                right_cursor[eco] + h,
                helper.ECOTYPE_COLORS[eco],
                alpha=0.70,
            )
            left_cursor[context] += h
            right_cursor[eco] += h

    for context, (y0, y1) in left_pos.items():
        ax.add_patch(Rectangle((x_left - 0.03, y0), 0.03, y1 - y0, color="#333333"))
        ax.text(
            x_left - 0.045,
            (y0 + y1) / 2,
            f"{helper.clean_context(context)} (n={context_sizes[contexts.index(context)]})",
            ha="right",
            va="center",
            fontsize=7.1,
        )
    short = {
        "NMF1": "basal/\ntumor",
        "NMF2": "lymphoid",
        "NMF3": "IFN/APC",
        "NMF4": "EMT/\nmyCAF",
    }
    for eco, (y0, y1) in right_pos.items():
        ax.add_patch(Rectangle((x_right, y0), 0.035, y1 - y0, color=helper.ECOTYPE_COLORS[eco]))
        ax.text(x_right + 0.045, (y0 + y1) / 2, f"{eco}\n{short[eco]}", ha="left", va="center", fontsize=7.8)
    ax.text(0.02, 0.98, "A", transform=ax.transAxes, fontsize=13, fontweight="bold", va="top")
    ax.set_title("Context-to-CAF-core ecotype architecture", loc="left", fontsize=10, fontweight="bold")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def plot_decoupling(ax: plt.Axes, helper) -> None:
    samples = pd.read_csv(TABLE_DIR / "spatial_ecotype_sample_summary.csv")
    ecotypes = ["NMF1", "NMF2", "NMF3", "NMF4"]
    labels = ["basal/\ntumor", "lymphoid", "IFN/APC", "EMT/\nmyCAF"]
    data = [samples.loc[samples["dominant_nmf_ecotype"].eq(e), "immune_decoupling_index"].dropna().values for e in ecotypes]
    parts = ax.violinplot(data, positions=np.arange(len(ecotypes)), showmedians=False, widths=0.75)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(helper.ECOTYPE_COLORS[ecotypes[i]])
        body.set_edgecolor("none")
        body.set_alpha(0.34)
    rng = np.random.default_rng(17)
    for i, vals in enumerate(data):
        ax.scatter(
            i + rng.uniform(-0.12, 0.12, size=len(vals)),
            vals,
            s=15,
            color=helper.ECOTYPE_COLORS[ecotypes[i]],
            alpha=0.72,
            edgecolor="white",
            linewidth=0.25,
        )
        if len(vals):
            ax.plot([i - 0.23, i + 0.23], [np.median(vals), np.median(vals)], color="black", lw=1.1)
    panel_label(ax, "B")
    ax.axhline(0, color="#777777", lw=0.8, ls=":")
    ax.set_xticks(np.arange(len(ecotypes)), labels)
    ax.set_ylabel("immune-decoupling index", fontsize=8)
    ax.set_title("Ecotypes stratify immune decoupling", loc="left", fontsize=10, fontweight="bold")
    clean_axes(ax)


def plot_candidate_axis_heatmap(ax: plt.Axes) -> None:
    axes = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_7A_B.csv")
    contexts = [
        "post_neoadjuvant_sections",
        "treatment_naive_primary",
        "primary_tumor",
        "liver_metastasis",
        "lymph_node_metastasis",
    ]
    context_labels = ["post-NACT", "treatment-\nnaive", "primary", "liver met", "LN met"]
    axis_order = ["SPP1-TAM/matrix", "TGF-beta/EMT invasive", "IFN/APC antigen", "B/plasma lymphoid", "T cell/checkpoint"]
    core = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_core_enrichment", aggfunc="median").reindex(index=axis_order, columns=contexts)
    interface = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_interface_enrichment", aggfunc="median").reindex(index=axis_order, columns=contexts)
    mat = np.concatenate([core.to_numpy(float), interface.to_numpy(float)], axis=1)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax.axvline(len(contexts) - 0.5, color="#333333", lw=0.9)
    ax.set_xticks(np.arange(mat.shape[1]), context_labels + context_labels, rotation=35, ha="right", fontsize=6.7)
    ax.set_yticks(np.arange(len(axis_order)), axis_order, fontsize=7.4)
    ax.text(0.24, 1.04, "CAF core", transform=ax.transAxes, ha="center", fontsize=8, fontweight="bold")
    ax.text(0.76, 1.04, "interface", transform=ax.transAxes, ha="center", fontsize=8, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=5.9)
    panel_label(ax, "C")
    ax.set_title("Candidate axes localize to CAF-core/interface regions", loc="left", fontsize=10, fontweight="bold", pad=20)
    cb = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cb.ax.tick_params(labelsize=7)
    cb.set_label("median enrichment", fontsize=7)


def plot_interface_maps(ax: plt.Axes) -> None:
    path = FIG_DIR / "extended_data_figure16_interface_compartment_maps.png"
    img = Image.open(path).convert("RGB")
    arr = np.asarray(img)
    mask = np.any(arr < 248, axis=2)
    if mask.any():
        ys, xs = np.where(mask)
        img = img.crop((max(0, xs.min() - 12), max(0, ys.min() - 12), min(img.width, xs.max() + 12), min(img.height, ys.max() + 12)))
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])
    panel_label(ax, "D")
    ax.set_title("Representative program-defined interface maps", loc="left", fontsize=10, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)


def save(fig: plt.Figure) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        path = OUT.with_suffix(f".{ext}")
        if ext == "png":
            fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
        else:
            fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    helper = load_ecotype_helpers()
    fig = plt.figure(figsize=(14.0, 9.2))
    gs = GridSpec(2, 3, figure=fig, height_ratios=[1.05, 1.05], width_ratios=[1.25, 1.0, 1.25], hspace=0.42, wspace=0.36)
    plot_flow(fig.add_subplot(gs[0, :2]), helper)
    plot_decoupling(fig.add_subplot(gs[0, 2]), helper)
    plot_candidate_axis_heatmap(fig.add_subplot(gs[1, :2]))
    plot_interface_maps(fig.add_subplot(gs[1, 2]))
    fig.suptitle("CAF-core ecotypes nominate immune-decoupled invasive-interface states", fontsize=14, fontweight="bold", y=0.985)
    fig.text(
        0.5,
        0.015,
        "Candidate main Figure 3 v2. Flow and compartments are expression-derived visual summaries; they do not imply temporal transitions, histologic annotations or causal ligand-receptor signaling.",
        ha="center",
        fontsize=8.4,
    )
    save(fig)
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
