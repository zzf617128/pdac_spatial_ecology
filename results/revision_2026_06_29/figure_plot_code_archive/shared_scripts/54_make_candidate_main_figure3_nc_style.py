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
from PIL import Image, ImageOps


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
OUT = FIG_DIR / "figure3_candidate_nc_style_ecotype_interface_story"
ACTIVE_OUT = FIG_DIR / "figure3_submission_ecotypes_mechanism_axes_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Fig_3_candidate_NC_style_panel_index.csv"


CONTEXTS = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
]
CONTEXT_LABELS = ["post-NACT", "treatment-\nnaive", "primary", "liver met", "LN met"]
AXIS_ORDER = ["SPP1-TAM/matrix", "TGF-beta/EMT invasive", "IFN/APC antigen", "B/plasma lymphoid", "T cell/checkpoint"]


def load_module(script_name: str, name: str):
    path = ROOT / "scripts" / script_name
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.05) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def plot_ecotype_loadings(ax: plt.Axes) -> None:
    nmf = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_5A.csv")
    programs = [
        "myCAF",
        "iCAF",
        "apCAF",
        "myeloid",
        "SPP1/TREM2 TAM",
        "TGF-beta",
        "EMT/invasion",
        "hypoxia",
        "basal-like",
        "tumor aggressive",
        "IFN/MHC",
        "immune core",
        "T cell",
        "B cell",
        "DC/APC",
        "plasma cell",
    ]
    mat = nmf.set_index("nmf_ecotype")[programs]
    mat = mat.div(mat.max(axis=1), axis=0)
    im = ax.imshow(mat.to_numpy(float), cmap="viridis", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(programs)), programs, rotation=50, ha="right", fontsize=6.1)
    ax.set_yticks(np.arange(len(mat)), mat.index, fontsize=7)
    ax.set_title("CAF-core ecotype programs", loc="left", fontsize=9.5, fontweight="bold")
    panel_label(ax, "A")
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cb.ax.tick_params(labelsize=6.5)


def plot_flow(ax: plt.Axes, flow) -> None:
    counts = pd.read_csv(TABLE_DIR / "spatial_ecotype_context_counts.csv")
    contexts = [c for c in CONTEXTS + ["normal_pancreas"] if c in set(counts["cohort_context"])]
    ecotypes = ["NMF1", "NMF2", "NMF3", "NMF4"]
    counts = counts[counts["cohort_context"].isin(contexts)].copy()
    context_sizes = [int(counts.loc[counts["cohort_context"].eq(c), "n_samples"].sum()) for c in contexts]
    ecotype_sizes = [int(counts.loc[counts["dominant_nmf_ecotype"].eq(e), "n_samples"].sum()) for e in ecotypes]
    left_pos = flow.stacked_positions(contexts, context_sizes, gap=0.018, y0=0.04, y1=0.88)
    right_pos = flow.stacked_positions(ecotypes, ecotype_sizes, gap=0.035, y0=0.04, y1=0.88)
    x_left, x_right = 0.18, 0.78
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
            flow.ribbon(
                ax,
                x_left + 0.025,
                x_right - 0.025,
                left_cursor[context],
                left_cursor[context] + h,
                right_cursor[eco],
                right_cursor[eco] + h,
                flow.ECOTYPE_COLORS[eco],
                alpha=0.68,
            )
            left_cursor[context] += h
            right_cursor[eco] += h
    for context, (y0, y1) in left_pos.items():
        ax.add_patch(Rectangle((x_left - 0.03, y0), 0.025, y1 - y0, color="#333333"))
        ax.text(
            x_left - 0.042,
            (y0 + y1) / 2,
            f"{flow.clean_context(context)}\n{context_sizes[contexts.index(context)]}",
            ha="right",
            va="center",
            fontsize=6.4,
        )
    short = {"NMF1": "basal/tumor", "NMF2": "lymphoid", "NMF3": "IFN/APC", "NMF4": "EMT/myCAF"}
    for eco, (y0, y1) in right_pos.items():
        ax.add_patch(Rectangle((x_right, y0), 0.028, y1 - y0, color=flow.ECOTYPE_COLORS[eco]))
        ax.text(x_right + 0.038, (y0 + y1) / 2, f"{eco}\n{short[eco]}", ha="left", va="center", fontsize=6.5)
    ax.set_title("Context-to-ecotype architecture", loc="left", fontsize=9.5, fontweight="bold")
    panel_label(ax, "B", x=0.0, y=0.98)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")


def plot_composition(ax: plt.Axes, flow) -> None:
    counts = pd.read_csv(TABLE_DIR / "spatial_ecotype_context_counts.csv")
    counts = counts[counts["cohort_context"].isin(CONTEXTS)].copy()
    short = {"NMF1": "basal/tumor", "NMF2": "lymphoid", "NMF3": "IFN/APC", "NMF4": "EMT/myCAF"}
    counts["short"] = counts["dominant_nmf_ecotype"].map(short)
    frac = counts.pivot_table(index="cohort_context", columns="short", values="n_samples", aggfunc="sum", fill_value=0).reindex(CONTEXTS)
    frac = frac.div(frac.sum(axis=1), axis=0)
    order = ["basal/tumor", "lymphoid", "IFN/APC", "EMT/myCAF"]
    colors = {
        "basal/tumor": flow.ECOTYPE_COLORS["NMF1"],
        "lymphoid": flow.ECOTYPE_COLORS["NMF2"],
        "IFN/APC": flow.ECOTYPE_COLORS["NMF3"],
        "EMT/myCAF": flow.ECOTYPE_COLORS["NMF4"],
    }
    bottom = np.zeros(len(frac))
    for label in order:
        vals = frac[label].to_numpy(float) if label in frac.columns else np.zeros(len(frac))
        ax.barh(np.arange(len(frac)), vals, left=bottom, color=colors[label], label=label, height=0.72)
        bottom += vals
    ax.set_yticks(np.arange(len(frac)), ["post-NACT", "treatment-naive", "primary", "liver met", "LN met"], fontsize=6.8)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction", fontsize=7)
    ax.set_title("Ecotype composition", loc="left", fontsize=9.5, fontweight="bold")
    # Ecotype colors are labeled in panel B; omitting the repeated legend prevents
    # overlap with the association panel in the compact NC-style layout.
    panel_label(ax, "C")
    clean_axes(ax)


def plot_decoupling(ax: plt.Axes, flow) -> None:
    samples = pd.read_csv(TABLE_DIR / "spatial_ecotype_sample_summary.csv")
    ecotypes = ["NMF1", "NMF2", "NMF3", "NMF4"]
    labels = ["basal/\ntumor", "lymphoid", "IFN/APC", "EMT/\nmyCAF"]
    data = [samples.loc[samples["dominant_nmf_ecotype"].eq(e), "immune_decoupling_index"].dropna().values for e in ecotypes]
    parts = ax.violinplot(data, positions=np.arange(len(ecotypes)), showmedians=False, widths=0.74)
    for i, body in enumerate(parts["bodies"]):
        body.set_facecolor(flow.ECOTYPE_COLORS[ecotypes[i]])
        body.set_edgecolor("none")
        body.set_alpha(0.34)
    rng = np.random.default_rng(17)
    for i, vals in enumerate(data):
        ax.scatter(i + rng.uniform(-0.12, 0.12, size=len(vals)), vals, s=12, color=flow.ECOTYPE_COLORS[ecotypes[i]], alpha=0.72, edgecolor="white", linewidth=0.22)
        if len(vals):
            ax.plot([i - 0.23, i + 0.23], [np.median(vals), np.median(vals)], color="black", lw=1.0)
    ax.axhline(0, color="#777777", lw=0.7, ls=":")
    ax.set_xticks(np.arange(len(ecotypes)), labels, fontsize=6.8)
    ax.set_ylabel("immune-decoupling index", fontsize=7)
    ax.set_title("Decoupling by ecotype", loc="left", fontsize=9.5, fontweight="bold")
    panel_label(ax, "D")
    clean_axes(ax)


def plot_axis_heatmap(ax: plt.Axes) -> None:
    axes = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_7A_B.csv")
    core = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_core_enrichment", aggfunc="median").reindex(index=AXIS_ORDER, columns=CONTEXTS)
    interface = axes.pivot_table(index="axis_label", columns="cohort_context", values="median_interface_enrichment", aggfunc="median").reindex(index=AXIS_ORDER, columns=CONTEXTS)
    mat = np.concatenate([core.to_numpy(float), interface.to_numpy(float)], axis=1)
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-1.2, vmax=1.2, aspect="auto")
    ax.axvline(len(CONTEXTS) - 0.5, color="#333333", lw=0.9)
    ax.set_xticks(np.arange(mat.shape[1]), CONTEXT_LABELS + CONTEXT_LABELS, rotation=35, ha="right", fontsize=6.3)
    ax.set_yticks(np.arange(len(AXIS_ORDER)), AXIS_ORDER, fontsize=7)
    ax.text(0.24, 1.04, "CAF core", transform=ax.transAxes, ha="center", fontsize=7.5, fontweight="bold")
    ax.text(0.76, 1.04, "interface", transform=ax.transAxes, ha="center", fontsize=7.5, fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=5.2)
    ax.set_title("Candidate-axis localization", loc="left", fontsize=9.5, fontweight="bold", pad=18)
    panel_label(ax, "E")
    cb = plt.colorbar(im, ax=ax, fraction=0.028, pad=0.015)
    cb.ax.tick_params(labelsize=6.5)
    cb.set_label("median enrichment", fontsize=6.5)


def plot_axis_correlation(ax: plt.Axes) -> None:
    corr = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_7C.csv")
    c = corr[corr["metric"].eq("core_enrichment")].set_index("axis_label").reindex(AXIS_ORDER)
    vals = c["rho_with_immune_decoupling_index"].astype(float)
    y = np.arange(len(vals))
    ax.barh(y, vals, color=["#4C78A8" if v > 0 else "#B279A2" for v in vals], height=0.72)
    ax.axvline(0, color="#333333", lw=0.8)
    ax.set_yticks(y, c.index, fontsize=6.6)
    ax.set_xlabel("rho", fontsize=7)
    ax.set_title("Association with immune decoupling", loc="left", fontsize=9.5, fontweight="bold")
    ax.set_xlim(-0.75, 0.55)
    for yy, val in zip(y, vals):
        ax.text(val + (0.03 if val >= 0 else -0.03), yy, f"{val:.2f}", va="center", ha="left" if val >= 0 else "right", fontsize=6.3)
    panel_label(ax, "F")
    clean_axes(ax)


def get_spatial_sample_rows() -> list[dict[str, str]]:
    manifest = pd.read_csv(TABLE_DIR / "gse272362_rds_overlay_manifest.csv")
    desired = [
        ("IU_PDA_T3", "primary"),
        ("IU_PDA_HM10", "liver met"),
        ("IU_PDA_LNM7", "LN met"),
    ]
    rows = []
    for sample_id, label in desired:
        row = manifest[manifest["sample_id"].eq(sample_id)].iloc[0].to_dict()
        row["label"] = label
        rows.append(row)
    return rows


def spatial_limits(spots: pd.DataFrame, scale: float, image: Image.Image) -> tuple[float, float, float, float]:
    x = spots["x_pixel"].to_numpy(float) * scale
    y = spots["y_pixel"].to_numpy(float) * scale
    xmin, xmax = np.nanquantile(x, [0.03, 0.97])
    ymin, ymax = np.nanquantile(y, [0.03, 0.97])
    dx = xmax - xmin
    dy = ymax - ymin
    return max(0, xmin - 0.10 * dx), min(image.width, xmax + 0.10 * dx), max(0, ymin - 0.10 * dy), min(image.height, ymax + 0.10 * dy)


def plot_compartment_map(ax: plt.Axes, helper, sample_row: dict[str, str], letter: str) -> None:
    spots = helper.classify_spot_compartments(helper.load_spots("GSE272362", sample_row["sample_id"]))
    image_path, scale_path = helper.image_paths_for_sample("GSE272362", sample_row["sample_id"], sample_row["geo_base"])
    image = Image.open(image_path).convert("RGB")
    scale = helper.read_scale(scale_path)
    xmin, xmax, ymin, ymax = spatial_limits(spots, scale, image)
    ax.imshow(ImageOps.grayscale(image).convert("RGB"), alpha=0.40)
    colors = {k: helper.PROGRAM_COLORS[k] for k in ["CAF core", "tumor-high", "interface", "other"]}
    for label, color in colors.items():
        sub = spots[spots["spatial_compartment"].eq(label)]
        if len(sub):
            ax.scatter(sub["x_pixel"].to_numpy(float) * scale, sub["y_pixel"].to_numpy(float) * scale, s=5.0, c=color, alpha=0.78, linewidths=0)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{sample_row['label']} compartment", fontsize=7.8, fontweight="bold")
    panel_label(ax, letter, x=-0.04, y=1.02)


def plot_program_map(ax: plt.Axes, helper, sample_row: dict[str, str], col: str, title: str, letter: str) -> None:
    spots = helper.classify_spot_compartments(helper.load_spots("GSE272362", sample_row["sample_id"]))
    image_path, scale_path = helper.image_paths_for_sample("GSE272362", sample_row["sample_id"], sample_row["geo_base"])
    image = Image.open(image_path).convert("RGB")
    scale = helper.read_scale(scale_path)
    xmin, xmax, ymin, ymax = spatial_limits(spots, scale, image)
    vals = spots[col].astype(float).to_numpy()
    ax.imshow(ImageOps.grayscale(image).convert("RGB"), alpha=0.35)
    ax.scatter(
        spots["x_pixel"].to_numpy(float) * scale,
        spots["y_pixel"].to_numpy(float) * scale,
        c=vals,
        cmap="magma",
        vmin=-1.5,
        vmax=1.5,
        s=5.0,
        linewidths=0,
        alpha=0.78,
    )
    interface = spots["spatial_compartment"].eq("interface")
    ax.scatter(
        spots.loc[interface, "x_pixel"].to_numpy(float) * scale,
        spots.loc[interface, "y_pixel"].to_numpy(float) * scale,
        facecolors="none",
        edgecolors="#55FFFF",
        s=15,
        linewidths=0.35,
        alpha=0.78,
    )
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymax, ymin)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{sample_row['label']} {title}", fontsize=7.8, fontweight="bold")
    panel_label(ax, letter, x=-0.04, y=1.02)


def save(fig: plt.Figure) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    for ext in ["pdf", "svg", "png"]:
        for base in [OUT, ACTIVE_OUT]:
            path = base.with_suffix(f".{ext}")
            if ext == "png":
                fig.savefig(path, dpi=360, bbox_inches="tight", facecolor="white")
            else:
                fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> int:
    flow = load_module("52_make_ecotype_flow_figure.py", "ecotype_flow")
    spatial = load_module("51_make_deep_reference_style_figures.py", "deep_reference_style")

    fig = plt.figure(figsize=(18.2, 16.2))
    gs = GridSpec(5, 1, figure=fig, height_ratios=[1.10, 1.05, 1.0, 1.0, 1.0], hspace=0.52)
    fig.subplots_adjust(top=0.94, bottom=0.075, left=0.045, right=0.985)

    top = gs[0, 0].subgridspec(1, 3, width_ratios=[1.5, 2.25, 1.0], wspace=0.38)
    mid = gs[1, 0].subgridspec(1, 3, width_ratios=[1.35, 2.55, 1.45], wspace=0.55)
    plot_ecotype_loadings(fig.add_subplot(top[0, 0]))
    plot_flow(fig.add_subplot(top[0, 1]), flow)
    plot_composition(fig.add_subplot(top[0, 2]), flow)
    plot_decoupling(fig.add_subplot(mid[0, 0]), flow)
    plot_axis_heatmap(fig.add_subplot(mid[0, 1]))
    plot_axis_correlation(fig.add_subplot(mid[0, 2]))

    sub = gs[2:, 0].subgridspec(3, 4, hspace=0.18, wspace=0.08)
    letters = list("GHIJKLMNOPQR")
    sample_rows = get_spatial_sample_rows()
    map_specs = [
        ("compartment", "compartment", ""),
        ("program", "z_spp1_tam", "SPP1/TAM"),
        ("program", "z_tgfb_pathway", "TGF-beta"),
        ("program", "score_tumor_aggressive", "tumor aggressive"),
    ]
    idx = 0
    for r, sample_row in enumerate(sample_rows):
        for c, spec in enumerate(map_specs):
            ax = fig.add_subplot(sub[r, c])
            if spec[0] == "compartment":
                plot_compartment_map(ax, spatial, sample_row, letters[idx])
            else:
                plot_program_map(ax, spatial, sample_row, spec[1], spec[2], letters[idx])
            idx += 1

    # Shared compact legends for bottom maps.
    legend_ax = fig.add_axes([0.08, 0.038, 0.38, 0.018])
    legend_ax.axis("off")
    comp_items = [("CAF core", "#D73027"), ("tumor-high", "#F0A202"), ("interface", "#7B68A6"), ("other", "#D8D8D8")]
    x = 0.0
    for label, color in comp_items:
        legend_ax.scatter([x], [0.5], s=24, color=color)
        legend_ax.text(x + 0.025, 0.5, label, va="center", fontsize=7)
        x += 0.18
    cax = fig.add_axes([0.58, 0.041, 0.18, 0.012])
    sm = mpl.cm.ScalarMappable(cmap="magma", norm=mpl.colors.Normalize(vmin=-1.5, vmax=1.5))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6)
    cb.set_label("relative program score; cyan outlines mark interface spots", fontsize=7)

    fig.suptitle("CAF-core ecotypes nominate immune-decoupled invasive-interface states", fontsize=15.2, fontweight="bold", y=0.982)
    save(fig)

    panel_index = pd.DataFrame(
        [
            ("A", "CAF-core ecotype loading heatmap"),
            ("B", "Context-to-CAF-core ecotype flow"),
            ("C", "Ecotype composition by context"),
            ("D", "Immune-decoupling index by dominant ecotype"),
            ("E", "Candidate-axis CAF-core/interface enrichment"),
            ("F", "Candidate-axis correlation with immune decoupling"),
            ("G-R", "Primary, liver-metastasis and lymph-node-metastasis interface spatial maps"),
        ],
        columns=["panel", "content"],
    )
    panel_index.to_csv(SOURCE_OUT, index=False)
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {ACTIVE_OUT.with_suffix('.pdf')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
