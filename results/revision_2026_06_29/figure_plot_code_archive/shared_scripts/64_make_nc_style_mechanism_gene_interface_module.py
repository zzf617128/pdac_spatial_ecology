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
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "extended_data_figure29_mechanism_gene_interface_module_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Extended_Data_Fig_29_mechanism_gene_interface_module.csv"
REPORT_OUT = REPORT_DIR / "extended_data_figure29_mechanism_gene_interface_module_notes.md"

CANDIDATES = ["matrix-integrin", "SPP1-CD44/integrin", "TGF-beta/TGFBR", "IL6-OSM/LIF-JAKSTAT"]
CONTEXTS = [
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
PROGRAM_FOR_CANDIDATE = {
    "matrix-integrin": "SPP1-TAM/matrix",
    "SPP1-CD44/integrin": "SPP1-TAM/matrix",
    "TGF-beta/TGFBR": "TGF-beta/EMT invasive",
    "IL6-OSM/LIF-JAKSTAT": "IFN/APC antigen",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def plot_evidence_matrix(ax: plt.Axes, mech: pd.DataFrame) -> pd.DataFrame:
    evidence_cols = [
        "evidence_caf_core",
        "evidence_interface_response",
        "evidence_directional_lr",
        "evidence_decoupling_association",
        "evidence_targeted_genes",
        "evidence_xenium_resolution",
        "evidence_tcga_bulk_context",
    ]
    labels = ["CAF core", "Interface", "Directional", "Decoupling", "Genes", "Xenium", "TCGA"]
    df = mech.set_index("candidate_axis").reindex(CANDIDATES).reset_index()
    arr = df[evidence_cols].to_numpy(float)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("evidence", ["#F1F2F4", "#A8C6D8", "#1F5F8B"])
    im = ax.imshow(arr, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["candidate_axis"], fontsize=7)
    ax.set_title("Triangulated evidence layers", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "A")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            token = "++" if arr[i, j] >= 0.95 else ("+" if arr[i, j] > 0 else "")
            ax.text(j, i, token, ha="center", va="center", fontsize=6, color="white" if arr[i, j] > 0.7 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.tick_params(labelsize=6)
    return df[["candidate_axis", *evidence_cols, "evidence_total_score", "priority_rank"]]


def plot_targeted_gene_heatmap(ax: plt.Axes, targeted: pd.DataFrame, metric: str, title: str, panel: str) -> pd.DataFrame:
    rows = []
    for cand in CANDIDATES:
        prog = PROGRAM_FOR_CANDIDATE[cand]
        sub = targeted[targeted["axis_label"].eq(prog) & targeted["cohort_context"].isin(CONTEXTS)].copy()
        for _, row in sub.iterrows():
            rows.append({"candidate_axis": cand, "cohort_context": row["cohort_context"], "program_axis": prog, "value": row[metric], "n_positive": row["n_core_positive"] if metric == "median_core_enrichment" else row["n_interface_positive"], "n_samples": row["n_samples"]})
    df = pd.DataFrame(rows)
    mat = df.pivot_table(index="candidate_axis", columns="cohort_context", values="value", aggfunc="median").reindex(index=CANDIDATES, columns=CONTEXTS)
    im = ax.imshow(mat.to_numpy(float), cmap="YlGnBu", vmin=-0.35 if metric.endswith("interface_enrichment") else -0.2, vmax=1.35, aspect="auto")
    ax.set_xticks(np.arange(len(CONTEXTS)))
    ax.set_xticklabels([CONTEXT_LABELS[c] for c in CONTEXTS], rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES, fontsize=7)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, panel)
    for i, cand in enumerate(CANDIDATES):
        for j, ctx in enumerate(CONTEXTS):
            val = mat.loc[cand, ctx]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if val > 0.75 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("median enrichment", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df.assign(metric=metric)


def plot_lr_metrics(ax: plt.Axes, ctx: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "ligand_core_enrichment",
        "receptor_interface_enrichment",
        "response_interface_enrichment",
        "directional_core_to_interface_score",
    ]
    labels = ["ligand core", "receptor iface", "response iface", "directional score"]
    df = ctx[ctx["axis"].isin(CANDIDATES) & ctx["metric"].isin(metrics)].copy()
    summary = df.groupby(["axis", "metric"], as_index=False).agg(value=("median_value", "median"), support=("n_positive", "sum"), n=("n_samples", "sum"))
    mat = summary.pivot_table(index="axis", columns="metric", values="value", aggfunc="median").reindex(index=CANDIDATES, columns=metrics)
    im = ax.imshow(mat.to_numpy(float), cmap="YlOrRd", vmin=-0.05, vmax=1.55, aspect="auto")
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES, fontsize=7)
    ax.set_title("Focused LR/interface metrics", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "D")
    for i, cand in enumerate(CANDIDATES):
        for j, metric in enumerate(metrics):
            val = mat.loc[cand, metric]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if val > 0.95 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("median value", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return summary


def plot_decoupling_correlations(ax: plt.Axes, corr: pd.DataFrame) -> pd.DataFrame:
    metric = "directional_core_to_interface_score"
    df = corr[corr["axis"].isin(CANDIDATES) & corr["metric"].eq(metric) & corr["target"].isin(["immune_decoupling_index", "stromal_tumor_core_coupling"])].copy()
    targets = ["immune_decoupling_index", "stromal_tumor_core_coupling"]
    colors = {"immune_decoupling_index": "#B23A48", "stromal_tumor_core_coupling": "#4C78A8"}
    offsets = {"immune_decoupling_index": -0.13, "stromal_tumor_core_coupling": 0.13}
    y_base = np.arange(len(CANDIDATES))
    ax.axvline(0, color="#333333", lw=0.65)
    for target in targets:
        sub = df[df["target"].eq(target)].set_index("axis").reindex(CANDIDATES)
        y = y_base + offsets[target]
        vals = sub["spearman_rho"].to_numpy(float)
        ax.hlines(y, 0, vals, color=colors[target], lw=1.6)
        ax.scatter(vals, y, color=colors[target], s=34, edgecolor="white", linewidth=0.5, label=target.replace("_", " "))
    ax.set_yticks(y_base)
    ax.set_yticklabels(CANDIDATES, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(-0.25, 0.85)
    ax.set_xlabel("Spearman rho", fontsize=7)
    ax.set_title("Directional score correlations", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax, axis="x")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    return df


def plot_priority(ax: plt.Axes, mech: pd.DataFrame) -> pd.DataFrame:
    df = mech.set_index("candidate_axis").reindex(CANDIDATES).reset_index()
    y = np.arange(len(df))
    colors = np.where(df["priority_rank"] <= 2, "#B23A48", "#8A8F98")
    ax.barh(y, df["evidence_total_score"], color=colors, height=0.55)
    ax.set_yticks(y)
    ax.set_yticklabels(df["candidate_axis"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 7.4)
    ax.set_xlabel("triangulated evidence score", fontsize=7)
    ax.set_title("Perturbation-priority ranking", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "F")
    clean_axes(ax)
    for i, row in df.iterrows():
        ax.text(row["evidence_total_score"] + 0.08, i, f"rank {int(row['priority_rank'])}", va="center", fontsize=6.5)
    return df[["candidate_axis", "evidence_total_score", "priority_rank", "xenium_median_delta_vs_random", "tcga_stromal_myeloid_rho"]]


def draw_mechanism_schema(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "G", x=-0.02, y=1.02)
    ax.text(0.02, 0.98, "Mechanism model for follow-up", fontsize=10, fontweight="bold", va="top")
    boxes = [
        (0.05, 0.50, 0.24, 0.22, "CAF-myeloid\ncore", "#2C7A51"),
        (0.38, 0.50, 0.24, 0.22, "tumor-stroma\ninterface", "#B23A48"),
        (0.71, 0.50, 0.24, 0.22, "immune/tumor\nprograms", "#4C78A8"),
    ]
    for x, y, w, h, label, color in boxes:
        ax.add_patch(Rectangle((x, y), w, h, facecolor="#F4F6F8", edgecolor=color, linewidth=1.2))
        ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=8, fontweight="bold", color="#222222")
    arrows = [
        ((0.29, 0.61), (0.38, 0.61), "matrix-integrin\nSPP1-CD44"),
        ((0.62, 0.61), (0.71, 0.61), "TGF-beta/EMT\nresponse"),
    ]
    for start, end, label in arrows:
        ax.add_patch(FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=10, lw=0.9, color="#555555"))
        ax.text((start[0] + end[0]) / 2, 0.75, label, ha="center", va="center", fontsize=6.5, color="#333333")
    ax.text(0.05, 0.25, "Use as candidate biology for perturbation, not observational causal proof.", fontsize=7, color="#555555")
    return pd.DataFrame(
        [
            {"item": "CAF-myeloid core", "role": "ligand-rich core"},
            {"item": "tumor-stroma interface", "role": "receptor/response interface"},
            {"item": "immune/tumor programs", "role": "organized downstream programs"},
        ]
    )


def write_source_data(*tables: pd.DataFrame) -> None:
    names = ["A_evidence", "B_targeted_core", "C_targeted_interface", "D_lr_metrics", "E_correlations", "F_priority", "G_schema"]
    rows: list[dict[str, object]] = []
    for name, table in zip(names, tables):
        for _, row in table.iterrows():
            item = row.get("candidate_axis", row.get("axis", row.get("item", row.get("program_axis", ""))))
            for col, value in row.items():
                if col in {"candidate_axis", "axis", "item", "program_axis", "rationale"}:
                    continue
                if isinstance(value, (int, float, np.number)) and np.isfinite(value):
                    rows.append({"panel": name, "item": item, "metric": col, "value": float(value)})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Extended Data Figure 29 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "NC-style mechanism gene/interface module. This figure deepens the perturbation-priority story by consolidating targeted-gene support, focused ligand/receptor/interface metrics, immune-decoupling associations and external context.\n\n"
        "## Panel contract\n\n"
        "- A: triangulated evidence layers for candidate axes.\n"
        "- B-C: targeted-gene CAF-core and interface enrichment across contexts.\n"
        "- D: focused ligand-core, receptor-interface, response-interface and directional metrics.\n"
        "- E: directional-score associations with immune decoupling and stromal-tumor coupling.\n"
        "- F: perturbation-priority ranking.\n"
        "- G: compact follow-up model.\n\n"
        "## Boundary\n\n"
        "This module ranks perturbation-ready candidates. It does not establish causal ligand-receptor signaling from observational spatial transcriptomics.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf')}`\n"
        f"- `{OUT.with_suffix('.svg')}`\n"
        f"- `{OUT.with_suffix('.png')}`\n"
        f"- `{SOURCE_OUT}`\n",
        encoding="utf-8",
    )


def main() -> None:
    mech = pd.read_csv(TABLE_DIR / "mechanism_triangulation_priority_matrix.csv")
    targeted = pd.read_csv(TABLE_DIR / "targeted_gene_axis_validation_summary.csv")
    ctx = pd.read_csv(TABLE_DIR / "gap3_focused_lr_interface_context_summary.csv")
    corr = pd.read_csv(TABLE_DIR / "gap3_focused_lr_interface_correlations.csv")

    fig = plt.figure(figsize=(14.8, 13.4), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[1.0, 1.0, 0.95, 0.58], hspace=0.92, wspace=0.88)
    fig.suptitle("Mechanism gene/interface module for perturbation-ready axes", fontsize=15, fontweight="bold", y=0.986)

    ax_a = fig.add_subplot(gs[0, 0:3])
    ax_b = fig.add_subplot(gs[0, 3:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:6])
    ax_e = fig.add_subplot(gs[2, 0:2])
    ax_f = fig.add_subplot(gs[2, 2:4])
    ax_g = fig.add_subplot(gs[2, 4:6])
    ax_note = fig.add_subplot(gs[3, :])

    a = plot_evidence_matrix(ax_a, mech)
    b = plot_targeted_gene_heatmap(ax_b, targeted, "median_core_enrichment", "Targeted-gene CAF-core enrichment", "B")
    c = plot_targeted_gene_heatmap(ax_c, targeted, "median_interface_enrichment", "Targeted-gene interface enrichment", "C")
    d = plot_lr_metrics(ax_d, ctx)
    e = plot_decoupling_correlations(ax_e, corr)
    f = plot_priority(ax_f, mech)
    g = draw_mechanism_schema(ax_g)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.72,
        "Interpretation: matrix-integrin and SPP1-CD44/integrin have the broadest triangulated support, including CAF-core ligand signal, interface response, directional core-to-interface structure, targeted-gene support, Xenium consistency and TCGA context. TGF-beta/TGFBR remains a strong invasive-interface follow-up axis.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.28,
        "Boundary: this module converts observational spatial evidence into a ranked perturbation agenda; it does not claim causal ligand-receptor signaling.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    write_source_data(a, b, c, d, e, f, g)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
