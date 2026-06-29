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

OUT = FIG_DIR / "supplementary_module4_mechanism_interface_priority"
SOURCE_OUT = SOURCE_DIR / "Source_Data_supplementary_module4_mechanism_interface_priority.csv"
REPORT_OUT = REPORT_DIR / "supplementary_module4_mechanism_interface_priority_notes.md"

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
EVIDENCE_COLS = [
    "evidence_caf_core",
    "evidence_interface_response",
    "evidence_directional_lr",
    "evidence_decoupling_association",
    "evidence_targeted_genes",
    "evidence_xenium_resolution",
    "evidence_tcga_bulk_context",
]
EVIDENCE_LABELS = ["CAF core", "Interface", "Directional", "Decoupling", "Genes", "Xenium", "TCGA"]


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def ordered_mech(mech: pd.DataFrame) -> pd.DataFrame:
    return mech.set_index("candidate_axis").reindex(CANDIDATES).reset_index()


def plot_evidence_matrix(ax: plt.Axes, mech: pd.DataFrame) -> pd.DataFrame:
    df = ordered_mech(mech)
    arr = df[EVIDENCE_COLS].to_numpy(float)
    cmap = mpl.colors.LinearSegmentedColormap.from_list("evidence", ["#F1F2F4", "#A8C6D8", "#1F5F8B"])
    im = ax.imshow(arr, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(EVIDENCE_LABELS)))
    ax.set_xticklabels(EVIDENCE_LABELS, rotation=32, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(df)))
    ax.set_yticklabels(df["candidate_axis"], fontsize=7.5)
    ax.set_title("Triangulated evidence layers", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "A")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            token = "++" if arr[i, j] >= 0.95 else ("+" if arr[i, j] > 0 else "")
            ax.text(j, i, token, ha="center", va="center", fontsize=6.3, color="white" if arr[i, j] > 0.7 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.034, pad=0.018)
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.tick_params(labelsize=6)
    return df[["candidate_axis", *EVIDENCE_COLS, "evidence_total_score", "priority_rank"]]


def plot_priority(ax: plt.Axes, mech: pd.DataFrame) -> pd.DataFrame:
    df = ordered_mech(mech)
    y = np.arange(len(df))
    colors = np.where(df["priority_rank"] <= 2, "#B23A48", "#8A8F98")
    ax.barh(y, df["evidence_total_score"], color=colors, height=0.56, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels(df["candidate_axis"], fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(0, 7.5)
    ax.set_xlabel("triangulated score / 7", fontsize=7.5)
    ax.set_title("Perturbation-priority rank", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "B")
    clean_axes(ax)
    for i, row in df.iterrows():
        ax.text(row["evidence_total_score"] + 0.10, i, f"rank {int(row['priority_rank'])}", va="center", fontsize=6.7)
    return df[["candidate_axis", "evidence_total_score", "priority_rank"]]


def targeted_gene_table(targeted: pd.DataFrame, metric: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for cand in CANDIDATES:
        prog = PROGRAM_FOR_CANDIDATE[cand]
        sub = targeted[targeted["axis_label"].eq(prog) & targeted["cohort_context"].isin(CONTEXTS)]
        for _, row in sub.iterrows():
            positive_col = "n_core_positive" if metric == "median_core_enrichment" else "n_interface_positive"
            rows.append(
                {
                    "candidate_axis": cand,
                    "program_axis": prog,
                    "cohort_context": row["cohort_context"],
                    "value": row[metric],
                    "n_positive": row[positive_col],
                    "n_samples": row["n_samples"],
                }
            )
    return pd.DataFrame(rows)


def plot_targeted_heatmap(
    ax: plt.Axes,
    targeted: pd.DataFrame,
    metric: str,
    title: str,
    panel: str,
    show_y_labels: bool = True,
) -> pd.DataFrame:
    df = targeted_gene_table(targeted, metric)
    mat = df.pivot_table(index="candidate_axis", columns="cohort_context", values="value", aggfunc="median").reindex(
        index=CANDIDATES, columns=CONTEXTS
    )
    cmap = plt.get_cmap("YlGnBu").copy()
    cmap.set_bad("#F4F4F4")
    im = ax.imshow(mat.to_numpy(float), cmap=cmap, vmin=-0.35 if "interface" in metric else -0.2, vmax=1.35, aspect="auto")
    ax.set_xticks(np.arange(len(CONTEXTS)))
    ax.set_xticklabels([CONTEXT_LABELS[c] for c in CONTEXTS], rotation=32, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES if show_y_labels else [], fontsize=7.5)
    ax.set_title(title, loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, panel)
    for i, cand in enumerate(CANDIDATES):
        for j, ctx in enumerate(CONTEXTS):
            val = mat.loc[cand, ctx]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.9, color="#FFFFFF" if val > 0.75 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return df.assign(metric=metric)


def plot_lr_metrics(ax: plt.Axes, ctx: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "ligand_core_enrichment",
        "receptor_interface_enrichment",
        "response_interface_enrichment",
        "directional_core_to_interface_score",
    ]
    labels = ["ligand\ncore", "receptor\ninterface", "response\ninterface", "directional\nscore"]
    df = ctx[ctx["axis"].isin(CANDIDATES) & ctx["metric"].isin(metrics)].copy()
    summary = df.groupby(["axis", "metric"], as_index=False).agg(
        value=("median_value", "median"),
        support=("n_positive", "sum"),
        n=("n_samples", "sum"),
    )
    mat = summary.pivot_table(index="axis", columns="metric", values="value", aggfunc="median").reindex(index=CANDIDATES, columns=metrics)
    im = ax.imshow(mat.to_numpy(float), cmap="YlOrRd", vmin=-0.05, vmax=1.55, aspect="auto")
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_yticks(np.arange(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES, fontsize=7.5)
    ax.set_title("Focused ligand-response interface metrics", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "E")
    for i, cand in enumerate(CANDIDATES):
        for j, metric in enumerate(metrics):
            val = mat.loc[cand, metric]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.9, color="#FFFFFF" if val > 0.95 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return summary


def plot_directional_correlations(ax: plt.Axes, corr: pd.DataFrame) -> pd.DataFrame:
    metric = "directional_core_to_interface_score"
    targets = ["immune_decoupling_index", "stromal_tumor_core_coupling"]
    df = corr[corr["axis"].isin(CANDIDATES) & corr["metric"].eq(metric) & corr["target"].isin(targets)].copy()
    colors = {"immune_decoupling_index": "#B23A48", "stromal_tumor_core_coupling": "#4C78A8"}
    offsets = {"immune_decoupling_index": -0.14, "stromal_tumor_core_coupling": 0.14}
    y_base = np.arange(len(CANDIDATES))
    ax.axvline(0, color="#333333", lw=0.65)
    for target in targets:
        sub = df[df["target"].eq(target)].set_index("axis").reindex(CANDIDATES)
        vals = sub["spearman_rho"].to_numpy(float)
        y = y_base + offsets[target]
        ax.hlines(y, 0, vals, color=colors[target], lw=1.6, zorder=2)
        sig = sub["p_value"].to_numpy(float) < 0.05
        ax.scatter(vals[sig], y[sig], color=colors[target], s=34, edgecolor="white", linewidth=0.55, zorder=3)
        ax.scatter(vals[~sig], y[~sig], facecolor="white", edgecolor=colors[target], s=34, linewidth=1.0, zorder=3)
    ax.set_yticks(y_base)
    ax.set_yticklabels(CANDIDATES, fontsize=7.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.25, 0.85)
    ax.set_xlabel("Spearman rho", fontsize=7.5)
    ax.set_title("Directional score association", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "F")
    clean_axes(ax)
    ax.text(0.02, -0.20, "filled: P < 0.05", transform=ax.transAxes, fontsize=6.5, color="#555555", va="top")
    ax.text(0.55, -0.20, "red: immune decoupling", transform=ax.transAxes, fontsize=6.5, color="#B23A48", va="top")
    ax.text(0.55, -0.33, "blue: stromal-tumor coupling", transform=ax.transAxes, fontsize=6.5, color="#4C78A8", va="top")
    return df


def plot_external_context(ax: plt.Axes, mech: pd.DataFrame) -> pd.DataFrame:
    cols = ["xenium_median_delta_vs_random", "tcga_stromal_myeloid_rho", "tcga_decoupling_like_rho"]
    labels = ["Xenium\ndelta", "TCGA stromal-\nmyeloid rho", "TCGA decoupling-\nlike rho"]
    df = ordered_mech(mech)
    mat = df[cols].to_numpy(float)
    vmax = max(0.95, np.nanmax(np.abs(mat)))
    im = ax.imshow(mat, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_yticks(np.arange(len(CANDIDATES)))
    ax.set_yticklabels(CANDIDATES, fontsize=7.5)
    ax.set_title("External-context anchors", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "G")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat[i, j]
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.9, color="#FFFFFF" if abs(val) > 0.55 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return df[["candidate_axis", *cols]]


def draw_perturbation_agenda(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "H", x=-0.03, y=1.02)
    ax.text(0.02, 0.97, "Mechanism claim ladder", fontsize=10.5, fontweight="bold", va="top")
    rows = [
        ("Observed", "CAF-myeloid cores align with interface programs", "#1F5F8B"),
        ("Prioritized", "matrix-integrin and SPP1-CD44/integrin lead ranking", "#B23A48"),
        ("Next test", "block ligand/receptor axes in spatial co-culture or organoid assays", "#2C7A51"),
    ]
    y0 = 0.75
    for i, (stage, claim, color) in enumerate(rows):
        y = y0 - i * 0.22
        ax.add_patch(Rectangle((0.04, y - 0.07), 0.22, 0.12, facecolor="#F4F6F8", edgecolor=color, linewidth=1.0))
        ax.text(0.15, y - 0.01, stage, ha="center", va="center", fontsize=7.2, fontweight="bold", color="#222222")
        ax.text(0.32, y - 0.01, claim, ha="left", va="center", fontsize=7.1, color="#333333")
        if i < len(rows) - 1:
            ax.add_patch(FancyArrowPatch((0.15, y - 0.095), (0.15, y - 0.15), arrowstyle="-|>", mutation_scale=8, lw=0.8, color="#777777"))
    ax.text(0.04, 0.05, "Boundary: ranked perturbation agenda, not causal proof from observation alone.", fontsize=6.8, color="#555555")
    return pd.DataFrame(
        [
            {"item": "Observed", "role": "spatial association and interface alignment"},
            {"item": "Prioritized", "role": "ranked candidate axes"},
            {"item": "Next test", "role": "experimental perturbation requirement"},
        ]
    )


def write_source_data(*tables: pd.DataFrame) -> None:
    names = [
        "A_evidence",
        "B_priority",
        "C_targeted_core",
        "D_targeted_interface",
        "E_lr_metrics",
        "F_directional_correlations",
        "G_external_context",
        "H_claim_ladder",
    ]
    rows: list[dict[str, object]] = []
    for name, table in zip(names, tables):
        for _, row in table.iterrows():
            item = row.get("candidate_axis", row.get("axis", row.get("item", row.get("program_axis", ""))))
            for col, value in row.items():
                if col in {"candidate_axis", "axis", "item", "program_axis", "role", "rationale"}:
                    continue
                if isinstance(value, (int, float, np.number)) and np.isfinite(value):
                    rows.append({"panel": name, "item": item, "metric": col, "value": float(value)})
                elif col in {"cohort_context", "metric", "target"} and isinstance(value, str):
                    rows.append({"panel": name, "item": item, "metric": col, "value": value})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Supplementary Module 4 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "Module 4 converts the observational spatial ecology results into a perturbation-ready mechanism agenda. "
        "It does not claim causal signaling; it ranks candidate axes using independent evidence layers.\n\n"
        "## Panel contract\n\n"
        "- A: triangulated evidence layers for four candidate axes.\n"
        "- B: total evidence score and perturbation-priority rank.\n"
        "- C-D: targeted-gene CAF-core and interface enrichment across spatial contexts.\n"
        "- E: focused ligand-core, receptor-interface, response-interface and directional metrics.\n"
        "- F: association of directional scores with immune decoupling and stromal-tumor coupling.\n"
        "- G: external-context anchors from Xenium-resolution and TCGA-level summaries.\n"
        "- H: claim ladder defining what is observed, prioritized and still requires perturbation.\n\n"
        "## Main conclusion\n\n"
        "matrix-integrin and SPP1-CD44/integrin are the leading perturbation-ready axes, while TGF-beta/TGFBR remains "
        "a strong invasive-interface follow-up axis.\n\n"
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

    fig = plt.figure(figsize=(15.2, 14.2), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[0.95, 1.0, 1.0, 0.78], hspace=0.96, wspace=1.12)
    fig.suptitle("Supplementary Module 4 | Mechanism, interface biology and perturbation-priority axes", fontsize=15, fontweight="bold", y=0.988)

    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 4:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:6])
    ax_e = fig.add_subplot(gs[2, 0:3])
    ax_f = fig.add_subplot(gs[2, 3:6])
    ax_g = fig.add_subplot(gs[3, 0:3])
    ax_h = fig.add_subplot(gs[3, 3:6])

    a = plot_evidence_matrix(ax_a, mech)
    b = plot_priority(ax_b, mech)
    c = plot_targeted_heatmap(ax_c, targeted, "median_core_enrichment", "Targeted-gene CAF-core enrichment", "C")
    d = plot_targeted_heatmap(
        ax_d,
        targeted,
        "median_interface_enrichment",
        "Targeted-gene interface enrichment",
        "D",
        show_y_labels=False,
    )
    e = plot_lr_metrics(ax_e, ctx)
    f = plot_directional_correlations(ax_f, corr)
    g = plot_external_context(ax_g, mech)
    h = draw_perturbation_agenda(ax_h)

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    write_source_data(a, b, c, d, e, f, g, h)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
