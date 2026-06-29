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
from matplotlib.ticker import NullFormatter, NullLocator


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "extended_data_figure24_review_risk_resolution_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Extended_Data_Fig_24_review_risk_resolution.csv"
REPORT_OUT = REPORT_DIR / "extended_data_figure24_review_risk_resolution_notes.md"


def clean_axes(ax: plt.Axes) -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)


def panel_label(ax: plt.Axes, label: str, x: float = -0.08, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def draw_evidence_spine(ax: plt.Axes) -> None:
    ax.axis("off")
    panel_label(ax, "A", x=-0.02, y=1.02)
    ax.text(0.02, 0.98, "Reviewer-risk resolution spine", fontsize=10, fontweight="bold", va="top")
    nodes = [
        ("Spatial\nspecificity", "1,000 random\ncores/sample", "#4C78A8"),
        ("Metastatic\ncontrast", "LN immune\nuncoupling", "#72B7B2"),
        ("Mechanism\npriority", "matrix-integrin\nSPP1-CD44", "#B23A48"),
        ("Clinical\ncontext", "TCGA bulk\nsurvival context", "#7B68A6"),
        ("TLS\nboundary", "16/198 stringent\nTLS-compatible", "#8A8F98"),
    ]
    xs = np.linspace(0.08, 0.92, len(nodes))
    for i, (title, subtitle, color) in enumerate(nodes):
        x = xs[i]
        rect = Rectangle((x - 0.075, 0.36), 0.15, 0.34, facecolor="#F4F6F8", edgecolor=color, linewidth=1.2)
        ax.add_patch(rect)
        ax.text(x, 0.62, title, ha="center", va="center", fontsize=8, fontweight="bold", color="#222222")
        ax.text(x, 0.45, subtitle, ha="center", va="center", fontsize=7, color="#444444")
        if i < len(nodes) - 1:
            arrow = FancyArrowPatch((x + 0.08, 0.53), (xs[i + 1] - 0.08, 0.53), arrowstyle="-|>", mutation_scale=9, lw=0.8, color="#666666")
            ax.add_patch(arrow)
    ax.text(
        0.02,
        0.17,
        "Goal: make each vulnerable claim testable, ranked or explicitly bounded.",
        fontsize=7,
        color="#444444",
    )


def draw_mechanism_matrix(ax: plt.Axes, mech: pd.DataFrame) -> None:
    evidence_cols = [
        "evidence_caf_core",
        "evidence_interface_response",
        "evidence_directional_lr",
        "evidence_decoupling_association",
        "evidence_targeted_genes",
        "evidence_xenium_resolution",
        "evidence_tcga_bulk_context",
    ]
    labels = ["Core", "Interface", "Directional", "Decoupling", "Genes", "Xenium", "TCGA"]
    arr = mech.sort_values("priority_rank")[evidence_cols].to_numpy(float)
    im = ax.imshow(arr, vmin=0, vmax=1, cmap=mpl.colors.LinearSegmentedColormap.from_list("ev", ["#F1F2F4", "#A7C7D9", "#1F5F8B"]), aspect="auto")
    ax.set_xticks(np.arange(len(labels)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=6.5)
    ax.set_yticks(np.arange(len(mech)))
    ax.set_yticklabels(mech.sort_values("priority_rank")["candidate_axis"], fontsize=7)
    panel_label(ax, "B")
    ax.set_title("Mechanism evidence layers", loc="left", fontsize=10, fontweight="bold")
    for i in range(arr.shape[0]):
        for j in range(arr.shape[1]):
            token = "++" if arr[i, j] >= 0.95 else ("+" if arr[i, j] > 0 else "")
            ax.text(j, i, token, ha="center", va="center", fontsize=6, color="white" if arr[i, j] > 0.7 else "#333333")
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_ticks([0, 0.5, 1])
    cbar.ax.tick_params(labelsize=6)


def draw_priority_lollipop(ax: plt.Axes, mech: pd.DataFrame) -> None:
    df = mech.sort_values("priority_rank").reset_index(drop=True)
    y = np.arange(len(df))
    colors = np.where(df["lr_p_with_immune_decoupling"] < 0.05, "#B23A48", "#8A8F98")
    ax.axvline(0, color="#333333", lw=0.6)
    ax.hlines(y, 0, df["lr_rho_with_immune_decoupling"], color=colors, lw=1.8)
    ax.scatter(df["lr_rho_with_immune_decoupling"], y, s=42, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["candidate_axis"], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("rho vs immune decoupling", fontsize=7)
    panel_label(ax, "C")
    ax.set_title("Prioritised interface axes", loc="left", fontsize=10, fontweight="bold")
    clean_axes(ax)
    for i, row in df.iterrows():
        label_x = max(row["lr_rho_with_immune_decoupling"] + 0.018, -0.03)
        ax.text(label_x, i, f"score {row['evidence_total_score']:.1f}", va="center", fontsize=6.5)


def draw_survival_forest(ax: plt.Axes, surv: pd.DataFrame) -> None:
    df = surv.sort_values("cox_hr_per_sd", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))
    colors = np.where(df["cox_fdr_bh"] < 0.10, "#B23A48", "#8A8F98")
    ax.axvline(1.0, color="#333333", lw=0.7)
    ax.hlines(y, df["cox_ci95_low"], df["cox_ci95_high"], color=colors, lw=1.5)
    ax.scatter(df["cox_hr_per_sd"], y, s=35, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.75, 1.85)
    ax.set_xticks([0.8, 1.0, 1.25, 1.6])
    ax.set_xticklabels(["0.8", "1.0", "1.25", "1.6"])
    ax.xaxis.set_minor_locator(NullLocator())
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels(df["axis_label"], fontsize=7)
    ax.set_xlabel("TCGA Cox HR per 1 s.d.", fontsize=7)
    panel_label(ax, "D")
    ax.set_title("Clinical context, not prediction", loc="left", fontsize=10, fontweight="bold")
    clean_axes(ax)
    for i, row in df.iterrows():
        ax.text(row["cox_ci95_high"] * 1.06, i, f"FDR {row['cox_fdr_bh']:.2g}", va="center", fontsize=6)


def draw_tls_gate(ax: plt.Axes, tls: pd.DataFrame) -> None:
    order = [
        "post_neoadjuvant_sections",
        "treatment_naive_primary",
        "primary_tumor",
        "liver_metastasis",
        "lymph_node_metastasis",
        "gse274557_primary",
        "gse274557_liver",
        "gse274557_lung",
        "gse274557_peritoneal",
    ]
    labels = {
        "post_neoadjuvant_sections": "post-NACT",
        "treatment_naive_primary": "treat-naive",
        "primary_tumor": "primary",
        "liver_metastasis": "liver met",
        "lymph_node_metastasis": "LN met",
        "gse274557_primary": "G274557 primary",
        "gse274557_liver": "G274557 liver",
        "gse274557_lung": "G274557 lung",
        "gse274557_peritoneal": "G274557 perit.",
    }
    df = tls[tls["cohort_context"].isin(order)].copy()
    df["ord"] = df["cohort_context"].map({v: i for i, v in enumerate(order)})
    df = df.sort_values("ord")
    y = np.arange(len(df))
    ax.barh(y, df["tls_compatible_fraction"], color="#7B68A6", height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels([labels[x] for x in df["cohort_context"]], fontsize=7)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Stringent TLS gate fraction", fontsize=7)
    panel_label(ax, "E")
    ax.set_title("TLS boundary test", loc="left", fontsize=10, fontweight="bold")
    clean_axes(ax)
    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(row["tls_compatible_fraction"] + 0.025, i, f"{int(row['tls_compatible_samples'])}/{int(row['n_samples'])}", va="center", fontsize=6.5)


def draw_claim_scope(ax: plt.Axes) -> None:
    ax.axis("off")
    panel_label(ax, "F", x=-0.02, y=1.02)
    ax.text(0.02, 0.98, "Claim scope after stress tests", fontsize=10, fontweight="bold", va="top")
    rows = [
        ("Spatial organizing core", "Supported", "#2C7A51"),
        ("Metastatic immune decoupling", "Supported", "#2C7A51"),
        ("Perturbation candidates", "Prioritized", "#4C78A8"),
        ("Clinical relevance", "Context only", "#7B68A6"),
        ("Mature TLS", "Not central", "#8A8F98"),
        ("Causal LR signaling", "Not proven", "#8A8F98"),
    ]
    y0 = 0.82
    for i, (claim, scope, color) in enumerate(rows):
        y = y0 - i * 0.12
        ax.add_patch(Rectangle((0.02, y - 0.04), 0.96, 0.075, facecolor="#F6F7F9", edgecolor="#DDDDDD", linewidth=0.5))
        ax.add_patch(Rectangle((0.02, y - 0.04), 0.025, 0.075, facecolor=color, edgecolor=color, linewidth=0))
        ax.text(0.065, y, claim, va="center", fontsize=7.5, color="#222222")
        ax.text(0.78, y, scope, va="center", ha="left", fontsize=7.5, fontweight="bold", color=color)
    ax.text(0.02, 0.05, "Use this as a boundary map for cover letter, rebuttal and figure legends.", fontsize=6.8, color="#444444")


def write_source(mech: pd.DataFrame, surv: pd.DataFrame, tls: pd.DataFrame) -> None:
    rows = []
    for _, row in mech.iterrows():
        rows.append({"panel": "B-C", "source_table": "mechanism_triangulation_priority_matrix.csv", "item": row["candidate_axis"], "value": row["evidence_total_score"]})
    for _, row in surv.iterrows():
        rows.append({"panel": "D", "source_table": "tcga_paad_survival_context_summary.csv", "item": row["axis_label"], "value": row["cox_hr_per_sd"]})
    for _, row in tls.iterrows():
        rows.append({"panel": "E", "source_table": "tls_maturity_stress_test_context_summary.csv", "item": row["cohort_context"], "value": row["tls_compatible_fraction"]})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    lines = [
        "# Extended Data Figure 24 Notes",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Figure role",
        "",
        "NC-style module figure that consolidates reviewer-risk resolution across mechanism, clinical-context and TLS-boundary analyses.",
        "",
        "## Panel contract",
        "",
        "- A: shows the evidence spine from spatial specificity to claim boundary.",
        "- B-C: show that candidate mechanisms are ranked rather than overclaimed as causal.",
        "- D: shows clinical context without turning the manuscript into a prognostic paper.",
        "- E: shows that mature TLS is a tested but non-central interpretation.",
        "- F: gives the permitted claim scope.",
        "",
        "## Outputs",
        "",
        f"- `{OUT.with_suffix('.pdf').relative_to(ROOT)}`",
        f"- `{OUT.with_suffix('.svg').relative_to(ROOT)}`",
        f"- `{OUT.with_suffix('.png').relative_to(ROOT)}`",
        f"- `{SOURCE_OUT.relative_to(ROOT)}`",
    ]
    REPORT_OUT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    mech = pd.read_csv(TABLE_DIR / "mechanism_triangulation_priority_matrix.csv")
    surv = pd.read_csv(TABLE_DIR / "tcga_paad_survival_context_summary.csv")
    tls = pd.read_csv(TABLE_DIR / "tls_maturity_stress_test_context_summary.csv")

    fig = plt.figure(figsize=(12.4, 9.0))
    gs = GridSpec(3, 3, figure=fig, height_ratios=[0.85, 1.05, 1.05], width_ratios=[1.25, 1.0, 1.0], wspace=0.55, hspace=0.68)
    draw_evidence_spine(fig.add_subplot(gs[0, :]))
    draw_mechanism_matrix(fig.add_subplot(gs[1, 0]), mech)
    draw_priority_lollipop(fig.add_subplot(gs[1, 1]), mech)
    draw_survival_forest(fig.add_subplot(gs[1, 2]), surv)
    draw_tls_gate(fig.add_subplot(gs[2, 0:2]), tls)
    draw_claim_scope(fig.add_subplot(gs[2, 2]))
    fig.suptitle("Reviewer-risk resolution module for the CAF-myeloid spatial niche model", fontsize=12, fontweight="bold", y=0.992)
    fig.text(
        0.01,
        0.012,
        "This synthesis figure consolidates existing source analyses; it defines claim scope and does not add causal, prognostic or mature-TLS claims.",
        fontsize=7,
        color="#444444",
    )
    for suffix in [".pdf", ".svg", ".png"]:
        fig.savefig(OUT.with_suffix(suffix), dpi=450 if suffix == ".png" else None, bbox_inches="tight")
    plt.close(fig)
    write_source(mech, surv, tls)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
