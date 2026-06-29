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
from matplotlib.patches import Circle, FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "supplementary_module6_spatial_architecture_mechanism_deepening"
SOURCE_OUT = SOURCE_DIR / "Source_Data_supplementary_module6_spatial_architecture_mechanism_deepening.csv"
REPORT_OUT = REPORT_DIR / "supplementary_module6_spatial_architecture_mechanism_deepening_notes.md"

CONTEXT_ORDER = ["post-NACT", "primary tumor", "liver metastasis", "lymph-node metastasis"]
CONTEXT_COLORS = {
    "post-NACT": "#4C78A8",
    "primary tumor": "#2C7A51",
    "liver metastasis": "#B23A48",
    "lymph-node metastasis": "#7B68A6",
}
PROGRAMS = ["SPP1/TAM", "TGFb/EMT", "tumor aggressive", "tumor epithelial"]
PROGRAM_COLORS = {
    "SPP1/TAM": "#B23A48",
    "TGFb/EMT": "#D88928",
    "tumor aggressive": "#7B68A6",
    "tumor epithelial": "#4C78A8",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.10, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "x") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def evidence_chain(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "A", x=-0.03)
    ax.text(0.02, 0.97, "Deep spatial architecture evidence chain", fontsize=10.5, fontweight="bold", va="top")
    nodes = [
        ("Anchor\nspecificity", "not a generic\nhigh-signal mask"),
        ("Cell-level\nneighborhoods", "Xenium CAF-myeloid\nproximity"),
        ("Core-to-interface\ngradient", "programs vary along\na tissue axis"),
        ("CAF-core\ngeometry", "architecture tunes\ncoupling strength"),
    ]
    xs = np.linspace(0.09, 0.91, len(nodes))
    y = 0.55
    colors = ["#4C78A8", "#B23A48", "#D88928", "#2C7A51"]
    for i, ((title, subtitle), x, color) in enumerate(zip(nodes, xs, colors)):
        ax.add_patch(Circle((x, y), 0.075, facecolor=color, edgecolor="white", linewidth=1.0))
        ax.text(x, y + 0.004, title, ha="center", va="center", fontsize=6.6, color="white", fontweight="bold")
        ax.text(x, y - 0.17, subtitle, ha="center", va="top", fontsize=5.9, color="#333333", linespacing=1.12)
        if i < len(nodes) - 1:
            ax.add_patch(
                FancyArrowPatch(
                    (x + 0.082, y),
                    (xs[i + 1] - 0.082, y),
                    arrowstyle="-|>",
                    mutation_scale=10,
                    linewidth=1.0,
                    color="#777777",
                )
            )
    ax.text(0.02, 0.08, "Claim level: spatial association and prioritization; perturbation is still required for causality.", fontsize=6.7, color="#555555")
    return pd.DataFrame(
        [
            {"step": title.replace("\n", " "), "interpretation": subtitle}
            for title, subtitle in nodes
        ]
    )


def anchor_specificity(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "alternative_biological_anchor_specificity_summary.csv")
    df = df[df["anchor"].eq("CAF-myeloid core")].copy()
    order = ["TGFb/EMT", "SPP1/TAM", "IFN/MHC", "tumor aggressive", "immune core", "tumor epithelial"]
    df["target_program"] = pd.Categorical(df["target_program"], categories=order, ordered=True)
    df = df.sort_values("target_program")
    y = np.arange(len(df))
    colors = mpl.colormaps["Blues"](0.25 + 0.65 * df["support_fraction"].to_numpy(float))
    ax.axvline(0, color="#333333", lw=0.65)
    ax.barh(y, df["median_delta"], color=colors, height=0.58, zorder=2)
    ax.scatter(df["median_observed_rho"], y, s=28, color="#B23A48", edgecolor="white", linewidth=0.45, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["target_program"].astype(str), fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlim(-0.56, 0.08)
    ax.set_xlabel("observed-minus-random spatial rho", fontsize=7.2, labelpad=4)
    ax.set_title("CAF-myeloid core specificity", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "B")
    clean_axes(ax)
    for i, row in df.reset_index(drop=True).iterrows():
        ax.text(0.075, i, f"{row['support_fraction']:.0%}", ha="right", va="center", fontsize=6.2, color="#333333")
    ax.text(0.02, -0.29, "red dot: observed rho; bar: observed - random", transform=ax.transAxes, fontsize=6.2, color="#555555")
    return df


def anchor_neighborhoods(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse274673_xenium_anchor_neighborhood_summary.csv")
    keep = ["CAF/matrix", "SPP1/TAM", "IFN/APC", "T/NK", "TGFb/EMT", "tumor epithelial"]
    df = df[df["neighbor_state"].isin(keep)].copy()
    anchors = ["CAF-SPP1/TAM domain", "CAF-APC domain"]
    y_base = np.arange(len(keep))
    offsets = {"CAF-SPP1/TAM domain": -0.18, "CAF-APC domain": 0.18}
    colors = {"CAF-SPP1/TAM domain": "#B23A48", "CAF-APC domain": "#4C78A8"}
    ax.axvline(0, color="#333333", lw=0.65)
    for anchor in anchors:
        sub = df[df["anchor"].eq(anchor)].set_index("neighbor_state").reindex(keep)
        y = y_base + offsets[anchor]
        ax.hlines(y, 0, sub["median_delta_fraction"], color=colors[anchor], lw=1.4, zorder=2)
        ax.scatter(
            sub["median_delta_fraction"],
            y,
            s=24 + 36 * sub["support_fraction"].fillna(0).to_numpy(float),
            color=colors[anchor],
            edgecolor="white",
            linewidth=0.45,
            label=anchor.replace(" domain", ""),
            zorder=3,
        )
    ax.set_yticks(y_base)
    ax.set_yticklabels(keep, fontsize=7.2)
    ax.invert_yaxis()
    ax.set_xlim(-0.22, 0.14)
    ax.set_xlabel("neighborhood fraction delta", fontsize=7.2)
    ax.set_title("Xenium CAF-domain neighborhoods", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "C")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=6.2, loc="lower right")
    return df


def xenium_network(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse274673_xenium_cell_state_adjacency_summary.csv")
    df = df[df["state_a"].ne(df["state_b"]) & (df["support_fraction"] >= 0.66)].copy()
    df = df.sort_values("median_log2_oe", ascending=False).head(8)
    states = ["CAF/matrix", "TGFb/EMT", "SPP1/TAM", "IFN/APC", "T/NK", "tumor epithelial"]
    angles = np.linspace(0, 2 * np.pi, len(states), endpoint=False) + np.pi / 8
    pos = {state: (0.5 + 0.34 * np.cos(a), 0.52 + 0.32 * np.sin(a)) for state, a in zip(states, angles)}
    ax.axis("off")
    panel_label(ax, "D", x=-0.03)
    ax.text(0.02, 0.97, "Cell-state adjacency network", fontsize=10.5, fontweight="bold", va="top")
    vmax = max(0.1, float(df["median_log2_oe"].max()))
    for _, row in df.iterrows():
        a, b = row["state_a"], row["state_b"]
        if a not in pos or b not in pos:
            continue
        x1, y1 = pos[a]
        x2, y2 = pos[b]
        color = "#B23A48" if row["median_log2_oe"] > 0 else "#8A8F98"
        ax.plot([x1, x2], [y1, y2], color=color, linewidth=0.6 + 2.2 * row["median_log2_oe"] / vmax, alpha=0.78, zorder=1)
    for state in states:
        x, y = pos[state]
        fill = "#B23A48" if state in {"SPP1/TAM", "IFN/APC"} else ("#2C7A51" if state == "CAF/matrix" else "#F4F6F8")
        text_color = "white" if state in {"CAF/matrix", "SPP1/TAM", "IFN/APC"} else "#222222"
        ax.add_patch(Circle((x, y), 0.072, facecolor=fill, edgecolor="#333333", linewidth=0.6, zorder=2))
        ax.text(x, y, state.replace("tumor ", "tumor\n"), ha="center", va="center", fontsize=6.0, color=text_color, fontweight="bold", zorder=3)
    ax.text(0.02, 0.04, "edges: positive observed/expected adjacency, support >= 4/6 sections", fontsize=6.2, color="#555555")
    return df


def transition_curves(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "core_to_interface_transition_context_summary.csv")
    df = df[df["context"].isin(CONTEXT_ORDER) & df["program"].isin(PROGRAMS)].copy()
    for program in PROGRAMS:
        sub = df[df["program"].eq(program)]
        curve = sub.groupby("axis_midpoint", as_index=False).agg(y=("median_mean_z", "median"))
        ax.plot(
            curve["axis_midpoint"],
            curve["y"],
            color=PROGRAM_COLORS[program],
            lw=1.8,
            marker="o",
            ms=3.2,
            label=program,
        )
    ax.axhline(0, color="#333333", lw=0.55)
    ax.set_xlim(0.04, 0.96)
    ax.set_xlabel("CAF-core to tumor-high axis", fontsize=7.2)
    ax.set_ylabel("median program z-score", fontsize=7.2)
    ax.set_title("Core-to-interface transition", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=6.2, loc="upper right", ncol=2)
    ax.text(0.02, -0.22, "pooled median across four main tissue contexts", transform=ax.transAxes, fontsize=6.2, color="#555555")
    return df


def geometry_context(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "caf_core_geometry_context_summary.csv")
    df = df[df["context"].isin(CONTEXT_ORDER)].copy()
    df["context"] = pd.Categorical(df["context"], categories=CONTEXT_ORDER, ordered=True)
    df = df.sort_values("context")
    x = np.arange(len(df))
    width = 0.34
    frag = df["median_fragmentation"].to_numpy(float)
    lcf = df["median_largest_component_fraction"].to_numpy(float)
    frag_scaled = (frag - np.nanmin(frag)) / (np.nanmax(frag) - np.nanmin(frag) + 1e-9)
    ax.bar(x - width / 2, frag_scaled, width=width, color="#8A8F98", label="fragmentation (scaled)", zorder=2)
    ax.bar(x + width / 2, lcf, width=width, color="#2C7A51", label="largest component fraction", zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(["post-\nNACT", "primary", "liver\nmet", "LN\nmet"], fontsize=7.0)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("scaled geometry metric", fontsize=7.2)
    ax.set_title("CAF-core geometry by context", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "F")
    clean_axes(ax, axis="y")
    ax.legend(frameon=False, fontsize=6.2, loc="upper left")
    return df.assign(median_fragmentation_scaled=frag_scaled)


def geometry_associations(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "caf_core_geometry_biological_correlations.csv")
    keep = [
        ("fragmentation_per_100_core_spots", "stromal_tumor_core_coupling"),
        ("largest_component_fraction", "stromal_tumor_core_coupling"),
        ("tumor_boundary_contacts_per_core_spot", "transition_rho__SPP1/TAM"),
        ("tumor_boundary_contacts_per_core_spot", "transition_rho__TGFb/EMT"),
        ("median_core_to_tumor_distance", "transition_rho__TGFb/EMT"),
        ("interface_core_fraction", "immune_decoupling_index"),
    ]
    parts = []
    for metric, readout in keep:
        parts.append(df[df["geometry_metric"].eq(metric) & df["biological_readout"].eq(readout)])
    out = pd.concat(parts, ignore_index=True)
    labels = [
        "fragmentation\nstromal-tumor",
        "dominant core\nstromal-tumor",
        "boundary contact\nSPP1/TAM",
        "boundary contact\nTGFb/EMT",
        "core distance\nTGFb/EMT",
        "interface fraction\nimmune decoupling",
    ]
    y = np.arange(len(out))
    colors = np.where(out["spearman_rho"] >= 0, "#B23A48", "#4C78A8")
    ax.axvline(0, color="#333333", lw=0.65)
    ax.hlines(y, 0, out["spearman_rho"], color=colors, lw=1.5, zorder=2)
    sig = out["p_value"].to_numpy(float) < 0.05
    ax.scatter(out.loc[sig, "spearman_rho"], y[sig], color=colors[sig], s=34, edgecolor="white", linewidth=0.45, zorder=3)
    ax.scatter(out.loc[~sig, "spearman_rho"], y[~sig], facecolor="white", edgecolor=colors[~sig], s=34, linewidth=1.0, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6.5)
    ax.invert_yaxis()
    ax.set_xlim(-0.72, 0.72)
    ax.set_xlabel("Spearman rho", fontsize=7.2)
    ax.set_title("Geometry-biological associations", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "G", x=-0.18)
    clean_axes(ax)
    ax.text(0.02, -0.19, "open dot: P >= 0.05", transform=ax.transAxes, fontsize=6.2, color="#555555")
    return out.assign(display_label=labels)


def claim_boundary(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "H", x=-0.04)
    ax.text(0.02, 0.98, "Claim boundary", fontsize=10.5, fontweight="bold", va="top")
    rows = [
        ("Strong", "specific spatial architecture"),
        ("Strong", "cell-level CAF-myeloid proximity"),
        ("Strong", "interface gradient organization"),
        ("Moderate", "geometry tunes coupling"),
        ("Open", "causal signaling direction"),
    ]
    colors = {"Strong": "#2C7A51", "Moderate": "#D88928", "Open": "#8A8F98"}
    for i, (level, text) in enumerate(rows):
        y = 0.82 - i * 0.15
        ax.add_patch(Rectangle((0.04, y - 0.045), 0.92, 0.085, facecolor="#F6F7F9", edgecolor="#DDDDDD", linewidth=0.5))
        ax.add_patch(Rectangle((0.04, y - 0.045), 0.08, 0.085, facecolor=colors[level], edgecolor=colors[level], linewidth=0))
        ax.text(0.08, y, level[0], ha="center", va="center", fontsize=7.0, color="white", fontweight="bold")
        ax.text(0.15, y, text, ha="left", va="center", fontsize=6.8, color="#333333")
    ax.text(0.04, 0.05, "Use as deepening evidence, not standalone causal proof.", fontsize=6.3, color="#555555")
    return pd.DataFrame(rows, columns=["claim_level", "claim"])


def write_source_data(*tables: pd.DataFrame) -> None:
    names = [
        "A_evidence_chain",
        "B_anchor_specificity",
        "C_anchor_neighborhoods",
        "D_xenium_adjacency_network",
        "E_core_to_interface_transition",
        "F_geometry_context",
        "G_geometry_associations",
        "H_claim_boundary",
    ]
    rows: list[dict[str, object]] = []
    for name, table in zip(names, tables):
        for _, row in table.iterrows():
            record = {"panel": name}
            for col, value in row.items():
                if isinstance(value, (np.floating, np.integer)):
                    record[col] = float(value)
                else:
                    record[col] = value
            rows.append(record)
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Supplementary Module 6 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "This synthesis module integrates Extended Data Figures 30-33 into a deep spatial-architecture evidence unit. "
        "It is designed to reduce figure-type monotony and make the mechanistic story more complete without overstating causality.\n\n"
        "## Core conclusion\n\n"
        "The CAF-myeloid niche is supported by anchor specificity, Xenium cell-level neighborhood enrichment, core-to-interface transition structure and CAF-core geometry associations. "
        "Together these results argue for a spatially organized tissue architecture rather than a generic expression-intensity artifact.\n\n"
        "## Panel contract\n\n"
        "- A: evidence chain and claim level.\n"
        "- B: specificity of the CAF-myeloid core against alternative biological anchors.\n"
        "- C: Xenium CAF-domain neighborhood enrichment/depletion.\n"
        "- D: Xenium cell-state adjacency network.\n"
        "- E: pooled core-to-interface program transition curves.\n"
        "- F: context-level CAF-core geometry summary.\n"
        "- G: geometry-biological association lollipop plot, including a negative boundary control.\n"
        "- H: claim boundary for manuscript language.\n\n"
        "## Boundary\n\n"
        "The module supports spatial association, architecture and perturbation-priority logic. It should not be written as direct causal proof without experimental perturbation.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf').relative_to(ROOT)}`\n"
        f"- `{OUT.with_suffix('.svg').relative_to(ROOT)}`\n"
        f"- `{OUT.with_suffix('.png').relative_to(ROOT)}`\n"
        f"- `{SOURCE_OUT.relative_to(ROOT)}`\n",
        encoding="utf-8",
    )


def main() -> None:
    fig = plt.figure(figsize=(15.8, 14.2), constrained_layout=False)
    gs = GridSpec(4, 8, figure=fig, height_ratios=[0.95, 1.06, 1.02, 0.56], hspace=0.92, wspace=0.92)
    fig.subplots_adjust(left=0.055, right=0.985, top=0.935, bottom=0.060)
    fig.suptitle("Supplementary Module 6 | Deep spatial architecture and mechanism-deepening evidence", fontsize=15.2, fontweight="bold", y=0.972)

    ax_a = fig.add_subplot(gs[0, 0:4])
    ax_b = fig.add_subplot(gs[0, 4:7])
    ax_h = fig.add_subplot(gs[0, 7:8])
    ax_c = fig.add_subplot(gs[1, 0:4])
    ax_d = fig.add_subplot(gs[1, 4:8])
    ax_e = fig.add_subplot(gs[2, 0:4])
    ax_f = fig.add_subplot(gs[2, 4:6])
    ax_g = fig.add_subplot(gs[2, 6:8])
    ax_note = fig.add_subplot(gs[3, :])

    a = evidence_chain(ax_a)
    b = anchor_specificity(ax_b)
    h = claim_boundary(ax_h)
    c = anchor_neighborhoods(ax_c)
    d = xenium_network(ax_d)
    e = transition_curves(ax_e)
    f = geometry_context(ax_f)
    g = geometry_associations(ax_g)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.70,
        "Interpretation: ED30-ED33 converge on a deeper model in which CAF-myeloid niches are spatially specific, reproduced at cell-neighborhood resolution, organized along a core-to-interface axis and modulated by CAF-core geometry.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.30,
        "Manuscript use: describe this as spatial architecture and perturbation-priority evidence. Reserve causal language for future blocking, organoid/co-culture or perturb-seq style validation.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

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
