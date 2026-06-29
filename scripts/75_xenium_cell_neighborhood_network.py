from __future__ import annotations

from collections import defaultdict
from pathlib import Path
import math

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
REPORT_DIR = PROJECT / "results" / "reports"
SOURCE_DIR = PROJECT / "results" / "source_data"
for d in [TABLE_DIR, FIG_DIR, REPORT_DIR, SOURCE_DIR]:
    d.mkdir(parents=True, exist_ok=True)

RNG = np.random.default_rng(20260628)
K_GRAPH = 8
K_NEIGHBORHOOD = 20
N_RANDOM_ANCHOR_SETS = 100
ANCHOR_FRACTION = 0.10

SCORE_COLS = {
    "CAF/matrix": "score_CAF_matrix",
    "SPP1/TAM": "score_SPP1_TAM",
    "IFN/APC": "score_IFN_APC",
    "T/NK": "score_T_NK",
    "TGFb/EMT": "score_TGFb_EMT",
    "tumor epithelial": "score_Tumor_epithelial",
}
STATE_ORDER = ["CAF/matrix", "SPP1/TAM", "IFN/APC", "T/NK", "TGFb/EMT", "tumor epithelial", "other"]
STATE_COLORS = {
    "CAF/matrix": "#2D6A8E",
    "SPP1/TAM": "#9C4E70",
    "IFN/APC": "#628F4E",
    "T/NK": "#4F7EC0",
    "TGFb/EMT": "#B57A3C",
    "tumor epithelial": "#6F5C99",
    "other": "#C8C8C8",
}
ANCHORS = {
    "CAF-APC domain": "anchor_CAF_APC",
    "CAF-SPP1/TAM domain": "anchor_CAF_SPP1TAM",
}


def assign_states(df: pd.DataFrame) -> pd.DataFrame:
    score_mat = df[list(SCORE_COLS.values())].to_numpy(float)
    labels = list(SCORE_COLS.keys())
    max_idx = np.nanargmax(score_mat, axis=1)
    max_val = score_mat[np.arange(len(df)), max_idx]
    out = df.copy()
    out["dominant_state"] = [labels[i] for i in max_idx]
    out.loc[~np.isfinite(max_val) | (max_val < 0), "dominant_state"] = "other"
    out["dominant_state_score"] = max_val
    return out


def graph_edges_for_sample(sample: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    xy = sample[["x_centroid", "y_centroid"]].to_numpy(float)
    labels = sample["dominant_state"].to_numpy(str)
    tree = cKDTree(xy)
    _, idx = tree.query(xy, k=K_GRAPH + 1)
    idx = idx[:, 1:]

    pair_counts: dict[tuple[str, str], int] = defaultdict(int)
    for i in range(len(sample)):
        a = labels[i]
        for j in idx[i]:
            if j <= i:
                continue
            b = labels[j]
            pair = tuple(sorted((a, b)))
            pair_counts[pair] += 1

    n_edges = sum(pair_counts.values())
    freqs = sample["dominant_state"].value_counts(normalize=True).reindex(STATE_ORDER, fill_value=0.0)
    rows = []
    for a in STATE_ORDER:
        for b in STATE_ORDER:
            if STATE_ORDER.index(a) > STATE_ORDER.index(b):
                continue
            observed = pair_counts.get(tuple(sorted((a, b))), 0)
            expected_fraction = freqs[a] * freqs[b] * (2 if a != b else 1)
            expected = n_edges * expected_fraction
            rows.append(
                {
                    "geo_accession": sample["geo_accession"].iloc[0],
                    "title": sample["title"].iloc[0],
                    "treatment": sample["treatment"].iloc[0],
                    "state_a": a,
                    "state_b": b,
                    "observed_edges": int(observed),
                    "expected_edges": float(expected),
                    "oe_ratio": float((observed + 1) / (expected + 1)),
                    "log2_oe": float(np.log2((observed + 1) / (expected + 1))),
                    "n_edges": int(n_edges),
                    "n_cells": int(len(sample)),
                }
            )
    edge_df = pd.DataFrame(rows)

    comp = (
        sample.groupby(["geo_accession", "title", "treatment", "dominant_state"], as_index=False)
        .agg(n_cells=("cell_id", "count"))
    )
    comp["fraction"] = comp["n_cells"] / len(sample)
    return edge_df, comp


def anchor_neighborhood_for_sample(sample: pd.DataFrame) -> pd.DataFrame:
    xy = sample[["x_centroid", "y_centroid"]].to_numpy(float)
    labels = sample["dominant_state"].to_numpy(str)
    tree = cKDTree(xy)
    rows = []
    for anchor_label, anchor_col in ANCHORS.items():
        vals = pd.to_numeric(sample[anchor_col], errors="coerce").to_numpy(float)
        n_anchor = max(50, int(math.ceil(ANCHOR_FRACTION * len(sample))))
        cutoff = np.nanquantile(vals, 1 - n_anchor / len(sample))
        anchor_idx = np.where(np.isfinite(vals) & (vals >= cutoff))[0]
        _, neigh_idx = tree.query(xy[anchor_idx], k=K_NEIGHBORHOOD + 1)
        observed_labels = labels[neigh_idx[:, 1:].ravel()]
        observed_counts = pd.Series(observed_labels).value_counts().reindex(STATE_ORDER, fill_value=0)

        random_fracs = []
        for _ in range(N_RANDOM_ANCHOR_SETS):
            rand_idx = RNG.choice(len(sample), size=len(anchor_idx), replace=False)
            _, rand_neigh_idx = tree.query(xy[rand_idx], k=K_NEIGHBORHOOD + 1)
            rand_labels = labels[rand_neigh_idx[:, 1:].ravel()]
            counts = pd.Series(rand_labels).value_counts().reindex(STATE_ORDER, fill_value=0)
            random_fracs.append((counts / counts.sum()).to_numpy(float))
        random_fracs = np.asarray(random_fracs)
        obs_frac = observed_counts / observed_counts.sum()
        for i, state in enumerate(STATE_ORDER):
            rows.append(
                {
                    "geo_accession": sample["geo_accession"].iloc[0],
                    "title": sample["title"].iloc[0],
                    "treatment": sample["treatment"].iloc[0],
                    "anchor": anchor_label,
                    "neighbor_state": state,
                    "n_anchor_cells": int(len(anchor_idx)),
                    "observed_fraction": float(obs_frac[state]),
                    "random_median_fraction": float(np.nanmedian(random_fracs[:, i])),
                    "delta_fraction": float(obs_frac[state] - np.nanmedian(random_fracs[:, i])),
                    "random_p05_fraction": float(np.nanpercentile(random_fracs[:, i], 5)),
                    "random_p95_fraction": float(np.nanpercentile(random_fracs[:, i], 95)),
                }
            )
    return pd.DataFrame(rows)


def select_roi(cells: pd.DataFrame) -> pd.DataFrame:
    sample = cells[cells["geo_accession"].eq("GSM8454446")].copy()
    if sample.empty:
        sample = cells.groupby("geo_accession").get_group(cells["geo_accession"].iloc[0]).copy()
    vals = sample["anchor_CAF_SPP1TAM"].to_numpy(float)
    xy_all = sample[["x_centroid", "y_centroid"]].to_numpy(float)
    xlo, xhi = np.nanquantile(xy_all[:, 0], [0.08, 0.92])
    ylo, yhi = np.nanquantile(xy_all[:, 1], [0.08, 0.92])
    central = np.where((xy_all[:, 0] >= xlo) & (xy_all[:, 0] <= xhi) & (xy_all[:, 1] >= ylo) & (xy_all[:, 1] <= yhi))[0]
    top_pool = central[np.argsort(vals[central])[-min(1200, len(central)) :]] if len(central) else np.argsort(vals)[-min(1200, len(sample)) :]
    xy = sample[["x_centroid", "y_centroid"]].to_numpy(float)
    state = sample["dominant_state"].to_numpy(str)
    best_idx = None
    best_score = -np.inf
    for idx in top_pool[:: max(1, len(top_pool) // 120)]:
        center = xy[idx]
        mask = (np.abs(xy[:, 0] - center[0]) <= 360.0) & (np.abs(xy[:, 1] - center[1]) <= 360.0)
        if mask.sum() < 500:
            continue
        states = pd.Series(state[mask]).value_counts(normalize=True)
        score = states.get("CAF/matrix", 0) + states.get("SPP1/TAM", 0) + states.get("IFN/APC", 0) + states.get("T/NK", 0) - states.get("tumor epithelial", 0)
        if score > best_score:
            best_score = score
            best_idx = int(idx)
    center = xy[best_idx if best_idx is not None else int(RNG.choice(top_pool))]
    radius = 320.0
    for r in [240.0, 320.0, 420.0, 560.0]:
        mask = (np.abs(xy[:, 0] - center[0]) <= r) & (np.abs(xy[:, 1] - center[1]) <= r)
        if 600 <= mask.sum() <= 5000:
            radius = r
            break
    roi = sample.loc[(np.abs(xy[:, 0] - center[0]) <= radius) & (np.abs(xy[:, 1] - center[1]) <= radius)].copy()
    return roi


def draw_network(ax, edge_summary: pd.DataFrame) -> None:
    pos = {
        "CAF/matrix": (0.10, 0.52),
        "SPP1/TAM": (0.38, 0.78),
        "IFN/APC": (0.78, 0.66),
        "T/NK": (0.88, 0.22),
        "TGFb/EMT": (0.36, 0.04),
        "tumor epithelial": (0.02, 0.12),
    }
    edges = edge_summary[
        (edge_summary["state_a"] != "other")
        & (edge_summary["state_b"] != "other")
        & (edge_summary["state_a"] != edge_summary["state_b"])
        & (edge_summary["median_log2_oe"] > 0.08)
    ].copy()
    edges = edges.sort_values("median_log2_oe", ascending=False).head(12)
    for _, row in edges.iterrows():
        a, b = row["state_a"], row["state_b"]
        xa, ya = pos[a]
        xb, yb = pos[b]
        lw = 0.6 + 2.8 * min(row["median_log2_oe"], 1.5) / 1.5
        ax.plot([xa, xb], [ya, yb], color="#4A4A4A", alpha=0.58, lw=lw, zorder=1)
    for state, (x, y) in pos.items():
        ax.scatter([x], [y], s=720, color=STATE_COLORS[state], edgecolor="white", linewidth=1.5, zorder=3)
        ax.text(x, y, state.replace(" ", "\n"), ha="center", va="center", fontsize=7.2, color="white", zorder=4)
    ax.set_xlim(-0.12, 1.02)
    ax.set_ylim(-0.12, 0.92)
    ax.axis("off")
    ax.set_title("B  Enriched cell-state adjacency network", loc="left", fontweight="bold")


def make_figure(cells: pd.DataFrame, edge_summary: pd.DataFrame, neigh_summary: pd.DataFrame, comp: pd.DataFrame) -> None:
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
            "pdf.fonttype": 42,
            "svg.fonttype": "none",
            "axes.linewidth": 0.8,
        }
    )
    fig = plt.figure(figsize=(14.2, 9.4))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.08, 1.0, 1.04], height_ratios=[1.0, 0.92], wspace=0.42, hspace=0.48)
    fig.suptitle("Xenium cell-neighborhood networks resolve CAF-domain spatial organization", fontsize=16, fontweight="bold")

    ax_a = fig.add_subplot(gs[0, 0])
    roi = select_roi(cells)
    for state in STATE_ORDER[::-1]:
        sub = roi[roi["dominant_state"].eq(state)]
        if sub.empty:
            continue
        size = 5 if state == "other" else 10
        alpha = 0.35 if state == "other" else 0.86
        ax_a.scatter(sub["x_centroid"], sub["y_centroid"], s=size, color=STATE_COLORS[state], label=state, alpha=alpha, linewidth=0)
    ax_a.set_aspect("equal")
    ax_a.invert_yaxis()
    ax_a.set_xticks([])
    ax_a.set_yticks([])
    ax_a.set_title("A  Representative cell-level neighborhood", loc="left", fontweight="bold")
    handles = [Line2D([0], [0], marker="o", color="w", label=s, markerfacecolor=STATE_COLORS[s], markersize=6) for s in STATE_ORDER if s != "other"]
    ax_a.legend(handles=handles, frameon=False, fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.03))

    ax_b = fig.add_subplot(gs[0, 1])
    draw_network(ax_b, edge_summary)

    ax_c = fig.add_subplot(gs[0, 2])
    mat = (
        edge_summary.pivot(index="state_a", columns="state_b", values="median_log2_oe")
        .reindex(index=STATE_ORDER[:-1], columns=STATE_ORDER[:-1])
    )
    im = ax_c.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.8, vmax=0.8)
    ax_c.set_xticks(range(len(mat.columns)))
    ax_c.set_xticklabels(mat.columns, rotation=45, ha="right", fontsize=8)
    ax_c.set_yticks(range(len(mat.index)))
    ax_c.set_yticklabels(mat.index, fontsize=8)
    ax_c.set_title("C  Observed/expected adjacency", loc="left", fontweight="bold")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iat[i, j]
            if pd.notna(val):
                ax_c.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7, color="white" if abs(val) > 0.55 else "#111")
    cb = fig.colorbar(im, ax=ax_c, fraction=0.046, pad=0.03)
    cb.set_label("median log2 observed/expected")

    ax_d = fig.add_subplot(gs[1, 0])
    comp_summary = (
        comp.groupby(["treatment", "dominant_state"], as_index=False)
        .agg(fraction=("fraction", "median"))
    )
    treatments = list(comp_summary["treatment"].drop_duplicates())
    bottom = np.zeros(len(treatments))
    x = np.arange(len(treatments))
    for state in STATE_ORDER[:-1]:
        vals = comp_summary[comp_summary["dominant_state"].eq(state)].set_index("treatment").reindex(treatments)["fraction"].fillna(0).to_numpy()
        ax_d.bar(x, vals, bottom=bottom, color=STATE_COLORS[state], label=state)
        bottom += vals
    ax_d.set_xticks(x)
    ax_d.set_xticklabels(treatments, rotation=20, ha="right", fontsize=8)
    ax_d.set_ylabel("median cell-state fraction")
    ax_d.set_title("D  Cell-state composition by treatment", loc="left", fontweight="bold")

    ax_e = fig.add_subplot(gs[1, 1])
    sub = neigh_summary[neigh_summary["anchor"].eq("CAF-SPP1/TAM domain") & neigh_summary["neighbor_state"].isin(STATE_ORDER[:-1])]
    sub = sub.set_index("neighbor_state").reindex(STATE_ORDER[:-1]).reset_index()
    ax_e.barh(np.arange(len(sub)), sub["median_delta_fraction"], color=[STATE_COLORS[s] for s in sub["neighbor_state"]])
    ax_e.axvline(0, color="#333", lw=0.8)
    ax_e.set_yticks(np.arange(len(sub)))
    ax_e.set_yticklabels(sub["neighbor_state"], fontsize=8)
    ax_e.invert_yaxis()
    ax_e.set_xlabel("delta versus random neighbors")
    ax_e.set_title("E  CAF-SPP1/TAM neighborhood composition", loc="left", fontweight="bold")

    ax_f = fig.add_subplot(gs[1, 2])
    sub = neigh_summary[neigh_summary["anchor"].eq("CAF-APC domain") & neigh_summary["neighbor_state"].isin(STATE_ORDER[:-1])]
    sub = sub.set_index("neighbor_state").reindex(STATE_ORDER[:-1]).reset_index()
    ax_f.barh(np.arange(len(sub)), sub["median_delta_fraction"], color=[STATE_COLORS[s] for s in sub["neighbor_state"]])
    ax_f.axvline(0, color="#333", lw=0.8)
    ax_f.set_yticks(np.arange(len(sub)))
    ax_f.set_yticklabels(sub["neighbor_state"], fontsize=8)
    ax_f.invert_yaxis()
    ax_f.set_xlabel("delta versus random neighbors")
    ax_f.set_title("F  CAF-APC neighborhood composition", loc="left", fontweight="bold")

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    out = FIG_DIR / "extended_data_figure31_xenium_cell_neighborhood_network"
    for ext in ["pdf", "svg", "png"]:
        fig.savefig(out.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_report(edge_summary: pd.DataFrame, neigh_summary: pd.DataFrame) -> None:
    top_edges = edge_summary[
        (edge_summary["state_a"] != edge_summary["state_b"])
        & ~edge_summary["state_a"].eq("other")
        & ~edge_summary["state_b"].eq("other")
    ].sort_values("median_log2_oe", ascending=False).head(8)
    lines = [
        "# GSE274673 Xenium Cell-Neighborhood Network",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Question",
        "",
        "This analysis tests whether cell-resolution Xenium data resolve a spatial neighborhood grammar around CAF-domain anchors, rather than only reproducing module-score distance gradients.",
        "",
        "## Method",
        "",
        f"Cells were assigned a dominant program state from CAF/matrix, SPP1/TAM, IFN/APC, T/NK, TGFb/EMT and tumor epithelial scores. A k={K_GRAPH} nearest-neighbor graph was built within each section. Pairwise state adjacency was summarized as observed/expected edge enrichment. CAF-domain neighborhood composition used k={K_NEIGHBORHOOD} nearest neighbors around top {int(ANCHOR_FRACTION*100)}% CAF-APC or CAF-SPP1/TAM anchor cells and {N_RANDOM_ANCHOR_SETS} matched random anchor sets per sample.",
        "",
        "## Strongest Enriched Cross-State Edges",
        "",
    ]
    for _, row in top_edges.iterrows():
        lines.append(f"- {row.state_a} - {row.state_b}: median log2 observed/expected {row.median_log2_oe:.2f}.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "The analysis adds a cell-neighborhood layer to the CAF-domain model. It should be interpreted as spatial adjacency enrichment in observational Xenium data, not as direct cell-cell signaling or perturbational causality.",
            "",
            "## Outputs",
            "",
            "- `results/tables/gse274673_xenium_cell_state_adjacency_per_sample.csv`",
            "- `results/tables/gse274673_xenium_cell_state_adjacency_summary.csv`",
            "- `results/tables/gse274673_xenium_anchor_neighborhood_composition.csv`",
            "- `results/tables/gse274673_xenium_anchor_neighborhood_summary.csv`",
            "- `results/source_data/Source_Data_Extended_Data_Fig_31_xenium_neighborhood_network.csv`",
            "- `results/figures/submission/extended_data_figure31_xenium_cell_neighborhood_network.pdf`",
        ]
    )
    (REPORT_DIR / "gse274673_xenium_cell_neighborhood_network_report.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    cells = pd.read_csv(TABLE_DIR / "gse274673_xenium_fixed_anchor_cell_scores.csv")
    cells = assign_states(cells)

    edge_frames = []
    comp_frames = []
    neigh_frames = []
    for accession, sample in cells.groupby("geo_accession", sort=False):
        edge_df, comp = graph_edges_for_sample(sample)
        neigh = anchor_neighborhood_for_sample(sample)
        edge_frames.append(edge_df)
        comp_frames.append(comp)
        neigh_frames.append(neigh)
        print(f"{accession}: cells={len(sample)} edges={edge_df['n_edges'].iloc[0]}")

    edges = pd.concat(edge_frames, ignore_index=True)
    comp = pd.concat(comp_frames, ignore_index=True)
    neigh = pd.concat(neigh_frames, ignore_index=True)

    edge_summary = (
        edges.groupby(["state_a", "state_b"], as_index=False)
        .agg(
            n_samples=("geo_accession", "nunique"),
            median_log2_oe=("log2_oe", "median"),
            median_oe_ratio=("oe_ratio", "median"),
            support_n=("log2_oe", lambda x: int(np.sum(np.asarray(x) > 0))),
        )
    )
    edge_summary["support_fraction"] = edge_summary["support_n"] / edge_summary["n_samples"]

    neigh_summary = (
        neigh.groupby(["anchor", "neighbor_state"], as_index=False)
        .agg(
            n_samples=("geo_accession", "nunique"),
            median_observed_fraction=("observed_fraction", "median"),
            median_random_fraction=("random_median_fraction", "median"),
            median_delta_fraction=("delta_fraction", "median"),
            support_n=("delta_fraction", lambda x: int(np.sum(np.asarray(x) > 0))),
        )
    )
    neigh_summary["support_fraction"] = neigh_summary["support_n"] / neigh_summary["n_samples"]

    edges.to_csv(TABLE_DIR / "gse274673_xenium_cell_state_adjacency_per_sample.csv", index=False)
    edge_summary.to_csv(TABLE_DIR / "gse274673_xenium_cell_state_adjacency_summary.csv", index=False)
    comp.to_csv(TABLE_DIR / "gse274673_xenium_cell_state_composition.csv", index=False)
    neigh.to_csv(TABLE_DIR / "gse274673_xenium_anchor_neighborhood_composition.csv", index=False)
    neigh_summary.to_csv(TABLE_DIR / "gse274673_xenium_anchor_neighborhood_summary.csv", index=False)

    edge_summary.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_31_xenium_neighborhood_network.csv", index=False)
    neigh_summary.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_31_anchor_neighborhood_summary.csv", index=False)
    comp.to_csv(SOURCE_DIR / "Source_Data_Extended_Data_Fig_31_cell_state_composition.csv", index=False)

    make_figure(cells, edge_summary, neigh_summary, comp)
    write_report(edge_summary, neigh_summary)
    print("Wrote Xenium cell-neighborhood network outputs.")


if __name__ == "__main__":
    main()
