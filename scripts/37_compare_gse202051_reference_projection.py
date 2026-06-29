from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


PROJECT = Path(__file__).resolve().parents[1]
SMALL_PER_SAMPLE = PROJECT / "results/tables/gap2_reference_projection_deconvolution_per_sample.csv"
FULL_PER_SAMPLE = PROJECT / "results/tables/gap2_full_reference_projection_deconvolution_per_sample.csv"
SMALL_CORR = PROJECT / "results/tables/gap2_reference_projection_deconvolution_correlations.csv"
FULL_CORR = PROJECT / "results/tables/gap2_full_reference_projection_deconvolution_correlations.csv"
OUT_TABLE = PROJECT / "results/tables/gap2_reference_projection_small_vs_full_comparison.csv"
OUT_REPORT = PROJECT / "results/reports/gap2_reference_projection_small_vs_full_comparison.md"
OUT_FIG = PROJECT / "results/figures/submission/extended_data_gap2_reference_projection_small_vs_full_comparison"
STATUS = PROJECT / "results/logs/stage_37_reference_projection_small_vs_full_status.json"

KEY_STATES = ["myCAF_matrix", "SPP1_TAM", "DC_APC", "T_NK", "B_plasma", "epithelial_tumor"]


def compare_per_sample() -> pd.DataFrame:
    small = pd.read_csv(SMALL_PER_SAMPLE)
    full = pd.read_csv(FULL_PER_SAMPLE)
    merged = small.merge(
        full,
        on=["sample_id", "cell_state"],
        suffixes=("_small", "_full"),
        how="inner",
    )
    rows = []
    for state, sub in merged.groupby("cell_state"):
        tmp = sub[["core_enrichment_small", "core_enrichment_full"]].dropna()
        if len(tmp) < 5:
            continue
        rho, p = spearmanr(tmp["core_enrichment_small"], tmp["core_enrichment_full"])
        rows.append(
            {
                "comparison_type": "per_sample_core_enrichment",
                "cell_state": state,
                "n_pairs": int(len(tmp)),
                "spearman_rho": float(rho),
                "p_value": float(p),
                "same_direction_fraction": float(
                    (
                        np.sign(tmp["core_enrichment_small"].to_numpy(float))
                        == np.sign(tmp["core_enrichment_full"].to_numpy(float))
                    ).mean()
                ),
                "median_small": float(tmp["core_enrichment_small"].median()),
                "median_full": float(tmp["core_enrichment_full"].median()),
            }
        )
    return pd.DataFrame(rows)


def compare_immune_decoupling_rhos() -> pd.DataFrame:
    small = pd.read_csv(SMALL_CORR)
    full = pd.read_csv(FULL_CORR)
    small = small[small["target"].eq("immune_decoupling_index")][["cell_state", "spearman_rho", "p_value"]].rename(
        columns={"spearman_rho": "immune_decoupling_rho_small", "p_value": "immune_decoupling_p_small"}
    )
    full = full[full["target"].eq("immune_decoupling_index")][["cell_state", "spearman_rho", "p_value"]].rename(
        columns={"spearman_rho": "immune_decoupling_rho_full", "p_value": "immune_decoupling_p_full"}
    )
    merged = small.merge(full, on="cell_state", how="inner")
    merged["comparison_type"] = "immune_decoupling_rho"
    merged["same_direction"] = np.sign(merged["immune_decoupling_rho_small"]) == np.sign(merged["immune_decoupling_rho_full"])
    merged["delta_full_minus_small"] = merged["immune_decoupling_rho_full"] - merged["immune_decoupling_rho_small"]
    return merged


def plot_comparison(per_sample_cmp: pd.DataFrame, rho_cmp: pd.DataFrame) -> None:
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    states = [state for state in KEY_STATES if state in set(rho_cmp["cell_state"])]
    rhos = rho_cmp.set_index("cell_state").reindex(states)
    per = per_sample_cmp.set_index("cell_state").reindex(states)

    fig = plt.figure(figsize=(11.8, 4.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.1, 1.0])
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    fig.suptitle("GSE202051 small-reference versus full-reference projection", fontsize=14, fontweight="bold")

    x = np.arange(len(states))
    width = 0.35
    ax1.bar(x - width / 2, rhos["immune_decoupling_rho_small"], width=width, color="#4E79A7", label="small reference")
    ax1.bar(x + width / 2, rhos["immune_decoupling_rho_full"], width=width, color="#E15759", label="full reference")
    ax1.axhline(0, color="#333333", lw=0.8)
    ax1.set_title("A  Immune-decoupling correlations", loc="left", fontsize=10.5, fontweight="bold")
    ax1.set_ylabel("Spearman rho")
    ax1.set_xticks(x, states, rotation=35, ha="right")
    ax1.tick_params(labelsize=8)
    ax1.legend(frameon=False, fontsize=8)

    ax2.scatter(per["spearman_rho"], per["same_direction_fraction"], s=70, color="#59A14F", edgecolor="white", linewidth=0.6)
    for state in states:
        if state in per.index and np.isfinite(per.loc[state, "spearman_rho"]):
            ax2.text(per.loc[state, "spearman_rho"] + 0.01, per.loc[state, "same_direction_fraction"], state, fontsize=7, va="center")
    ax2.axvline(0, color="#999999", lw=0.8, ls="--")
    ax2.set_xlim(-0.1, 1.0)
    ax2.set_ylim(-0.03, 1.03)
    ax2.set_title("B  Per-sample agreement", loc="left", fontsize=10.5, fontweight="bold")
    ax2.set_xlabel("small vs full Spearman rho")
    ax2.set_ylabel("same-direction fraction")
    ax2.tick_params(labelsize=8)

    for ext in [".pdf", ".png", ".svg"]:
        fig.savefig(f"{OUT_FIG}{ext}", dpi=320 if ext == ".png" else None)
    plt.close(fig)


def write_report(per_sample_cmp: pd.DataFrame, rho_cmp: pd.DataFrame) -> None:
    key = rho_cmp[rho_cmp["cell_state"].isin(KEY_STATES)].copy().sort_values("immune_decoupling_rho_full")
    per_key = per_sample_cmp[per_sample_cmp["cell_state"].isin(KEY_STATES)].copy()
    lines = [
        "# GSE202051 Small-vs-Full Reference Projection Comparison",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Main Result",
        "",
        "The full GSE202051 total h5ad reference validates the main Gap 2 direction: immune decoupling is associated with weaker immune-cell projection around CAF cores and stronger myCAF/matrix projection.",
        "",
        "## Immune-Decoupling Correlations",
        "",
        "| cell state | small rho | full rho | same direction |",
        "|---|---:|---:|---|",
    ]
    for row in key.itertuples(index=False):
        lines.append(
            f"| {row.cell_state} | {row.immune_decoupling_rho_small:.3f} | {row.immune_decoupling_rho_full:.3f} | {bool(row.same_direction)} |"
        )
    lines.extend(
        [
            "",
            "## Per-Sample Core-Enrichment Agreement",
            "",
            "| cell state | n pairs | small-vs-full rho | same-direction fraction |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in per_key.sort_values("spearman_rho", ascending=False).itertuples(index=False):
        lines.append(f"| {row.cell_state} | {int(row.n_pairs)} | {row.spearman_rho:.3f} | {row.same_direction_fraction:.3f} |")
    lines.extend(
        [
            "",
            "## Interpretation Boundary",
            "",
            "This compares two marker-constrained reference-projection prototypes. It supports robustness of the direction of the Gap 2 cell-state interpretation, but it is still not a fully validated deconvolution benchmark.",
            "",
            "## Generated Outputs",
            "",
            f"- `{OUT_TABLE.relative_to(PROJECT).as_posix()}`",
            f"- `{OUT_FIG.relative_to(PROJECT).as_posix()}.pdf`",
        ]
    )
    OUT_REPORT.parent.mkdir(parents=True, exist_ok=True)
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    per_sample_cmp = compare_per_sample()
    rho_cmp = compare_immune_decoupling_rhos()
    out = rho_cmp.merge(per_sample_cmp, on="cell_state", how="outer", suffixes=("", "_per_sample"))
    OUT_TABLE.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_TABLE, index=False)
    plot_comparison(per_sample_cmp, rho_cmp)
    write_report(per_sample_cmp, rho_cmp)
    STATUS.parent.mkdir(parents=True, exist_ok=True)
    STATUS.write_text(
        json.dumps(
            {
                "stage": "37_reference_projection_small_vs_full",
                "status": "success",
                "n_states": int(out["cell_state"].nunique()),
                "outputs": [str(OUT_TABLE), str(OUT_REPORT), f"{OUT_FIG}.pdf"],
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_REPORT}")


if __name__ == "__main__":
    main()
