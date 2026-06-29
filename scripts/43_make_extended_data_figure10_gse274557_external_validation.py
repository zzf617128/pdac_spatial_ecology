from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
TABLE_DIR = PROJECT / "results" / "tables"
FIG_DIR = PROJECT / "results" / "figures" / "submission"
FIG_DIR.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "extended_data_figure10_gse274557_external_validation"

TISSUE_ORDER = ["Primary PDAC", "Liver metastasis", "Lung metastasis", "Peritoneal metastasis"]
PROGRAM_ORDER = ["IFN/MHC", "immune_core", "tumor_aggressive", "SPP1_TAM", "TGFb_EMT"]
PROGRAM_LABELS = {
    "IFN/MHC": "IFN/MHC",
    "immune_core": "immune core",
    "tumor_aggressive": "tumor aggressive",
    "SPP1_TAM": "SPP1/TAM",
    "TGFb_EMT": "TGF-beta/EMT",
}
COLORS = {
    "Primary PDAC": "#4C78A8",
    "Liver metastasis": "#59A14F",
    "Lung metastasis": "#E15759",
    "Peritoneal metastasis": "#B279A2",
}


def main() -> None:
    context = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_context_summary.csv")
    gradients = pd.read_csv(TABLE_DIR / "gse274557_full_caf_core_gradients.csv")

    counts = (
        gradients[["sample_id", "tissue"]]
        .drop_duplicates()
        .groupby("tissue")
        .size()
        .reindex(TISSUE_ORDER)
        .fillna(0)
        .astype(int)
    )

    plot_df = context[context["target_program"].isin(PROGRAM_ORDER)].copy()
    matrix = plot_df.pivot(index="target_program", columns="tissue", values="median_delta_vs_random").reindex(
        index=PROGRAM_ORDER, columns=TISSUE_ORDER
    )
    labels = (
        plot_df.assign(
            label=lambda x: x["median_delta_vs_random"].map(lambda v: f"{v:.2f}")
            + "\n"
            + x["n_more_negative"].astype(int).astype(str)
            + "/"
            + x["n_samples"].astype(int).astype(str)
        )
        .pivot(index="target_program", columns="tissue", values="label")
        .reindex(index=PROGRAM_ORDER, columns=TISSUE_ORDER)
    )

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

    fig = plt.figure(figsize=(11.4, 7.2), constrained_layout=False)
    gs = fig.add_gridspec(
        2,
        2,
        width_ratios=[0.95, 1.55],
        height_ratios=[1.0, 0.72],
        left=0.07,
        right=0.96,
        top=0.88,
        bottom=0.12,
        wspace=0.42,
        hspace=0.62,
    )
    fig.suptitle(
        "Extended Data Fig. 10 | Independent GSE274557 metastatic PDAC validation",
        fontsize=15,
        fontweight="bold",
        y=0.96,
    )

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[:, 1])
    ax_c = fig.add_subplot(gs[1, 0])

    x = np.arange(len(counts))
    ax_a.bar(x, counts.values, color=[COLORS[t] for t in counts.index], edgecolor="#333333", linewidth=0.6)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(["Primary", "Liver", "Lung", "Peritoneal"], rotation=25, ha="right")
    ax_a.set_ylabel("samples")
    ax_a.set_title("A  External cohort composition", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_a.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_a.set_axisbelow(True)
    for xi, val in zip(x, counts.values):
        ax_a.text(xi, val + 0.4, str(val), ha="center", va="bottom", fontsize=8)

    im = ax_b.imshow(matrix.to_numpy(float), cmap="RdBu_r", vmin=-0.50, vmax=0.20, aspect="auto")
    ax_b.set_title("B  CAF-core gradients versus 1,000 random cores", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_b.set_xticks(np.arange(len(TISSUE_ORDER)))
    ax_b.set_xticklabels(TISSUE_ORDER, rotation=35, ha="right")
    ax_b.set_yticks(np.arange(len(PROGRAM_ORDER)))
    ax_b.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER])
    ax_b.tick_params(length=0)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            value = matrix.iat[i, j]
            text_color = "white" if pd.notna(value) and value <= -0.32 else "#111111"
            ax_b.text(j, i, labels.iat[i, j], ha="center", va="center", fontsize=7.4, color=text_color)
    cb = fig.colorbar(im, ax=ax_b, fraction=0.046, pad=0.025)
    cb.set_label("median observed-minus-random rho")

    tumor = gradients[gradients["target_program"].eq("tumor_aggressive")].copy()
    rng = np.random.default_rng(20260628)
    for i, tissue in enumerate(TISSUE_ORDER):
        vals = tumor.loc[tumor["tissue"].eq(tissue), "delta_vs_random_median"].dropna().to_numpy(float)
        jitter = rng.uniform(-0.12, 0.12, size=len(vals))
        ax_c.scatter(np.full(len(vals), i) + jitter, vals, s=20, color=COLORS[tissue], edgecolor="#333333", linewidth=0.4, alpha=0.82)
        if len(vals):
            ax_c.plot([i - 0.24, i + 0.24], [np.median(vals), np.median(vals)], color="#111111", lw=1.4)
    ax_c.axhline(0, color="#333333", lw=0.8, ls="--")
    ax_c.set_xticks(np.arange(len(TISSUE_ORDER)))
    ax_c.set_xticklabels(["Primary", "Liver", "Lung", "Peritoneal"], rotation=35, ha="right")
    ax_c.set_ylabel("delta vs random")
    ax_c.set_title("C  Tumor-aggressive variation", loc="left", fontsize=10, fontweight="bold", pad=8)
    ax_c.grid(axis="y", color="#dddddd", linewidth=0.6)
    ax_c.set_axisbelow(True)

    for ext in ["pdf", "png", "svg"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300)
    plt.close(fig)
    print(f"Wrote {OUT}.pdf/.png/.svg")


if __name__ == "__main__":
    main()
