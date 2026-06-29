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
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
SOURCE_DIR = ROOT / "results" / "source_data"
FIG_DIR = ROOT / "results" / "figures" / "submission"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "extended_data_figure27_metastatic_immune_decoupling_module_nc_style"
SOURCE_OUT = SOURCE_DIR / "Source_Data_Extended_Data_Fig_27_metastatic_immune_decoupling_module.csv"
REPORT_OUT = REPORT_DIR / "extended_data_figure27_metastatic_immune_decoupling_module_notes.md"

SITE_ORDER = ["primary_tumor", "liver_metastasis", "lymph_node_metastasis"]
SITE_LABELS = {
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
}
PROGRAM_ORDER = ["mycaf", "spp1_tam", "tgfb", "tumor_aggressive", "ifn_mhc", "immune_core", "dc_apc", "t_cell", "b_cell"]
PROGRAM_LABELS = {
    "mycaf": "myCAF",
    "spp1_tam": "SPP1/TAM",
    "tgfb": "TGF-beta",
    "tumor_aggressive": "Tumor aggr.",
    "ifn_mhc": "IFN/MHC",
    "immune_core": "Immune core",
    "dc_apc": "DC/APC",
    "t_cell": "T cell",
    "b_cell": "B cell",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def plot_site_scale(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse272362_rds_sample_specimen_summary.csv")
    df = df[df["specimen_type"].isin(SITE_ORDER)].copy()
    summary = df.groupby("specimen_type", as_index=False).agg(n_samples=("sample_id", "nunique"), n_spots=("n_spots", "sum"))
    summary["ord"] = summary["specimen_type"].map({v: i for i, v in enumerate(SITE_ORDER)})
    summary = summary.sort_values("ord")
    x = np.arange(len(summary))
    ax.bar(x, summary["n_spots"] / 1000, color=["#4C78A8", "#2C7A51", "#B23A48"], width=0.58)
    ax.set_xticks(x)
    ax.set_xticklabels([SITE_LABELS[s] for s in summary["specimen_type"]], fontsize=7)
    ax.set_ylabel("spots (thousand)", fontsize=7)
    ax.set_title("GSE272362 metastatic-site scale", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "A")
    clean_axes(ax)
    for i, row in summary.reset_index(drop=True).iterrows():
        ax.text(i, row["n_spots"] / 1000 + 1.2, f"n={int(row['n_samples'])}", ha="center", va="bottom", fontsize=6.5)
    return summary


def plot_random_support(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(SOURCE_DIR / "Source_Data_Fig_2.csv")
    df = df[df["panel"].eq("B") & df["group"].isin(SITE_ORDER)].copy()
    keep_metrics = ["ifn_mhc delta_vs_null", "immune_core delta_vs_null", "tumor_aggressive delta_vs_null"]
    df = df[df["metric"].isin(keep_metrics)].copy()
    df["program"] = df["metric"].str.replace(" delta_vs_null", "", regex=False)
    label_map = {"ifn_mhc": "IFN/MHC", "immune_core": "Immune core", "tumor_aggressive": "Tumor aggr."}
    df["program_label"] = df["program"].map(label_map)
    mat = df.pivot_table(index="program_label", columns="group", values="value", aggfunc="median").reindex(
        index=["IFN/MHC", "Immune core", "Tumor aggr."], columns=SITE_ORDER
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.45, vmax=0.45, aspect="auto")
    ax.set_xticks(np.arange(len(SITE_ORDER)))
    ax.set_xticklabels([SITE_LABELS[x] for x in SITE_ORDER], rotation=25, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=7)
    ax.set_title("Random-core support by site", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "B")
    for i, program in enumerate(mat.index):
        for j, site in enumerate(SITE_ORDER):
            val = mat.loc[program, site]
            row = df[(df["program_label"].eq(program)) & (df["group"].eq(site))]
            support = row["support"].iloc[0] if len(row) else ""
            ax.text(j, i, f"{val:.2f}\n{support}", ha="center", va="center", fontsize=6, color="#FFFFFF" if abs(val) > 0.27 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("delta vs random", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_subprogram_gradients(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse272362_caf_core_subprogram_gradient_summary.csv")
    df = df[df["specimen_type"].isin(SITE_ORDER) & df["program"].isin(PROGRAM_ORDER)].copy()
    mat = df.pivot_table(index="program", columns="specimen_type", values="median_rho", aggfunc="median").reindex(
        index=PROGRAM_ORDER, columns=SITE_ORDER
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.55, vmax=0.30, aspect="auto")
    ax.set_xticks(np.arange(len(SITE_ORDER)))
    ax.set_xticklabels([SITE_LABELS[x] for x in SITE_ORDER], fontsize=7)
    ax.set_yticks(np.arange(len(PROGRAM_ORDER)))
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], fontsize=7)
    ax.set_title("CAF-core subprogram gradients", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "C")
    for i, program in enumerate(PROGRAM_ORDER):
        for j, site in enumerate(SITE_ORDER):
            val = mat.loc[program, site]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > 0.32 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("rho to CAF-core distance", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_decoupling(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "immune_decoupling_context_summary.csv")
    df = df[df["cohort_context"].isin(SITE_ORDER)].copy()
    df["ord"] = df["cohort_context"].map({v: i for i, v in enumerate(SITE_ORDER)})
    df = df.sort_values("ord")
    x = np.arange(len(df))
    ax.plot(x, df["median_stromal_tumor_core_coupling"], marker="o", color="#B23A48", lw=1.8, label="stromal-tumor")
    ax.plot(x, df["median_immune_core_coupling"], marker="o", color="#2C7A51", lw=1.8, label="immune-core")
    ax.plot(x, df["median_immune_decoupling_index"], marker="o", color="#4C78A8", lw=1.8, label="decoupling index")
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_xticks(x)
    ax.set_xticklabels([SITE_LABELS[s] for s in df["cohort_context"]], fontsize=7)
    ax.set_ylabel("median coupling/index", fontsize=7)
    ax.set_title("Selective immune decoupling", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "D")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=6, loc="upper left")
    return df


def plot_matched_delta_summary(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse272362_patient_matched_site_delta_summary.csv")
    df = df[df["metric"].eq("rho_distance_to_caf_core") & df["comparison_specimen_type"].isin(["liver_metastasis", "lymph_node_metastasis"])].copy()
    df = df[df["program"].isin(PROGRAM_ORDER)].copy()
    mat = df.pivot_table(index="program", columns="comparison_specimen_type", values="median_delta_vs_primary", aggfunc="median").reindex(
        index=PROGRAM_ORDER, columns=["liver_metastasis", "lymph_node_metastasis"]
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.30, vmax=0.30, aspect="auto")
    ax.set_xticks(np.arange(2))
    ax.set_xticklabels(["liver-primary", "LN-primary"], fontsize=7)
    ax.set_yticks(np.arange(len(PROGRAM_ORDER)))
    ax.set_yticklabels([PROGRAM_LABELS[p] for p in PROGRAM_ORDER], fontsize=7)
    ax.set_title("Patient-matched site deltas", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "E")
    for i, program in enumerate(PROGRAM_ORDER):
        for j, site in enumerate(["liver_metastasis", "lymph_node_metastasis"]):
            val = mat.loc[program, site]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8, color="#FFFFFF" if abs(val) > 0.18 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.05, pad=0.02)
    cbar.set_label("delta rho vs primary", fontsize=6)
    cbar.ax.tick_params(labelsize=6)
    return df


def plot_patient_delta_points(ax: plt.Axes) -> pd.DataFrame:
    df = pd.read_csv(TABLE_DIR / "gse272362_patient_matched_site_deltas.csv")
    keep = ["immune_core", "ifn_mhc", "tumor_aggressive"]
    df = df[
        df["metric"].eq("rho_distance_to_caf_core")
        & df["comparison_specimen_type"].isin(["liver_metastasis", "lymph_node_metastasis"])
        & df["program"].isin(keep)
    ].copy()
    colors = {"immune_core": "#2C7A51", "ifn_mhc": "#4C78A8", "tumor_aggressive": "#B23A48"}
    labels = {"immune_core": "Immune core", "ifn_mhc": "IFN/MHC", "tumor_aggressive": "Tumor aggr."}
    site_offsets = {"liver_metastasis": -0.13, "lymph_node_metastasis": 0.13}
    for i, program in enumerate(keep):
        sub = df[df["program"].eq(program)]
        for site in ["liver_metastasis", "lymph_node_metastasis"]:
            vals = sub[sub["comparison_specimen_type"].eq(site)]["delta_vs_primary"].to_numpy(float)
            x = np.full_like(vals, i + site_offsets[site], dtype=float)
            ax.scatter(x, vals, s=22, color=colors[program], alpha=0.65, edgecolor="white", linewidth=0.4)
            if len(vals):
                ax.scatter(i + site_offsets[site], np.median(vals), s=62, color=colors[program], edgecolor="#222222", linewidth=0.5, zorder=4)
    ax.axhline(0, color="#333333", lw=0.65)
    ax.set_xticks(np.arange(len(keep)))
    ax.set_xticklabels([labels[p] for p in keep], fontsize=7)
    ax.set_ylabel("sample delta rho vs primary", fontsize=7)
    ax.set_title("Matched-patient delta distribution", loc="left", fontsize=10, fontweight="bold")
    panel_label(ax, "F")
    clean_axes(ax)
    ax.text(0.02, 0.96, "left: liver; right: LN", transform=ax.transAxes, fontsize=6.5, va="top", color="#444444")
    return df


def plot_claim_scope(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "G", x=-0.02, y=1.02)
    ax.text(0.02, 0.98, "Interpretation boundary", fontsize=10, fontweight="bold", va="top")
    rows = [
        ("CAF-core organization", "primary + liver + LN tested", "#2C7A51"),
        ("LN immune uncoupling", "strong lead, n=5", "#B23A48"),
        ("Tumor/stromal coupling", "retained in LN", "#4C78A8"),
        ("Clinical subtype claim", "not established", "#8A8F98"),
        ("Mechanism of LN remodeling", "not causal", "#8A8F98"),
    ]
    for i, (claim, scope, color) in enumerate(rows):
        y = 0.82 - i * 0.14
        ax.add_patch(Rectangle((0.02, y - 0.045), 0.96, 0.085, facecolor="#F6F7F9", edgecolor="#DDDDDD", linewidth=0.5))
        ax.add_patch(Rectangle((0.02, y - 0.045), 0.025, 0.085, facecolor=color, edgecolor=color, linewidth=0))
        ax.text(0.06, y, claim, va="center", fontsize=7.2)
        ax.text(0.58, y, scope, va="center", fontsize=7.2, color=color, fontweight="bold")
    return pd.DataFrame(rows, columns=["claim", "scope", "color"])


def write_source_data(*dfs: pd.DataFrame) -> None:
    labels = ["A_site_scale", "B_random_support", "C_subprogram_gradients", "D_decoupling", "E_matched_delta_summary", "F_matched_delta_points", "G_claim_scope"]
    rows: list[dict[str, object]] = []
    for label, df in zip(labels, dfs):
        for _, row in df.iterrows():
            rows.append({"panel": label, **{k: row[k] for k in df.columns if k in row.index and k != "ord"}})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Extended Data Figure 27 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "NC-style metastatic immune-decoupling module focused on GSE272362 primary, liver-metastasis and lymph-node-metastasis spatial contrasts.\n\n"
        "## Panel contract\n\n"
        "- A: site-level sample and spot scale.\n"
        "- B: random-core support for IFN/MHC, immune-core and tumor-aggressive programs by site.\n"
        "- C: CAF-core subprogram distance gradients across primary, liver and LN sites.\n"
        "- D: stromal-tumor coupling, immune-core coupling and immune-decoupling index by site.\n"
        "- E-F: patient-matched primary-to-metastasis deltas for selected programs.\n"
        "- G: interpretation boundary for the LN result.\n\n"
        "## Boundary\n\n"
        "The lymph-node result is a spatial decoupling contrast from five LN metastasis samples. It supports a strong metastatic-site biology lead, not a definitive clinical subtype or causal mechanism.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf')}`\n"
        f"- `{OUT.with_suffix('.svg')}`\n"
        f"- `{OUT.with_suffix('.png')}`\n"
        f"- `{SOURCE_OUT}`\n",
        encoding="utf-8",
    )


def main() -> None:
    fig = plt.figure(figsize=(14.8, 13.2), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[0.82, 1.12, 1.05, 0.62], hspace=0.95, wspace=0.95)
    fig.suptitle("Metastatic-site immune decoupling module", fontsize=15, fontweight="bold", y=0.986)

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:4])
    ax_d = fig.add_subplot(gs[0, 4:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_e = fig.add_subplot(gs[1, 3:6])
    ax_f = fig.add_subplot(gs[2, 0:4])
    ax_g = fig.add_subplot(gs[2, 4:6])
    ax_note = fig.add_subplot(gs[3, :])

    a = plot_site_scale(ax_a)
    b = plot_random_support(ax_b)
    d = plot_decoupling(ax_d)
    c = plot_subprogram_gradients(ax_c)
    e = plot_matched_delta_summary(ax_e)
    f = plot_patient_delta_points(ax_f)
    g = plot_claim_scope(ax_g)

    ax_note.axis("off")
    ax_note.text(
        0.01,
        0.72,
        "Interpretation: GSE272362 supports a metastatic-site contrast in which primary tumors and liver metastases preserve CAF-core-centered immune and tumor-aggressive organization, "
        "whereas lymph-node metastases retain stromal-tumor coupling but selectively weaken immune/IFN coupling to the CAF core.",
        fontsize=8.1,
        color="#333333",
        wrap=True,
    )
    ax_note.text(
        0.01,
        0.28,
        "Boundary: the LN finding is a biologically sharp spatial contrast from five LN metastasis samples and should be framed as a metastatic immune-remodeling lead, not as a clinical subtype or causal mechanism.",
        fontsize=7.5,
        color="#555555",
        wrap=True,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    write_source_data(a, b, c, d, e, f, g)
    write_report()
    print(f"Wrote {OUT.with_suffix('.pdf')}")
    print(f"Wrote {SOURCE_OUT}")
    print(f"Wrote {REPORT_OUT}")


if __name__ == "__main__":
    main()
