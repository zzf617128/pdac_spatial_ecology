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
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT = FIG_DIR / "supplementary_module5_pathology_tcga_tls_boundaries"
SOURCE_OUT = SOURCE_DIR / "Source_Data_supplementary_module5_pathology_tcga_tls_boundaries.csv"
REPORT_OUT = REPORT_DIR / "supplementary_module5_pathology_tcga_tls_boundaries_notes.md"

TARGET_ORDER = ["CAF-myeloid", "tumor aggressive", "IFN/MHC", "immune core"]
TARGET_LABELS = ["CAF-myeloid", "Tumor aggressive", "IFN/MHC", "Immune core"]
TARGET_COLORS = {
    "CAF-myeloid": "#8C2D04",
    "tumor aggressive": "#B23A48",
    "IFN/MHC": "#3B6FB6",
    "immune core": "#2C7A51",
}
FEATURE_ORDER = [
    "red OD",
    "stain density",
    "blue OD",
    "green OD",
    "brightness",
    "red-blue",
    "purple fraction",
    "edge density",
]
TCGA_AXES = [
    "myCAF_matrix",
    "SPP1_TAM",
    "TGFb_EMT",
    "matrix_integrin",
    "DC_APC",
    "T_NK",
    "B_plasma",
    "bulk_decoupling_like_index",
]
TCGA_LABELS = ["myCAF", "SPP1/TAM", "TGF-beta", "Matrix-int", "DC/APC", "T/NK", "B/plasma", "Decoupling"]
TLS_CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treat-naive",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "gse274557_primary": "GSE274557 primary",
    "gse274557_liver": "GSE274557 liver",
    "gse274557_lung": "GSE274557 lung",
    "gse274557_peritoneal": "GSE274557 peritoneal",
}


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(x, y, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="bottom")


def clean_axes(ax: plt.Axes, axis: str = "y") -> None:
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(labelsize=7)
    ax.grid(axis=axis, color="#E7E7E7", linewidth=0.55, zorder=0)


def load_cv_summary(cv: pd.DataFrame) -> pd.DataFrame:
    observed = cv[cv["metric_scope"].eq("observed_all")].copy()
    null = (
        cv[cv["metric_scope"].eq("permuted_all")]
        .groupby(["target", "target_label"], as_index=False)
        .agg(
            null_median=("spearman_rho", "median"),
            null_q05=("spearman_rho", lambda s: s.quantile(0.05)),
            null_q95=("spearman_rho", lambda s: s.quantile(0.95)),
        )
    )
    out = observed.merge(null, on=["target", "target_label"], how="left")
    out["target_label"] = pd.Categorical(out["target_label"], categories=TARGET_ORDER, ordered=True)
    return out.sort_values("target_label").reset_index(drop=True)


def plot_he_cv(ax: plt.Axes, cv: pd.DataFrame) -> pd.DataFrame:
    summary = load_cv_summary(cv)
    x = np.arange(len(summary))
    colors = [TARGET_COLORS[str(label)] for label in summary["target_label"]]
    ax.bar(x, summary["spearman_rho"], color=colors, width=0.62, zorder=3)
    ax.errorbar(
        x,
        summary["null_median"],
        yerr=[
            summary["null_median"] - summary["null_q05"],
            summary["null_q95"] - summary["null_median"],
        ],
        fmt="o",
        color="#222222",
        markersize=4,
        capsize=3,
        lw=0.9,
        zorder=4,
    )
    for i, row in summary.iterrows():
        ax.text(i, row["spearman_rho"] + 0.018, f"{row['spearman_rho']:.2f}", ha="center", fontsize=6.8)
    ax.axhline(0, color="#333333", lw=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(TARGET_LABELS, rotation=28, ha="right", fontsize=7)
    ax.set_ylim(-0.05, 0.52)
    ax.set_ylabel("held-out Spearman rho", fontsize=7.5)
    ax.set_title("H&E patch bridge exceeds shuffled target", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "A")
    clean_axes(ax)
    ax.text(0.02, 0.94, "dot: shuffled-target null", transform=ax.transAxes, fontsize=6.6, color="#555555", va="top")
    return summary


def plot_he_feature_heatmap(ax: plt.Axes, corr: pd.DataFrame) -> pd.DataFrame:
    sub = corr[corr["feature_label"].isin(FEATURE_ORDER)].copy()
    mat = (
        sub.pivot_table(index="feature_label", columns="target_label", values="median_rho", aggfunc="median")
        .reindex(index=FEATURE_ORDER, columns=TARGET_ORDER)
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-0.42, vmax=0.42, aspect="auto")
    ax.set_xticks(np.arange(len(TARGET_ORDER)))
    ax.set_xticklabels(TARGET_LABELS, rotation=30, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(FEATURE_ORDER)))
    ax.set_yticklabels(FEATURE_ORDER, fontsize=7)
    ax.set_title("H&E feature-direction map", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "B")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.8)
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return sub


def plot_tcga_correlation(ax: plt.Axes, corr: pd.DataFrame) -> pd.DataFrame:
    sub = corr[corr["axis_x"].isin(TCGA_AXES) & corr["axis_y"].isin(TCGA_AXES)].copy()
    mat = sub.pivot_table(index="axis_y", columns="axis_x", values="spearman_rho", aggfunc="median").reindex(
        index=TCGA_AXES, columns=TCGA_AXES
    )
    im = ax.imshow(mat.to_numpy(float), cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(len(TCGA_AXES)))
    ax.set_xticklabels(TCGA_LABELS, rotation=35, ha="right", fontsize=6.8)
    ax.set_yticks(np.arange(len(TCGA_AXES)))
    ax.set_yticklabels(TCGA_LABELS, fontsize=6.8)
    ax.set_title("TCGA PAAD bulk-axis context", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "C")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            val = mat.iloc[i, j]
            if np.isfinite(val):
                ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5.0, color="#FFFFFF" if abs(val) > 0.65 else "#222222")
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return sub


def plot_tcga_survival(ax: plt.Axes, surv: pd.DataFrame) -> pd.DataFrame:
    df = surv.sort_values("cox_hr_per_sd", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))
    ax.axvline(1.0, color="#333333", lw=0.7)
    colors = np.where(df["cox_p"] < 0.05, "#B23A48", "#8A8F98")
    ax.hlines(y, df["cox_ci95_low"], df["cox_ci95_high"], color=colors, lw=1.6)
    ax.scatter(df["cox_hr_per_sd"], y, s=36, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax.set_xscale("log")
    ax.set_xlim(0.78, 1.82)
    ax.set_xticks([0.8, 1.0, 1.25, 1.5])
    ax.set_xticklabels(["0.8", "1.0", "1.25", "1.5"])
    ax.xaxis.set_minor_locator(mpl.ticker.NullLocator())
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_yticks(y)
    ax.set_yticklabels(df["axis_label"], fontsize=7)
    ax.set_xlabel("Cox HR per 1 s.d.", fontsize=7.5)
    ax.set_title("Exploratory TCGA survival context", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "D")
    clean_axes(ax)
    ax.text(0.98, 0.95, "univariable, non-spatial", transform=ax.transAxes, fontsize=6.5, color="#555555", ha="right", va="top")
    return df


def plot_tls_compatibility(ax: plt.Axes, tls: pd.DataFrame) -> pd.DataFrame:
    df = tls.copy()
    df["context_label"] = df["cohort_context"].map(TLS_CONTEXT_LABELS).fillna(df["cohort_context"])
    df = df.sort_values("tls_compatible_fraction", ascending=True).reset_index(drop=True)
    y = np.arange(len(df))
    ax.barh(y - 0.16, df["tls_loose_compatible_fraction"], height=0.28, color="#A8C6D8", label="loose")
    ax.barh(y + 0.16, df["tls_compatible_fraction"], height=0.28, color="#1F5F8B", label="stringent")
    ax.set_yticks(y)
    ax.set_yticklabels(df["context_label"], fontsize=6.7)
    ax.set_xlim(0, 1.02)
    ax.set_xlabel("compatible fraction", fontsize=7.5)
    ax.set_title("TLS-maturity stress test", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "E")
    clean_axes(ax)
    ax.legend(frameon=False, fontsize=6.8, loc="lower right")
    return df


def plot_tls_compartment_heatmap(ax: plt.Axes, tls: pd.DataFrame) -> pd.DataFrame:
    metrics = [
        "median_tls_chemokine",
        "median_b_cell",
        "median_t_cell",
        "median_plasma_cell",
        "median_dc_apc",
        "median_fdc_gc_like",
        "median_caf_myeloid",
    ]
    labels = ["TLS chem", "B", "T", "plasma", "DC/APC", "FDC/GC", "CAF-myeloid"]
    df = tls.copy()
    df["context_label"] = df["cohort_context"].map(TLS_CONTEXT_LABELS).fillna(df["cohort_context"])
    mat = df.set_index("context_label")[metrics]
    im = ax.imshow(mat.to_numpy(float), cmap="YlGnBu", vmin=-0.25, vmax=2.25, aspect="auto")
    ax.set_xticks(np.arange(len(metrics)))
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=6.7)
    ax.set_yticks(np.arange(len(mat.index)))
    ax.set_yticklabels(mat.index, fontsize=6.7)
    ax.set_title("Immune-hub markers do not equal mature TLS", loc="left", fontsize=10.5, fontweight="bold")
    panel_label(ax, "F")
    cbar = plt.colorbar(im, ax=ax, fraction=0.030, pad=0.014)
    cbar.ax.tick_params(labelsize=6)
    return df[["cohort_context", *metrics, "median_caf_minus_tls_maturity"]]


def draw_claim_boundary(ax: plt.Axes) -> pd.DataFrame:
    ax.axis("off")
    panel_label(ax, "G", x=-0.03, y=1.02)
    ax.text(0.02, 0.97, "Boundary map for reviewer-risk claims", fontsize=10.5, fontweight="bold", va="top")
    rows = [
        ("Pathology bridge", "supported", "patch features track spatial programs"),
        ("TCGA context", "context only", "bulk, non-spatial external layer"),
        ("TLS framing", "bounded", "strict mature TLS support is limited"),
        ("Survival/prognosis", "not claimed", "exploratory univariable TCGA only"),
        ("Causal signaling", "not claimed", "requires perturbation"),
    ]
    colors = {"supported": "#2C7A51", "context only": "#4C78A8", "bounded": "#C77C2D", "not claimed": "#8A8F98"}
    y0 = 0.78
    for i, (domain, status, note) in enumerate(rows):
        y = y0 - i * 0.14
        ax.add_patch(Rectangle((0.04, y - 0.05), 0.28, 0.09, facecolor="#F4F6F8", edgecolor=colors[status], linewidth=1.0))
        ax.text(0.18, y - 0.005, domain, ha="center", va="center", fontsize=6.7, fontweight="bold")
        ax.text(0.37, y - 0.005, status, ha="left", va="center", fontsize=6.8, color=colors[status], fontweight="bold")
        ax.text(0.57, y - 0.005, note, ha="left", va="center", fontsize=6.6, color="#333333")
    ax.text(0.04, 0.04, "Module role: strengthens translational plausibility while narrowing unsupported claims.", fontsize=6.8, color="#555555")
    return pd.DataFrame(rows, columns=["item", "status", "note"])


def write_source_data(*tables: pd.DataFrame) -> None:
    names = [
        "A_he_cv",
        "B_he_features",
        "C_tcga_bulk_context",
        "D_tcga_survival",
        "E_tls_compatibility",
        "F_tls_compartment",
        "G_claim_boundary",
    ]
    rows: list[dict[str, object]] = []
    for name, table in zip(names, tables):
        for _, row in table.iterrows():
            item = row.get("target_label", row.get("axis_label", row.get("cohort_context", row.get("item", row.get("axis_x", "")))))
            for col, value in row.items():
                if col in {"target", "target_label", "axis", "axis_label", "cohort_context", "item", "note", "status"}:
                    continue
                if isinstance(value, (int, float, np.number)) and np.isfinite(value):
                    rows.append({"panel": name, "item": item, "metric": col, "value": float(value)})
                elif col in {"feature", "feature_label", "axis_x", "axis_y", "context_label"} and isinstance(value, str):
                    rows.append({"panel": name, "item": item, "metric": col, "value": value})
    pd.DataFrame(rows).to_csv(SOURCE_OUT, index=False)


def write_report() -> None:
    REPORT_OUT.write_text(
        "# Supplementary Module 5 Notes\n\n"
        "Last updated: 2026-06-28\n\n"
        "## Figure role\n\n"
        "Module 5 is a reviewer-risk and translational-context module. It shows that H&E patch features provide an "
        "exploratory pathology bridge, TCGA PAAD provides non-spatial bulk context, and TLS stress testing argues "
        "against reframing the story as mature TLS biology.\n\n"
        "## Panel contract\n\n"
        "- A-B: H&E feature-model performance and feature directions.\n"
        "- C-D: TCGA PAAD bulk-axis correlation and exploratory survival context.\n"
        "- E-F: TLS compatibility and immune-compartment stress testing.\n"
        "- G: claim-boundary map.\n\n"
        "## Boundary\n\n"
        "This module supports translational plausibility and claim control. It does not establish clinical-grade "
        "pathology prediction, spatially localized TCGA biology, mature TLS biology, prognosis, therapy response or causal signaling.\n\n"
        "## Outputs\n\n"
        f"- `{OUT.with_suffix('.pdf')}`\n"
        f"- `{OUT.with_suffix('.svg')}`\n"
        f"- `{OUT.with_suffix('.png')}`\n"
        f"- `{SOURCE_OUT}`\n",
        encoding="utf-8",
    )


def main() -> None:
    he_cv = pd.read_csv(TABLE_DIR / "mvp_he_patch_grouped_cv_metrics.csv")
    he_corr = pd.read_csv(TABLE_DIR / "mvp_he_patch_feature_correlation_summary.csv")
    tcga_corr = pd.read_csv(TABLE_DIR / "tcga_paad_bulk_context_axis_correlations.csv")
    tcga_surv = pd.read_csv(TABLE_DIR / "tcga_paad_survival_context_summary.csv")
    tls = pd.read_csv(TABLE_DIR / "tls_maturity_stress_test_context_summary.csv")

    fig = plt.figure(figsize=(15.2, 14.1), constrained_layout=False)
    gs = GridSpec(4, 6, figure=fig, height_ratios=[0.95, 1.05, 1.0, 0.82], hspace=0.98, wspace=1.05)
    fig.suptitle("Supplementary Module 5 | Pathology bridge, TCGA context and TLS/claim boundaries", fontsize=15, fontweight="bold", y=0.988)

    ax_a = fig.add_subplot(gs[0, 0:2])
    ax_b = fig.add_subplot(gs[0, 2:6])
    ax_c = fig.add_subplot(gs[1, 0:3])
    ax_d = fig.add_subplot(gs[1, 3:6])
    ax_e = fig.add_subplot(gs[2, 0:3])
    ax_f = fig.add_subplot(gs[2, 3:6])
    ax_g = fig.add_subplot(gs[3, :])

    a = plot_he_cv(ax_a, he_cv)
    b = plot_he_feature_heatmap(ax_b, he_corr)
    c = plot_tcga_correlation(ax_c, tcga_corr)
    d = plot_tcga_survival(ax_d, tcga_surv)
    e = plot_tls_compatibility(ax_e, tls)
    f = plot_tls_compartment_heatmap(ax_f, tls)
    g = draw_claim_boundary(ax_g)

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
