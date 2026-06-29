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


ROOT = Path(__file__).resolve().parents[1]
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [TABLE_DIR, FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT_PER_SAMPLE = TABLE_DIR / "tls_maturity_stress_test_per_sample.csv"
OUT_CONTEXT = TABLE_DIR / "tls_maturity_stress_test_context_summary.csv"
OUT_SOURCE = SOURCE_DIR / "Source_Data_Extended_Data_Fig_23_TLS_maturity_stress_test.csv"
OUT_REPORT = REPORT_DIR / "gap3_tls_maturity_resolution_plan.md"
OUT_FIG = FIG_DIR / "extended_data_figure23_tls_maturity_stress_test"
STRINGENT_GATE = 0.25

CORE_SCORE = "score_caf_myeloid_barrier"
MODULES = [
    ("tls_chemokine", "TLS chemokine"),
    ("b_cell", "B cell"),
    ("t_cell", "T cell"),
    ("plasma_cell", "Plasma cell"),
    ("dc_apc", "DC/APC"),
    ("fdc_gc_like", "FDC/GC-like"),
    ("immune_hub_core", "Immune hub core"),
    ("immune_hub_maturity", "Immune maturity-like"),
    ("caf_myeloid_barrier", "CAF-myeloid"),
]

CONTEXT_MEDIAN_COLS = {
    "tls_chemokine": "median_tls_chemokine",
    "b_cell": "median_b_cell",
    "t_cell": "median_t_cell",
    "plasma_cell": "median_plasma_cell",
    "dc_apc": "median_dc_apc",
    "fdc_gc_like": "median_fdc_gc_like",
    "immune_hub_core": "median_immune_hub_core",
    "immune_hub_maturity": "median_immune_hub_maturity",
    "caf_myeloid_barrier": "median_caf_myeloid",
}

CONTEXT_ORDER = [
    "post_neoadjuvant_sections",
    "treatment_naive_primary",
    "external_anchor",
    "primary_tumor",
    "liver_metastasis",
    "lymph_node_metastasis",
    "normal_pancreas",
    "gse274557_primary",
    "gse274557_liver",
    "gse274557_lung",
    "gse274557_peritoneal",
]

CONTEXT_LABELS = {
    "post_neoadjuvant_sections": "post-NACT",
    "treatment_naive_primary": "treatment-\nnaive",
    "external_anchor": "external\nanchor",
    "primary_tumor": "primary",
    "liver_metastasis": "liver met",
    "lymph_node_metastasis": "LN met",
    "normal_pancreas": "normal",
    "gse274557_primary": "GSE274557\nprimary",
    "gse274557_liver": "GSE274557\nliver",
    "gse274557_lung": "GSE274557\nlung",
    "gse274557_peritoneal": "GSE274557\nperitoneal",
}


def zscore(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    sd = values.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=values.index, dtype=float)
    return (values - values.mean()) / sd


def mvp_context(row: pd.Series) -> str:
    dataset = str(row.get("dataset_id", ""))
    if dataset == "GSE282302":
        return "post_neoadjuvant_sections"
    if dataset == "GSE274103":
        return "treatment_naive_primary"
    if dataset == "GSE235315":
        return "external_anchor"
    return dataset or "unknown"


def gse274557_context(row: pd.Series) -> str:
    tissue = str(row.get("tissue", "")).lower()
    if "primary" in tissue or "pdac" in tissue and "met" not in tissue:
        return "gse274557_primary"
    if "liver" in tissue:
        return "gse274557_liver"
    if "lung" in tissue:
        return "gse274557_lung"
    if "peritone" in tissue:
        return "gse274557_peritoneal"
    return "gse274557_other"


def read_required(path: Path, extra_cols: list[str], source: str) -> pd.DataFrame:
    score_cols = [CORE_SCORE] + [f"score_{m}" for m, _ in MODULES]
    z_cols = [f"z_{m}" for m, _ in MODULES if m not in {"immune_hub_core", "immune_hub_maturity", "caf_myeloid_barrier"}]
    usecols = [c for c in ["dataset_id", "sample_id", "patient_id", *extra_cols, *score_cols, *z_cols] if c]
    df = pd.read_csv(path, usecols=lambda c: c in set(usecols), dtype={"patient_id": "string"}, low_memory=False)
    df["source_layer"] = source
    return df


def load_all() -> pd.DataFrame:
    frames: list[pd.DataFrame] = []

    mvp = read_required(
        TABLE_DIR / "mvp_spot_level_scores_with_edge_qc.csv",
        ["edge_or_background_risk"],
        "discovery_support",
    )
    if "edge_or_background_risk" in mvp.columns:
        mvp = mvp[~mvp["edge_or_background_risk"].astype(bool)].copy()
    mvp["cohort_context"] = mvp.apply(mvp_context, axis=1)
    frames.append(mvp)

    gse272 = read_required(
        TABLE_DIR / "gse272362_rds_spot_level_scores.csv",
        ["specimen_type"],
        "gse272362",
    )
    gse272["cohort_context"] = gse272["specimen_type"].astype(str)
    frames.append(gse272)

    gse274557 = read_required(
        TABLE_DIR / "gse274557_full_spot_scores.csv",
        ["tissue"],
        "gse274557",
    )
    gse274557["cohort_context"] = gse274557.apply(gse274557_context, axis=1)
    frames.append(gse274557)

    df = pd.concat(frames, ignore_index=True, sort=False)
    return df[df["cohort_context"].isin(CONTEXT_ORDER)].copy()


def summarize_sample(group: pd.DataFrame) -> dict[str, object]:
    sample_id = str(group["sample_id"].iloc[0])
    context = str(group["cohort_context"].iloc[0])
    source_layer = str(group["source_layer"].iloc[0])
    core_threshold = float(group[CORE_SCORE].quantile(0.90))
    is_core = group[CORE_SCORE] >= core_threshold
    out: dict[str, object] = {
        "source_layer": source_layer,
        "cohort_context": context,
        "sample_id": sample_id,
        "n_spots": int(len(group)),
        "n_caf_core_spots": int(is_core.sum()),
        "caf_core_threshold": core_threshold,
    }
    support_flags: dict[str, bool] = {}
    for module, label in MODULES:
        col = f"z_{module}"
        if col not in group.columns:
            score_col = f"score_{module}"
            values = zscore(group[score_col])
        else:
            values = pd.to_numeric(group[col], errors="coerce")
        core_mean = float(values[is_core].mean())
        noncore_mean = float(values[~is_core].mean())
        enrich = core_mean - noncore_mean
        out[f"{module}_core_enrichment"] = enrich
        out[f"{module}_core_mean_z"] = core_mean
        support_flags[module] = bool(enrich > 0)

    lymphoid_max = max(
        float(out["b_cell_core_enrichment"]),
        float(out["t_cell_core_enrichment"]),
        float(out["plasma_cell_core_enrichment"]),
    )
    lymphoid_support = support_flags["b_cell"] or support_flags["t_cell"] or support_flags["plasma_cell"]
    mature_support_loose = support_flags["tls_chemokine"] and lymphoid_support and support_flags["fdc_gc_like"]
    mature_support_stringent = (
        float(out["tls_chemokine_core_enrichment"]) > STRINGENT_GATE
        and lymphoid_max > STRINGENT_GATE
        and float(out["fdc_gc_like_core_enrichment"]) > STRINGENT_GATE
        and float(out["immune_hub_maturity_core_enrichment"]) > STRINGENT_GATE
    )
    compartment_score = (
        float(out["tls_chemokine_core_enrichment"]) + lymphoid_max
        + float(out["fdc_gc_like_core_enrichment"])
    ) / 3.0
    out["tls_lymphoid_any_positive"] = lymphoid_support
    out["tls_maturity_compatible_loose"] = mature_support_loose
    out["tls_maturity_compatible"] = mature_support_stringent
    out["tls_lymphoid_max_core_enrichment"] = lymphoid_max
    out["tls_three_compartment_score"] = compartment_score
    out["caf_minus_tls_maturity"] = float(out["caf_myeloid_barrier_core_enrichment"]) - float(out["immune_hub_maturity_core_enrichment"])
    return out


def build_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    df = load_all()
    per_sample = pd.DataFrame(
        [summarize_sample(g) for _, g in df.groupby(["source_layer", "cohort_context", "sample_id"], sort=False)]
    )
    context = (
        per_sample.groupby("cohort_context")
        .agg(
            n_samples=("sample_id", "nunique"),
            tls_compatible_samples=("tls_maturity_compatible", "sum"),
            tls_compatible_fraction=("tls_maturity_compatible", "mean"),
            tls_loose_compatible_samples=("tls_maturity_compatible_loose", "sum"),
            tls_loose_compatible_fraction=("tls_maturity_compatible_loose", "mean"),
            median_tls_chemokine=("tls_chemokine_core_enrichment", "median"),
            median_b_cell=("b_cell_core_enrichment", "median"),
            median_t_cell=("t_cell_core_enrichment", "median"),
            median_plasma_cell=("plasma_cell_core_enrichment", "median"),
            median_dc_apc=("dc_apc_core_enrichment", "median"),
            median_fdc_gc_like=("fdc_gc_like_core_enrichment", "median"),
            median_immune_hub_core=("immune_hub_core_core_enrichment", "median"),
            median_immune_hub_maturity=("immune_hub_maturity_core_enrichment", "median"),
            median_caf_myeloid=("caf_myeloid_barrier_core_enrichment", "median"),
            median_three_compartment=("tls_three_compartment_score", "median"),
            median_caf_minus_tls_maturity=("caf_minus_tls_maturity", "median"),
        )
        .reset_index()
    )
    context["context_order"] = context["cohort_context"].map({c: i for i, c in enumerate(CONTEXT_ORDER)})
    context = context.sort_values("context_order").drop(columns="context_order")
    return per_sample, context


def plot(per_sample: pd.DataFrame, context: pd.DataFrame) -> None:
    heat_modules = [
        ("tls_chemokine", "TLS\nchemokine"),
        ("b_cell", "B cell"),
        ("t_cell", "T cell"),
        ("plasma_cell", "Plasma"),
        ("fdc_gc_like", "FDC/GC-\nlike"),
        ("immune_hub_core", "Immune\ncore"),
        ("immune_hub_maturity", "Immune\nmaturity"),
        ("caf_myeloid_barrier", "CAF-\nmyeloid"),
    ]
    contexts = context["cohort_context"].tolist()
    matrix = []
    for ctx in contexts:
        row = context[context["cohort_context"] == ctx].iloc[0]
        matrix.append([row[CONTEXT_MEDIAN_COLS[m]] for m, _ in heat_modules])
    matrix = np.array(matrix, dtype=float)

    fig = plt.figure(figsize=(11.6, 7.6))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.55, 1.0], height_ratios=[1.05, 0.95], wspace=0.34, hspace=0.42)

    ax_a = fig.add_subplot(gs[0, 0])
    vmax = max(1.0, float(np.nanmax(np.abs(matrix))))
    im = ax_a.imshow(matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax_a.set_xticks(np.arange(len(heat_modules)))
    ax_a.set_xticklabels([label for _, label in heat_modules], fontsize=7, rotation=35, ha="right")
    ax_a.set_yticks(np.arange(len(contexts)))
    ax_a.set_yticklabels([CONTEXT_LABELS.get(c, c) for c in contexts], fontsize=7)
    ax_a.set_title("A  CAF-core enrichment of TLS-maturity modules", loc="left", fontsize=10, fontweight="bold")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax_a.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=5.5, color="white" if abs(matrix[i, j]) > vmax * 0.55 else "#222222")
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.028, pad=0.02)
    cbar.ax.tick_params(labelsize=7)
    cbar.set_label("Median core enrichment", fontsize=7)

    ax_b = fig.add_subplot(gs[0, 1])
    y = np.arange(len(context))
    frac = context["tls_compatible_fraction"].to_numpy(float)
    ax_b.barh(y, frac, color="#7B68A6", height=0.58)
    ax_b.set_yticks(y)
    ax_b.set_yticklabels([CONTEXT_LABELS.get(c, c) for c in contexts], fontsize=7)
    ax_b.invert_yaxis()
    ax_b.set_xlim(0, 1)
    ax_b.set_xlabel("Samples passing stringent TLS-maturity gate", fontsize=8)
    ax_b.set_title("B  Stringent three-compartment TLS gate", loc="left", fontsize=10, fontweight="bold")
    ax_b.spines[["top", "right"]].set_visible(False)
    ax_b.tick_params(labelsize=7)
    for i, (_, row) in enumerate(context.iterrows()):
        ax_b.text(row["tls_compatible_fraction"] + 0.025, i, f"{int(row['tls_compatible_samples'])}/{int(row['n_samples'])}", va="center", fontsize=7)

    ax_c = fig.add_subplot(gs[1, 0])
    order = [c for c in CONTEXT_ORDER if c in per_sample["cohort_context"].unique()]
    box_data = [per_sample.loc[per_sample["cohort_context"] == c, "tls_three_compartment_score"].dropna().to_numpy(float) for c in order]
    bp = ax_c.boxplot(box_data, patch_artist=True, showfliers=False, widths=0.55)
    for patch in bp["boxes"]:
        patch.set(facecolor="#D8DEE8", edgecolor="#4A5568", linewidth=0.8)
    for element in ["whiskers", "caps", "medians"]:
        for line in bp[element]:
            line.set(color="#4A5568", linewidth=0.8)
    rng = np.random.default_rng(7)
    for i, data in enumerate(box_data, start=1):
        if len(data):
            x = rng.normal(i, 0.045, size=len(data))
            ax_c.scatter(x, data, s=9, alpha=0.35, color="#2F6F8F", linewidth=0)
    ax_c.axhline(0, color="#333333", lw=0.6)
    ax_c.set_xticks(np.arange(1, len(order) + 1))
    ax_c.set_xticklabels([CONTEXT_LABELS.get(c, c) for c in order], fontsize=7, rotation=35, ha="right")
    ax_c.set_ylabel("TLS three-compartment score", fontsize=8)
    ax_c.set_title("C  Chemokine + lymphoid + FDC/GC-like support", loc="left", fontsize=10, fontweight="bold")
    ax_c.tick_params(labelsize=7)
    ax_c.spines[["top", "right"]].set_visible(False)

    ax_d = fig.add_subplot(gs[1, 1])
    ax_d.axvline(0, color="#333333", lw=0.6)
    delta = context["median_caf_minus_tls_maturity"].to_numpy(float)
    colors = np.where(delta >= 0, "#B23A48", "#4C78A8")
    ax_d.barh(y, delta, color=colors, height=0.58)
    ax_d.set_yticks(y)
    ax_d.set_yticklabels([CONTEXT_LABELS.get(c, c) for c in contexts], fontsize=7)
    ax_d.invert_yaxis()
    ax_d.set_xlabel("CAF-myeloid minus immune-maturity enrichment", fontsize=8)
    ax_d.set_title("D  CAF niche versus TLS-maturity balance", loc="left", fontsize=10, fontweight="bold")
    ax_d.tick_params(labelsize=7)
    ax_d.spines[["top", "right"]].set_visible(False)

    fig.suptitle("TLS-maturity stress test for the CAF-myeloid spatial niche", fontsize=12, fontweight="bold", y=0.985)
    fig.text(
        0.01,
        0.012,
        f"Stringent TLS-maturity compatibility requires CAF-core enrichment > {STRINGENT_GATE:.2f} for TLS chemokine, lymphoid/plasma, FDC/GC-like and immune-maturity modules; module scores are expression-derived and not histologic TLS annotation.",
        fontsize=7,
        color="#444444",
    )
    for suffix in [".pdf", ".svg", ".png"]:
        fig.savefig(OUT_FIG.with_suffix(suffix), dpi=450 if suffix == ".png" else None, bbox_inches="tight")
    plt.close(fig)


def write_report(per_sample: pd.DataFrame, context: pd.DataFrame) -> None:
    total = int(per_sample["sample_id"].nunique())
    compatible = int(per_sample["tls_maturity_compatible"].sum())
    loose = int(per_sample["tls_maturity_compatible_loose"].sum())
    tumor_context = context[context["cohort_context"] != "normal_pancreas"].copy()
    best = tumor_context.sort_values("tls_compatible_fraction", ascending=False).iloc[0]
    lines = [
        "# Gap 3: TLS Maturity Resolution Plan",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Reviewer-facing problem",
        "",
        "The project began near a TLS framing, but the current spatial evidence is stronger for a CAF-myeloid stromal niche than for mature tertiary lymphoid structures. A purely defensive statement is weak; the issue should be resolved by testing whether CAF-core regions satisfy a mature TLS-compatible expression gate.",
        "",
        "## Resolution strategy",
        "",
        f"We tested CAF-core enrichment of TLS chemokine, B-cell, T-cell, plasma-cell, DC/APC, FDC/germinal-center-like, immune-hub-core, immune-maturity-like and CAF-myeloid modules across discovery/support Visium, GSE272362 primary/metastatic Visium and GSE274557 external metastatic Visium layers. A sample was considered stringently TLS-maturity compatible only if CAF cores showed enrichment > {STRINGENT_GATE:.2f} for TLS chemokine, at least one lymphoid/plasma module, the FDC/GC-like module and the immune-maturity-like module. A loose positive-only gate was retained for audit.",
        "",
        "## Result",
        "",
        f"- Stringent TLS-maturity compatible samples: {compatible}/{total}.",
        f"- Loose positive-only TLS-compatible samples: {loose}/{total}.",
        f"- Highest tumor context-level stringent compatibility: {best['cohort_context']} ({int(best['tls_compatible_samples'])}/{int(best['n_samples'])}).",
        "- The analysis supports inflammatory/immune CAF-core organization but does not support reframing the central claim as mature TLS biology without histologic TLS annotation or FDC/germinal-center validation.",
        "",
        "## Claim after resolution",
        "",
        "The manuscript should describe a CAF-myeloid spatial niche with immune-hub features, not a mature TLS program. TLS language can be used only as a tested alternative that was not sufficiently supported by the current data.",
        "",
        "## Outputs",
        "",
        f"- `{OUT_PER_SAMPLE.relative_to(ROOT)}`",
        f"- `{OUT_CONTEXT.relative_to(ROOT)}`",
        f"- `{OUT_SOURCE.relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.pdf').relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.svg').relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.png').relative_to(ROOT)}`",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    per_sample, context = build_tables()
    per_sample.to_csv(OUT_PER_SAMPLE, index=False)
    context.to_csv(OUT_CONTEXT, index=False)
    context.to_csv(OUT_SOURCE, index=False)
    plot(per_sample, context)
    write_report(per_sample, context)
    print(f"Wrote {OUT_PER_SAMPLE}")
    print(f"Wrote {OUT_CONTEXT}")
    print(f"Wrote {OUT_REPORT}")
    print(f"Wrote {OUT_FIG.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
