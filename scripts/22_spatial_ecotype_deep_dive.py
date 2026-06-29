from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.cluster import KMeans
from sklearn.decomposition import NMF
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE = "22_spatial_ecotype_deep_dive"
RANDOM_SEED = 20260624

PROGRAMS = [
    ("score_mycaf", "myCAF"),
    ("score_icaf", "iCAF"),
    ("score_apcaf", "apCAF"),
    ("score_pan_caf", "pan-CAF"),
    ("score_myeloid", "myeloid"),
    ("score_spp1_tam", "SPP1/TREM2 TAM"),
    ("score_tgfb_pathway", "TGF-beta"),
    ("score_emt_invasion", "EMT/invasion"),
    ("score_hypoxia", "hypoxia"),
    ("score_pdac_basal_like", "basal-like"),
    ("score_pdac_classical_like", "classical-like"),
    ("score_tumor_epithelial", "tumor epithelial"),
    ("score_tumor_aggressive", "tumor aggressive"),
    ("score_ifn_antigen_presentation", "IFN/MHC"),
    ("score_immune_hub_core", "immune core"),
    ("score_t_cell", "T cell"),
    ("score_b_cell", "B cell"),
    ("score_dc_apc", "DC/APC"),
    ("score_plasma_cell", "plasma cell"),
    ("score_neural_schwann", "neural/Schwann"),
]
CORE_PROGRAMS = [
    "myCAF",
    "iCAF",
    "apCAF",
    "myeloid",
    "SPP1/TREM2 TAM",
    "TGF-beta",
    "EMT/invasion",
    "hypoxia",
    "basal-like",
    "tumor aggressive",
    "IFN/MHC",
    "immune core",
    "T cell",
    "B cell",
    "DC/APC",
    "plasma cell",
]
STROMAL_TUMOR_LABELS = ["myCAF", "myeloid", "SPP1/TREM2 TAM", "TGF-beta", "EMT/invasion", "tumor aggressive"]
IMMUNE_LABELS = ["IFN/MHC", "immune core", "T cell", "B cell", "DC/APC"]
INTERFACE_LABELS = [
    "tumor aggressive",
    "EMT/invasion",
    "TGF-beta",
    "SPP1/TREM2 TAM",
    "IFN/MHC",
    "immune core",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_status(stage: str, status: str, payload: dict) -> None:
    base = {
        "stage": stage,
        "status": status,
        "timestamp_utc": now_iso(),
        "n_errors": 0,
        "critical_errors": [],
        "noncritical_warnings": [],
        "next_manual_check": [],
    }
    base.update(payload)
    path = PROJECT_ROOT / f"results/logs/stage_{stage}_status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def zscore(values: pd.Series) -> pd.Series:
    arr = values.to_numpy(float)
    sd = np.nanstd(arr)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.zeros(len(values)), index=values.index)
    return pd.Series((arr - np.nanmean(arr)) / sd, index=values.index)


def nearest_distance(points: np.ndarray, centers: np.ndarray) -> np.ndarray:
    if len(centers) == 0:
        return np.full(len(points), np.nan)
    nn = NearestNeighbors(n_neighbors=1)
    nn.fit(centers)
    return nn.kneighbors(points, return_distance=True)[0][:, 0]


def median_neighbor_distance(points: np.ndarray) -> float:
    if len(points) < 3:
        return 1.0
    nn = NearestNeighbors(n_neighbors=2)
    nn.fit(points)
    d = nn.kneighbors(points, return_distance=True)[0][:, 1]
    med = float(np.nanmedian(d))
    return med if np.isfinite(med) and med > 0 else 1.0


def load_spots() -> pd.DataFrame:
    cols = ["dataset_id", "sample_id", "patient_id", "specimen_type", "barcode", "x_pixel", "y_pixel"]
    score_cols = sorted(set([col for col, _ in PROGRAMS] + ["score_caf_myeloid_barrier"]))

    mvp = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv")
    mvp = mvp[~mvp["edge_or_background_risk"].astype(str).str.lower().eq("true")].copy()
    keep_mvp = [c for c in cols + score_cols if c in mvp.columns]
    mvp = mvp[keep_mvp].copy()
    mvp["cohort_context"] = np.where(
        mvp["dataset_id"].eq("GSE274103"),
        "treatment_naive_primary",
        "post_neoadjuvant_sections",
    )

    gse = pd.read_csv(PROJECT_ROOT / "results/tables/gse272362_rds_spot_level_scores.csv")
    keep_gse = [c for c in cols + score_cols if c in gse.columns]
    gse = gse[keep_gse].copy()
    gse["cohort_context"] = gse["specimen_type"].astype(str)

    return pd.concat([mvp, gse], ignore_index=True)


def summarize_sample(sample: pd.DataFrame) -> tuple[dict, list[dict]]:
    sample = sample.copy()
    points = sample[["x_pixel", "y_pixel"]].to_numpy(float)
    med_nn = median_neighbor_distance(points)

    for col, _ in PROGRAMS:
        sample[f"z_{col}"] = zscore(sample[col])

    caf_threshold = float(np.nanpercentile(sample["score_caf_myeloid_barrier"], 90))
    tumor_threshold = float(np.nanpercentile(sample["score_tumor_epithelial"], 80))
    caf_core = sample["score_caf_myeloid_barrier"].to_numpy(float) >= caf_threshold
    tumor_high = sample["score_tumor_epithelial"].to_numpy(float) >= tumor_threshold
    d_caf = nearest_distance(points, points[caf_core]) / med_nn
    d_tumor = nearest_distance(points, points[tumor_high]) / med_nn
    caf_near = d_caf <= 2.0
    tumor_near = d_tumor <= 2.0
    interface = caf_near & tumor_near
    caf_only = caf_near & ~tumor_near
    tumor_only = tumor_near & ~caf_near
    other = ~(interface | caf_only | tumor_only)

    first = sample.iloc[0]
    row: dict = {
        "dataset_id": first["dataset_id"],
        "sample_id": first["sample_id"],
        "patient_id": first["patient_id"],
        "specimen_type": first["specimen_type"],
        "cohort_context": first["cohort_context"],
        "n_spots": int(len(sample)),
        "n_caf_core": int(caf_core.sum()),
        "n_tumor_high": int(tumor_high.sum()),
        "n_interface": int(interface.sum()),
        "n_caf_only_near": int(caf_only.sum()),
        "n_tumor_only_near": int(tumor_only.sum()),
        "n_other": int(other.sum()),
        "fraction_interface": float(interface.mean()),
        "median_neighbor_distance_px": med_nn,
    }
    long_rows: list[dict] = []
    for col, label in PROGRAMS:
        z_col = f"z_{col}"
        core_mean = float(sample.loc[caf_core, z_col].mean())
        noncore_mean = float(sample.loc[~caf_core, z_col].mean())
        interface_mean = float(sample.loc[interface, z_col].mean()) if interface.sum() >= 20 else np.nan
        noninterface_mean = float(sample.loc[~interface, z_col].mean()) if (~interface).sum() >= 20 else np.nan
        caf_only_mean = float(sample.loc[caf_only, z_col].mean()) if caf_only.sum() >= 20 else np.nan
        tumor_only_mean = float(sample.loc[tumor_only, z_col].mean()) if tumor_only.sum() >= 20 else np.nan
        other_mean = float(sample.loc[other, z_col].mean()) if other.sum() >= 20 else np.nan
        row[f"core_mean_z__{label}"] = core_mean
        row[f"core_enrichment__{label}"] = core_mean - noncore_mean
        row[f"interface_enrichment__{label}"] = interface_mean - noninterface_mean
        row[f"interface_vs_caf_only__{label}"] = interface_mean - caf_only_mean
        row[f"interface_vs_tumor_only__{label}"] = interface_mean - tumor_only_mean
        long_rows.append(
            {
                "dataset_id": first["dataset_id"],
                "sample_id": first["sample_id"],
                "patient_id": first["patient_id"],
                "specimen_type": first["specimen_type"],
                "cohort_context": first["cohort_context"],
                "program_label": label,
                "core_mean_z": core_mean,
                "core_enrichment": core_mean - noncore_mean,
                "interface_mean_z": interface_mean,
                "interface_enrichment": interface_mean - noninterface_mean,
                "interface_vs_caf_only": interface_mean - caf_only_mean,
                "interface_vs_tumor_only": interface_mean - tumor_only_mean,
            }
        )

    stromal = np.nanmean([row[f"core_enrichment__{label}"] for label in STROMAL_TUMOR_LABELS])
    immune = np.nanmean([row[f"core_enrichment__{label}"] for label in IMMUNE_LABELS])
    row["stromal_tumor_core_coupling"] = float(stromal)
    row["immune_core_coupling"] = float(immune)
    row["immune_decoupling_index"] = float(stromal - immune)
    row["interface_tumor_stroma_score"] = float(
        np.nanmean([row[f"interface_enrichment__{label}"] for label in ["tumor aggressive", "EMT/invasion", "TGF-beta"]])
    )
    row["interface_inflammatory_score"] = float(
        np.nanmean([row[f"interface_enrichment__{label}"] for label in ["IFN/MHC", "immune core"]])
    )
    return row, long_rows


def assign_ecotypes(sample_summary: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    feature_cols = [f"core_enrichment__{label}" for label in CORE_PROGRAMS]
    matrix = sample_summary[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    scaled = StandardScaler().fit_transform(matrix)
    kmeans = KMeans(n_clusters=4, random_state=RANDOM_SEED, n_init=100)
    sample_summary = sample_summary.copy()
    sample_summary["caf_core_ecotype_cluster"] = kmeans.fit_predict(scaled) + 1

    nonnegative = matrix.copy()
    for col in nonnegative.columns:
        nonnegative[col] = nonnegative[col] - nonnegative[col].min() + 0.01
    nmf = NMF(n_components=4, init="nndsvda", random_state=RANDOM_SEED, max_iter=2000)
    weights = nmf.fit_transform(nonnegative.to_numpy(float))
    loadings = nmf.components_
    nmf_weight = pd.DataFrame(
        weights,
        columns=[f"nmf_ecotype_{i + 1}_weight" for i in range(weights.shape[1])],
    )
    nmf_weight.insert(0, "sample_id", sample_summary["sample_id"].to_numpy())
    nmf_loading = pd.DataFrame(loadings, columns=CORE_PROGRAMS)
    nmf_loading.insert(0, "nmf_ecotype", [f"NMF{i + 1}" for i in range(loadings.shape[0])])
    nmf_labels: list[str] = []
    for _, row in nmf_loading.iterrows():
        top = row.drop(labels=["nmf_ecotype"]).sort_values(ascending=False).head(3).index.tolist()
        nmf_labels.append(" / ".join(top))
    nmf_loading["nmf_ecotype_label"] = nmf_labels
    dominant = weights.argmax(axis=1)
    sample_summary["dominant_nmf_ecotype"] = [f"NMF{i + 1}" for i in dominant]
    sample_summary["dominant_nmf_ecotype_label"] = [nmf_labels[i] for i in dominant]

    labels: dict[int, str] = {}
    for cluster_id, group in sample_summary.groupby("caf_core_ecotype_cluster"):
        med = group[feature_cols].median().rename(lambda c: c.replace("core_enrichment__", ""))
        top = med.sort_values(ascending=False).head(3).index.tolist()
        label = " / ".join(top)
        labels[int(cluster_id)] = label
    sample_summary["caf_core_ecotype_label"] = sample_summary["caf_core_ecotype_cluster"].map(labels)
    return sample_summary, nmf_loading, nmf_weight


def summarize_groups(sample_summary: pd.DataFrame, long_programs: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group_cols = ["cohort_context", "dominant_nmf_ecotype", "dominant_nmf_ecotype_label"]
    ecotype_counts = (
        sample_summary.groupby(group_cols, dropna=False)
        .agg(n_samples=("sample_id", "nunique"), median_decoupling=("immune_decoupling_index", "median"))
        .reset_index()
    )
    decoupling = (
        sample_summary.groupby("cohort_context", dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_stromal_tumor_core_coupling=("stromal_tumor_core_coupling", "median"),
            median_immune_core_coupling=("immune_core_coupling", "median"),
            median_immune_decoupling_index=("immune_decoupling_index", "median"),
            n_decoupled=("immune_decoupling_index", lambda s: int((s > 0).sum())),
        )
        .reset_index()
    )
    interface = (
        long_programs[long_programs["program_label"].isin(INTERFACE_LABELS)]
        .groupby(["cohort_context", "program_label"], dropna=False)
        .agg(
            n_samples=("sample_id", "nunique"),
            median_interface_enrichment=("interface_enrichment", "median"),
            n_positive_interface=("interface_enrichment", lambda s: int((s > 0).sum())),
            median_interface_vs_caf_only=("interface_vs_caf_only", "median"),
            median_interface_vs_tumor_only=("interface_vs_tumor_only", "median"),
        )
        .reset_index()
    )
    return ecotype_counts, decoupling, interface


def bh_adjust(pvalues: list[float]) -> list[float]:
    p = np.asarray(pvalues, dtype=float)
    order = np.argsort(p)
    ranked = p[order]
    q = np.empty_like(ranked)
    m = len(p)
    prev = 1.0
    for i in range(m - 1, -1, -1):
        val = ranked[i] * m / (i + 1)
        prev = min(prev, val)
        q[i] = prev
    out = np.empty_like(q)
    out[order] = q
    return out.tolist()


def pairwise_decoupling_tests(sample_summary: pd.DataFrame) -> pd.DataFrame:
    target_context = "lymph_node_metastasis"
    rows: list[dict] = []
    target = sample_summary.loc[sample_summary["cohort_context"].eq(target_context), "immune_decoupling_index"].dropna()
    for context, group in sample_summary.groupby("cohort_context"):
        if context == target_context:
            continue
        values = group["immune_decoupling_index"].dropna()
        if len(target) < 2 or len(values) < 2:
            continue
        test = mannwhitneyu(target, values, alternative="greater")
        rows.append(
            {
                "comparison": f"{target_context}_greater_than_{context}",
                "reference_context": target_context,
                "comparison_context": context,
                "n_reference": int(len(target)),
                "n_comparison": int(len(values)),
                "median_reference": float(target.median()),
                "median_comparison": float(values.median()),
                "median_difference": float(target.median() - values.median()),
                "test": "Mann-Whitney U, one-sided greater",
                "p_value": float(test.pvalue),
            }
        )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["q_value"] = bh_adjust(out["p_value"].tolist())
    return out


def make_figure(
    sample_summary: pd.DataFrame,
    nmf_loading: pd.DataFrame,
    decoupling_summary: pd.DataFrame,
    interface_summary: pd.DataFrame,
    output_base: Path,
) -> None:
    context_order = [
        "post_neoadjuvant_sections",
        "treatment_naive_primary",
        "primary_tumor",
        "liver_metastasis",
        "lymph_node_metastasis",
        "normal_pancreas",
    ]
    context_labels = {
        "post_neoadjuvant_sections": "post-NACT",
        "treatment_naive_primary": "treatment-naive",
        "primary_tumor": "primary",
        "liver_metastasis": "liver met",
        "lymph_node_metastasis": "LN met",
        "normal_pancreas": "normal",
    }

    fig = plt.figure(figsize=(13.5, 9.2), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, width_ratios=[1.16, 1.0], height_ratios=[1.0, 1.0])

    ax0 = fig.add_subplot(gs[0, 0])
    heat = nmf_loading.set_index("nmf_ecotype")[CORE_PROGRAMS].to_numpy(float)
    heat = heat / np.maximum(heat.max(axis=1, keepdims=True), 1e-9)
    im = ax0.imshow(heat, cmap="viridis", aspect="auto", vmin=0, vmax=1)
    ax0.set_xticks(np.arange(len(CORE_PROGRAMS)), CORE_PROGRAMS, rotation=45, ha="right", fontsize=8)
    ax0.set_yticks(np.arange(len(nmf_loading)), nmf_loading["nmf_ecotype"])
    ax0.set_title("CAF-core ecotype programs (NMF loadings)")
    fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.02)

    ax1 = fig.add_subplot(gs[0, 1])
    counts = (
        sample_summary.groupby(["cohort_context", "dominant_nmf_ecotype"], dropna=False)["sample_id"]
        .nunique()
        .reset_index(name="n_samples")
    )
    pivot = counts.pivot(index="cohort_context", columns="dominant_nmf_ecotype", values="n_samples").fillna(0)
    pivot = pivot.reindex([c for c in context_order if c in pivot.index])
    frac = pivot.div(pivot.sum(axis=1), axis=0)
    bottom = np.zeros(len(frac))
    colors = ["#4C78A8", "#F58518", "#54A24B", "#B279A2"]
    nmf_label_map = dict(zip(nmf_loading["nmf_ecotype"], nmf_loading["nmf_ecotype_label"]))
    for i, col in enumerate(sorted(frac.columns)):
        label = f"{col}: {nmf_label_map.get(col, '')}"
        ax1.bar(np.arange(len(frac)), frac[col], bottom=bottom, color=colors[i % len(colors)], label=label)
        bottom += frac[col].to_numpy(float)
    ax1.set_xticks(np.arange(len(frac)), [context_labels.get(c, c) for c in frac.index], rotation=35, ha="right")
    ax1.set_ylabel("fraction of samples")
    ax1.set_ylim(0, 1.02)
    ax1.set_title("CAF-core ecotype composition")
    ax1.legend(frameon=False, fontsize=7, ncol=1, loc="upper left", bbox_to_anchor=(0.0, -0.30))

    ax2 = fig.add_subplot(gs[1, 0])
    data = [
        sample_summary.loc[sample_summary["cohort_context"].eq(c), "immune_decoupling_index"].dropna().to_numpy(float)
        for c in context_order
        if c in set(sample_summary["cohort_context"])
    ]
    labels = [context_labels.get(c, c) for c in context_order if c in set(sample_summary["cohort_context"])]
    bp = ax2.boxplot(data, patch_artist=True, showfliers=False)
    for patch in bp["boxes"]:
        patch.set_facecolor("#D8DEE9")
        patch.set_edgecolor("#333333")
    for i, values in enumerate(data, start=1):
        jitter = np.linspace(-0.16, 0.16, len(values)) if len(values) > 1 else np.array([0.0])
        ax2.scatter(np.full(len(values), i) + jitter, values, s=18, color="#4C78A8", alpha=0.62, linewidths=0)
    ax2.axhline(0, color="#555555", linewidth=0.8)
    ax2.set_xticks(np.arange(1, len(labels) + 1), labels, rotation=35, ha="right")
    ax2.set_ylabel("stromal-tumor coupling minus immune coupling")
    ax2.set_title("Immune-decoupling index")

    ax3 = fig.add_subplot(gs[1, 1])
    sub = interface_summary[interface_summary["cohort_context"].isin(context_order)].copy()
    sub["program_label"] = pd.Categorical(sub["program_label"], categories=INTERFACE_LABELS, ordered=True)
    collapsed = (
        sub.groupby("program_label", observed=False)
        .agg(median_interface_enrichment=("median_interface_enrichment", "median"))
        .reset_index()
    )
    vals = collapsed["median_interface_enrichment"].to_numpy(float)
    ax3.barh(np.arange(len(collapsed)), vals, color=["#B279A2" if v < 0 else "#4C78A8" for v in vals])
    ax3.axvline(0, color="#555555", linewidth=0.8)
    ax3.set_yticks(np.arange(len(collapsed)), collapsed["program_label"])
    ax3.set_xlabel("median interface enrichment z")
    ax3.set_title("Tumor-stroma interface programs")

    fig.suptitle("CAF-core spatial ecotypes, immune decoupling and invasive interfaces", fontsize=14)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=220)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def write_report(
    sample_summary: pd.DataFrame,
    nmf_loading: pd.DataFrame,
    ecotype_counts: pd.DataFrame,
    decoupling_summary: pd.DataFrame,
    interface_summary: pd.DataFrame,
    decoupling_tests: pd.DataFrame,
) -> None:
    lines = [
        "# Stage 22 Spatial Ecotype Deep Dive",
        "",
        "## Purpose",
        "",
        "This analysis deepens the CAF-myeloid niche story by asking whether CAF cores separate into spatial ecotypes, whether stromal-tumor coupling decouples from immune coupling across tissue contexts, and whether aggressive programs concentrate at tumor-stroma interfaces.",
        "",
        "## CAF-Core Ecotypes",
        "",
    ]
    for _, row in nmf_loading.iterrows():
        top = row.drop(labels=["nmf_ecotype", "nmf_ecotype_label"]).sort_values(ascending=False).head(5)
        lines.append(f"- {row['nmf_ecotype']}: " + ", ".join([f"{idx} ({val:.2f})" for idx, val in top.items()]))
    lines.extend(["", "## Ecotype Composition by Context", ""])
    for _, row in ecotype_counts.sort_values(["cohort_context", "dominant_nmf_ecotype"]).iterrows():
        lines.append(
            f"- {row['cohort_context']} {row['dominant_nmf_ecotype']}: "
            f"{int(row['n_samples'])} samples; median decoupling {row['median_decoupling']:.3f}; "
            f"label {row['dominant_nmf_ecotype_label']}."
        )
    lines.extend(["", "## Immune-Decoupling Index", ""])
    for _, row in decoupling_summary.sort_values("median_immune_decoupling_index", ascending=False).iterrows():
        lines.append(
            f"- {row['cohort_context']}: median decoupling {row['median_immune_decoupling_index']:.3f}; "
            f"stromal-tumor coupling {row['median_stromal_tumor_core_coupling']:.3f}; "
            f"immune coupling {row['median_immune_core_coupling']:.3f}; "
            f"{int(row['n_decoupled'])}/{int(row['n_samples'])} samples > 0."
        )
    if not decoupling_tests.empty:
        lines.extend(["", "Pairwise one-sided tests comparing lymph-node metastasis decoupling with other contexts:", ""])
        for _, row in decoupling_tests.sort_values("q_value").iterrows():
            lines.append(
                f"- LN metastasis > {row['comparison_context']}: median delta {row['median_difference']:.3f}, "
                f"p={row['p_value']:.3g}, q={row['q_value']:.3g}."
            )
    lines.extend(["", "## Tumor-Stroma Interface Programs", ""])
    overall_interface = (
        interface_summary.groupby("program_label", observed=False)
        .agg(median_interface_enrichment=("median_interface_enrichment", "median"))
        .sort_values("median_interface_enrichment", ascending=False)
        .reset_index()
    )
    for _, row in overall_interface.iterrows():
        lines.append(f"- {row['program_label']}: median interface enrichment {row['median_interface_enrichment']:.3f}.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "This analysis supports a stronger conceptual framing: CAF-myeloid cores are not a single structure but a family of spatial ecotypes. The immune-decoupling index converts the lymph-node observation into a measurable axis, and the interface analysis nominates tumor-stroma contact zones as candidate sites where CAF-myeloid, TGF-beta/EMT and tumor-aggressive programs intersect. These remain spatial association analyses and should not be written as causal mechanism without perturbation or orthogonal validation.",
        ]
    )
    path = PROJECT_ROOT / "results/reports/spatial_ecotype_deep_dive_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    spots = load_spots()
    rows: list[dict] = []
    long_rows: list[dict] = []
    warnings: list[str] = []
    for sample_id, sample in spots.groupby("sample_id", sort=True):
        if len(sample) < 200:
            warnings.append(f"Skipped {sample_id}: fewer than 200 spots after filtering")
            continue
        try:
            row, long = summarize_sample(sample)
            rows.append(row)
            long_rows.extend(long)
            print(f"Summarized CAF-core ecotype features: {sample_id} ({len(sample)} spots)")
        except Exception as exc:
            warnings.append(f"{sample_id}: {exc}")

    sample_summary = pd.DataFrame(rows)
    long_programs = pd.DataFrame(long_rows)
    sample_summary, nmf_loading, nmf_weight = assign_ecotypes(sample_summary)
    ecotype_counts, decoupling_summary, interface_summary = summarize_groups(sample_summary, long_programs)
    decoupling_tests = pairwise_decoupling_tests(sample_summary)

    tables = PROJECT_ROOT / "results/tables"
    tables.mkdir(parents=True, exist_ok=True)
    sample_summary.to_csv(tables / "spatial_ecotype_sample_summary.csv", index=False)
    long_programs.to_csv(tables / "spatial_ecotype_program_long.csv", index=False)
    nmf_loading.to_csv(tables / "spatial_ecotype_nmf_loadings.csv", index=False)
    nmf_weight.to_csv(tables / "spatial_ecotype_nmf_sample_weights.csv", index=False)
    ecotype_counts.to_csv(tables / "spatial_ecotype_context_counts.csv", index=False)
    decoupling_summary.to_csv(tables / "immune_decoupling_context_summary.csv", index=False)
    decoupling_tests.to_csv(tables / "immune_decoupling_pairwise_tests.csv", index=False)
    interface_summary.to_csv(tables / "tumor_stroma_interface_program_summary.csv", index=False)

    source = PROJECT_ROOT / "results/source_data"
    source.mkdir(parents=True, exist_ok=True)
    nmf_loading.to_csv(source / "Source_Data_Fig_5A.csv", index=False)
    ecotype_counts.to_csv(source / "Source_Data_Fig_5B.csv", index=False)
    decoupling_summary.to_csv(source / "Source_Data_Fig_5C.csv", index=False)
    interface_summary.to_csv(source / "Source_Data_Fig_5D.csv", index=False)
    decoupling_tests.to_csv(source / "Source_Data_Extended_Immune_Decoupling_Tests.csv", index=False)

    make_figure(
        sample_summary,
        nmf_loading,
        decoupling_summary,
        interface_summary,
        PROJECT_ROOT / "results/figures/main/figure5_spatial_ecotype_deep_dive",
    )
    write_report(sample_summary, nmf_loading, ecotype_counts, decoupling_summary, interface_summary, decoupling_tests)

    write_status(
        STAGE,
        "success" if not warnings else "partial_success",
        {
            "n_samples_processed": int(sample_summary["sample_id"].nunique()),
            "n_program_rows": int(len(long_programs)),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": warnings,
            "outputs": [
                "results/tables/spatial_ecotype_sample_summary.csv",
                "results/tables/immune_decoupling_context_summary.csv",
                "results/tables/immune_decoupling_pairwise_tests.csv",
                "results/tables/tumor_stroma_interface_program_summary.csv",
                "results/figures/main/figure5_spatial_ecotype_deep_dive.pdf",
                "results/reports/spatial_ecotype_deep_dive_report.md",
            ],
            "next_manual_check": [
                "Inspect Figure 5 for whether ecotype labels are biologically interpretable.",
                "Decide whether interface analysis should be promoted to main text or kept as extended data.",
            ],
        },
    )
    print("Stage 22 spatial ecotype deep dive complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
