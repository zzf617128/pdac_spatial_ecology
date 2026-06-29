from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.stats import spearmanr
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(__file__).resolve().parents[1]
STAGE = "21_mvp_he_morphology_screen"
TARGETS = [
    ("score_caf_myeloid_barrier", "CAF-myeloid"),
    ("score_tumor_aggressive", "tumor aggressive"),
    ("z_ifn_antigen_presentation", "IFN/MHC"),
    ("score_immune_hub_core", "immune core"),
]
FEATURES = [
    ("tissue_fraction", "tissue fraction"),
    ("white_fraction", "white/background"),
    ("mean_gray", "brightness"),
    ("std_gray", "gray texture"),
    ("mean_saturation", "saturation"),
    ("mean_od_sum", "stain density"),
    ("mean_od_red", "red OD"),
    ("mean_od_green", "green OD"),
    ("mean_od_blue", "blue OD"),
    ("mean_red_minus_blue", "red-blue"),
    ("mean_blue_minus_red", "blue-red"),
    ("pink_fraction", "pink fraction"),
    ("purple_fraction", "purple fraction"),
    ("mean_edge", "edge density"),
]
RANDOM_SEED = 20260624


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


def read_scale(path: Path) -> tuple[float, float]:
    if not path.exists():
        return 1.0, 55.0
    payload = json.loads(path.read_text(encoding="utf-8"))
    return float(payload.get("tissue_hires_scalef", 1.0)), float(payload.get("spot_diameter_fullres", 55.0))


def integral_image(channel: np.ndarray) -> np.ndarray:
    return np.pad(channel.cumsum(axis=0).cumsum(axis=1), ((1, 0), (1, 0)), mode="constant")


def rect_mean(ii: np.ndarray, x0: np.ndarray, y0: np.ndarray, x1: np.ndarray, y1: np.ndarray, area: np.ndarray) -> np.ndarray:
    values = ii[y1, x1] - ii[y0, x1] - ii[y1, x0] + ii[y0, x0]
    return values / np.maximum(area, 1)


def image_channels(image_path: Path) -> dict[str, np.ndarray]:
    image = Image.open(image_path).convert("RGB")
    rgb = np.asarray(image, dtype=np.float32) / 255.0
    r = rgb[:, :, 0]
    g = rgb[:, :, 1]
    b = rgb[:, :, 2]
    gray = (r + g + b) / 3.0
    maxc = np.maximum.reduce([r, g, b])
    minc = np.minimum.reduce([r, g, b])
    saturation = (maxc - minc) / np.maximum(maxc, 1e-4)
    od_r = -np.log(np.clip(r, 1 / 255, 1.0))
    od_g = -np.log(np.clip(g, 1 / 255, 1.0))
    od_b = -np.log(np.clip(b, 1 / 255, 1.0))
    od_sum = od_r + od_g + od_b
    white = ((gray > 0.92) & (saturation < 0.12)).astype(np.float32)
    tissue = ((gray < 0.93) & (saturation > 0.03)).astype(np.float32)
    pink = ((r - b > 0.035) & (r - g > -0.02) & (saturation > 0.08) & (gray < 0.94)).astype(np.float32)
    purple = ((b - r > 0.005) & (saturation > 0.07) & (gray < 0.88)).astype(np.float32)
    gy, gx = np.gradient(gray)
    edge = np.sqrt(gx * gx + gy * gy)
    return {
        "tissue_fraction": tissue,
        "white_fraction": white,
        "mean_gray": gray,
        "gray_sq": gray * gray,
        "mean_saturation": saturation,
        "mean_od_sum": od_sum,
        "mean_od_red": od_r,
        "mean_od_green": od_g,
        "mean_od_blue": od_b,
        "mean_red_minus_blue": r - b,
        "mean_blue_minus_red": b - r,
        "pink_fraction": pink,
        "purple_fraction": purple,
        "mean_edge": edge,
    }


def extract_patch_features(sample_spots: pd.DataFrame, image_path: Path, scalefactors_path: Path) -> pd.DataFrame:
    scale, spot_diameter_fullres = read_scale(scalefactors_path)
    channels = image_channels(image_path)
    height, width = next(iter(channels.values())).shape
    patch_radius = max(10, int(round(0.50 * spot_diameter_fullres * scale)))
    x = np.rint(sample_spots["x_pixel"].to_numpy(float) * scale).astype(int)
    y = np.rint(sample_spots["y_pixel"].to_numpy(float) * scale).astype(int)
    x0 = np.clip(x - patch_radius, 0, width)
    x1 = np.clip(x + patch_radius + 1, 0, width)
    y0 = np.clip(y - patch_radius, 0, height)
    y1 = np.clip(y + patch_radius + 1, 0, height)
    area = (x1 - x0) * (y1 - y0)

    out = sample_spots[
        [
            "dataset_id",
            "sample_id",
            "patient_id",
            "barcode",
            "x_pixel",
            "y_pixel",
            "edge_or_background_risk",
            "local_white_background_fraction",
        ]
        + [target for target, _ in TARGETS]
    ].copy()
    out["he_patch_radius_hires_px"] = patch_radius
    gray_mean = None
    gray_sq_mean = None
    for key, channel in channels.items():
        mean_values = rect_mean(integral_image(channel), x0, y0, x1, y1, area)
        if key == "gray_sq":
            gray_sq_mean = mean_values
            continue
        out[key] = mean_values
        if key == "mean_gray":
            gray_mean = mean_values
    if gray_mean is not None and gray_sq_mean is not None:
        out["std_gray"] = np.sqrt(np.maximum(gray_sq_mean - gray_mean * gray_mean, 0.0))
    return out


def within_sample_zscore(df: pd.DataFrame, columns: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in columns:
        mean = out.groupby("sample_id")[col].transform("mean")
        std = out.groupby("sample_id")[col].transform("std").replace(0, np.nan)
        out[f"{col}__within_sample_z"] = ((out[col] - mean) / std).fillna(0.0)
    return out


def signed_spearman(x: pd.Series, y: pd.Series) -> float:
    mask = np.isfinite(x.to_numpy(float)) & np.isfinite(y.to_numpy(float))
    if mask.sum() < 50:
        return np.nan
    rho = spearmanr(x.to_numpy(float)[mask], y.to_numpy(float)[mask]).statistic
    return float(rho) if np.isfinite(rho) else np.nan


def summarize_feature_correlations(feature_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows: list[dict] = []
    eligible = feature_df[feature_df["analysis_eligible"]].copy()
    for sample_id, group in eligible.groupby("sample_id"):
        if len(group) < 200:
            continue
        first = group.iloc[0]
        for target_col, target_label in TARGETS:
            for feature_col, feature_label in FEATURES:
                rows.append(
                    {
                        "dataset_id": first["dataset_id"],
                        "sample_id": sample_id,
                        "patient_id": first["patient_id"],
                        "target": target_col,
                        "target_label": target_label,
                        "feature": feature_col,
                        "feature_label": feature_label,
                        "n_spots": int(len(group)),
                        "spearman_rho": signed_spearman(group[feature_col], group[target_col]),
                    }
                )
    per_sample = pd.DataFrame(rows)
    summary_rows: list[dict] = []
    for (target_col, target_label, feature_col, feature_label), group in per_sample.groupby(
        ["target", "target_label", "feature", "feature_label"], dropna=False
    ):
        rho = group["spearman_rho"].dropna()
        summary_rows.append(
            {
                "target": target_col,
                "target_label": target_label,
                "feature": feature_col,
                "feature_label": feature_label,
                "n_samples": int(rho.shape[0]),
                "median_rho": float(rho.median()) if not rho.empty else np.nan,
                "mean_abs_rho": float(rho.abs().mean()) if not rho.empty else np.nan,
                "n_positive": int((rho > 0).sum()),
                "n_negative": int((rho < 0).sum()),
                "fraction_same_sign_as_median": float(((np.sign(rho) == np.sign(rho.median())) & (rho != 0)).mean())
                if not rho.empty
                else np.nan,
            }
        )
    return per_sample, pd.DataFrame(summary_rows)


def evaluate_grouped_ridge(feature_df: pd.DataFrame, n_permutations: int = 20) -> pd.DataFrame:
    eligible = feature_df[feature_df["analysis_eligible"]].copy()
    feature_cols = [f"{feature}__within_sample_z" for feature, _ in FEATURES]
    groups = eligible["dataset_id"].astype(str) + ":" + eligible["patient_id"].astype(str)
    n_splits = min(5, groups.nunique())
    splitter = GroupKFold(n_splits=n_splits)
    rng = np.random.default_rng(RANDOM_SEED)
    rows: list[dict] = []
    x = eligible[feature_cols].to_numpy(float)
    for target_col, target_label in TARGETS:
        y_col = f"{target_col}__within_sample_z"
        y = eligible[y_col].to_numpy(float)
        fold_pred = np.zeros_like(y, dtype=float)
        for fold_idx, (train_idx, test_idx) in enumerate(splitter.split(x, y, groups=groups), start=1):
            model = make_pipeline(StandardScaler(), Ridge(alpha=5.0, random_state=RANDOM_SEED))
            model.fit(x[train_idx], y[train_idx])
            fold_pred[test_idx] = model.predict(x[test_idx])
            rows.append(
                {
                    "target": target_col,
                    "target_label": target_label,
                    "metric_scope": "observed_fold",
                    "fold": fold_idx,
                    "permutation": -1,
                    "n_train": int(len(train_idx)),
                    "n_test": int(len(test_idx)),
                    "spearman_rho": signed_spearman(pd.Series(y[test_idx]), pd.Series(fold_pred[test_idx])),
                    "r2": float(r2_score(y[test_idx], fold_pred[test_idx])),
                }
            )
        rows.append(
            {
                "target": target_col,
                "target_label": target_label,
                "metric_scope": "observed_all",
                "fold": 0,
                "permutation": -1,
                "n_train": int(len(y)),
                "n_test": int(len(y)),
                "spearman_rho": signed_spearman(pd.Series(y), pd.Series(fold_pred)),
                "r2": float(r2_score(y, fold_pred)),
            }
        )

        for perm in range(n_permutations):
            shuffled = eligible[["sample_id", y_col]].copy()
            shuffled[y_col] = shuffled.groupby("sample_id")[y_col].transform(
                lambda s: rng.permutation(s.to_numpy(float))
            )
            y_perm = shuffled[y_col].to_numpy(float)
            perm_pred = np.zeros_like(y_perm, dtype=float)
            for train_idx, test_idx in splitter.split(x, y_perm, groups=groups):
                model = make_pipeline(StandardScaler(), Ridge(alpha=5.0, random_state=RANDOM_SEED))
                model.fit(x[train_idx], y_perm[train_idx])
                perm_pred[test_idx] = model.predict(x[test_idx])
            rows.append(
                {
                    "target": target_col,
                    "target_label": target_label,
                    "metric_scope": "permuted_all",
                    "fold": 0,
                    "permutation": perm,
                    "n_train": int(len(y_perm)),
                    "n_test": int(len(y_perm)),
                    "spearman_rho": signed_spearman(pd.Series(y_perm), pd.Series(perm_pred)),
                    "r2": float(r2_score(y_perm, perm_pred)),
                }
            )
    return pd.DataFrame(rows)


def make_figure(corr_summary: pd.DataFrame, model_metrics: pd.DataFrame, output_base: Path) -> None:
    feature_labels = [label for _, label in FEATURES]
    target_labels = [label for _, label in TARGETS]
    matrix = np.full((len(feature_labels), len(target_labels)), np.nan)
    for i, (_, feature_label) in enumerate(FEATURES):
        for j, (_, target_label) in enumerate(TARGETS):
            row = corr_summary[
                (corr_summary["feature_label"] == feature_label) & (corr_summary["target_label"] == target_label)
            ]
            if not row.empty:
                matrix[i, j] = float(row.iloc[0]["median_rho"])

    observed = model_metrics[model_metrics["metric_scope"] == "observed_all"].copy()
    permuted = model_metrics[model_metrics["metric_scope"] == "permuted_all"].copy()
    perm_summary = (
        permuted.groupby(["target", "target_label"], as_index=False)["spearman_rho"]
        .agg(perm_median="median", perm_q05=lambda s: s.quantile(0.05), perm_q95=lambda s: s.quantile(0.95))
        .merge(observed[["target", "target_label", "spearman_rho", "r2"]], on=["target", "target_label"], how="left")
    )
    perm_summary["target_label"] = pd.Categorical(perm_summary["target_label"], categories=target_labels, ordered=True)
    perm_summary = perm_summary.sort_values("target_label")

    fig = plt.figure(figsize=(11.2, 7.2), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.25, 1.0])
    ax0 = fig.add_subplot(gs[0, 0])
    im = ax0.imshow(matrix, cmap="coolwarm", vmin=-0.18, vmax=0.18, aspect="auto")
    ax0.set_xticks(np.arange(len(target_labels)), target_labels, rotation=35, ha="right")
    ax0.set_yticks(np.arange(len(feature_labels)), feature_labels)
    ax0.set_title("Median within-sample Spearman rho")
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            ax0.text(j, i, f"{matrix[i, j]:.2f}", ha="center", va="center", fontsize=7.2)
    fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.02)

    ax1 = fig.add_subplot(gs[0, 1])
    x = np.arange(len(perm_summary))
    ax1.bar(x, perm_summary["spearman_rho"], color="#4C78A8", width=0.58, label="observed")
    ax1.errorbar(
        x,
        perm_summary["perm_median"],
        yerr=[
            perm_summary["perm_median"] - perm_summary["perm_q05"],
            perm_summary["perm_q95"] - perm_summary["perm_median"],
        ],
        fmt="o",
        color="#333333",
        capsize=3,
        label="within-sample target shuffle",
    )
    ax1.axhline(0, color="#555555", linewidth=0.8)
    ax1.set_xticks(x, perm_summary["target_label"], rotation=35, ha="right")
    ax1.set_ylabel("Held-out grouped CV Spearman rho")
    ax1.set_title("Multifeature H&E predictability")
    ax1.legend(frameon=False, fontsize=8)
    fig.suptitle("H&E patch morphology carries limited but testable spatial ecology signal", fontsize=13)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=220)
    fig.savefig(output_base.with_suffix(".pdf"))
    plt.close(fig)


def write_report(corr_summary: pd.DataFrame, model_metrics: pd.DataFrame, feature_df: pd.DataFrame) -> None:
    observed = model_metrics[model_metrics["metric_scope"] == "observed_all"].copy()
    permuted = model_metrics[model_metrics["metric_scope"] == "permuted_all"].copy()
    perm_summary = permuted.groupby(["target", "target_label"], as_index=False)["spearman_rho"].agg(
        perm_median="median", perm_q95=lambda s: s.quantile(0.95)
    )
    observed = observed.merge(perm_summary, on=["target", "target_label"], how="left")
    top_corr = (
        corr_summary.assign(abs_median_rho=lambda d: d["median_rho"].abs())
        .sort_values(["target_label", "abs_median_rho"], ascending=[True, False])
        .groupby("target_label")
        .head(3)
    )
    lines = [
        "# Stage 21 H&E morphology screen",
        "",
        "## Purpose",
        "",
        "This analysis asks whether simple H&E patch features around Visium spots carry spatial signal for the CAF-myeloid, tumor-aggressive, IFN/MHC, and immune-core programs in the MVP cohorts (GSE282302 and GSE274103). It is intended as a cautious pathology bridge, not as evidence that morphology alone mechanistically defines the transcriptomic niche.",
        "",
        "## Inputs and QC",
        "",
        f"- Eligible spot-level rows after edge/background and tissue-patch filters: {int(feature_df['analysis_eligible'].sum()):,}.",
        f"- Samples with extracted H&E features: {feature_df['sample_id'].nunique():,}.",
        "- Model evaluation used within-sample z-scored targets and features, with grouped cross-validation by dataset-patient group.",
        "- Negative control permuted each target within sample before the same grouped CV procedure.",
        "",
        "## Multifeature predictability",
        "",
    ]
    for _, row in observed.sort_values("target_label").iterrows():
        above_null = row["spearman_rho"] - row["perm_median"]
        lines.append(
            f"- {row['target_label']}: observed CV Spearman rho {row['spearman_rho']:.3f}, R2 {row['r2']:.3f}, "
            f"null median {row['perm_median']:.3f}, delta {above_null:.3f}."
        )
    lines.extend(["", "## Strongest univariate feature associations", ""])
    for _, row in top_corr.iterrows():
        lines.append(
            f"- {row['target_label']} vs {row['feature_label']}: median within-sample rho {row['median_rho']:.3f} "
            f"across {int(row['n_samples'])} samples."
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "A reproducible but modest signal would be useful as a figure or extended-data panel because it links the spatial transcriptomic niche to pathologist-visible tissue context. Claims should remain bounded: these color and texture features are proxies for tissue composition and staining patterns, and do not replace expert histopathology or image-model validation.",
        ]
    )
    path = PROJECT_ROOT / "results/reports/he_morphology_screen_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    manifest = pd.read_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv")
    manifest = manifest[manifest["dataset_id"].isin(["GSE282302", "GSE274103"])].copy()
    manifest = manifest[manifest["include_primary"].astype(str).str.lower().eq("true")].copy()
    spot = pd.read_csv(PROJECT_ROOT / "results/tables/mvp_spot_level_scores_with_edge_qc.csv")
    spot = spot[spot["dataset_id"].isin(["GSE282302", "GSE274103"])].copy()

    rows: list[pd.DataFrame] = []
    warnings: list[str] = []
    for _, manifest_row in manifest.iterrows():
        dataset_id = str(manifest_row["dataset_id"])
        sample_id = str(manifest_row["sample_id"])
        sample_spots = spot[(spot["dataset_id"] == dataset_id) & (spot["sample_id"] == sample_id)].copy()
        if sample_spots.empty:
            warnings.append(f"No scored spots for {dataset_id} {sample_id}")
            continue
        image_path = Path(str(manifest_row["image_path"]))
        scalefactors_path = Path(str(manifest_row["scalefactors_path"]))
        if not image_path.exists():
            warnings.append(f"Missing image for {dataset_id} {sample_id}: {image_path}")
            continue
        rows.append(extract_patch_features(sample_spots, image_path, scalefactors_path))
        print(f"Extracted H&E patch features: {dataset_id} {sample_id} ({len(sample_spots)} spots)")

    feature_df = pd.concat(rows, ignore_index=True)
    feature_df["edge_or_background_risk"] = feature_df["edge_or_background_risk"].astype(str).str.lower().eq("true")
    feature_df["analysis_eligible"] = (
        ~feature_df["edge_or_background_risk"]
        & (feature_df["local_white_background_fraction"].fillna(1.0) <= 0.20)
        & (feature_df["tissue_fraction"].fillna(0.0) >= 0.35)
    )
    feature_df = within_sample_zscore(feature_df, [col for col, _ in FEATURES] + [col for col, _ in TARGETS])

    tables_dir = PROJECT_ROOT / "results/tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    feature_columns = [
        "dataset_id",
        "sample_id",
        "patient_id",
        "barcode",
        "x_pixel",
        "y_pixel",
        "analysis_eligible",
        "he_patch_radius_hires_px",
    ] + [col for col, _ in TARGETS] + [col for col, _ in FEATURES]
    feature_df[feature_columns].to_csv(tables_dir / "mvp_he_patch_morphology_features.csv", index=False)

    per_sample_corr, corr_summary = summarize_feature_correlations(feature_df)
    model_metrics = evaluate_grouped_ridge(feature_df, n_permutations=20)
    per_sample_corr.to_csv(tables_dir / "mvp_he_patch_feature_correlations_per_sample.csv", index=False)
    corr_summary.to_csv(tables_dir / "mvp_he_patch_feature_correlation_summary.csv", index=False)
    model_metrics.to_csv(tables_dir / "mvp_he_patch_grouped_cv_metrics.csv", index=False)

    source_dir = PROJECT_ROOT / "results/source_data"
    source_dir.mkdir(parents=True, exist_ok=True)
    corr_summary.to_csv(source_dir / "Source_Data_Fig_4A.csv", index=False)
    model_metrics.to_csv(source_dir / "Source_Data_Fig_4B.csv", index=False)

    make_figure(corr_summary, model_metrics, PROJECT_ROOT / "results/figures/main/figure4_he_morphology_screen")
    write_report(corr_summary, model_metrics, feature_df)

    write_status(
        STAGE,
        "success" if not warnings else "partial_success",
        {
            "n_samples_processed": int(feature_df["sample_id"].nunique()),
            "n_spots_with_features": int(len(feature_df)),
            "n_spots_analysis_eligible": int(feature_df["analysis_eligible"].sum()),
            "n_errors": 0,
            "critical_errors": [],
            "noncritical_warnings": warnings,
            "outputs": [
                "results/tables/mvp_he_patch_morphology_features.csv",
                "results/tables/mvp_he_patch_feature_correlation_summary.csv",
                "results/tables/mvp_he_patch_grouped_cv_metrics.csv",
                "results/figures/main/figure4_he_morphology_screen.pdf",
                "results/reports/he_morphology_screen_report.md",
            ],
            "next_manual_check": [
                "Inspect whether Figure 4 should remain exploratory or move to Extended Data.",
                "Do not claim histology-only prediction without independent image-model validation.",
            ],
        },
    )
    print("Stage 21 H&E morphology screen complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
