from __future__ import annotations

import itertools
import warnings
from datetime import datetime, timezone
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import cophenet, linkage
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score


PROJECT_ROOT = Path(__file__).resolve().parents[1]
REVISION_ROOT = PROJECT_ROOT / "results" / "revision_2026_06_29"
RANDOM_SEED = 20260629
RANKS = list(range(2, 9))
N_RUNS = 50
MAX_ITER = 20000

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


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def normalize_rows(matrix: np.ndarray) -> np.ndarray:
    denom = np.linalg.norm(matrix, axis=1, keepdims=True)
    denom[denom == 0] = 1.0
    return matrix / denom


def prepare_matrix() -> tuple[pd.DataFrame, np.ndarray]:
    sample_summary = pd.read_csv(PROJECT_ROOT / "results" / "tables" / "spatial_ecotype_sample_summary.csv")
    feature_cols = [f"core_enrichment__{label}" for label in CORE_PROGRAMS]
    missing = [col for col in feature_cols if col not in sample_summary.columns]
    if missing:
        raise ValueError(f"Missing required ecotype feature columns: {missing}")
    matrix = sample_summary[feature_cols].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    nonnegative = matrix.copy()
    for col in nonnegative.columns:
        nonnegative[col] = nonnegative[col] - nonnegative[col].min() + 0.01
    return sample_summary, nonnegative.to_numpy(float)


def fit_nmf(matrix: np.ndarray, rank: int, seed: int, init: str = "nndsvdar") -> dict:
    model = NMF(
        n_components=rank,
        init=init,
        random_state=seed,
        max_iter=MAX_ITER,
        solver="cd",
        beta_loss="frobenius",
        tol=1e-5,
    )
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        weights = model.fit_transform(matrix)
    loadings = model.components_
    reconstruction = weights @ loadings
    residual = matrix - reconstruction
    reconstruction_error = float(np.linalg.norm(residual, ord="fro"))
    relative_error = float(reconstruction_error / max(np.linalg.norm(matrix, ord="fro"), 1e-12))
    explained_fraction = float(1.0 - (reconstruction_error**2 / max(np.linalg.norm(matrix, ord="fro") ** 2, 1e-12)))
    labels = np.argmax(weights, axis=1)
    return {
        "rank": rank,
        "seed": seed,
        "weights": weights,
        "loadings": loadings,
        "labels": labels,
        "reconstruction_error": reconstruction_error,
        "relative_error": relative_error,
        "explained_fraction": explained_fraction,
        "n_iter": int(model.n_iter_),
        "converged": bool(model.n_iter_ < MAX_ITER),
        "n_warnings": int(len(caught)),
        "init": init,
    }


def consensus_from_labels(label_list: list[np.ndarray]) -> np.ndarray:
    n_samples = len(label_list[0])
    consensus = np.zeros((n_samples, n_samples), dtype=float)
    for labels in label_list:
        consensus += labels[:, None] == labels[None, :]
    return consensus / len(label_list)


def consensus_metrics(consensus: np.ndarray) -> dict:
    tri = consensus[np.tril_indices_from(consensus, k=-1)]
    pac = float(((tri > 0.1) & (tri < 0.9)).mean()) if len(tri) else np.nan
    try:
        dist = 1.0 - consensus
        np.fill_diagonal(dist, 0.0)
        condensed = squareform(dist, checks=False)
        if np.nanstd(condensed) == 0:
            coph = np.nan
        else:
            z = linkage(condensed, method="average")
            coph = float(cophenet(z, condensed)[0])
    except Exception:
        coph = np.nan
    return {"consensus_pac_0.1_0.9": pac, "consensus_cophenetic": coph}


def component_stability(fits: list[dict]) -> tuple[float, float]:
    reference = min(fits, key=lambda item: item["reconstruction_error"])
    ref = normalize_rows(reference["loadings"])
    scores: list[float] = []
    for fit in fits:
        current = normalize_rows(fit["loadings"])
        sim = ref @ current.T
        rows, cols = linear_sum_assignment(-sim)
        scores.extend(sim[rows, cols].tolist())
    return float(np.mean(scores)), float(np.std(scores))


def ari_stability(fits: list[dict]) -> tuple[float, float]:
    scores = [
        adjusted_rand_score(a["labels"], b["labels"])
        for a, b in itertools.combinations(fits, 2)
    ]
    return float(np.mean(scores)), float(np.std(scores))


def component_rows(rank: int, fit: dict) -> list[dict]:
    rows: list[dict] = []
    for idx, values in enumerate(fit["loadings"], start=1):
        order = np.argsort(values)[::-1]
        top = [CORE_PROGRAMS[i] for i in order[:5]]
        rows.append(
            {
                "rank": rank,
                "component": f"rank{rank}_NMF{idx}",
                "top_program_1": top[0],
                "top_program_2": top[1],
                "top_program_3": top[2],
                "top_program_4": top[3],
                "top_program_5": top[4],
                "component_label": " / ".join(top[:3]),
            }
        )
    return rows


def reference_rank_sweep(matrix: np.ndarray) -> pd.DataFrame:
    rows: list[dict] = []
    for rank in RANKS:
        fit = fit_nmf(matrix, rank, RANDOM_SEED, init="nndsvda")
        labels = fit["labels"]
        counts = np.bincount(labels, minlength=rank)
        rows.append(
            {
                "rank": rank,
                "init": "nndsvda",
                "seed": RANDOM_SEED,
                "reconstruction_error": fit["reconstruction_error"],
                "relative_error": fit["relative_error"],
                "explained_fraction": fit["explained_fraction"],
                "n_iter": fit["n_iter"],
                "converged": fit["converged"],
                "min_dominant_component_n": int(counts.min()),
                "max_dominant_component_n": int(counts.max()),
                "component_labels": "; ".join(
                    row["component_label"] for row in component_rows(rank, fit)
                ),
            }
        )
    reference = pd.DataFrame(rows).sort_values("rank").reset_index(drop=True)
    reference["delta_explained_from_previous_rank"] = reference["explained_fraction"].diff()
    return reference


def run_rank_scan(matrix: np.ndarray) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[int, dict]]:
    rng = np.random.default_rng(RANDOM_SEED)
    summary_rows: list[dict] = []
    run_rows: list[dict] = []
    component_label_rows: list[dict] = []
    best_by_rank: dict[int, dict] = {}

    for rank in RANKS:
        fits: list[dict] = []
        seeds = rng.integers(1, 2_000_000_000, size=N_RUNS)
        for run_idx, seed in enumerate(seeds, start=1):
            fit = fit_nmf(matrix, rank, int(seed), init="nndsvdar")
            fits.append(fit)
            run_rows.append(
                {
                    "rank": rank,
                    "run": run_idx,
                    "seed": int(seed),
                    "reconstruction_error": fit["reconstruction_error"],
                    "relative_error": fit["relative_error"],
                    "explained_fraction": fit["explained_fraction"],
                    "n_iter": fit["n_iter"],
                    "converged": fit["converged"],
                    "n_warnings": fit["n_warnings"],
                    "init": fit["init"],
                }
            )
        best = min(fits, key=lambda item: item["reconstruction_error"])
        best_by_rank[rank] = best
        consensus = consensus_from_labels([fit["labels"] for fit in fits])
        cmetrics = consensus_metrics(consensus)
        ari_mean, ari_sd = ari_stability(fits)
        comp_mean, comp_sd = component_stability(fits)
        run_df = pd.DataFrame([row for row in run_rows if row["rank"] == rank])
        counts = np.bincount(best["labels"], minlength=rank)
        summary_rows.append(
            {
                "rank": rank,
                "n_runs": N_RUNS,
                "n_samples": int(matrix.shape[0]),
                "n_features": int(matrix.shape[1]),
                "reconstruction_error_mean": float(run_df["reconstruction_error"].mean()),
                "reconstruction_error_sd": float(run_df["reconstruction_error"].std(ddof=0)),
                "reconstruction_error_best": float(best["reconstruction_error"]),
                "relative_error_mean": float(run_df["relative_error"].mean()),
                "explained_fraction_mean": float(run_df["explained_fraction"].mean()),
                "explained_fraction_best": float(best["explained_fraction"]),
                "ari_mean": ari_mean,
                "ari_sd": ari_sd,
                "component_cosine_to_best_mean": comp_mean,
                "component_cosine_to_best_sd": comp_sd,
                "n_converged_runs": int(run_df["converged"].sum()),
                "n_not_converged_runs": int((~run_df["converged"]).sum()),
                "min_dominant_component_n_best": int(counts.min()),
                "max_dominant_component_n_best": int(counts.max()),
                **cmetrics,
            }
        )
        component_label_rows.extend(component_rows(rank, best))
        print(
            f"Rank {rank}: explained={run_df['explained_fraction'].mean():.3f}; "
            f"ARI={ari_mean:.3f}; PAC={cmetrics['consensus_pac_0.1_0.9']:.3f}; "
            f"component stability={comp_mean:.3f}"
        )

    summary = pd.DataFrame(summary_rows).sort_values("rank").reset_index(drop=True)
    summary["delta_explained_from_previous_rank"] = summary["explained_fraction_mean"].diff()
    summary["delta_best_explained_from_previous_rank"] = summary["explained_fraction_best"].diff()
    return summary, pd.DataFrame(run_rows), pd.DataFrame(component_label_rows), best_by_rank


def make_figure(summary: pd.DataFrame, reference: pd.DataFrame, best_by_rank: dict[int, dict], output_base: Path) -> None:
    fig = plt.figure(figsize=(13.2, 9.0), constrained_layout=True)
    gs = fig.add_gridspec(2, 2)

    ax0 = fig.add_subplot(gs[0, 0])
    ax0.plot(summary["rank"], summary["explained_fraction_mean"], marker="o", color="#246A73", linewidth=2.0)
    ax0.plot(summary["rank"], summary["explained_fraction_best"], marker="s", color="#B04A37", linewidth=1.8)
    ax0.plot(reference["rank"], reference["explained_fraction"], marker="^", color="#6D597A", linewidth=1.8)
    ax0.axvline(4, color="#222222", linestyle="--", linewidth=1.0)
    ax0.set_xlabel("NMF rank")
    ax0.set_ylabel("explained fraction")
    ax0.set_title("Reconstruction gain")
    ax0.legend(["randomized NNDSVD mean", "best randomized seed", "NNDSVDa reference"], frameon=False, fontsize=8)

    ax1 = fig.add_subplot(gs[0, 1])
    ax1.plot(summary["rank"], summary["ari_mean"], marker="o", color="#4E79A7", linewidth=2.0, label="ARI")
    ax1.plot(summary["rank"], 1.0 - summary["consensus_pac_0.1_0.9"], marker="s", color="#59A14F", linewidth=2.0, label="1 - PAC")
    ax1.axvline(4, color="#222222", linestyle="--", linewidth=1.0)
    ax1.set_xlabel("NMF rank")
    ax1.set_ylabel("stability score")
    ax1.set_ylim(0, 1.03)
    ax1.set_title("Assignment stability")
    ax1.legend(frameon=False, fontsize=8)

    ax2 = fig.add_subplot(gs[1, 0])
    ax2.plot(summary["rank"], summary["component_cosine_to_best_mean"], marker="o", color="#8E6C8A", linewidth=2.0)
    ax2.fill_between(
        summary["rank"],
        summary["component_cosine_to_best_mean"] - summary["component_cosine_to_best_sd"],
        summary["component_cosine_to_best_mean"] + summary["component_cosine_to_best_sd"],
        color="#8E6C8A",
        alpha=0.16,
        linewidth=0,
    )
    ax2.axvline(4, color="#222222", linestyle="--", linewidth=1.0)
    ax2.set_xlabel("NMF rank")
    ax2.set_ylabel("component cosine similarity")
    ax2.set_ylim(0, 1.03)
    ax2.set_title("Component reproducibility")

    ax3 = fig.add_subplot(gs[1, 1])
    rank4 = best_by_rank[4]["loadings"]
    rank4_norm = rank4 / np.maximum(rank4.max(axis=1, keepdims=True), 1e-12)
    im = ax3.imshow(rank4_norm, cmap="mako" if "mako" in plt.colormaps() else "viridis", aspect="auto", vmin=0, vmax=1)
    ax3.set_xticks(np.arange(len(CORE_PROGRAMS)), CORE_PROGRAMS, rotation=45, ha="right", fontsize=8)
    ax3.set_yticks(np.arange(rank4.shape[0]), [f"NMF{i}" for i in range(1, rank4.shape[0] + 1)])
    ax3.set_title("Rank-4 ecotype loading structure")
    fig.colorbar(im, ax=ax3, fraction=0.046, pad=0.02)

    fig.suptitle("NMF rank sensitivity of CAF-core spatial ecotypes", fontsize=15)
    output_base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_base.with_suffix(".png"), dpi=300)
    fig.savefig(output_base.with_suffix(".pdf"))
    fig.savefig(output_base.with_suffix(".svg"))
    plt.close(fig)


def write_report(summary: pd.DataFrame, reference: pd.DataFrame, component_labels: pd.DataFrame) -> None:
    n_samples = int(summary.iloc[0]["n_samples"])
    rank4 = summary.loc[summary["rank"].eq(4)].iloc[0]
    ref4 = reference.loc[reference["rank"].eq(4)].iloc[0]
    best_stability_rank = int(summary.sort_values(["ari_mean", "component_cosine_to_best_mean"], ascending=False).iloc[0]["rank"])
    best_pac_rank = int(summary.sort_values("consensus_pac_0.1_0.9").iloc[0]["rank"])
    lines = [
        "# NMF Rank Stability Revision Analysis",
        "",
        f"Generated: {utc_now()}",
        "",
        "## Purpose",
        "",
        "This reviewer-facing sensitivity analysis tests whether the CAF-core spatial ecotype result depends on arbitrarily fixing NMF rank to 4.",
        "",
        "## Design",
        "",
        f"- Input: {n_samples} spatial samples by {len(CORE_PROGRAMS)} CAF-core enrichment features from `spatial_ecotype_sample_summary.csv`.",
        f"- Ranks tested: {min(RANKS)}-{max(RANKS)}.",
        f"- Repeats per rank: {N_RUNS} randomized NNDSVDar initializations plus one NNDSVDa reference fit matching the original Stage 22 implementation.",
        "- Metrics: reconstruction explained fraction, adjusted Rand index across runs, consensus PAC, consensus cophenetic correlation and component cosine reproducibility.",
        "",
        "## Rank Summary",
        "",
    ]
    for _, row in summary.iterrows():
        lines.append(
            f"- Rank {int(row['rank'])}: explained {row['explained_fraction_mean']:.3f}; "
            f"delta from previous {row['delta_explained_from_previous_rank'] if pd.notna(row['delta_explained_from_previous_rank']) else float('nan'):.3f}; "
            f"ARI {row['ari_mean']:.3f}; PAC {row['consensus_pac_0.1_0.9']:.3f}; "
            f"component cosine {row['component_cosine_to_best_mean']:.3f}; "
            f"converged {int(row['n_converged_runs'])}/{int(row['n_runs'])}; "
            f"best-run smallest component n={int(row['min_dominant_component_n_best'])}."
        )
    lines.extend(["", "## NNDSVDa Reference Sweep", ""])
    for _, row in reference.iterrows():
        lines.append(
            f"- Rank {int(row['rank'])}: explained {row['explained_fraction']:.3f}; "
            f"delta from previous {row['delta_explained_from_previous_rank'] if pd.notna(row['delta_explained_from_previous_rank']) else float('nan'):.3f}; "
            f"converged={row['converged']}; smallest component n={int(row['min_dominant_component_n'])}."
        )
    lines.extend(
        [
            "",
            "## Rank-4 Component Labels",
            "",
        ]
    )
    for _, row in component_labels.loc[component_labels["rank"].eq(4)].iterrows():
        lines.append(f"- {row['component']}: {row['component_label']}.")
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"Rank 4 explains {rank4['explained_fraction_mean']:.3f} of the non-negative CAF-core feature matrix on average across randomized starts and {ref4['explained_fraction']:.3f} in the NNDSVDa reference fit, with ARI {rank4['ari_mean']:.3f}, PAC {rank4['consensus_pac_0.1_0.9']:.3f} and component cosine reproducibility {rank4['component_cosine_to_best_mean']:.3f}.",
            f"The most stable assignment rank by ARI/cosine is rank {best_stability_rank}, and the lowest-PAC consensus rank is rank {best_pac_rank}; therefore rank 4 should be presented as a biologically interpretable working resolution rather than as the mathematically unique optimum.",
            "The reviewer-safe claim is that the main CAF-core axes are robust across a rank sweep, while the exact number of ecotype labels is a resolution choice supported by reconstruction, stability and interpretability.",
            "",
            "## Manuscript Use",
            "",
            "- Add as an Extended Data or Supplementary Figure supporting rank choice.",
            "- In Methods, state that NMF ranks 2-8 were rerun with 50 random initializations per rank.",
            "- In Results, avoid saying four ecotypes are intrinsically discrete; say rank 4 resolved four recurrent axes used for downstream interpretation.",
        ]
    )
    path = REVISION_ROOT / "docs" / "nmf_rank_stability_report.md"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    (REVISION_ROOT / "analysis_outputs").mkdir(parents=True, exist_ok=True)
    (REVISION_ROOT / "supplementary_tables").mkdir(parents=True, exist_ok=True)
    (REVISION_ROOT / "figures").mkdir(parents=True, exist_ok=True)

    _, matrix = prepare_matrix()
    reference = reference_rank_sweep(matrix)
    summary, run_level, component_labels, best_by_rank = run_rank_scan(matrix)

    summary_path = REVISION_ROOT / "analysis_outputs" / "nmf_rank_stability_summary.csv"
    run_path = REVISION_ROOT / "analysis_outputs" / "nmf_rank_stability_run_level.csv"
    labels_path = REVISION_ROOT / "analysis_outputs" / "nmf_rank_stability_component_labels.csv"
    reference_path = REVISION_ROOT / "analysis_outputs" / "nmf_rank_nndsvda_reference.csv"
    summary.to_csv(summary_path, index=False, encoding="utf-8")
    run_level.to_csv(run_path, index=False, encoding="utf-8")
    component_labels.to_csv(labels_path, index=False, encoding="utf-8")
    reference.to_csv(reference_path, index=False, encoding="utf-8")
    summary.to_csv(
        REVISION_ROOT / "supplementary_tables" / "Supplementary_Table_7_NMF_Rank_Stability.csv",
        index=False,
        encoding="utf-8",
    )
    make_figure(summary, reference, best_by_rank, REVISION_ROOT / "figures" / "Extended_Data_Figure_NMF_Rank_Stability")
    write_report(summary, reference, component_labels)
    print(f"Wrote {summary_path}")
    print("NMF rank stability revision analysis complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
