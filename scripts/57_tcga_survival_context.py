from __future__ import annotations

from pathlib import Path
from statistics import NormalDist

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
from scipy.stats import chi2
from statsmodels.duration.hazard_regression import PHReg


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external" / "TCGA_PAAD"
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR = ROOT / "results" / "figures" / "submission"
SOURCE_DIR = ROOT / "results" / "source_data"
REPORT_DIR = ROOT / "results" / "reports"

for directory in [TABLE_DIR, FIG_DIR, SOURCE_DIR, REPORT_DIR]:
    directory.mkdir(parents=True, exist_ok=True)

OUT_TABLE = TABLE_DIR / "tcga_paad_survival_context_summary.csv"
OUT_SOURCE = SOURCE_DIR / "Source_Data_Extended_Data_Fig_22_TCGA_survival_context.csv"
OUT_REPORT = REPORT_DIR / "gap2_clinical_outcome_resolution_plan.md"
OUT_FIG = FIG_DIR / "extended_data_figure22_tcga_survival_context"

AXES = [
    ("matrix_integrin", "Matrix-integrin"),
    ("SPP1_TAM", "SPP1/TAM"),
    ("TGFb_EMT", "TGF-beta/EMT"),
    ("myCAF_matrix", "myCAF/matrix"),
    ("bulk_decoupling_like_index", "Bulk decoupling-like"),
    ("DC_APC", "DC/APC"),
    ("T_NK", "T/NK"),
    ("B_plasma", "B/plasma"),
]


def load_survival_frame() -> pd.DataFrame:
    scores = pd.read_csv(TABLE_DIR / "tcga_paad_bulk_context_scores.csv")
    surv = pd.read_csv(DATA_DIR / "TCGA-PAAD.survival.tsv.gz", sep="\t")
    surv = surv.rename(columns={"_PATIENT": "patient_id"})
    keep = surv[["patient_id", "OS.time", "OS"]].dropna().copy()
    keep["OS.time"] = pd.to_numeric(keep["OS.time"], errors="coerce")
    keep["OS"] = pd.to_numeric(keep["OS"], errors="coerce")
    keep = keep.dropna(subset=["OS.time", "OS"])
    keep = keep[keep["OS.time"] > 0]
    merged = scores.merge(keep, on="patient_id", how="inner")
    return merged


def normal_2sided_p(z: float) -> float:
    return 2.0 * (1.0 - NormalDist().cdf(abs(z)))


def fit_univariable_cox(df: pd.DataFrame, axis: str) -> dict[str, float]:
    tmp = df[[axis, "OS.time", "OS"]].dropna().copy()
    tmp[axis] = (tmp[axis] - tmp[axis].mean()) / tmp[axis].std(ddof=0)
    model = PHReg(tmp["OS.time"], tmp[[axis]], status=tmp["OS"])
    fit = model.fit(disp=False)
    beta = float(fit.params[0])
    se = float(fit.bse[0])
    hr = float(np.exp(beta))
    lo = float(np.exp(beta - 1.96 * se))
    hi = float(np.exp(beta + 1.96 * se))
    p = normal_2sided_p(beta / se) if se > 0 else np.nan
    return {"cox_hr_per_sd": hr, "cox_ci95_low": lo, "cox_ci95_high": hi, "cox_p": p, "n": int(len(tmp)), "events": int(tmp["OS"].sum())}


def logrank_p(time: np.ndarray, event: np.ndarray, group_high: np.ndarray) -> float:
    order = np.argsort(time)
    time = time[order]
    event = event[order].astype(int)
    group_high = group_high[order].astype(bool)
    event_times = np.unique(time[event == 1])
    observed = 0.0
    expected = 0.0
    variance = 0.0
    for t in event_times:
        at_risk = time >= t
        n = int(at_risk.sum())
        n_high = int((at_risk & group_high).sum())
        d = int(((time == t) & (event == 1)).sum())
        d_high = int(((time == t) & (event == 1) & group_high).sum())
        if n <= 1:
            continue
        observed += d_high
        expected += d * (n_high / n)
        variance += (n_high * (n - n_high) * d * (n - d)) / (n * n * (n - 1))
    if variance <= 0:
        return np.nan
    stat = (observed - expected) ** 2 / variance
    return float(chi2.sf(stat, 1))


def km_curve(time: np.ndarray, event: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    order = np.argsort(time)
    time = time[order]
    event = event[order].astype(int)
    event_times = np.unique(time[event == 1])
    xs = [0.0]
    ys = [1.0]
    survival = 1.0
    for t in event_times:
        at_risk = np.sum(time >= t)
        deaths = np.sum((time == t) & (event == 1))
        if at_risk > 0:
            survival *= 1.0 - deaths / at_risk
        xs.extend([float(t), float(t)])
        ys.extend([ys[-1], survival])
    if len(time):
        xs.append(float(np.max(time)))
        ys.append(ys[-1])
    return np.array(xs), np.array(ys)


def build_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for axis, label in AXES:
        cox = fit_univariable_cox(df, axis)
        tmp = df[[axis, "OS.time", "OS"]].dropna().copy()
        median = float(tmp[axis].median())
        tmp["group"] = np.where(tmp[axis] >= median, "high", "low")
        p_lr = logrank_p(tmp["OS.time"].to_numpy(float), tmp["OS"].to_numpy(int), (tmp["group"] == "high").to_numpy(bool))
        rows.append(
            {
                "axis": axis,
                "axis_label": label,
                **cox,
                "median_split_cutpoint": median,
                "n_high": int((tmp["group"] == "high").sum()),
                "events_high": int(tmp.loc[tmp["group"] == "high", "OS"].sum()),
                "n_low": int((tmp["group"] == "low").sum()),
                "events_low": int(tmp.loc[tmp["group"] == "low", "OS"].sum()),
                "logrank_p_median_split": p_lr,
            }
        )
    out = pd.DataFrame(rows)
    out["cox_fdr_bh"] = bh_fdr(out["cox_p"].to_numpy(float))
    out["logrank_fdr_bh"] = bh_fdr(out["logrank_p_median_split"].to_numpy(float))
    return out.sort_values("cox_p").reset_index(drop=True)


def bh_fdr(pvalues: np.ndarray) -> np.ndarray:
    pvalues = np.asarray(pvalues, dtype=float)
    out = np.full_like(pvalues, np.nan)
    ok = ~np.isnan(pvalues)
    p = pvalues[ok]
    order = np.argsort(p)
    ranked = p[order]
    m = len(ranked)
    adjusted = np.empty(m)
    running = 1.0
    for i in range(m - 1, -1, -1):
        running = min(running, ranked[i] * m / (i + 1))
        adjusted[i] = running
    restored = np.empty(m)
    restored[order] = adjusted
    out[ok] = np.minimum(restored, 1.0)
    return out


def plot_km(ax: plt.Axes, df: pd.DataFrame, axis: str, label: str, summary: pd.DataFrame) -> None:
    tmp = df[[axis, "OS.time", "OS"]].dropna().copy()
    cut = float(tmp[axis].median())
    high = tmp[axis] >= cut
    colors = {"high": "#B23A48", "low": "#4C78A8"}
    for group, mask in [("high", high), ("low", ~high)]:
        xs, ys = km_curve(tmp.loc[mask, "OS.time"].to_numpy(float), tmp.loc[mask, "OS"].to_numpy(int))
        ax.step(xs / 30.44, ys, where="post", color=colors[group], lw=1.8, label=f"{label} {group}")
    p = float(summary.loc[summary["axis"] == axis, "logrank_p_median_split"].iloc[0])
    ax.set_ylim(0, 1.03)
    ax.set_xlim(left=0)
    ax.set_xlabel("Overall survival (months)", fontsize=8)
    ax.set_ylabel("Survival probability", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.legend(frameon=False, fontsize=7, loc="lower left")
    ax.text(0.98, 0.08, f"log-rank p={p:.2g}", transform=ax.transAxes, ha="right", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)


def plot(summary: pd.DataFrame, df: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(10.8, 7.2))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.1, 1.0], height_ratios=[1, 1], wspace=0.36, hspace=0.42)

    ax_a = fig.add_subplot(gs[:, 0])
    s = summary.sort_values("cox_hr_per_sd", ascending=True).reset_index(drop=True)
    y = np.arange(len(s))
    ax_a.axvline(1.0, color="#333333", lw=0.7)
    colors = np.where(s["cox_p"] < 0.05, "#B23A48", "#8A8F98")
    ax_a.hlines(y, s["cox_ci95_low"], s["cox_ci95_high"], color=colors, lw=1.6)
    ax_a.scatter(s["cox_hr_per_sd"], y, s=40, color=colors, edgecolor="white", linewidth=0.5, zorder=3)
    ax_a.set_xscale("log")
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(s["axis_label"], fontsize=8)
    ax_a.set_xlabel("Univariable Cox HR per 1 s.d. increase", fontsize=8)
    ax_a.set_title("A  TCGA PAAD bulk-axis survival context", loc="left", fontsize=10, fontweight="bold")
    ax_a.tick_params(labelsize=7)
    ax_a.spines[["top", "right"]].set_visible(False)
    for i, row in s.iterrows():
        ax_a.text(row["cox_ci95_high"] * 1.06, i, f"p={row['cox_p']:.2g}", va="center", fontsize=7)

    ax_b = fig.add_subplot(gs[0, 1])
    plot_km(ax_b, df, "matrix_integrin", "Matrix-integrin", summary)
    ax_b.set_title("B  Matrix-integrin median split", loc="left", fontsize=10, fontweight="bold")

    ax_c = fig.add_subplot(gs[1, 1])
    plot_km(ax_c, df, "bulk_decoupling_like_index", "Decoupling-like", summary)
    ax_c.set_title("C  Bulk decoupling-like median split", loc="left", fontsize=10, fontweight="bold")

    fig.suptitle("TCGA PAAD clinical context for nominated stromal-myeloid axes", fontsize=12, fontweight="bold", y=0.985)
    fig.text(
        0.01,
        0.012,
        "Bulk RNA-seq survival context is non-spatial and exploratory; it is not used as clinical validation of CAF-core localization.",
        fontsize=7,
        color="#444444",
    )
    for suffix in [".pdf", ".svg", ".png"]:
        fig.savefig(OUT_FIG.with_suffix(suffix), dpi=450 if suffix == ".png" else None, bbox_inches="tight")
    plt.close(fig)


def write_report(summary: pd.DataFrame) -> None:
    best = summary.sort_values("cox_p").iloc[0]
    n_nominal = int((summary["cox_p"] < 0.05).sum())
    n_fdr = int((summary["cox_fdr_bh"] < 0.10).sum())
    lines = [
        "# Gap 2: Clinical Outcome Resolution Plan",
        "",
        "Last updated: 2026-06-28",
        "",
        "## Reviewer-facing problem",
        "",
        "The spatial cohorts do not provide enough harmonized patient-level outcome metadata to support a clinical response, prognosis or survival claim.",
        "",
        "## Resolution strategy",
        "",
        "We used TCGA PAAD bulk RNA-seq as an exploratory non-spatial clinical-context layer for the nominated stromal-myeloid, matrix and immune axes. Univariable Cox models used standardized axis scores; median-split Kaplan-Meier/log-rank analyses were used only as a visual stress test.",
        "",
        "## Result",
        "",
        f"- Strongest nominal Cox association: **{best['axis_label']}** HR {best['cox_hr_per_sd']:.2f} per 1 s.d., p = {best['cox_p']:.2g}, FDR = {best['cox_fdr_bh']:.2g}.",
        f"- Nominal Cox p < 0.05 axes: {n_nominal}/{len(summary)}.",
        f"- FDR < 0.10 axes: {n_fdr}/{len(summary)}.",
        "",
        "## Claim after resolution",
        "",
        "This analysis does not turn the study into a clinical prognostic paper. It does, however, documents that clinical-outcome support was actively tested rather than merely avoided. Any survival language should remain exploratory and non-spatial unless an independent outcome-annotated spatial cohort is added.",
        "",
        "## Outputs",
        "",
        f"- `{OUT_TABLE.relative_to(ROOT)}`",
        f"- `{OUT_SOURCE.relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.pdf').relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.svg').relative_to(ROOT)}`",
        f"- `{OUT_FIG.with_suffix('.png').relative_to(ROOT)}`",
    ]
    OUT_REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    df = load_survival_frame()
    summary = build_summary(df)
    summary.to_csv(OUT_TABLE, index=False)
    summary.to_csv(OUT_SOURCE, index=False)
    plot(summary, df)
    write_report(summary)
    print(f"Wrote {OUT_TABLE}")
    print(f"Wrote {OUT_REPORT}")
    print(f"Wrote {OUT_FIG.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
