from __future__ import annotations

import gzip
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
REF_DIR = PROJECT / "data" / "external" / "GSE111672"
TABLE_DIR = PROJECT / "results" / "tables"
REPORT_DIR = PROJECT / "results" / "reports"

FILES = [
    ("PDAC-A", "GSE111672_PDAC-A-indrop-filtered-expMat.txt.gz"),
    ("PDAC-B", "GSE111672_PDAC-B-indrop-filtered-expMat.txt.gz"),
]

MARKER_SETS = {
    "myCAF_matrix": ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "ACTA2", "TAGLN"],
    "SPP1_TAM": ["SPP1", "APOE", "TREM2", "CD68", "CD163", "C1QA", "C1QB", "LST1"],
    "DC_APC": ["HLA-DRA", "HLA-DRB1", "CD74", "LAMP3", "CCR7", "CLEC10A", "FCER1A"],
    "T_NK": ["CD3D", "CD3E", "CD4", "CD8A", "IL7R", "NKG7", "GZMB"],
    "epithelial_tumor": ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT17", "MSLN", "CEACAM6"],
    "endothelial": ["PECAM1", "VWF", "EMCN", "KDR"],
    "normal_acinar": ["PRSS1", "CPA1", "REG1A"],
}


def read_header(path: Path) -> list[str]:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\n").split("\t")[1:]


def normalize_counts(matrix: pd.DataFrame) -> pd.DataFrame:
    counts = matrix.astype(float)
    library = counts.sum(axis=0).replace(0, np.nan)
    cpm = counts.div(library, axis=1) * 1_000_000
    return np.log1p(cpm.fillna(0))


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    REPORT_DIR.mkdir(parents=True, exist_ok=True)

    label_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []

    for sample, filename in FILES:
        path = REF_DIR / filename
        labels = read_header(path)
        counts_by_label = Counter(labels)
        for label, n_cells in sorted(counts_by_label.items()):
            label_rows.append(
                {
                    "reference_dataset": "GSE111672",
                    "sample": sample,
                    "cell_label": label,
                    "n_cells": n_cells,
                }
            )

        matrix = pd.read_csv(path, sep="\t", compression="gzip", index_col=0)
        matrix.index = matrix.index.astype(str).str.upper()
        log_cpm = normalize_counts(matrix)

        label_to_cols: dict[str, list[str]] = defaultdict(list)
        for col, label in zip(log_cpm.columns, labels):
            label_to_cols[label].append(col)

        for label, cols in sorted(label_to_cols.items()):
            label_expr = log_cpm[cols]
            for marker_set, genes in MARKER_SETS.items():
                present = [gene for gene in genes if gene in label_expr.index]
                if not present:
                    score = np.nan
                else:
                    score = float(label_expr.loc[present].mean(axis=0).mean())
                marker_rows.append(
                    {
                        "reference_dataset": "GSE111672",
                        "sample": sample,
                        "cell_label": label,
                        "marker_set": marker_set,
                        "n_cells": len(cols),
                        "n_genes_present": len(present),
                        "mean_log1p_cpm": score,
                    }
                )

    label_df = pd.DataFrame(label_rows)
    marker_df = pd.DataFrame(marker_rows)

    label_path = TABLE_DIR / "gse111672_reference_cell_label_summary.csv"
    marker_path = TABLE_DIR / "gse111672_reference_marker_means.csv"
    report_path = REPORT_DIR / "gse111672_reference_download_summary.md"

    label_df.to_csv(label_path, index=False)
    marker_df.to_csv(marker_path, index=False)

    total_cells = int(label_df["n_cells"].sum())
    n_labels = int(label_df["cell_label"].nunique())
    top_labels = (
        label_df.groupby("cell_label", as_index=False)["n_cells"].sum()
        .sort_values("n_cells", ascending=False)
        .head(12)
    )

    lines = [
        "# GSE111672 Reference Download Summary",
        "",
        "Last updated: 2026-06-25",
        "",
        "## Purpose",
        "",
        "This small PDAC scRNA/ST reference was downloaded as a low-cost first candidate for formal cell-state attribution or deconvolution prototyping.",
        "",
        "## Downloaded Files",
        "",
    ]
    for _, filename in FILES:
        lines.append(f"- `data/external/GSE111672/{filename}`")
    lines.extend(
        [
            "- `data/external/GSE111672/GSE111672_README.rtf`",
            "- `data/external/GSE111672/filelist.txt`",
            "- `data/external/GSE111672/GSE111672_1000L3_barcodes.txt.gz`",
            "- `data/external/GSE111672/GSE111672_LP_slide.batch_10001.version2.txt.gz`",
            "",
            "## Reference Label Summary",
            "",
            f"- Total filtered inDrop cells in downloaded PDAC-A/PDAC-B matrices: {total_cells}.",
            f"- Unique header-derived cell labels: {n_labels}.",
            "",
            "| cell label | n cells |",
            "|---|---:|",
        ]
    )
    for row in top_labels.itertuples(index=False):
        lines.append(f"| {row.cell_label} | {int(row.n_cells)} |")
    lines.extend(
        [
            "",
            "## Generated Outputs",
            "",
            f"- `{label_path.relative_to(PROJECT).as_posix()}`",
            f"- `{marker_path.relative_to(PROJECT).as_posix()}`",
            "",
            "## Interpretation Boundary",
            "",
            "The downloaded processed matrices provide a practical reference prototype with header-derived cell labels. They are not yet harmonized into a full deconvolution model for the current Visium cohorts.",
        ]
    )
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {label_path}")
    print(f"Wrote {marker_path}")
    print(f"Wrote {report_path}")


if __name__ == "__main__":
    main()
