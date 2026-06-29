from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pandas as pd


PROJECT = Path(__file__).resolve().parents[1]
DATA = PROJECT / "data" / "external"
OUT_TABLE = PROJECT / "results" / "tables"
OUT_REPORT = PROJECT / "results" / "reports"
OUT_TABLE.mkdir(parents=True, exist_ok=True)
OUT_REPORT.mkdir(parents=True, exist_ok=True)


def parse_geo_series_matrix(path: Path) -> pd.DataFrame:
    sample_ids: list[str] | None = None
    rows: dict[str, list[str]] = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            parts = [p.strip().strip('"') for p in line.rstrip("\n").split("\t")]
            key = parts[0].replace("!Sample_", "")
            values = parts[1:]
            if key == "geo_accession":
                sample_ids = values
            rows.setdefault(key, []).append(values)

    if sample_ids is None:
        return pd.DataFrame()

    out = pd.DataFrame({"geo_accession": sample_ids})
    for key, entries in rows.items():
        if key == "geo_accession":
            continue
        if key == "characteristics_ch1":
            for values in entries:
                parsed = []
                for value in values:
                    if ": " in value:
                        parsed.append(value.split(": ", 1)[1])
                    else:
                        parsed.append(value)
                label = values[0].split(": ", 1)[0].replace(" ", "_").lower() if values else "characteristic"
                out[f"characteristic__{label}"] = parsed
        elif len(entries) == 1:
            out[key] = entries[0]
        else:
            for idx, values in enumerate(entries, start=1):
                out[f"{key}_{idx}"] = values
    return out


def summarize_gzip_csv(path: Path, nrows: int = 5) -> dict:
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        header = handle.readline().rstrip("\n").split(",")
    preview = pd.read_csv(path, compression="gzip", nrows=nrows)
    rows = 0
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        next(handle, None)
        for rows, _ in enumerate(handle, start=1):
            pass
    return {
        "file": str(path.relative_to(PROJECT)),
        "size_mb": round(path.stat().st_size / 1024 / 1024, 3),
        "n_columns": len(header),
        "n_rows": rows,
        "columns_first_30": header[:30],
        "preview_columns": list(preview.columns),
    }


def parse_visium_filelist(path: Path) -> pd.DataFrame:
    records = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            parts = re.split(r"\s+", line)
            filename = parts[-1]
            sample = filename.split("_")[0]
            records.append(
                {
                    "filename": filename,
                    "sample_prefix": sample,
                    "kind": infer_file_kind(filename),
                }
            )
    return pd.DataFrame(records)


def infer_file_kind(filename: str) -> str:
    lower = filename.lower()
    if lower.endswith(".h5"):
        return "expression_h5"
    if "tissue_positions" in lower:
        return "tissue_positions"
    if "scalefactors" in lower:
        return "scalefactors"
    if "filtered_feature_bc_matrix" in lower or "matrix.mtx" in lower:
        return "matrix"
    if lower.endswith((".png", ".jpg", ".jpeg", ".tif", ".tiff")):
        return "image"
    return "other"


def main() -> None:
    gse274 = DATA / "GSE274557"
    gse277 = DATA / "GSE277782"

    gse274_meta = parse_geo_series_matrix(gse274 / "GSE274557_series_matrix.txt.gz")
    gse277_meta = parse_geo_series_matrix(gse277 / "GSE277782_series_matrix.txt.gz")
    filelist = parse_visium_filelist(gse274 / "filelist.txt")

    gse274_meta.to_csv(OUT_TABLE / "gse274557_series_metadata.csv", index=False)
    gse277_meta.to_csv(OUT_TABLE / "gse277782_series_metadata.csv", index=False)
    filelist.to_csv(OUT_TABLE / "gse274557_raw_filelist_summary.csv", index=False)

    cosmx_files = [
        gse277 / "GSE277782_CosMx_SCT_data.csv.gz",
        gse277 / "GSE277782_Meta.data_CoxMx.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide1_exprMat_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide1_metadata_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide1_fov_positions_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide2_exprMat_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide2_metadata_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide2_fov_positions_file.csv.gz",
    ]
    cosmx_summary = [summarize_gzip_csv(path) for path in cosmx_files if path.exists()]
    pd.DataFrame(cosmx_summary).to_csv(OUT_TABLE / "gse277782_cosmx_file_feasibility.csv", index=False)

    # More detailed column audit for metadata tables.
    meta_audit = []
    for path in [
        gse277 / "GSE277782_Meta.data_CoxMx.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide1_metadata_file.csv.gz",
        gse277 / "GSE277782_Run5961_PDAC_Slide2_metadata_file.csv.gz",
    ]:
        df = pd.read_csv(path, compression="gzip", nrows=200)
        meta_audit.append(
            {
                "file": str(path.relative_to(PROJECT)),
                "columns": "; ".join(df.columns.astype(str).tolist()),
                "sample_values": json.dumps(
                    {col: df[col].dropna().astype(str).head(5).tolist() for col in df.columns[:20]},
                    ensure_ascii=False,
                ),
            }
        )
    pd.DataFrame(meta_audit).to_csv(OUT_TABLE / "gse277782_metadata_column_audit.csv", index=False)

    file_counts = (
        filelist.groupby(["sample_prefix", "kind"]).size().reset_index(name="n_files").sort_values(["sample_prefix", "kind"])
    )
    file_counts.to_csv(OUT_TABLE / "gse274557_raw_file_counts_by_sample.csv", index=False)

    gse274_site_cols = [c for c in gse274_meta.columns if "site" in c.lower() or "tissue" in c.lower() or "organ" in c.lower()]
    gse277_site_cols = [c for c in gse277_meta.columns if "site" in c.lower() or "tissue" in c.lower() or "organ" in c.lower()]

    report = []
    report.append("# Nature 2025 Metastatic PDAC Data Feasibility Audit\n")
    report.append("Last updated: 2026-06-27\n")
    report.append("## Downloaded lightweight files\n")
    report.append("- GSE274557 series matrix and RAW `filelist.txt`.\n")
    report.append("- GSE277782 series matrix, CosMx SCT data, cell metadata, slide-level expression matrices, FOV positions and cell metadata.\n")
    report.append("## GSE274557 Visium metadata\n")
    report.append(f"- GEO samples parsed: {len(gse274_meta)}.\n")
    report.append(f"- Candidate tissue/site columns: {', '.join(gse274_site_cols) if gse274_site_cols else 'none found in series matrix'}.\n")
    report.append(f"- RAW file entries: {len(filelist)}.\n")
    report.append("- RAW file types by sample are saved to `results/tables/gse274557_raw_file_counts_by_sample.csv`.\n")
    report.append("## GSE277782 CosMx metadata\n")
    report.append(f"- GEO samples parsed: {len(gse277_meta)}.\n")
    report.append(f"- Candidate tissue/site columns: {', '.join(gse277_site_cols) if gse277_site_cols else 'none found in series matrix'}.\n")
    report.append("- CosMx file dimensions are saved to `results/tables/gse277782_cosmx_file_feasibility.csv`.\n")
    report.append("- Metadata column audit is saved to `results/tables/gse277782_metadata_column_audit.csv`.\n")
    report.append("## Initial feasibility interpretation\n")
    report.append("- GSE274557 is the best candidate for broad external metastatic validation, but requires downloading or selectively extracting the 6.3 GB RAW tar.\n")
    report.append("- GSE277782 is immediately usable for high-resolution cell-level marker/module scoring if metadata include sample or slide/site labels.\n")
    report.append("- Next step: inspect GSE274557 RAW file names to identify whether Visium expression H5, tissue positions and images can be selectively downloaded/extracted, then run a minimal CAF-core scoring pilot on one or two samples.\n")
    (OUT_REPORT / "nature2025_metastatic_pdac_data_feasibility.md").write_text("\n".join(report), encoding="utf-8")

    print("Wrote feasibility outputs:")
    for path in [
        OUT_TABLE / "gse274557_series_metadata.csv",
        OUT_TABLE / "gse277782_series_metadata.csv",
        OUT_TABLE / "gse274557_raw_filelist_summary.csv",
        OUT_TABLE / "gse274557_raw_file_counts_by_sample.csv",
        OUT_TABLE / "gse277782_cosmx_file_feasibility.csv",
        OUT_TABLE / "gse277782_metadata_column_audit.csv",
        OUT_REPORT / "nature2025_metastatic_pdac_data_feasibility.md",
    ]:
        print(path)


if __name__ == "__main__":
    main()
