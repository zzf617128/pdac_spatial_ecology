from __future__ import annotations

import csv
import gzip
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parents[1]
ACCESSIONS = ["GSE282302", "GSE274103", "GSE272362"]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def normalize_key(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return value.strip("_")


def clean_value(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == '"' and value[-1] == '"':
        value = value[1:-1]
    return value


def parse_series_matrix(path: Path) -> pd.DataFrame:
    sample_ids: list[str] = []
    sample_lines: list[tuple[str, list[str]]] = []
    rows_by_sample: dict[str, dict] = {}
    repeated_counts: dict[str, int] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for parts in reader:
            if not parts or not parts[0].startswith("!Sample_"):
                continue
            key = parts[0].removeprefix("!Sample_")
            values = [clean_value(value) for value in parts[1:]]
            if key == "geo_accession":
                sample_ids = values
            sample_lines.append((key, values))

    if not sample_ids:
        return pd.DataFrame()
    rows_by_sample = {sample_id: {"gsm_id": sample_id} for sample_id in sample_ids}

    for key, values in sample_lines:
            if key == "geo_accession":
                continue
            if len(values) != len(sample_ids):
                continue
            if key == "characteristics_ch1":
                repeated_counts[key] = repeated_counts.get(key, 0) + 1
                for sample_id, value in zip(sample_ids, values):
                    row = rows_by_sample[sample_id]
                    row.setdefault("characteristics_ch1_all", [])
                    row["characteristics_ch1_all"].append(value)
                    if ":" in value:
                        char_key, char_value = value.split(":", 1)
                        row[f"char_{normalize_key(char_key)}"] = char_value.strip()
                    else:
                        row[f"characteristics_ch1_{repeated_counts[key]}"] = value
                continue

            repeated_counts[key] = repeated_counts.get(key, 0) + 1
            out_key = normalize_key(key)
            if repeated_counts[key] > 1:
                out_key = f"{out_key}_{repeated_counts[key]}"
            for sample_id, value in zip(sample_ids, values):
                rows_by_sample[sample_id][out_key] = value

    rows = []
    for row in rows_by_sample.values():
        if "characteristics_ch1_all" in row:
            row["characteristics_ch1_all"] = " | ".join(row["characteristics_ch1_all"])
        rows.append(row)
    return pd.DataFrame(rows)


def infer_patient_section(sample_id: str, title: str, dataset_id: str) -> dict[str, str]:
    text = f"{sample_id} {title}"
    out = {
        "patient_id_inferred": "",
        "section_id_inferred": "",
        "roi_id_inferred": "",
        "metadata_inference_notes": "",
    }
    if dataset_id == "GSE282302":
        # Titles/files look like C1_D10_ROI1_s1. Treat Dxx as patient/case until better metadata is found.
        match = re.search(r"(C\d+)_(D\d+)_(ROI\d+)(?:_(s\d+))?", text)
        if match:
            cohort, donor, roi, section = match.groups()
            out["patient_id_inferred"] = f"{cohort}_{donor}"
            out["section_id_inferred"] = section or "section_unknown"
            out["roi_id_inferred"] = roi
            out["metadata_inference_notes"] = "patient_id_inferred_from_Cx_Dx_pattern_not_clinically_verified"
        return out
    if dataset_id == "GSE274103":
        match = re.search(r"PDAC-p(\d+)", text)
        if match:
            out["patient_id_inferred"] = f"PDAC_p{match.group(1)}"
            out["section_id_inferred"] = "single_section"
            out["roi_id_inferred"] = "roi_unknown"
            out["metadata_inference_notes"] = "patient_id_inferred_from_PDAC-p_sample_name"
        return out
    if dataset_id == "GSE272362":
        subject_match = re.search(r"_S(\d+)_", text, flags=re.I)
        match = re.search(r"Sample[_ ]?(\d+)", text, flags=re.I)
        if subject_match or match:
            out["patient_id_inferred"] = f"S{subject_match.group(1)}" if subject_match else f"patient_unknown_sample_{match.group(1)}"
            out["section_id_inferred"] = f"Sample_{match.group(1)}"
            out["roi_id_inferred"] = "roi_unknown"
            out["metadata_inference_notes"] = "GSE272362_patient_id_inferred_from_Sxx_filename_pattern_not_verified"
        return out
    return out


def normalize_specimen(row: pd.Series, dataset_id: str) -> str:
    candidates = [
        str(row.get("char_tissue_specimen_type", "")),
        str(row.get("source_name_ch1", "")),
        str(row.get("char_tissue", "")),
        str(row.get("title", "")),
    ]
    text = " ".join(candidates).lower()
    if "normal pancreas" in text:
        return "normal_pancreas"
    if "primary tumor" in text or "primary" in text:
        return "primary_tumor"
    if "hepatic" in text or "liver" in text:
        return "liver_metastasis"
    if "lymph" in text:
        return "lymph_node_metastasis"
    if dataset_id == "GSE274103":
        return "treatment_naive_pdac"
    if dataset_id == "GSE282302":
        return "pdac_residual_or_primary_metadata_required"
    return "metadata_required"


def treatment_context(row: pd.Series, dataset_id: str) -> str:
    text = " ".join(str(row.get(col, "")) for col in row.index).lower()
    if dataset_id == "GSE274103":
        return "treatment_naive"
    if dataset_id == "GSE282302":
        if "folfirinox" in text or "neoadjuvant" in text:
            return "post_neoadjuvant_or_folfirinox_context_from_series"
        return "metadata_required_gse282302_treatment_context_ambiguous"
    return "not_applicable_or_metadata_required"


def update_manifest(sample_meta: pd.DataFrame) -> None:
    manifest_path = PROJECT_ROOT / "metadata/dataset_manifest_curated.csv"
    manifest = pd.read_csv(manifest_path)
    meta_lookup = sample_meta.set_index(["dataset_id", "gsm_id"])

    for idx, row in manifest.iterrows():
        key = (row["dataset_id"], row["gsm_id"])
        if key not in meta_lookup.index:
            continue
        meta = meta_lookup.loc[key]
        local_inference = infer_patient_section(str(row["sample_id"]), str(row["sample_id"]), str(row["dataset_id"]))
        patient_id = meta.get("patient_id_inferred", "") or local_inference.get("patient_id_inferred", "") or row["patient_id"]
        section_id = meta.get("section_id_inferred", "") or local_inference.get("section_id_inferred", "") or row["section_id"]
        manifest.loc[idx, "patient_id"] = patient_id
        manifest.loc[idx, "section_id"] = section_id
        manifest.loc[idx, "specimen_type"] = meta.get("specimen_type_curated", row["specimen_type"]) or row["specimen_type"]
        manifest.loc[idx, "treatment_context"] = meta.get("treatment_context_curated", row["treatment_context"]) or row["treatment_context"]
        note = str(row.get("notes", ""))
        infer_note = str(meta.get("metadata_inference_notes", "") or local_inference.get("metadata_inference_notes", ""))
        if infer_note and infer_note not in note:
            manifest.loc[idx, "notes"] = (note + ";" + infer_note).strip(";")

    manifest.to_csv(manifest_path, index=False)

    mapping_cols = ["dataset_id", "sample_id", "gsm_id", "patient_id", "section_id", "specimen_type"]
    manifest[mapping_cols].to_csv(PROJECT_ROOT / "metadata/patient_sample_mapping.csv", index=False)


def write_status(stage: str, status: str, payload: dict) -> None:
    base = {
        "stage": stage,
        "status": status,
        "timestamp_utc": now_iso(),
        "n_datasets_processed": 0,
        "n_samples_processed": 0,
        "n_errors": 0,
        "critical_errors": [],
        "noncritical_warnings": [],
        "next_manual_check": [],
    }
    base.update(payload)
    path = PROJECT_ROOT / f"results/logs/stage_{stage}_status.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(base, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    frames = []
    errors: list[str] = []
    for accession in ACCESSIONS:
        path = PROJECT_ROOT / "metadata/geo" / f"{accession}_series_matrix.txt.gz"
        try:
            frame = parse_series_matrix(path)
            frame.insert(0, "dataset_id", accession)
            inferred = [
                infer_patient_section(str(row["gsm_id"]), str(row.get("title", "")), accession)
                for _, row in frame.iterrows()
            ]
            infer_df = pd.DataFrame(inferred)
            frame = pd.concat([frame.reset_index(drop=True), infer_df], axis=1)
            frame["specimen_type_curated"] = [normalize_specimen(row, accession) for _, row in frame.iterrows()]
            frame["treatment_context_curated"] = [treatment_context(row, accession) for _, row in frame.iterrows()]
            frames.append(frame)
        except Exception as exc:
            errors.append(f"{accession}: {exc}")

    if not frames:
        raise RuntimeError("No GEO metadata parsed.")
    sample_meta = pd.concat(frames, ignore_index=True, sort=False)
    sample_meta.to_csv(PROJECT_ROOT / "metadata/geo_sample_metadata.csv", index=False)
    update_manifest(sample_meta)

    summary = (
        sample_meta.groupby(["dataset_id", "specimen_type_curated"], dropna=False)
        .size()
        .reset_index(name="n_samples")
    )
    summary.to_csv(PROJECT_ROOT / "metadata/geo_metadata_summary.csv", index=False)

    write_status(
        "06_geo_metadata",
        "success" if not errors else "partial_success",
        {
            "n_datasets_processed": int(sample_meta["dataset_id"].nunique()),
            "n_samples_processed": len(sample_meta),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors
            + [
                "GSE282302 patient IDs are inferred from Cx_Dx pattern and require manual verification.",
                "GSE272362 patient matching is not resolved from series matrix alone.",
            ],
            "next_manual_check": [
                "Inspect metadata/geo_sample_metadata.csv.",
                "Verify GSE282302 Cx_Dx patient inference against article/supplement if available.",
                "Do not make patient-level GSE282302 claims until inference is confirmed.",
            ],
        },
    )
    print(f"Parsed GEO metadata for {len(sample_meta)} samples")
    print("Wrote metadata/geo_sample_metadata.csv")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
