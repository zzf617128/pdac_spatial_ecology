from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MVP_ACCESSIONS = ["GSE272362", "GSE282302", "GSE274103"]


@dataclass
class SampleFiles:
    dataset_id: str
    accession: str
    sample_id: str
    gsm_id: str
    file_names: list[str] = field(default_factory=list)
    expression_file: str = ""
    raw_expression_file: str = ""
    coordinates_file: str = ""
    image_file: str = ""
    scalefactors_file: str = ""
    spatial_bundle_file: str = ""
    notes: list[str] = field(default_factory=list)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_filelist(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if not row.get("Name") or row.get("#Archive/File") != "File":
                continue
            rows.append(row)
    return rows


def strip_compression(name: str) -> str:
    for suffix in [".tar.gz", ".csv.gz", ".json.gz", ".jpg.gz", ".png.gz"]:
        if name.endswith(suffix):
            return name[: -len(suffix)] + suffix[:-3]
    return name


def sample_key_from_name(file_name: str) -> tuple[str, str]:
    match = re.match(r"^(GSM\d+)_(.+)$", file_name)
    if not match:
        return "", ""
    gsm_id, rest = match.groups()
    known_suffixes = [
        "_filtered_feature_bc_matrix.h5",
        "_raw_feature_bc_matrix.h5",
        "_spatial.tar.gz",
        "_spatial.tar",
        "_spatial_enrichment.csv.gz",
        "_spatial_enrichment.csv",
        "_tissue_positions_list.csv.gz",
        "_tissue_positions_list.csv",
        "_tissue_positions.csv.gz",
        "_tissue_positions.csv",
        "_scalefactors_json.json.gz",
        "_scalefactors_json.json",
        "_tissue_hires_image.png.gz",
        "_tissue_hires_image.png",
        "_tissue_lowres_image.png.gz",
        "_tissue_lowres_image.png",
        "_detected_tissue_image.jpg.gz",
        "_detected_tissue_image.jpg",
        "_aligned_fiducials.jpg.gz",
        "_aligned_fiducials.jpg",
    ]
    sample_part = rest
    for suffix in known_suffixes:
        if rest.endswith(suffix):
            sample_part = rest[: -len(suffix)]
            break
    return gsm_id, f"{gsm_id}_{sample_part}"


def choose_file(current: str, candidate: str, priority: list[str]) -> str:
    if not current:
        return candidate
    current_rank = next((i for i, token in enumerate(priority) if token in current), len(priority))
    candidate_rank = next((i for i, token in enumerate(priority) if token in candidate), len(priority))
    return candidate if candidate_rank < current_rank else current


def accession_metadata(accession: str) -> dict[str, str | bool]:
    if accession == "GSE272362":
        return {
            "cohort_role": "main_primary_metastasis",
            "treatment_context": "not_applicable_mixed_primary_metastatic",
            "platform": "spatial_transcriptomics_author_processed",
            "include_primary": True,
        }
    if accession == "GSE282302":
        return {
            "cohort_role": "main_residual_large_he_st",
            "treatment_context": "metadata_required_post_neoadjuvant_or_primary_untreated_ambiguity",
            "platform": "10x Visium",
            "include_primary": True,
        }
    if accession == "GSE274103":
        return {
            "cohort_role": "main_treatment_naive",
            "treatment_context": "treatment_naive",
            "platform": "10x Visium",
            "include_primary": True,
        }
    return {
        "cohort_role": "unclassified",
        "treatment_context": "unknown",
        "platform": "unknown",
        "include_primary": False,
    }


def build_samples(accession: str) -> dict[str, SampleFiles]:
    filelist_path = PROJECT_ROOT / "data/raw" / accession / "filelist.txt"
    if not filelist_path.exists():
        raise FileNotFoundError(f"Missing filelist: {filelist_path}")

    samples: dict[str, SampleFiles] = {}
    for row in read_filelist(filelist_path):
        name = row["Name"]
        gsm_id, sample_id = sample_key_from_name(name)
        if not gsm_id:
            continue
        sample = samples.setdefault(
            sample_id,
            SampleFiles(dataset_id=accession, accession=accession, sample_id=sample_id, gsm_id=gsm_id),
        )
        sample.file_names.append(name)

        if "filtered_feature_bc_matrix.h5" in name:
            sample.expression_file = choose_file(sample.expression_file, name, ["filtered_feature"])
        elif "raw_feature_bc_matrix.h5" in name:
            sample.raw_expression_file = name
            if not sample.expression_file:
                sample.expression_file = name
        elif "spatial_enrichment.csv" in name:
            sample.expression_file = choose_file(sample.expression_file, name, ["spatial_enrichment"])
            sample.notes.append("expression_file_is_spatial_enrichment_csv_not_raw_count_matrix")

        if "tissue_positions" in name:
            sample.coordinates_file = choose_file(
                sample.coordinates_file, name, ["tissue_positions.csv", "tissue_positions_list.csv"]
            )
        if "scalefactors_json" in name:
            sample.scalefactors_file = name
        if "spatial.tar" in name:
            sample.spatial_bundle_file = name
            if not sample.coordinates_file:
                sample.coordinates_file = name
            if not sample.image_file:
                sample.image_file = name
            if not sample.scalefactors_file:
                sample.scalefactors_file = name
        if any(token in name for token in ["tissue_hires_image", "tissue_lowres_image", "detected_tissue_image"]):
            sample.image_file = choose_file(
                sample.image_file,
                name,
                ["tissue_hires_image", "tissue_lowres_image", "detected_tissue_image"],
            )
    return samples


def expected_local_path(accession: str, file_name: str) -> str:
    if not file_name:
        return ""
    return str(PROJECT_ROOT / "data/raw" / accession / strip_compression(file_name))


def spatial_bundle_dir(accession: str, file_name: str) -> Path:
    if file_name.endswith("_spatial.tar.gz"):
        return PROJECT_ROOT / "data/raw" / accession / file_name[: -len(".tar.gz")]
    if file_name.endswith("_spatial.tar"):
        return PROJECT_ROOT / "data/raw" / accession / file_name[: -len(".tar")]
    return PROJECT_ROOT / "data/raw" / accession


def expected_coordinate_path(sample: SampleFiles) -> str:
    if sample.coordinates_file and not sample.coordinates_file.endswith(("_spatial.tar.gz", "_spatial.tar")):
        return expected_local_path(sample.accession, sample.coordinates_file)
    if sample.spatial_bundle_file:
        return str(spatial_bundle_dir(sample.accession, sample.spatial_bundle_file) / "tissue_positions.csv")
    return ""


def expected_image_path(sample: SampleFiles) -> str:
    if sample.image_file and not sample.image_file.endswith(("_spatial.tar.gz", "_spatial.tar")):
        return expected_local_path(sample.accession, sample.image_file)
    if sample.spatial_bundle_file:
        return str(spatial_bundle_dir(sample.accession, sample.spatial_bundle_file) / "tissue_hires_image.png")
    return ""


def expected_scalefactors_path(sample: SampleFiles) -> str:
    if sample.scalefactors_file and not sample.scalefactors_file.endswith(("_spatial.tar.gz", "_spatial.tar")):
        return expected_local_path(sample.accession, sample.scalefactors_file)
    if sample.spatial_bundle_file:
        return str(spatial_bundle_dir(sample.accession, sample.spatial_bundle_file) / "scalefactors_json.json")
    return ""


def write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


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


def write_manifest_report(rows: list[dict], errors: list[str]) -> None:
    report_path = PROJECT_ROOT / "docs/dataset_audit_report.md"
    by_dataset: dict[str, list[dict]] = {}
    for row in rows:
        by_dataset.setdefault(row["dataset_id"], []).append(row)

    lines = [
        "# Dataset Audit Report",
        "",
        f"Last updated UTC: {now_iso()}",
        "",
        "## MVP Manifest Summary",
        "",
        "| dataset | n_samples | expression_expected | coordinates_expected | image_expected | notes |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for dataset_id, dataset_rows in by_dataset.items():
        n_expression = sum(row["has_expression_expected"] == "true" for row in dataset_rows)
        n_coordinates = sum(row["has_coordinates_expected"] == "true" for row in dataset_rows)
        n_image = sum(row["has_image_expected"] == "true" for row in dataset_rows)
        dataset_notes = sorted({row["notes"] for row in dataset_rows if row["notes"]})
        note_text = "; ".join(dataset_notes[:3]) if dataset_notes else "ok"
        lines.append(
            f"| {dataset_id} | {len(dataset_rows)} | {n_expression} | {n_coordinates} | {n_image} | {note_text} |"
        )

    if errors:
        lines.extend(["", "## Errors", ""])
        lines.extend([f"- {error}" for error in errors])

    lines.extend(
        [
            "",
            "## Manual Checks Before Claims",
            "",
            "- Confirm GSE282302 treatment context at sample level.",
            "- Confirm GSE272362 `spatial_enrichment.csv` content after downloading/extracting.",
            "- Do not make H&E claims for samples lacking paired tissue images after extraction.",
            "- Do not make patient-level claims until patient IDs are curated from GEO/sample metadata.",
            "",
        ]
    )
    report_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an expected sample manifest from GEO filelists.")
    parser.add_argument("--accessions", nargs="+", default=MVP_ACCESSIONS)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    rows: list[dict] = []
    errors: list[str] = []

    for accession in args.accessions:
        metadata = accession_metadata(accession)
        try:
            samples = build_samples(accession)
        except Exception as exc:
            errors.append(f"{accession}: {exc}")
            continue

        for sample in sorted(samples.values(), key=lambda item: item.sample_id):
            notes = sorted(set(sample.notes))
            include_reason = "mvp_expected_complete" if (
                sample.expression_file and sample.coordinates_file and sample.image_file
            ) else ""
            exclude_reason = "" if include_reason else "missing_expected_expression_coordinates_or_image"
            rows.append(
                {
                    "dataset_id": sample.dataset_id,
                    "accession": sample.accession,
                    "cohort_role": metadata["cohort_role"],
                    "sample_id": sample.sample_id,
                    "gsm_id": sample.gsm_id,
                    "patient_id": f"patient_id_unknown_{sample.sample_id}",
                    "section_id": sample.sample_id,
                    "specimen_type": "metadata_required",
                    "disease": "PDAC",
                    "treatment_context": metadata["treatment_context"],
                    "platform": metadata["platform"],
                    "expression_path": expected_local_path(accession, sample.expression_file),
                    "raw_expression_path": expected_local_path(accession, sample.raw_expression_file),
                    "spatial_path": expected_local_path(accession, sample.spatial_bundle_file),
                    "image_path": expected_image_path(sample),
                    "coordinates_path": expected_coordinate_path(sample),
                    "scalefactors_path": expected_scalefactors_path(sample),
                    "metadata_source": "GEO_filelist",
                    "has_expression_expected": str(bool(sample.expression_file)).lower(),
                    "has_coordinates_expected": str(bool(sample.coordinates_file)).lower(),
                    "has_image_expected": str(bool(sample.image_file)).lower(),
                    "public_no_auth": "true",
                    "include_primary": str(bool(metadata["include_primary"])).lower(),
                    "include_reason": include_reason,
                    "exclude_reason": exclude_reason,
                    "notes": ";".join(notes),
                }
            )

    fieldnames = [
        "dataset_id",
        "accession",
        "cohort_role",
        "sample_id",
        "gsm_id",
        "patient_id",
        "section_id",
        "specimen_type",
        "disease",
        "treatment_context",
        "platform",
        "expression_path",
        "raw_expression_path",
        "spatial_path",
        "image_path",
        "coordinates_path",
        "scalefactors_path",
        "metadata_source",
        "has_expression_expected",
        "has_coordinates_expected",
        "has_image_expected",
        "public_no_auth",
        "include_primary",
        "include_reason",
        "exclude_reason",
        "notes",
    ]
    write_csv(PROJECT_ROOT / "metadata/dataset_manifest_raw.csv", rows, fieldnames)
    write_csv(PROJECT_ROOT / "metadata/dataset_manifest_curated.csv", rows, fieldnames)

    mapping_fields = ["dataset_id", "sample_id", "gsm_id", "patient_id", "section_id", "specimen_type"]
    write_csv(
        PROJECT_ROOT / "metadata/patient_sample_mapping.csv",
        [{key: row[key] for key in mapping_fields} for row in rows],
        mapping_fields,
    )
    write_manifest_report(rows, errors)

    status = "success" if not errors else "partial_success"
    write_status(
        "02",
        status,
        {
            "n_datasets_processed": len(args.accessions),
            "n_samples_processed": len(rows),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Inspect metadata/dataset_manifest_curated.csv.",
                "Resolve patient_id and specimen_type from GEO metadata before patient-level claims.",
                "Download and inspect one pilot sample from each dataset before full extraction.",
            ],
        },
    )
    print(f"Wrote manifest with {len(rows)} expected samples")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
