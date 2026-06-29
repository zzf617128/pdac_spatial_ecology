from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]

REQUIRED_DIRS = [
    "config",
    "data/raw/GSE272362",
    "data/raw/GSE282302",
    "data/raw/GSE274103",
    "data/interim",
    "data/processed/anndata",
    "data/processed/anndata_qc",
    "data/processed/anndata_scored",
    "data/processed/patches",
    "data/processed/embeddings",
    "data/external",
    "metadata",
    "src/pdac_spatial_ecology",
    "scripts",
    "notebooks",
    "results/tables",
    "results/figures/mvp",
    "results/figures/qc",
    "results/models",
    "results/logs",
    "results/reports",
    "docs",
]

REQUIRED_FILES = [
    "README.md",
    "environment.yml",
    "pyproject.toml",
    "config/datasets.yaml",
    "config/analysis_params.yaml",
    "config/signatures.yaml",
    "docs/analysis_decision_log.md",
    "docs/mvp_decision_report.md",
]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main() -> int:
    created_dirs: list[str] = []
    missing_files: list[str] = []

    for rel_dir in REQUIRED_DIRS:
        path = PROJECT_ROOT / rel_dir
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
            created_dirs.append(rel_dir)

    for rel_file in REQUIRED_FILES:
        if not (PROJECT_ROOT / rel_file).exists():
            missing_files.append(rel_file)

    status = {
        "stage": "00",
        "status": "success" if not missing_files else "partial_success",
        "timestamp_utc": now_iso(),
        "project_initialized": True,
        "project_root": str(PROJECT_ROOT),
        "python_executable": sys.executable,
        "python_version": sys.version,
        "platform": platform.platform(),
        "n_datasets_processed": 0,
        "n_samples_processed": 0,
        "n_errors": len(missing_files),
        "critical_errors": [],
        "noncritical_warnings": [f"Missing required file: {item}" for item in missing_files],
        "created_dirs": created_dirs,
        "next_manual_check": [
            "Run scripts/01_download_geo.py --mvp --audit-only before downloading large files."
        ],
    }
    if missing_files:
        status["critical_errors"].append("Project scaffold is incomplete.")

    write_json(PROJECT_ROOT / "results/logs/pipeline_status.json", status)
    write_json(PROJECT_ROOT / "results/logs/stage_00_status.json", status)

    print(f"Initialized project at {PROJECT_ROOT}")
    print(f"Status: {status['status']}")
    if missing_files:
        print("Missing files:")
        for item in missing_files:
            print(f"  - {item}")
    return 0 if not missing_files else 1


if __name__ == "__main__":
    raise SystemExit(main())

