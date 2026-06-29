from __future__ import annotations

import argparse
import gzip
import json
import shutil
import sys
import tarfile
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MVP_ACCESSIONS = ["GSE272362", "GSE274103"]


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


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


def decompress_gzip(path: Path, overwrite: bool = False) -> Path:
    if path.suffix != ".gz":
        raise ValueError(f"Not a .gz file: {path}")
    output = path.with_suffix("")
    if output.exists() and not overwrite:
        return output
    with gzip.open(path, "rb") as src, output.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return output


def extract_spatial_bundle(path: Path, overwrite: bool = False) -> Path:
    if not path.name.endswith("_spatial.tar.gz"):
        raise ValueError(f"Not a spatial bundle: {path}")
    output_dir = path.parent / path.name[: -len(".tar.gz")]
    output_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            target = output_dir / Path(member.name).name
            if target.exists() and not overwrite:
                continue
            src = archive.extractfile(member)
            if src is None:
                continue
            with target.open("wb") as dst:
                shutil.copyfileobj(src, dst)
    return output_dir


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract downloaded GEO archives for MVP datasets.")
    parser.add_argument("--accessions", nargs="+", default=MVP_ACCESSIONS)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    errors: list[str] = []
    n_gzip = 0
    n_bundles = 0

    for accession in args.accessions:
        raw_dir = PROJECT_ROOT / "data/raw" / accession
        if not raw_dir.exists():
            errors.append(f"Missing raw dir: {raw_dir}")
            continue

        for gz_path in sorted(raw_dir.glob("*.gz")):
            if gz_path.name.endswith("_spatial.tar.gz"):
                try:
                    extract_spatial_bundle(gz_path, overwrite=args.overwrite)
                    n_bundles += 1
                except Exception as exc:
                    errors.append(f"{gz_path}: {exc}")
            else:
                try:
                    decompress_gzip(gz_path, overwrite=args.overwrite)
                    n_gzip += 1
                except Exception as exc:
                    errors.append(f"{gz_path}: {exc}")

    status = "success" if not errors else "partial_success"
    write_status(
        "02_extract",
        status,
        {
            "n_datasets_processed": len(args.accessions),
            "n_gzip_decompressed": n_gzip,
            "n_spatial_bundles_extracted": n_bundles,
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "next_manual_check": [
                "Re-run scripts/02_build_manifest.py after extraction.",
                "Inspect one image and one coordinate file per dataset.",
            ],
        },
    )
    print(f"Decompressed {n_gzip} gzip files")
    print(f"Extracted {n_bundles} spatial bundles")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())

