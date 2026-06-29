from __future__ import annotations

import argparse
import csv
import hashlib
import html.parser
import json
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
MVP_ACCESSIONS = ["GSE272362", "GSE282302", "GSE274103"]
DEFAULT_ACCESSIONS = [
    "GSE272362",
    "GSE282302",
    "GSE274103",
    "GSE235315",
    "GSE202740",
    "GSE111672",
]
USER_AGENT = "pdac-spatial-ecology-mvp/0.1"


@dataclass
class RemoteFile:
    accession: str
    file_name: str
    url: str
    file_size: str | None = None


class LinkParser(html.parser.HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.links: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        if tag.lower() != "a":
            return
        href = dict(attrs).get("href")
        if href:
            self.links.append(href)


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def accession_to_suppl_url(accession: str) -> str:
    digits = "".join(ch for ch in accession if ch.isdigit())
    if not accession.startswith("GSE") or len(digits) < 3:
        raise ValueError(f"Unsupported GEO accession: {accession}")
    prefix = f"GSE{digits[:3]}nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}/{accession}/suppl/"


def request_url(url: str, method: str = "GET", timeout: int = 60) -> urllib.response.addinfourl:
    request = urllib.request.Request(url, method=method, headers={"User-Agent": USER_AGENT})
    return urllib.request.urlopen(request, timeout=timeout)


def list_remote_files(accession: str) -> list[RemoteFile]:
    base_url = accession_to_suppl_url(accession)
    with request_url(base_url) as response:
        html = response.read().decode("utf-8", errors="replace")

    parser = LinkParser()
    parser.feed(html)

    files: list[RemoteFile] = []
    for href in parser.links:
        if href in {"../", "/geo/series/"} or href.endswith("/") or href == "index.html":
            continue
        file_url = urllib.parse.urljoin(base_url, href)
        file_name = urllib.parse.unquote(Path(urllib.parse.urlparse(file_url).path).name)
        if not file_name:
            continue
        files.append(RemoteFile(accession=accession, file_name=file_name, url=file_url))

    for remote_file in files:
        remote_file.file_size = get_remote_size(remote_file.url)
    return files


def get_remote_size(url: str) -> str | None:
    try:
        with request_url(url, method="HEAD", timeout=30) as response:
            value = response.headers.get("Content-Length")
            return value if value else None
    except Exception:
        return None


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download_file(remote_file: RemoteFile, dest_dir: Path, overwrite: bool = False) -> tuple[str, str]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    local_path = dest_dir / remote_file.file_name
    if local_path.exists() and not overwrite:
        return "exists", sha256_file(local_path)

    tmp_path = local_path.with_suffix(local_path.suffix + ".part")
    with request_url(remote_file.url, timeout=300) as response, tmp_path.open("wb") as output:
        shutil.copyfileobj(response, output)
    tmp_path.replace(local_path)
    return "downloaded", sha256_file(local_path)


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


def write_audit_report(accessions: list[str], rows: list[dict], errors: list[str], audit_only: bool) -> None:
    report_path = PROJECT_ROOT / "docs/dataset_audit_report.md"
    lines = [
        "# Dataset Audit Report",
        "",
        f"Last updated UTC: {now_iso()}",
        "",
        "## Scope",
        "",
        "First-pass MVP accessions:",
        "",
    ]
    lines.extend([f"- {accession}" for accession in accessions])
    lines.extend(
        [
            "",
            "Mode:",
            "",
            f"- audit_only: {str(audit_only).lower()}",
            "",
            "## GEO Supplementary File Audit",
            "",
            "| accession | n_remote_files | n_downloaded_or_existing | notes |",
            "|---|---:|---:|---|",
        ]
    )

    by_accession: dict[str, list[dict]] = {accession: [] for accession in accessions}
    for row in rows:
        by_accession.setdefault(row["accession"], []).append(row)

    for accession in accessions:
        accession_rows = by_accession.get(accession, [])
        n_remote = sum(1 for row in accession_rows if row["download_status"] != "list_failed")
        n_local = sum(1 for row in accession_rows if row["download_status"] in {"downloaded", "exists"})
        notes = "ok" if accession_rows else "no rows"
        lines.append(f"| {accession} | {n_remote} | {n_local} | {notes} |")

    if errors:
        lines.extend(["", "## Errors", ""])
        lines.extend([f"- {error}" for error in errors])

    lines.extend(
        [
            "",
            "## Bias-Control Notes",
            "",
            "- This audit does not prove analyzability; it only checks remote supplementary file visibility.",
            "- Samples without expression, coordinates, and paired image files must be excluded from H&E claims.",
            "- GSE282302 treatment context must be resolved at sample level before any treatment-related claim.",
            "",
        ]
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text("\n".join(lines), encoding="utf-8")


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Audit or download GEO supplementary files.")
    parser.add_argument("--accessions", nargs="+", default=None)
    parser.add_argument("--mvp", action="store_true", help="Use first-pass MVP accessions.")
    parser.add_argument("--audit-only", action="store_true", help="List remote files without downloading.")
    parser.add_argument("--download", action="store_true", help="Download remote files.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files.")
    parser.add_argument("--sleep", type=float, default=0.2, help="Seconds to sleep between accessions.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv or sys.argv[1:])
    if args.mvp:
        accessions = MVP_ACCESSIONS
    elif args.accessions:
        accessions = args.accessions
    else:
        accessions = DEFAULT_ACCESSIONS

    audit_only = args.audit_only or not args.download
    log_path = PROJECT_ROOT / "results/logs/download_log.tsv"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    rows: list[dict] = []
    errors: list[str] = []

    for accession in accessions:
        raw_dir = PROJECT_ROOT / "data/raw" / accession
        raw_dir.mkdir(parents=True, exist_ok=True)
        try:
            remote_files = list_remote_files(accession)
            if not remote_files:
                message = f"No supplementary files listed for {accession}"
                errors.append(message)
                rows.append(
                    {
                        "accession": accession,
                        "file_name": "",
                        "file_size": "",
                        "remote_url": accession_to_suppl_url(accession),
                        "local_path": "",
                        "download_status": "no_remote_files",
                        "checksum": "",
                        "error_message": message,
                    }
                )
                (raw_dir / "MANUAL_DOWNLOAD_REQUIRED.txt").write_text(message + "\n", encoding="utf-8")
                continue

            for remote_file in remote_files:
                status = "remote_listed"
                checksum = ""
                error_message = ""
                local_path = raw_dir / remote_file.file_name
                if not audit_only:
                    try:
                        status, checksum = download_file(remote_file, raw_dir, overwrite=args.overwrite)
                    except Exception as exc:
                        status = "download_failed"
                        error_message = str(exc)
                        errors.append(f"{accession}: {remote_file.file_name}: {exc}")
                rows.append(
                    {
                        "accession": accession,
                        "file_name": remote_file.file_name,
                        "file_size": remote_file.file_size or "",
                        "remote_url": remote_file.url,
                        "local_path": str(local_path) if local_path.exists() else "",
                        "download_status": status,
                        "checksum": checksum,
                        "error_message": error_message,
                    }
                )
        except Exception as exc:
            message = f"{accession}: failed to list supplementary files: {exc}"
            errors.append(message)
            rows.append(
                {
                    "accession": accession,
                    "file_name": "",
                    "file_size": "",
                    "remote_url": "",
                    "local_path": "",
                    "download_status": "list_failed",
                    "checksum": "",
                    "error_message": str(exc),
                }
            )
            (raw_dir / "MANUAL_DOWNLOAD_REQUIRED.txt").write_text(message + "\n", encoding="utf-8")
        time.sleep(args.sleep)

    fieldnames = [
        "accession",
        "file_name",
        "file_size",
        "remote_url",
        "local_path",
        "download_status",
        "checksum",
        "error_message",
    ]
    with log_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    write_audit_report(accessions, rows, errors, audit_only=audit_only)

    status = "success" if not errors else "partial_success"
    write_status(
        "01",
        status,
        {
            "n_datasets_processed": len(accessions),
            "n_remote_files": sum(1 for row in rows if row["download_status"] != "list_failed"),
            "n_errors": len(errors),
            "critical_errors": [],
            "noncritical_warnings": errors,
            "audit_only": audit_only,
            "download_log": str(log_path),
            "next_manual_check": [
                "Inspect docs/dataset_audit_report.md.",
                "Confirm which files are expression matrices, coordinates, scalefactors, and H&E images.",
                "Run with --download only after the audit looks reasonable.",
            ],
        },
    )

    print(f"Audited {len(accessions)} accessions")
    print(f"Wrote {log_path}")
    print(f"Status: {status}")
    return 0 if not errors else 1


if __name__ == "__main__":
    raise SystemExit(main())
