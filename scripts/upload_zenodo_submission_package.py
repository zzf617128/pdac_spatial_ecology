from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_METADATA = ROOT / "docs" / "zenodo_metadata_v2026.07.01-ed10-submission.json"
DEFAULT_FILE = ROOT / "results" / "PDAC_spatial_ecology_source_data.zip"


def request_json(method: str, url: str, token: str, payload: dict | None = None) -> dict:
    data = None
    headers = {
        "Authorization": f"Bearer {token}",
        "Accept": "application/json",
    }
    if payload is not None:
        data = json.dumps(payload).encode("utf-8")
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = resp.read().decode("utf-8")
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise RuntimeError(f"Zenodo API error {exc.code} for {method} {url}: {detail}") from exc


def upload_file(bucket_url: str, token: str, path: Path) -> dict:
    if not path.exists():
        raise FileNotFoundError(path)
    target = bucket_url.rstrip("/") + "/" + urllib.parse.quote(path.name)
    cmd = [
        "curl.exe",
        "--fail-with-body",
        "--retry",
        "3",
        "--retry-delay",
        "10",
        "--connect-timeout",
        "60",
        "--max-time",
        "0",
        "-X",
        "PUT",
        "-H",
        f"Authorization: Bearer {token}",
        "-H",
        "Content-Type: application/octet-stream",
        "-T",
        str(path),
        target,
    ]
    result = subprocess.run(cmd, text=True, capture_output=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Zenodo file upload failed with curl exit code "
            f"{result.returncode}: {result.stderr.strip()} {result.stdout.strip()}"
        )
    output = result.stdout.strip()
    if not output:
        return {}
    try:
        return json.loads(output)
    except json.JSONDecodeError:
        return {"raw_response": output}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Create a Zenodo draft deposition and upload the source-data package."
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--file", type=Path, action="append", help="File to upload. Can be supplied multiple times.")
    parser.add_argument("--deposition-id", help="Reuse an existing Zenodo draft deposition instead of creating a new one.")
    parser.add_argument("--sandbox", action="store_true", help="Use sandbox.zenodo.org instead of production Zenodo.")
    parser.add_argument("--publish", action="store_true", help="Publish the deposition after upload. Default is draft only.")
    args = parser.parse_args()

    token = os.environ.get("ZENODO_TOKEN")
    if not token:
        print("ERROR: Set a fresh Zenodo token in the ZENODO_TOKEN environment variable.", file=sys.stderr)
        return 2

    metadata_path = args.metadata.resolve()
    payload = json.loads(metadata_path.read_text(encoding="utf-8"))
    api_root = "https://sandbox.zenodo.org/api" if args.sandbox else "https://zenodo.org/api"

    if args.deposition_id:
        deposition = request_json("GET", f"{api_root}/deposit/depositions/{args.deposition_id}", token)
        print(f"Using existing Zenodo draft deposition: {args.deposition_id}")
    else:
        deposition = request_json("POST", f"{api_root}/deposit/depositions", token, payload={})
        print(f"Created Zenodo draft deposition: {deposition['id']}")
    dep_id = deposition["id"]
    bucket_url = deposition["links"]["bucket"]

    files = args.file if args.file else [DEFAULT_FILE]
    for raw_file in files:
        path = raw_file.resolve()
        print(f"Uploading {path.name} ({path.stat().st_size} bytes)...")
        upload_file(bucket_url, token, path)

    request_json("PUT", f"{api_root}/deposit/depositions/{dep_id}", token, payload=payload)
    print("Updated Zenodo metadata.")

    record = request_json("GET", f"{api_root}/deposit/depositions/{dep_id}", token)
    html = record.get("links", {}).get("html", "")
    doi = record.get("metadata", {}).get("prereserve_doi", {}).get("doi", "")
    print(f"Draft URL: {html}")
    if doi:
        print(f"Reserved DOI: {doi}")

    if args.publish:
        print("Publishing deposition...")
        published = request_json("POST", f"{api_root}/deposit/depositions/{dep_id}/actions/publish", token)
        print(f"Published URL: {published.get('links', {}).get('html', html)}")
        print(f"DOI: {published.get('doi', doi)}")
    else:
        print("Draft left unpublished for manual review. Re-run with --publish only after checking the Zenodo draft.")

    time.sleep(0.1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
