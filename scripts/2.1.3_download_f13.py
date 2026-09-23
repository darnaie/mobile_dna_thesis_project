#!/usr/bin/env python3
"""
download_f13.py  —  Section 2.1.3
====================================
input : none
output: BASE/files/f13_dataset.txt

Stages the proGenomes3 taxonomy/specI lookup table
(GCA -> specI -> taxid -> species/strain, tab-separated) as
BASE/files/f13_dataset.txt.

HONEST NOTE: section_2_1["2.1.3"].download_url in config.yaml is
INVENTED -- there is no confirmed real source for this table yet (ask
your supervisor). Real download logic IS implemented below (fetch ->
gunzip if needed -> validate column count -> write out), so once you have
the real URL you can just update download_url and flip placeholder to
false. Until then, leave placeholder=true and point path_to_output at
your existing copy, or run with placeholder=false to see it fail loudly
against the made-up URL (which is expected and fine -- better than
silently producing nothing).

Run standalone:
    python 2.1.3_download_f13.py --config path/to/config.yaml
    python 2.1.3_download_f13.py --url https://... --out-dir ./out
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import sys
import urllib.request
import urllib.error
from pathlib import Path

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.1.3"


def download_file(url: str, dest: Path, log, timeout: int = 120) -> None:
    log.info(f"Downloading {url} -> {dest} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "thesis_rep/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as response, open(dest, "wb") as out_fh:
        shutil.copyfileobj(response, out_fh)


def validate_table(path: Path, log) -> int:
    """Sanity check: every non-empty, non-comment line should have >= 4
    tab-separated columns (GCA, specI, taxid, species...). Returns line count."""
    n_lines = 0
    n_bad = 0
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            n_lines += 1
            if len(line.split("\t")) < 4:
                n_bad += 1
    if n_bad:
        log.warning(f"{n_bad}/{n_lines} lines have fewer than 4 tab-separated columns "
                    f"-- the downloaded file may not be in the expected format.")
    return n_lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--url", default=None, help="Override the download URL")
    ap.add_argument("--out-dir", default=None, help="Override BASE/files")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = Path(args.out_dir) if args.out_dir else paths["files_dir"]
    files_dir.mkdir(parents=True, exist_ok=True)
    section_cfg = cfg["section_2_1"][SECTION_KEY]

    log = get_logger("2.1.3_download_f13", paths["logs_dir"])

    out_path = files_dir / "f13_dataset.txt"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Taxonomy/specI dataset staging")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        n_lines = validate_table(out_path, log)
        log.info(f"Staged {out_path} ({n_lines:,} data lines)")
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_1_3_download_f13"):
        return

    url = args.url or section_cfg.get("download_url")
    if not url:
        log.error("placeholder=false but no download_url is set in config.yaml "
                   "(section_2_1['2.1.3'].download_url) -- nothing to download.")
        raise SystemExit(1)

    log.warning(
        "download_url is an INVENTED placeholder (no confirmed real source yet) "
        "-- this download is expected to fail until you have the real URL. "
        "See the script docstring."
    )

    tmp_path = files_dir / "f13_dataset.download.tmp"
    try:
        download_file(url, tmp_path, log)
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, OSError) as exc:
        log.error(f"Download failed: {exc}")
        tmp_path.unlink(missing_ok=True)
        raise SystemExit(1)

    # gunzip if the response was gzip-compressed
    with open(tmp_path, "rb") as f:
        magic = f.read(2)
    if magic == b"\x1f\x8b":
        log.info("Response is gzip-compressed -- decompressing ...")
        with gzip.open(tmp_path, "rb") as fin, open(out_path, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        tmp_path.unlink()
    else:
        tmp_path.rename(out_path)

    n_lines = validate_table(out_path, log)
    log.info(f"Wrote {out_path} ({n_lines:,} data lines)")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()