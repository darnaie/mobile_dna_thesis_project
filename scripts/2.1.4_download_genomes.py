#!/usr/bin/env python3
"""
download_genomes.py  —  Section 2.1.4
========================================
input : BASE/files/non_nested_is_tn.gff3 (Section 2.1.1)
output: BASE/genome_seqs/  (a directory of <taxid.sampleid>.fasta files)

  1. Reads BASE/files/non_nested_is_tn.gff3 and collects the unique
     (taxid.sampleid) genome identifiers referenced by its contig column
     (column 1).
  2. For each unique genome, downloads it from proGenomes3
     (dumpSequence.cgi) with `wget -c`.
  3. Immediately decompresses and renames it to `<taxid.sampleid>.fasta`
     (using the identifier already known from the GFF3, not whatever the
     downloaded FASTA's own header happens to say — see the comment in
     download_and_prepare() for why) in the same step.
  4. Records progress in a manifest TSV (BASE/files/<manifest_file>) so
     reruns are idempotent — already-downloaded genomes are skipped.

If placeholder=true in config.yaml, this instead SYMLINKS
BASE/genome_seqs to the existing directory named by path_to_output (no
copying -- genome directories are far too large for that) and does no
downloading at all.

Run standalone:
    python 2.1.4_download_genomes.py --config path/to/config.yaml
    python 2.1.4_download_genomes.py --gff /path/to/non_nested_is_tn.gff3 --out-dir /path/to/genome_seqs
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.1.4"


def read_manifest(manifest_path: Path) -> dict:
    manifest = {}
    if not manifest_path.exists():
        return manifest
    with open(manifest_path) as fh:
        next(fh, None)  # header
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                manifest[parts[0]] = {"status": parts[1], "filename": parts[2]}
    return manifest


def append_manifest(manifest_path: Path, full_code: str, status: str, filename: str) -> None:
    write_header = not manifest_path.exists()
    with open(manifest_path, "a") as fh:
        if write_header:
            fh.write("full_code\tstatus\tfilename\n")
        fh.write(f"{full_code}\t{status}\t{filename}\n")


def collect_unique_genomes(gff_path: Path) -> set[tuple[str, str]]:
    """Returns a set of (full_code, first_code), e.g. ('573.SAMN03280285', '573')."""
    genomes = set()
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            contig_id = line.split("\t")[0]
            parts = contig_id.split(".")
            if len(parts) < 2:
                continue
            full_code = f"{parts[0]}.{parts[1]}"
            first_code = parts[0]
            genomes.add((full_code, first_code))
    return genomes


def download_and_prepare(
    full_code: str, first_code: str, base_url: str, genome_dir: Path,
    timeout: int, log,
) -> tuple[str, str]:
    """Returns (status, final_filename_or_empty)."""
    url = f"{base_url}?p={full_code}&t=c&a={first_code}"
    tmp_path = genome_dir / f"{full_code}.download.tmp"

    try:
        subprocess.run(
            ["wget", "-c", "-O", str(tmp_path), url],
            check=False, timeout=timeout,
            stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
        )
    except Exception as exc:
        log.warning(f"  [FAIL download] {full_code}: {exc}")
        tmp_path.unlink(missing_ok=True)
        return "failed_download", ""

    if not tmp_path.exists() or tmp_path.stat().st_size == 0:
        log.warning(f"  [FAIL download] {full_code}: empty/missing response")
        tmp_path.unlink(missing_ok=True)
        return "failed_download", ""

    try:
        with open(tmp_path, "rb") as f:
            magic = f.read(2)
        if magic != b"\x1f\x8b":
            log.warning(f"  [FAIL not-gzip] {full_code}")
            tmp_path.unlink(missing_ok=True)
            return "failed_not_gzip", ""

        decompressed = genome_dir / f"{full_code}.decompressed.tmp"
        with gzip.open(tmp_path, "rb") as fin, open(decompressed, "wb") as fout:
            shutil.copyfileobj(fin, fout)
        tmp_path.unlink(missing_ok=True)

        with open(decompressed) as f:
            header = f.readline().strip()
        if not header.startswith(">"):
            decompressed.unlink(missing_ok=True)
            log.warning(f"  [FAIL bad-fasta] {full_code}")
            return "failed_bad_fasta", ""

        # Name the file by the taxid.sampleid we already know (full_code),
        # NOT by whatever the downloaded FASTA's own header happens to say.
        # stage1d_recombinase_seq.py looks genomes up by the first two
        # dot-separated tokens of the FILENAME (taxid.sampleid) to match
        # against the GFF3 contig column — keying off full_code guarantees
        # the two scripts agree, regardless of what convention proGenomes3
        # uses inside the header itself.
        final_path = genome_dir / f"{full_code}.fasta"

        if final_path.exists():
            decompressed.unlink(missing_ok=True)
            return "skipped_exists", final_path.name

        decompressed.rename(final_path)
        return "downloaded", final_path.name

    except Exception as exc:
        log.warning(f"  [FAIL processing] {full_code}: {exc}")
        tmp_path.unlink(missing_ok=True)
        return "failed_processing", ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--gff", default=None, help="Override BASE/files/non_nested_is_tn.gff3")
    ap.add_argument("--out-dir", default=None, help="Override BASE/genome_seqs")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    genome_dir = Path(args.out_dir) if args.out_dir else paths["genome_dir"]
    section_cfg = cfg["section_2_1"][SECTION_KEY]
    dl_cfg = section_cfg

    log = get_logger("2.1.4_download_genomes", paths["logs_dir"])

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Genome download + decompress/rename")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, genome_dir, log):
        n_fasta = len(list(genome_dir.glob("*.fasta"))) if genome_dir.exists() else 0
        log.info(f"{genome_dir} now has {n_fasta:,} .fasta files (via symlink)")
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_1_4_download_genomes"):
        return

    t0 = time.time()
    genome_dir.mkdir(parents=True, exist_ok=True)

    gff_path = Path(args.gff) if args.gff else files_dir / "non_nested_is_tn.gff3"

    log.info(f"Source GFF3     : {gff_path}")
    log.info(f"Genome dir      : {genome_dir}")

    if not gff_path.exists():
        log.error(f"{gff_path} not found — run extract_mges.py (Section 2.1.1) first.")
        raise SystemExit(1)

    genomes = collect_unique_genomes(gff_path)
    log.info(f"Unique genomes referenced in {gff_path.name}: {len(genomes):,}")

    manifest_path = files_dir / dl_cfg["manifest_file"]
    manifest = read_manifest(manifest_path)
    log.info(f"Manifest has {len(manifest):,} prior entries ({manifest_path})")

    stats = Counter()

    for full_code, first_code in sorted(genomes):
        prior = manifest.get(full_code)
        if prior and prior["status"] in ("downloaded", "skipped_exists"):
            final = genome_dir / prior["filename"]
            if prior["filename"] and final.exists():
                stats["already_done"] += 1
                continue

        status, filename = download_and_prepare(
            full_code, first_code, dl_cfg["base_url"], genome_dir,
            dl_cfg["wget_timeout_seconds"], log,
        )
        append_manifest(manifest_path, full_code, status, filename)
        stats[status] += 1

        n_done = sum(stats.values())
        if n_done % 200 == 0:
            log.info(f"  {n_done}/{len(genomes)} genomes handled so far ({dict(stats)})")

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Genomes referenced by GFF3   : {len(genomes):,}")
    log.info(f"Already present (skipped)    : {stats['already_done']:,}")
    log.info(f"Newly downloaded             : {stats['downloaded']:,}")
    log.info(f"Already-existing on disk     : {stats['skipped_exists']:,}")
    log.info(f"Failed (download)            : {stats['failed_download']:,}")
    log.info(f"Failed (not gzip)            : {stats['failed_not_gzip']:,}")
    log.info(f"Failed (bad fasta header)    : {stats['failed_bad_fasta']:,}")
    log.info(f"Failed (other processing)    : {stats['failed_processing']:,}")
    log.info(f"Manifest                     : {manifest_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()