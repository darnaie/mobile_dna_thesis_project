#!/usr/bin/env python3
"""
2.1.1_extract_mges.py  —  Section 2.1.1
====================================
input : paths.progenomes3_dir
output: files/master_mges.gff3
        files/non_nested_is_tn.gff3
        files/non_nested_is_tn.fasta[.gz]

Walks paths.progenomes3_dir for every per-genome .gff3 and .ffn.gz file.

No longer configurable to nested/merged variants (that choice has been
removed per the thesis-methods writeup, which only ever uses non-nested
is_tn islands) -- this script always produces exactly the three files
above, nothing else.

  1. master_mges.gff3
       EVERY mobile_genetic_element line, all MGE types, unfiltered.
       (Transposon isolation happens below, not here.)
  2. non_nested_is_tn.gff3
       mobile_genetic_element lines whose 'mge=' field contains an is_tn
       entry AND whose mge_type is non-nested.
  3. non_nested_is_tn.fasta[.gz]
       the corresponding *gene-level* nucleotide sequences for (2), read
       out of the matching .ffn.gz files.

BUG FIX (see chat -- this was producing an empty .fasta): the .ffn.gz
header matching now goes through utils.parse_ffn_header_mge_fields(),
which lower-cases the WHOLE header before searching for 'mge=' /
'mge_type=', mirroring the original (known-working) pipeline script's
`header.lower()` approach. The previous version of this script used a
structured key=value parser that did not normalize case, so any case
variation in the real .ffn.gz headers (e.g. 'Mge=IS_TN:1') silently
matched nothing and produced an empty output file with no error. This
version also logs (a) a sample of raw headers actually seen, and (b) a
hard warning if it finds MGE rows in the .gff3 files but extracts zero
sequences from the .ffn.gz files, so a repeat of this bug would be loud
instead of silent.

BUG FIX #2 (see chat -- placeholder key mismatch): the config keys this
script asked resolve_placeholder_multi() for (path_to_output,
path_to_output_non_nested_gff, path_to_output_non_nested_fasta) didn't
match the more descriptive names actually used in config.yaml
(path_to_output_master_mges, path_to_output_non_nested_is_tn_gff,
path_to_output_non_nested_is_tn_fasta) -- so placeholder mode failed with
"3 of 3 required path(s) are not set" even when all three WERE set, just
under different key names. resolve_placeholder_multi() itself (in
utils.py) is fine -- it's fully generic, takes whatever output_map you
give it, and correctly requires every key in that map to resolve to a
real path when placeholder=true; the mismatch was entirely in what keys
THIS script was asking for. Fixed by using the same key names config.yaml
already has, rather than changing utils.py or asking you to rename your
config.

Counting note: a single mobile_genetic_element row can bundle more than
one recombinase, e.g. 'mge=is_tn:2,ce:1' is ONE row but represents 2
is_tn + 1 ce recombinase instances. The summary log reports BOTH the
number of matching rows/islands AND the summed instance count -- do not
read the row count alone as "how many is_tn there are".

Run:
    python 2.1.1_extract_mges.py [--config path/to/config.yaml]
"""

from __future__ import annotations

import argparse
import gzip
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from utils import (
    load_config,
    get_paths,
    get_logger,
    parse_gff_attributes,
    parse_mge_field,
    parse_ffn_header_mge_fields,
    resolve_placeholder_multi,
    maybe_submit_and_exit,
)

SECTION_KEY = "2.1.1"
MAX_SAMPLE_HEADERS = 10


# ──────────────────────────────────────────────────────────────────────────
# Worker (runs in a subprocess — keep it self-contained / picklable)
# ──────────────────────────────────────────────────────────────────────────

def process_file(path: str, mge_filter: str) -> dict:
    result = {
        "gff_all": [],
        "gff_non_nested": [],
        "ffn_non_nested": [],
        "sample_headers": [],
        "stats": Counter(),
    }

    if path.endswith(".gff3"):
        result["stats"]["gff3_files"] += 1
        with open(path) as fh:
            for line in fh:
                if "\tmobile_genetic_element\t" not in line:
                    continue
                result["gff_all"].append(line)
                result["stats"]["mge_rows_total"] += 1

                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    continue
                attrs = parse_gff_attributes(cols[8])
                mge_type = attrs.get("mge_type", "")
                counts = parse_mge_field(attrs.get("mge", ""))

                if mge_filter not in counts or mge_type != "non-nested":
                    continue

                result["gff_non_nested"].append(line)
                result["stats"]["is_tn_rows_non_nested"] += 1
                result["stats"]["is_tn_instances_non_nested"] += counts[mge_filter]

    elif path.endswith(".ffn.gz"):
        result["stats"]["ffn_files"] += 1
        with gzip.open(path, "rt") as fh:
            header, seq = None, []
            for raw in fh:
                raw = raw.rstrip("\n")
                if raw.startswith(">"):
                    if header is not None:
                        _consume_ffn_record(header, "".join(seq), mge_filter, result)
                    header, seq = raw, []
                else:
                    seq.append(raw)
            if header is not None:
                _consume_ffn_record(header, "".join(seq), mge_filter, result)

    return result


def _consume_ffn_record(header: str, seq: str, mge_filter: str, result: dict) -> None:
    result["stats"]["ffn_records_total"] += 1
    if len(result["sample_headers"]) < 2:
        result["sample_headers"].append(header)

    counts, mge_type = parse_ffn_header_mge_fields(header)
    if mge_filter not in counts or mge_type != "non-nested":
        return

    result["ffn_non_nested"].append(f"{header}\n{seq}\n")
    result["stats"]["is_tn_seqs_non_nested"] += 1


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None, help="Path to config.yaml")
    ap.add_argument("--pg3-dir", default=None, help="Override paths.progenomes3_dir")
    ap.add_argument("--out-dir", default=None, help="Override BASE/files")
    ap.add_argument("--threads", type=int, default=None, help="Worker processes (default: CPU count)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = Path(args.out_dir) if args.out_dir else paths["files_dir"]
    files_dir.mkdir(parents=True, exist_ok=True)
    pg3_dir = Path(args.pg3_dir) if args.pg3_dir else paths["progenomes3_dir"]
    n_threads = args.threads or os.cpu_count() or 4
    section_cfg = cfg["section_2_1"][SECTION_KEY]
    mge_filter = section_cfg["mge_type_filter"]
    gzip_out = section_cfg["gzip_fasta_output"]

    log = get_logger("2.1.1_extract_mges", paths["logs_dir"])

    master_gff_out = files_dir / "master_mges.gff3"
    non_nested_gff_out = files_dir / "non_nested_is_tn.gff3"
    ext = "fasta.gz" if gzip_out else "fasta"
    non_nested_fasta_out = files_dir / f"non_nested_is_tn.{ext}"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — MGE aggregation + is_tn isolation (non-nested only)")
    log.info("=" * 70)

    # BUG FIX (see chat / module docstring): these three key names now
    # match what config.yaml actually calls them
    # (path_to_output_master_mges / path_to_output_non_nested_is_tn_gff /
    # path_to_output_non_nested_is_tn_fasta), instead of the shorter,
    # generic names this used to ask for (path_to_output /
    # path_to_output_non_nested_gff / path_to_output_non_nested_fasta),
    # which didn't exist under those names in config.yaml at all and made
    # placeholder mode fail even when every path WAS set.
    if resolve_placeholder_multi(section_cfg, {
        "path_to_output_master_mges": master_gff_out,
        "path_to_output_non_nested_is_tn_gff": non_nested_gff_out,
        "path_to_output_non_nested_is_tn_fasta": non_nested_fasta_out,
    }, log):
        log.info(f"Section {SECTION_KEY} satisfied entirely via PLACEHOLDER (all three outputs staged).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_1_1_extract_mges"):
        return

    log.info(f"proGenomes3 dir : {pg3_dir}")
    log.info(f"Output dir      : {files_dir}")
    log.info(f"mge_type_filter : {mge_filter}")

    if not pg3_dir.exists():
        log.error(f"progenomes3_dir does not exist: {pg3_dir}")
        raise SystemExit(1)

    t0 = time.time()
    files = []
    for root, _dirs, fnames in os.walk(pg3_dir):
        for fname in fnames:
            if fname.endswith(".gff3") or fname.endswith(".ffn.gz"):
                files.append(os.path.join(root, fname))

    log.info(f"Discovered {len(files)} files (.gff3 + .ffn.gz) under proGenomes3 dir")
    if not files:
        log.error("No input files found — check paths.progenomes3_dir in config.yaml")
        raise SystemExit(1)

    gff_all_out = open(master_gff_out, "w")
    gff_nn_out = open(non_nested_gff_out, "w")
    fasta_opener = gzip.open if gzip_out else open
    fasta_mode = "wt" if gzip_out else "w"
    ffn_nn_out = fasta_opener(non_nested_fasta_out, fasta_mode)

    total_stats = Counter()
    sample_headers: list[str] = []

    with ProcessPoolExecutor(max_workers=n_threads) as exe:
        futures = {exe.submit(process_file, f, mge_filter): f for f in files}

        for i, fut in enumerate(as_completed(futures), 1):
            src = futures[fut]
            try:
                res = fut.result()
            except Exception as exc:
                log.warning(f"Worker failed for {src}: {exc}")
                continue

            for line in res["gff_all"]:
                gff_all_out.write(line)
            for line in res["gff_non_nested"]:
                gff_nn_out.write(line)
            for rec in res["ffn_non_nested"]:
                ffn_nn_out.write(rec)

            if len(sample_headers) < MAX_SAMPLE_HEADERS:
                sample_headers.extend(res["sample_headers"])

            total_stats.update(res["stats"])

            if i % 500 == 0 or i == len(files):
                log.info(f"  Processed {i}/{len(files)} files")

    gff_all_out.close()
    gff_nn_out.close()
    ffn_nn_out.close()

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Genomes processed (gff3 files)  : {total_stats['gff3_files']:,}")
    log.info(f".ffn.gz files processed         : {total_stats['ffn_files']:,}")
    log.info(f"Total MGE island rows (unfiltered, all types) : {total_stats['mge_rows_total']:,}")
    log.info(f"non_nested is_tn islands (rows)         : {total_stats['is_tn_rows_non_nested']:,}")
    log.info(f"non_nested is_tn instances (summed)     : {total_stats['is_tn_instances_non_nested']:,}")
    log.info(f"Total .ffn.gz gene records scanned      : {total_stats['ffn_records_total']:,}")
    log.info(f"non_nested is_tn sequences extracted    : {total_stats['is_tn_seqs_non_nested']:,}")

    if total_stats["is_tn_rows_non_nested"] > 0 and total_stats["is_tn_seqs_non_nested"] == 0:
        log.warning(
            "*** non_nested_is_tn.fasta is EMPTY despite finding non-nested is_tn "
            "rows in the GFF3 files. This is exactly the bug this rewrite targeted "
            "-- if you still see this, the real .ffn.gz header format differs from "
            "what parse_ffn_header_mge_fields() expects. Sample raw headers seen:"
        )
        for h in sample_headers:
            log.warning(f"    {h}")
    elif sample_headers:
        log.info("Sample .ffn.gz headers seen (for reference):")
        for h in sample_headers[:3]:
            log.info(f"    {h}")

    log.info("-" * 70)
    log.info(f"master_mges.gff3        : {master_gff_out}")
    log.info(f"non_nested_is_tn.gff3   : {non_nested_gff_out}")
    log.info(f"non_nested_is_tn.{ext:<6}: {non_nested_fasta_out}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()