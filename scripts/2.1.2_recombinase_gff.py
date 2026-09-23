#!/usr/bin/env python3
"""
recombinase_gff.py  —  Section 2.1.2
=======================================
input : BASE/files/master_mges.gff3 (for island typing), paths.progenomes3_dir
        (for the actual gene rows)
output: 1. BASE/files/recomb.gff3                        (master, any MGE type)
        2. BASE/files/non_nested_is_tn_recombinase.gff3

GFF-only now -- sequence extraction is a separate step (Section 2.1.5,
recombinase_seq.py), since it needs genome FASTAs that aren't downloaded
until Section 2.1.4.

  Pass 1 -- read BASE/files/master_mges.gff3 ONCE (not a re-walk of every
            raw per-genome file) to record, per MGE island ID: whether its
            'mge=' field contains is_tn, and its mge_type (nested /
            non-nested).
  Pass 2 -- walk paths.progenomes3_dir's raw per-genome GFF3 files (still
            needed, since gene-level rows aren't present in
            master_mges.gff3) for every gene row carrying a
            'recombinase=' attribute, and look up its Parent= island from
            Pass 1 to classify it.

Run standalone:
    python 2.1.2_recombinase_gff.py --config path/to/config.yaml
    python 2.1.2_recombinase_gff.py --pg3-dir /path/to/pg3 --master-gff master_mges.gff3 --out-dir ./out
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from utils import (
    get_logger,
    get_paths,
    load_config,
    parse_gff_attributes,
    parse_mge_field,
    resolve_placeholder_multi,
    maybe_submit_and_exit,
)

SECTION_KEY = "2.1.2"

# ──────────────────────────────────────────────────────────────────────────
# Pass 1 — loaded ONCE per worker process via the pool initializer, not
# re-walked per file.
# ──────────────────────────────────────────────────────────────────────────

_worker_island_info: dict[str, tuple[bool, str]] = {}


def load_island_info(master_gff_path: str, mge_filter: str) -> dict[str, tuple[bool, str]]:
    island_info: dict[str, tuple[bool, str]] = {}
    with open(master_gff_path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.split("\t")
            if len(cols) < 9 or cols[2] != "mobile_genetic_element":
                continue
            attrs = parse_gff_attributes(cols[8])
            mge_id = attrs.get("ID", "")
            if not mge_id:
                continue
            counts = parse_mge_field(attrs.get("mge", ""))
            has_filter = mge_filter in counts
            mge_type = attrs.get("mge_type", "")
            island_info[mge_id] = (has_filter, mge_type)
    return island_info


def _init_worker(master_gff_path: str, mge_filter: str) -> None:
    global _worker_island_info
    _worker_island_info = load_island_info(master_gff_path, mge_filter)


# ──────────────────────────────────────────────────────────────────────────
# Pass 2 — recombinase gene rows
# ──────────────────────────────────────────────────────────────────────────

def process_gff3(path: str) -> dict:
    result = {
        "total_rows": [],
        "non_nested_is_tn_rows": [],
        "stats": Counter(),
    }
    result["stats"]["gff3_files"] += 1

    with open(path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = parse_gff_attributes(cols[8])
            if "recombinase" not in attrs:
                continue
            parent = attrs.get("Parent", "")
            info = _worker_island_info.get(parent)
            if info is None:
                result["stats"]["orphan_recombinase_rows"] += 1
                continue

            has_filter, mge_type = info
            row = raw if raw.endswith("\n") else raw + "\n"

            result["total_rows"].append(row)
            result["stats"]["total_recombinase_rows"] += 1

            if has_filter and mge_type == "non-nested":
                result["non_nested_is_tn_rows"].append(row)
                result["stats"]["non_nested_is_tn_recombinase_rows"] += 1

    return result


# ──────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--pg3-dir", default=None, help="Override paths.progenomes3_dir")
    ap.add_argument("--master-gff", default=None, help="Override BASE/files/master_mges.gff3")
    ap.add_argument("--out-dir", default=None, help="Override BASE/files")
    ap.add_argument("--threads", type=int, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = Path(args.out_dir) if args.out_dir else paths["files_dir"]
    files_dir.mkdir(parents=True, exist_ok=True)
    pg3_dir = Path(args.pg3_dir) if args.pg3_dir else paths["progenomes3_dir"]
    n_threads = args.threads or os.cpu_count() or 4
    section_cfg = cfg["section_2_1"][SECTION_KEY]
    mge_filter = cfg["section_2_1"]["2.1.1"]["mge_type_filter"]

    log = get_logger("2.1.2_recombinase_gff", paths["logs_dir"])

    total_gff_out = files_dir / "recomb.gff3"
    nn_gff_out = files_dir / "non_nested_is_tn_recombinase.gff3"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Recombinase GFF extraction (master + non-nested-is_tn only)")
    log.info("=" * 70)

    if resolve_placeholder_multi(section_cfg, {
        "path_to_output": total_gff_out,
        "path_to_output_non_nested_gff": nn_gff_out,
    }, log):
        log.info(f"Section {SECTION_KEY} satisfied entirely via PLACEHOLDER (both outputs staged).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_1_2_recombinase_gff"):
        return

    master_gff_path = Path(args.master_gff) if args.master_gff else files_dir / "master_mges.gff3"
    if not master_gff_path.exists():
        log.error(f"{master_gff_path} not found -- run extract_mges.py (Section 2.1.1) first.")
        raise SystemExit(1)

    t0 = time.time()

    files = []
    for root, _dirs, fnames in os.walk(pg3_dir):
        for fname in fnames:
            if fname.endswith(".gff3"):
                files.append(os.path.join(root, fname))
    log.info(f"Discovered {len(files)} GFF3 files under {pg3_dir}")
    if not files:
        log.error("No GFF3 files found -- check paths.progenomes3_dir in config.yaml")
        raise SystemExit(1)

    total_out = open(total_gff_out, "w")
    total_out.write("##gff-version 3\n")
    nn_out = open(nn_gff_out, "w")
    nn_out.write("##gff-version 3\n")

    total_stats = Counter()

    with ProcessPoolExecutor(
        max_workers=n_threads, initializer=_init_worker,
        initargs=(str(master_gff_path), mge_filter),
    ) as pool:
        futures = {pool.submit(process_gff3, p): p for p in files}
        for i, fut in enumerate(as_completed(futures), 1):
            src = futures[fut]
            try:
                res = fut.result()
            except Exception as exc:
                log.warning(f"Worker crashed for {src}: {exc}")
                continue

            for row in res["total_rows"]:
                total_out.write(row)
            for row in res["non_nested_is_tn_rows"]:
                nn_out.write(row)

            total_stats.update(res["stats"])

            if i % 500 == 0 or i == len(files):
                log.info(f"  {i:>6}/{len(files)} files | "
                         f"{total_stats['total_recombinase_rows']:,} total recombinase rows so far")

    total_out.close()
    nn_out.close()

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Genomes processed (gff3 files) : {total_stats['gff3_files']:,}")
    log.info(f"Total recombinase gene rows (any parent MGE type) : "
             f"{total_stats['total_recombinase_rows']:,}")
    log.info(f"Orphan recombinase rows (Parent not found)        : "
             f"{total_stats['orphan_recombinase_rows']:,}")
    log.info(f"non-nested is_tn recombinase rows                 : "
             f"{total_stats['non_nested_is_tn_recombinase_rows']:,}")
    log.info("-" * 70)
    log.info(f"recomb.gff3                          : {total_gff_out}")
    log.info(f"non_nested_is_tn_recombinase.gff3     : {nn_gff_out}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()