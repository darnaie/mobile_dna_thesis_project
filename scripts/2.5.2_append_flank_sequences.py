#!/usr/bin/env python3
"""
2.5.2_append_flank_sequences.py  —  Section 2.5.2
======================================================
input : BASE/files/c_mge_recombinase_table.tsv (2.4.2),
        BASE/flanks_split_acc_to_duplicates/ (2.5.1)
output: BASE/files/d_mge_recombinase_table.tsv

For every row of c_mge_recombinase_table.tsv, writes TWO adjacent output
rows -- one for its LEFT flank, one for its RIGHT flank -- each carrying
all of that row's original columns plus:

BUG FIX: this script's own COLUMNS passthrough list previously didn't
include "cluster_size" or "duplicate_group_size" (both added upstream,
in Sections 2.3.2/2.4.2, after this script was first written) -- those
columns were silently being dropped here. They're now carried through
unchanged.

  col 17  side            : "left" or "right"
  col 18  flank_sequence  : the actual nucleotide sequence, read off the
                             matching *_left.fasta / *_right.fasta file
                             2.5.1 wrote for that mge_id.

The file is located using mge_id + side together with rec_type/cluster_number/
duplicate_group_number/singleton from THIS SAME row -- via
utils.flank_output_dir(), the exact same directory-naming logic 2.5.1
used when it WROTE the files, so a lookup here always lands in the same
place 2.5.1 put it (no separate/duplicated path logic to drift out of sync).

Run standalone:
    python 2.5.2_append_flank_sequences.py --config path/to/config.yaml
    python 2.5.2_append_flank_sequences.py --table c_table.tsv --flanks-dir ./flanks --out d_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, flank_output_dir,
    resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.5.2"

COLUMNS = [
    "mge_id", "mge_start", "mge_end", "n_genes", "mgeR",
    "recombinase_id", "recombinase_start", "recombinase_end", "recombinase_type",
    "is_representative", "cluster_number", "specI", "cluster_size",
    "singleton", "duplicate_group_number", "duplicate_group_size",
]


def read_flank_sequence(path: Path) -> str | None:
    """Returns None if the file doesn't exist -- caller decides how to
    represent that in the table (this script uses "NA")."""
    if not path.exists():
        return None
    seq_lines = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            seq_lines.append(line.strip())
    return "".join(seq_lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--table", default=None, help="Override BASE/files/c_mge_recombinase_table.tsv")
    ap.add_argument("--flanks-dir", default=None, help="Override BASE/flanks_split_acc_to_duplicates")
    ap.add_argument("--out", default=None, help="Override BASE/files/d_mge_recombinase_table.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_5"][SECTION_KEY]

    log = get_logger("2.5.2_append_flank_sequences", paths["logs_dir"])

    table_in = Path(args.table) if args.table else files_dir / "c_mge_recombinase_table.tsv"
    flanks_dir = Path(args.flanks_dir) if args.flanks_dir else paths["base_dir"] / "flanks_split_acc_to_duplicates"
    out_path = Path(args.out) if args.out else files_dir / "d_mge_recombinase_table.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Append flank sequences (doubles every row: left + right)")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_5_2_append_flank_sequences"):
        return

    if not table_in.exists():
        log.error(f"{table_in} not found -- run 2.4.2_append_duplicate_info.py first.")
        raise SystemExit(1)
    if not flanks_dir.exists():
        log.error(f"{flanks_dir} not found -- run 2.5.1_extract_flanks.py first.")
        raise SystemExit(1)

    t0 = time.time()

    log.info(f"Loading input table: {table_in}")
    with open(table_in, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    log.info(f"  Rows loaded: {len(rows):,}")

    n_found = n_missing = 0
    out_rows = []

    for row in rows:
        mge_id = row.get("mge_id", "")
        rec_type = row.get("recombinase_type", "")
        cluster_number = row.get("cluster_number", "NA")
        duplicate_group_number = row.get("duplicate_group_number", "NA")
        singleton = row.get("singleton", "yes")

        mge_dir = flank_output_dir(flanks_dir, rec_type, cluster_number, duplicate_group_number, singleton)

        for side in ("left", "right"):
            flank_path = mge_dir / f"{mge_id}_{side}.fasta"
            seq = read_flank_sequence(flank_path)
            if seq is None:
                n_missing += 1
                seq = "NA"
            else:
                n_found += 1

            new_row = {k: row.get(k, "") for k in COLUMNS}
            new_row["side"] = side
            new_row["flank_sequence"] = seq
            out_rows.append(new_row)

    log.info(f"Writing output table: {out_path}")
    out_columns = COLUMNS + ["side", "flank_sequence"]
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Input rows                     : {len(rows):,}")
    log.info(f"Output rows (2x: left + right) : {len(out_rows):,}")
    log.info(f"  flank sequence found         : {n_found:,}")
    log.info(f"  flank sequence MISSING       : {n_missing:,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()