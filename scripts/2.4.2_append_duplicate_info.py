#!/usr/bin/env python3
"""
append_duplicate_info.py  —  Section 2.4.2
==============================================
input : BASE/files/b_mge_recombinase_table.tsv (Section 2.3.2),
        BASE/vclust_deduplication/*.fasta.duplicates.txt (Section 2.4.1)
output: BASE/files/c_mge_recombinase_table.tsv

Appends two columns to every row of b_mge_recombinase_table.tsv, keyed on
mge_id (column 1) against the vclust duplicates file (whose tokens are
mge_ids, since Section 2.4.1 deduplicated whole is_tn elements, not
individual recombinase genes):

BUG FIX: this script's own COLUMNS passthrough list previously didn't
include "cluster_size" (added to b_mge_recombinase_table.tsv by Section
2.3.2 after this script was first written) -- that column was silently
being dropped here. It's now carried through unchanged, so
c_mge_recombinase_table.tsv's column order is: the 9 base columns, then
is_representative, cluster_number, specI, cluster_size (all from 2.3.2),
then the two/three new columns below (from this script).

  col 14  singleton            : "yes" if this row's mge_id never appears
                                  as a token on ANY line of the duplicates
                                  file (no duplicate copy was ever found
                                  for it), "no" otherwise.
  col 15  duplicate_group_number : every line of the duplicates file is
                                  one duplicate group, numbered 1..N in
                                  file order. If singleton is "yes", this
                                  is "NA".
  col 16  duplicate_group_size   : how many DISTINCT mge_ids in
                                  b_mge_recombinase_table.tsv (the actual
                                  dataset, restricted the same way as
                                  2.3.2's cluster_size) share this row's
                                  duplicate_group_number -- counts unique
                                  MGEs, not raw rows (a multi-recombinase
                                  MGE contributes >1 row for the same
                                  group). "NA" if singleton is "yes".

Run standalone:
    python append_duplicate_info.py --config path/to/config.yaml
    python append_duplicate_info.py --table b_table.tsv --duplicates dups.txt --out c_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, parse_duplicate_groups_file,
    resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.4.2"

COLUMNS = [
    "mge_id", "mge_start", "mge_end", "n_genes", "mgeR",
    "recombinase_id", "recombinase_start", "recombinase_end", "recombinase_type",
    "is_representative", "cluster_number", "specI", "cluster_size",
]


def find_duplicates_file(base_dir: Path) -> Path | None:
    vclust_dir = base_dir / "vclust_deduplication"
    if not vclust_dir.exists():
        return None
    candidates = sorted(vclust_dir.glob("*.duplicates.txt"))
    return candidates[0] if candidates else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--table", default=None, help="Override BASE/files/b_mge_recombinase_table.tsv")
    ap.add_argument("--duplicates", default=None, help="Override the vclust *.duplicates.txt path")
    ap.add_argument("--out", default=None, help="Override BASE/files/c_mge_recombinase_table.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_4"][SECTION_KEY]

    log = get_logger("2.4.2_append_duplicate_info", paths["logs_dir"])

    table_in = Path(args.table) if args.table else files_dir / "b_mge_recombinase_table.tsv"
    out_path = Path(args.out) if args.out else files_dir / "c_mge_recombinase_table.tsv"
    duplicates_file = Path(args.duplicates) if args.duplicates else find_duplicates_file(paths["base_dir"])

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Append duplicate-group info")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_4_2_append_duplicate_info"):
        return

    if not table_in.exists():
        log.error(f"{table_in} not found -- run 2.3.2_append_clustering_info.py (Section 2.3.2) first.")
        raise SystemExit(1)
    if duplicates_file is None or not duplicates_file.exists():
        log.error(
            f"No *.duplicates.txt found under {paths['base_dir'] / 'vclust_deduplication'} "
            f"-- run 2.4.1_vclust_dedup.sh (Section 2.4.1) first, or pass --duplicates explicitly."
        )
        raise SystemExit(1)

    t0 = time.time()

    log.info(f"Loading duplicate groups: {duplicates_file}")
    groups = parse_duplicate_groups_file(duplicates_file)
    log.info(f"  {len(groups):,} duplicate groups")

    # mge_id -> 1-based group number (file order)
    mge_id_to_group_num: dict = {}
    for i, group in enumerate(groups, start=1):
        for mge_id in group:
            mge_id_to_group_num[mge_id] = i

    log.info(f"Loading input table: {table_in}")
    with open(table_in, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    log.info(f"  Rows loaded: {len(rows):,}")

    n_singleton = n_non_singleton = 0
    row_group_nums = []
    group_num_to_mge_ids: dict = defaultdict(set)
    for row in rows:
        mge_id = row.get("mge_id", "")
        group_num = mge_id_to_group_num.get(mge_id)
        row_group_nums.append(group_num)
        if group_num is not None:
            group_num_to_mge_ids[group_num].add(mge_id)

    enriched_rows = []
    for row, group_num in zip(rows, row_group_nums):
        if group_num is None:
            singleton, dup_group, dup_group_size = "yes", "NA", "NA"
            n_singleton += 1
        else:
            singleton, dup_group = "no", str(group_num)
            # Distinct mge_ids sharing this group_num, restricted to what's
            # actually IN b_mge_recombinase_table.tsv -- NOT the raw
            # duplicates file's group size, which may include mge_ids never
            # even considered here (e.g. non-representative genomes). Counts
            # unique MGEs, not raw rows, since one MGE with >1 recombinase
            # contributes >1 row here for the SAME group.
            dup_group_size = str(len(group_num_to_mge_ids[group_num]))
            n_non_singleton += 1

        new_row = {k: row.get(k, "") for k in COLUMNS}
        new_row["singleton"] = singleton
        new_row["duplicate_group_number"] = dup_group
        new_row["duplicate_group_size"] = dup_group_size
        enriched_rows.append(new_row)

    out_columns = COLUMNS + ["singleton", "duplicate_group_number", "duplicate_group_size"]
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(enriched_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Rows written        : {len(enriched_rows):,}")
    log.info(f"  singleton = yes   : {n_singleton:,}")
    log.info(f"  singleton = no    : {n_non_singleton:,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()