#!/usr/bin/env python3
"""
append_clustering_info.py  —  Section 2.3.2
================================================
input : BASE/files/a_mge_recombinase_table.tsv (Section 2.2.4),
        BASE/files/recombinase_cluster_<ani>/recombinase_cluster_cluster.tsv
        (Section 2.3.1), BASE/files/f13_dataset.txt (Section 2.1.3)
output: BASE/files/b_mge_recombinase_table.tsv

Appends three columns to every row of a_mge_recombinase_table.tsv:
  col 10  is_representative  : "+|{rep_id}" (this row's recombinase IS a
                                cluster representative), "-|{rep_id}"
                                (it's a member of rep_id's cluster), or
                                "NA|NA" (no cluster match found at all)
  col 11  cluster_number     : integer per representative, assigned
                                *within each recombinase_type* (1-based,
                                ordered by first appearance in the table);
                                "NA" if col 10 is "NA|NA"
  col 12  specI               : looked up via the mge_id's GCA accession
                                against BASE/files/f13_dataset.txt; "NA"
                                if no taxonomy entry
  col 13  cluster_size        : how many rows of a_mge_recombinase_table.tsv
                                (the representative-genome dataset -- NOT
                                mmseqs' raw cluster membership, which can
                                include recombinases never even in that
                                table) share this row's cluster_number
                                within its recombinase_type; "NA" if col 10
                                is "NA|NA" (no cluster to size)

Adapted from the original Phase6.1_table_rearrangement.py, with two
changes to match how this repo's tables actually look:
  * a_mge_recombinase_table.tsv HAS a header row (unlike the old script's
    input) -- read via csv.DictReader instead of raw positional columns.
  * specI is looked up directly via GCA (utils.gca_from_mge_id +
    utils.load_taxonomy_map), not by parsing a "specI:" token out of an
    old-format taxonomy results file.

Run standalone:
    python append_clustering_info.py --config path/to/config.yaml
    python append_clustering_info.py --table a_table.tsv --cluster-tsv cluster.tsv --taxonomy f13.txt --out b_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter, OrderedDict, defaultdict
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, load_taxonomy_map, gca_from_mge_id,
    parse_mmseqs_cluster_tsv, resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.3.2"

COLUMNS = [
    "mge_id", "mge_start", "mge_end", "n_genes", "mgeR",
    "recombinase_id", "recombinase_start", "recombinase_end", "recombinase_type",
]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--table", default=None, help="Override BASE/files/a_mge_recombinase_table.tsv")
    ap.add_argument("--cluster-tsv", default=None,
                     help="Override BASE/files/recombinase_cluster_<ani>/recombinase_cluster_cluster.tsv")
    ap.add_argument("--taxonomy", default=None, help="Override BASE/files/f13_dataset.txt")
    ap.add_argument("--out", default=None, help="Override BASE/files/b_mge_recombinase_table.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_3"][SECTION_KEY]

    log = get_logger("2.3.2_append_clustering_info", paths["logs_dir"])

    table_in = Path(args.table) if args.table else files_dir / "a_mge_recombinase_table.tsv"
    taxonomy_file = Path(args.taxonomy) if args.taxonomy else files_dir / "f13_dataset.txt"
    out_path = Path(args.out) if args.out else files_dir / "b_mge_recombinase_table.tsv"

    if args.cluster_tsv:
        cluster_tsv = Path(args.cluster_tsv)
    else:
        ani_threshold = cfg["section_2_3"]["2.3.1"]["ani_threshold"]
        cluster_tsv = files_dir / f"recombinase_cluster_{ani_threshold}" / "recombinase_cluster_cluster.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Append clustering info + specI")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_3_2_append_clustering_info"):
        return

    for p, label in (
        (table_in, "a_mge_recombinase_table.tsv (Section 2.2.4)"),
        (cluster_tsv, "recombinase cluster.tsv (Section 2.3.1)"),
        (taxonomy_file, "f13 taxonomy table (Section 2.1.3)"),
    ):
        if not p.exists():
            log.error(f"{p} not found -- {label} must exist first.")
            raise SystemExit(1)

    t0 = time.time()

    log.info(f"Loading taxonomy: {taxonomy_file}")
    taxonomy_map = load_taxonomy_map(taxonomy_file)
    log.info(f"  {len(taxonomy_map):,} GCA -> specI entries")

    log.info(f"Loading cluster file: {cluster_tsv}")
    representatives, member_to_rep = parse_mmseqs_cluster_tsv(cluster_tsv)
    log.info(f"  Unique cluster representatives : {len(representatives):,}")
    log.info(f"  Total member-to-rep mappings   : {len(member_to_rep):,}")

    log.info(f"Loading input table: {table_in}")
    with open(table_in, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    log.info(f"  Rows loaded: {len(rows):,}")

    log.info("Grouping by recombinase type ...")
    type_to_rows: dict = defaultdict(list)
    for row in rows:
        type_to_rows[row.get("recombinase_type", "UNKNOWN")].append(row)

    log.info("Enriching rows ...")
    n_is_rep = n_not_rep = n_no_clust = n_no_tax = 0
    enriched_rows = []

    for rtype in sorted(type_to_rows):
        rrows = type_to_rows[rtype]

        # Assign cluster numbers within this type, ordered by first occurrence
        rep_to_clust_num: "OrderedDict[str, int]" = OrderedDict()
        cluster_counter = 1
        for row in rrows:
            rec_id = row.get("recombinase_id", "")
            if rec_id in representatives:
                rep_id = rec_id
            elif rec_id in member_to_rep:
                rep_id = member_to_rep[rec_id]
            else:
                rep_id = None
            if rep_id is not None and rep_id not in rep_to_clust_num:
                rep_to_clust_num[rep_id] = cluster_counter
                cluster_counter += 1

        # Resolve every row's is_representative/cluster_number/specI first,
        # tallying how many rows land in each cluster_number as we go --
        # this is cluster_size: how many members of a cluster are actually
        # PRESENT IN a_mge_recombinase_table.tsv (representative-genome
        # dataset), not the raw mmseqs cluster's full membership (which
        # can include non-representative-genome recombinases never even
        # considered here).
        row_results = []
        cluster_number_counts: Counter = Counter()
        for row in rrows:
            rec_id = row.get("recombinase_id", "")
            mge_id = row.get("mge_id", "")

            if rec_id in representatives:
                rep_id = rec_id
                col10 = f"+|{rep_id}"
                n_is_rep += 1
            elif rec_id in member_to_rep:
                rep_id = member_to_rep[rec_id]
                col10 = f"-|{rep_id}"
                n_not_rep += 1
            else:
                rep_id = None
                col10 = "NA|NA"
                n_no_clust += 1

            col11 = str(rep_to_clust_num[rep_id]) if rep_id and rep_id in rep_to_clust_num else "NA"
            if col11 != "NA":
                cluster_number_counts[col11] += 1

            gca = gca_from_mge_id(mge_id)
            info = taxonomy_map.get(gca) if gca else None
            if info:
                col12 = info[0]
            else:
                col12 = "NA"
                n_no_tax += 1

            row_results.append((row, col10, col11, col12))

        for row, col10, col11, col12 in row_results:
            new_row = {k: row.get(k, "") for k in COLUMNS}
            new_row["is_representative"] = col10
            new_row["cluster_number"] = col11
            new_row["specI"] = col12
            new_row["cluster_size"] = str(cluster_number_counts[col11]) if col11 != "NA" else "NA"
            enriched_rows.append(new_row)

    log.info(f"Writing output table: {out_path}")
    out_columns = COLUMNS + ["is_representative", "cluster_number", "specI", "cluster_size"]
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(enriched_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Rows written                          : {len(enriched_rows):,}")
    log.info(f"  marked representative (+)           : {n_is_rep:,}")
    log.info(f"  marked non-representative (-)       : {n_not_rep:,}")
    log.info(f"  no cluster match (NA)                : {n_no_clust:,}")
    log.info(f"  no taxonomy entry                    : {n_no_tax:,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()