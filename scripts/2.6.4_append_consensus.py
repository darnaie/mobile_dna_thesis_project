#!/usr/bin/env python3
"""
2.6.4_append_consensus.py  —  Section 2.6.4
================================================
input : BASE/files/d_mge_recombinase_table.tsv (Section 2.5.2),
        BASE/files/cluster_consensus_results.tsv (Section 2.6.2)
output: BASE/files/e_mge_recombinase_table.tsv

Every row of d_mge_recombinase_table.tsv is per-MGE, per-side (left/right)
-- but consensus is built at CLUSTER level (2.6.2), not per-MGE. So every
MGE row belonging to the same (rec_type, cluster_number, side) gets the
SAME consensus_method / consensus_sequence appended -- i.e. this is a
many-to-one lookup (many MGE rows -> one cluster's consensus result),
not a 1:1 join.

If consensus_source was "both" in Section 2.6.2 (both stacking AND
alignment results computed for the same cluster/side), the STACKING
result is preferred here (matches build_consensus's own stated default
method) -- both are never ambiguous to pick between silently otherwise,
since a table needs exactly one value per column per row.

BUG FIX (see chat): the lookup this script joins consensus results
through is keyed on (rec_type, cluster_number, side) ONLY -- it has no
concept of duplicate_group_number or singleton status. Since a single
cluster can (and often does) contain BOTH duplicate_group_* members AND
singleton MGEs side by side, EVERY MGE row in that cluster was
previously getting stamped with the same cluster-level consensus,
including singleton rows -- even though a singleton was never part of
building it (2.6.1 already excludes singleton/ folders from cluster-level
pooling; the consensus itself was always built correctly, only the
join-back here wasn't respecting who actually contributed to it). Fixed:
every row with singleton == "yes" now gets consensus_method,
consensus_start, consensus_length, and consensus_sequence forced to "NA"
unconditionally, regardless of what its cluster's lookup entry says.

Appends four columns (was two -- consensus_start and consensus_length
are NEW, see chat):

  col 19  consensus_method    : e.g. "stacking_A_run2" or "alignment_C"/
                                 "alignment_D_run2"/"alignment_E_dyn_run2"
                                 (see 2.6.2's own docstring for what each
                                 alignment-mode label means) -- exactly
                                 the `method` value 2.6.2 recorded for
                                 that cluster/side; "NA" if this row is a
                                 singleton (see BUG FIX above) or if no
                                 consensus result exists for its cluster
                                 at all.
  col 20  consensus_start     : where the reported consensus starts, in
                                 the original (untrimmed) alignment/
                                 stack's own column numbering -- straight
                                 passthrough of 2.6.2's own consensus_start
                                 column (see its docstring); "NA" under
                                 the same conditions as consensus_method.
  col 21  consensus_length    : how many nt long the consensus is --
                                 passthrough of 2.6.2's consensus_length
                                 (this existed in 2.6.2's own output
                                 already, but wasn't being carried through
                                 to this table until now); "NA" under the
                                 same conditions.
  col 22  consensus_sequence  : the actual consensus text; "NA" under the
                                 same conditions.

Note this reads FROM d_mge_recombinase_table.tsv but writes to
e_mge_recombinase_table.tsv, NOT overwriting d -- consistent with the
established a -> b -> c -> d -> e table-lineage convention elsewhere in
this pipeline (2.2.4 through 2.5.2), even though "append the consensus
sequence to ... d_mge_recombinase_table.tsv" could be read literally as
in-place -- flag if you actually wanted d overwritten instead.

Run standalone:
    python 2.6.4_append_consensus.py --config path/to/config.yaml
    python 2.6.4_append_consensus.py --table d_table.tsv --results cluster_consensus_results.tsv --out e_table.tsv
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.6.4"

COLUMNS = [
    "mge_id", "mge_start", "mge_end", "n_genes", "mgeR",
    "recombinase_id", "recombinase_start", "recombinase_end", "recombinase_type",
    "is_representative", "cluster_number", "specI", "cluster_size",
    "singleton", "duplicate_group_number", "duplicate_group_size",
    "side", "flank_sequence",
]


def load_consensus_lookup(results_path: Path, log) -> dict:
    """(rec_type, cluster_number, side) -> {source: (method, start, length, consensus_sequence)}.
    "stacking" preferred at lookup time if multiple sources exist for the
    same key -- see module docstring."""
    lookup: dict = {}
    with open(results_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            key = (row["rec_type"], row["cluster_number"], row["side"])
            source = row["consensus_source"]
            lookup.setdefault(key, {})[source] = (
                row["method"], row.get("consensus_start", "NA"),
                row.get("consensus_length", "NA"), row["consensus_sequence"],
            )

    n_both = sum(1 for v in lookup.values() if len(v) > 1)
    if n_both:
        log.info(f"{n_both:,} cluster/side keys have BOTH stacking and alignment results -- "
                 f"preferring stacking for each.")
    return lookup


def pick_consensus(lookup: dict, rec_type: str, cluster_number: str, side: str) -> tuple[str, str, str, str]:
    """Returns (method, consensus_start, consensus_length, consensus_sequence),
    all "NA" if there's no result for this (rec_type, cluster_number, side)
    at all. Callers MUST ALSO check singleton status separately -- see
    main(): this lookup is cluster-level only and has no idea which MGEs
    within that cluster actually contributed to the consensus (only
    duplicate-group members do -- see module docstring's BUG FIX note)."""
    entry = lookup.get((rec_type, cluster_number, side))
    if not entry:
        return "NA", "NA", "NA", "NA"
    if "stacking" in entry:
        return entry["stacking"]
    # only alignment available
    return next(iter(entry.values()))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--table", default=None, help="Override BASE/files/d_mge_recombinase_table.tsv")
    ap.add_argument("--results", default=None, help="Override BASE/files/cluster_consensus_results.tsv")
    ap.add_argument("--out", default=None, help="Override BASE/files/e_mge_recombinase_table.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_6"][SECTION_KEY]

    log = get_logger("2.6.4_append_consensus", paths["logs_dir"])

    table_in = Path(args.table) if args.table else files_dir / "d_mge_recombinase_table.tsv"
    results_path = Path(args.results) if args.results else files_dir / "cluster_consensus_results.tsv"
    out_path = Path(args.out) if args.out else files_dir / "e_mge_recombinase_table.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Append cluster consensus onto the per-MGE table")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_6_4_append_consensus"):
        return

    if not table_in.exists():
        log.error(f"{table_in} not found -- run 2.5.2_append_flank_sequences.py first.")
        raise SystemExit(1)
    if not results_path.exists():
        log.error(f"{results_path} not found -- run 2.6.3_build_consensus.py first.")
        raise SystemExit(1)

    t0 = time.time()

    log.info(f"Loading consensus results: {results_path}")
    lookup = load_consensus_lookup(results_path, log)
    log.info(f"  {len(lookup):,} (rec_type, cluster, side) keys have a consensus result")

    log.info(f"Loading input table: {table_in}")
    with open(table_in, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        rows = list(reader)
    log.info(f"  Rows loaded: {len(rows):,}")

    n_found = n_missing = n_singleton_suppressed = 0
    out_rows = []
    for row in rows:
        is_singleton = row.get("singleton", "").strip().lower() == "yes"

        if is_singleton:
            # BUG FIX (see chat / module docstring): the lookup below is
            # keyed on (rec_type, cluster_number, side) only -- it has no
            # idea whether THIS row's specific MGE was one of the
            # duplicate-group members that actually went into building
            # the cluster's consensus, or a singleton that just happens
            # to share the same cluster. A singleton was never part of
            # any alignment/stacking this consensus came from (2.6.1
            # excludes singleton/ folders from cluster-level pooling), so
            # it never gets a value here, regardless of what the cluster
            # itself has.
            method, cons_start, cons_length, seq = "NA", "NA", "NA", "NA"
            n_singleton_suppressed += 1
        else:
            method, cons_start, cons_length, seq = pick_consensus(
                lookup, row.get("recombinase_type", ""), row.get("cluster_number", "NA"), row.get("side", "")
            )

        if method != "NA":
            n_found += 1
        else:
            n_missing += 1

        new_row = {k: row.get(k, "") for k in COLUMNS}
        new_row["consensus_method"] = method
        new_row["consensus_start"] = cons_start
        new_row["consensus_length"] = cons_length
        new_row["consensus_sequence"] = seq
        out_rows.append(new_row)

    log.info(f"Writing output table: {out_path}")
    out_columns = COLUMNS + ["consensus_method", "consensus_start", "consensus_length", "consensus_sequence"]
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Rows written                  : {len(out_rows):,}")
    log.info(f"  consensus found             : {n_found:,}")
    log.info(f"  consensus NOT found (NA)    : {n_missing:,}")
    log.info(f"    of which singleton rows (correctly forced to NA): {n_singleton_suppressed:,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()