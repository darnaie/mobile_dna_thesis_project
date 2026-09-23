#!/usr/bin/env python3
"""
build_recombinase_table.py  —  Section 2.2.4
================================================
input : BASE/files/non_nested_is_tn.gff3 (2.1.1, MGE-island-level data:
        start/end/n_genes/mgeR), BASE/files/non_nested_is_tn_recombinase.gff3
        (2.1.2, recombinase gene-level data), BASE/files/
        is_tn_recombinase_non_nested.fasta (2.1.5, already-extracted
        sequences for ALL is_tn recombinases -- REUSED here, not
        re-extracted from raw genome FASTAs), BASE/genome_seqs_repr/
        (2.2.3, restricts everything below to representative/deduplicated
        genomes only)
output: BASE/files/a_mge_recombinase_table.tsv,
        BASE/files/is_tn_recombinase_non_nested_repr.fasta

Column layout of the TSV:
  1  mge_id
  2  mge_start
  3  mge_end
  4  n_genes
  5  mgeR
  6  recombinase_id
  7  recombinase_start
  8  recombinase_end
  9  recombinase_type

Rules (adapted from the original phase4b table-builder):
  * Only MGEs on a REPRESENTATIVE genome (present in BASE/genome_seqs_repr/)
    are considered at all.
  * If an MGE has >1 recombinase, ALL are kept, each its own row.
  * Sanity check: recombinase_start >= mge_start AND recombinase_end <= mge_end.
    Rows failing this are excluded (counted in the summary).
  * Only MGEs with >=1 recombinase passing sanity AND with a sequence
    available in the 2.1.5 fasta appear in the table and output fasta.

NOTE: the output fasta is named is_tn_recombinase_non_nested_repr.fasta
(the "_repr" suffix marks it as the representative-genome-only subset) --
it does NOT overwrite 2.1.5's is_tn_recombinase_non_nested.fasta (which
covers ALL is_tn recombinases, not just those on representative genomes).
The two files are deliberately distinct: 2.1.5's is the full set, this
one is what everything from Section 2.3 onward should actually use.

Run standalone:
    python 2.2.4_build_recombinase_table.py --config path/to/config.yaml
    python 2.2.4_build_recombinase_table.py --mge-gff ... --recomb-gff ... --recomb-fasta ... --repr-dir ... --out-dir ./out
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, parse_gff_attributes, read_fasta,
    taxid_sample_from_mge_id, resolve_placeholder_multi, maybe_submit_and_exit,
)

SECTION_KEY = "2.2.4"


def load_mge_data(gff_path: Path) -> dict:
    mge_data = {}
    with open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "mobile_genetic_element":
                continue
            attrs = parse_gff_attributes(cols[8])
            mge_id = attrs.get("ID", "")
            if not mge_id:
                continue
            mge_data[mge_id] = {
                "start": int(cols[3]), "end": int(cols[4]),
                "n_genes": attrs.get("n_genes", ""), "mgeR": attrs.get("mgeR", ""),
            }
    return mge_data


def load_recomb_map(gff_path: Path) -> dict:
    recomb_map = defaultdict(list)
    with open(gff_path) as fh:
        for raw in fh:
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue
            attrs = parse_gff_attributes(cols[8])
            parent = attrs.get("Parent", "")
            recomb_id = attrs.get("ID", "")
            recomb_t = attrs.get("recombinase", "")
            if not parent or not recomb_id:
                continue
            try:
                r_start, r_end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            recomb_map[parent].append((recomb_id, r_start, r_end, recomb_t))
    return recomb_map


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--mge-gff", default=None, help="Override BASE/files/non_nested_is_tn.gff3")
    ap.add_argument("--recomb-gff", default=None, help="Override BASE/files/non_nested_is_tn_recombinase.gff3")
    ap.add_argument("--recomb-fasta", default=None, help="Override BASE/files/is_tn_recombinase_non_nested.fasta (input, from 2.1.5)")
    ap.add_argument("--repr-dir", default=None, help="Override BASE/genome_seqs_repr")
    ap.add_argument("--out-dir", default=None, help="Override BASE/files")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    out_dir = Path(args.out_dir) if args.out_dir else files_dir
    mge_gff = Path(args.mge_gff) if args.mge_gff else files_dir / "non_nested_is_tn.gff3"
    recomb_gff = Path(args.recomb_gff) if args.recomb_gff else files_dir / "non_nested_is_tn_recombinase.gff3"
    recomb_fasta_in = Path(args.recomb_fasta) if args.recomb_fasta else files_dir / "is_tn_recombinase_non_nested.fasta"
    repr_dir = Path(args.repr_dir) if args.repr_dir else paths["genome_seqs_repr"]
    section_cfg = cfg["section_2_2"][SECTION_KEY]

    log = get_logger("2.2.4_build_recombinase_table", paths["logs_dir"])

    table_out = out_dir / "a_mge_recombinase_table.tsv"
    # Renamed vs. 2.1.5's own is_tn_recombinase_non_nested.fasta -- 2.1.5
    # extracts sequences for ALL is_tn recombinases; this file is the
    # narrower, representative-genome-only subset, so it now gets its own
    # name instead of silently overwriting 2.1.5's output.
    fasta_out = out_dir / "is_tn_recombinase_non_nested_repr.fasta"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Build MGE/recombinase table (representative genomes only)")
    log.info("=" * 70)

    if resolve_placeholder_multi(section_cfg, {
        "path_to_output_a_rec_table": table_out,
        "path_to_output_is_tn_rec_fasta": fasta_out,
    }, log):
        log.info(f"Section {SECTION_KEY} satisfied entirely via PLACEHOLDER (both outputs staged).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_2_4_build_recombinase_table"):
        return

    for p, label in (
        (mge_gff, "MGE island GFF3 (Section 2.1.1)"),
        (recomb_gff, "recombinase GFF3 (Section 2.1.2)"),
        (recomb_fasta_in, "recombinase FASTA (Section 2.1.5)"),
        (repr_dir, "representative genomes (Section 2.2.3)"),
    ):
        if not p.exists():
            log.error(f"{p} not found -- {label} must exist first.")
            raise SystemExit(1)

    t0 = time.time()

    log.info(f"[1] Loading MGE island data from {mge_gff}")
    mge_data = load_mge_data(mge_gff)
    log.info(f"    Loaded {len(mge_data):,} MGEs")

    log.info(f"[2] Loading recombinase gene rows from {recomb_gff}")
    recomb_map = load_recomb_map(recomb_gff)
    log.info(f"    Recombinase parents found: {len(recomb_map):,}")

    log.info(f"[3] Loading recombinase sequences from {recomb_fasta_in}")
    recomb_seqs = read_fasta(recomb_fasta_in)
    log.info(f"    Loaded {len(recomb_seqs):,} sequences")

    log.info(f"[4] Loading representative genome set from {repr_dir}")
    repr_taxid_samples = {
        f.stem for f in repr_dir.iterdir() if f.is_file() and f.suffix.lower() in (".fasta", ".fna")
    }
    log.info(f"    {len(repr_taxid_samples):,} representative genomes")

    log.info("[5] Building table ...")
    stats = Counter()
    table_rows = []
    seen_fasta = set()

    for mge_id, mge in mge_data.items():
        taxid_sample = taxid_sample_from_mge_id(mge_id)
        if taxid_sample not in repr_taxid_samples:
            stats["not_representative"] += 1
            continue

        candidates = recomb_map.get(mge_id, [])
        if not candidates:
            stats["mge_no_match"] += 1
            continue
        stats["mge_matched_1" if len(candidates) == 1 else "mge_matched_multi"] += 1

        mge_contributed = False
        for recomb_id, r_start, r_end, recomb_t in candidates:
            if not (mge["start"] <= r_start and mge["end"] >= r_end):
                stats["rows_failed_sanity"] += 1
                continue
            if recomb_id not in recomb_seqs:
                stats["rows_no_seq"] += 1
                continue

            table_rows.append([
                mge_id, str(mge["start"]), str(mge["end"]), mge["n_genes"], mge["mgeR"],
                recomb_id, str(r_start), str(r_end), recomb_t,
            ])
            stats["rows_total"] += 1
            seen_fasta.add(recomb_id)
            mge_contributed = True

        if mge_contributed:
            stats["mge_written"] += 1

    log.info(f"    Table rows written: {stats['rows_total']:,}")

    log.info(f"[6] Writing table to {table_out}")
    header = "\t".join([
        "mge_id", "mge_start", "mge_end", "n_genes", "mgeR",
        "recombinase_id", "recombinase_start", "recombinase_end", "recombinase_type",
    ])
    with open(table_out, "w") as tf:
        tf.write(header + "\n")
        for row in table_rows:
            tf.write("\t".join(row) + "\n")

    log.info(f"Writing FASTA to {fasta_out}")
    written = set()
    with open(fasta_out, "w") as ff:
        for row in table_rows:
            rid = row[5]
            if rid not in written:
                ff.write(f">{rid}\n{recomb_seqs[rid]}\n")
                written.add(rid)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Total MGEs (non-nested is_tn)         : {len(mge_data):,}")
    log.info(f"  not on a representative genome      : {stats['not_representative']:,}")
    log.info(f"  no recombinase match                : {stats['mge_no_match']:,}")
    log.info(f"  matched, exactly 1 recombinase       : {stats['mge_matched_1']:,}")
    log.info(f"  matched, >1 recombinase (all kept)   : {stats['mge_matched_multi']:,}")
    log.info(f"MGEs with >=1 valid row in table        : {stats['mge_written']:,}")
    log.info(f"Table rows total                        : {stats['rows_total']:,}")
    log.info(f"  excluded (failed sanity check)        : {stats['rows_failed_sanity']:,}")
    log.info(f"  excluded (no sequence in 2.1.5 fasta) : {stats['rows_no_seq']:,}")
    log.info(f"Unique recombinase sequences written    : {len(written):,}")
    log.info(f"-> {table_out}")
    log.info(f"-> {fasta_out}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()