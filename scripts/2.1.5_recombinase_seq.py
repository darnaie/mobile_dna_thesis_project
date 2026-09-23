#!/usr/bin/env python3
"""
recombinase_seq.py  —  Section 2.1.5
=======================================
input : BASE/files/non_nested_is_tn_recombinase.gff3 (Section 2.1.2),
        BASE/genome_seqs/ (Section 2.1.4)
output: BASE/files/is_tn_recombinase_non_nested.fasta

Extracts each non-nested is_tn recombinase gene's raw nucleotide sequence
from its genome FASTA. Simplified from the earlier multi-variant version
of this script -- there's only ever one output now (non-nested), matching
the rest of Section 2.1.

Design notes carried over as-is from the original version of this script
(please double check both of these against your intended methods text):

  * NO reverse-complementing is applied, regardless of strand. The strand
    is preserved in the output header for reference only.
  * Genome files in genome_seqs/ are matched by the first two
    dot-separated tokens of their filename, lower-cased (e.g.
    '573.samn03280285'), against the same key derived from the GFF3
    contig ID (column 1).

Run standalone:
    python 2.1.5_recombinase_seq.py --config path/to/config.yaml
    python 2.1.5_recombinase_seq.py --gff /path/to/gff3 --genome-dir /path/to/genomes --out /path/to/out.fasta
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, parse_gff_attributes, read_fasta,
    taxid_sample_from_contig_id, resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.1.5"
MAX_GENOME_CACHE = 50


class GenomeCache:
    """Small LRU-style cache so we don't re-parse the same genome FASTA
    over and over when its genes are scattered through the GFF3."""

    def __init__(self, max_size: int = MAX_GENOME_CACHE):
        self.max_size = max_size
        self._cache: dict[str, dict[str, str]] = {}
        self._order: list[str] = []

    def get(self, path: Path) -> dict[str, str]:
        key = str(path)
        if key not in self._cache:
            if len(self._order) >= self.max_size:
                del self._cache[self._order.pop(0)]
            self._cache[key] = read_fasta(path)
            self._order.append(key)
        return self._cache[key]


def build_genome_index(genome_dir: Path) -> dict[str, Path]:
    index = {}
    for f in genome_dir.iterdir():
        if f.is_file() and f.suffix.lower() in (".fasta", ".fna"):
            index[taxid_sample_from_contig_id(f.name)] = f
    return index


def extract_sequences(gff_path: Path, out_path: Path, genome_index: dict) -> Counter:
    stats = Counter()
    cache = GenomeCache()
    with open(out_path, "w") as out_fh, open(gff_path) as gff:
        for lineno, raw in enumerate(gff, start=1):
            if raw.startswith("#") or not raw.strip():
                continue
            cols = raw.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "gene":
                continue

            contig_id = cols[0]
            strand = cols[6]
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                stats["bad_coords"] += 1
                continue

            attrs = parse_gff_attributes(cols[8])
            gene_id = attrs.get("ID", f"gene_{lineno}")

            genome_key = taxid_sample_from_contig_id(contig_id)
            genome_path = genome_index.get(genome_key)
            if genome_path is None:
                stats["no_genome"] += 1
                continue

            genome = cache.get(genome_path)
            if contig_id not in genome:
                stats["no_contig"] += 1
                continue

            seq = genome[contig_id][start - 1: end]  # no revcomp -- see docstring
            out_fh.write(f">{gene_id} {contig_id}:{start}-{end}({strand})\n{seq}\n")
            stats["extracted"] += 1

    return stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--gff", default=None, help="Override BASE/files/non_nested_is_tn_recombinase.gff3")
    ap.add_argument("--genome-dir", default=None, help="Override BASE/genome_seqs")
    ap.add_argument("--out", default=None, help="Override BASE/files/is_tn_recombinase_non_nested.fasta")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    genome_dir = Path(args.genome_dir) if args.genome_dir else paths["genome_dir"]
    section_cfg = cfg["section_2_1"][SECTION_KEY]

    log = get_logger("2.1.5_recombinase_seq", paths["logs_dir"])

    gff_path = Path(args.gff) if args.gff else files_dir / "non_nested_is_tn_recombinase.gff3"
    out_path = Path(args.out) if args.out else files_dir / "is_tn_recombinase_non_nested.fasta"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Recombinase sequence extraction (non-nested is_tn only)")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_1_5_recombinase_seq"):
        return

    log.info(f"Source GFF3 : {gff_path}")
    log.info(f"Genome dir  : {genome_dir}")
    log.info("NOTE: no reverse-complementing is applied to minus-strand genes "
             "(strand kept in header only) — confirm this is intended.")

    if not gff_path.exists():
        log.error(f"{gff_path} not found -- run recombinase_gff.py (Section 2.1.2) first.")
        raise SystemExit(1)
    if not genome_dir.exists() or not any(genome_dir.iterdir()):
        log.error(
            f"{genome_dir} does not exist or is empty -- run download_genomes.py "
            f"(Section 2.1.4) first, or point --genome-dir at a populated directory."
        )
        raise SystemExit(1)

    t0 = time.time()

    log.info(f"Indexing genomes in {genome_dir} ...")
    genome_index = build_genome_index(genome_dir)
    log.info(f"Indexed {len(genome_index):,} genome FASTA files")

    stats = extract_sequences(gff_path, out_path, genome_index)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Extracted            : {stats['extracted']:,}")
    log.info(f"  skipped (no_genome) : {stats['no_genome']:,}")
    log.info(f"  skipped (no_contig) : {stats['no_contig']:,}")
    log.info(f"  skipped (bad_coords): {stats['bad_coords']:,}")
    log.info(f"Genomes indexed      : {len(genome_index):,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()