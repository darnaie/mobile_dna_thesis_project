#!/usr/bin/env python3
"""
divide_genomes_by_speci.py  —  Section 2.2.1
================================================
input : BASE/genome_seqs/ (Section 2.1.4), BASE/files/f13_dataset.txt
        (Section 2.1.3), BASE/files/master_mges.gff3 (Section 2.1.1)
output: BASE/genome_seqs_div_by_speci/specI_vX_NNNNN/<genome>.fasta

Organizes genome FASTA files into one subfolder per specI, e.g.:

    genome_seqs_div_by_speci/
      specI_v4_00001/
        1001989.SAMN02469619.fasta
        1005433.SAMN01112588.fasta
      specI_v4_00002/
        ...

Genome files are named by taxid.sampleid (e.g. '573.SAMN03280285.fasta'),
but the f13 taxonomy table is keyed on GCA accession -- there's no direct
taxid.sampleid -> specI mapping anywhere on its own. This script bridges
the two by scanning master_mges.gff3's MGE IDs (which embed BOTH the GCA
and the taxid.sampleid together, e.g.
'MGE_GCA_001031825.1_573.SAMN03280285.KQ088489:...') to build a
taxid.sampleid -> GCA lookup, then GCA -> specI via the f13 table.

Files are SYMLINKED into their specI subfolder, not copied (genome
directories are far too large to duplicate on disk just to reorganize them).

Run standalone:
    python 2.2.1_divide_genomes_by_speci.py --config path/to/config.yaml
    python 2.2.1_divide_genomes_by_speci.py --genome-dir ./genomes --taxonomy f13.txt --gca-source master_mges.gff3 --out-dir ./div_by_speci
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, load_taxonomy_map,
    build_taxid_sample_to_gca_map, resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.2.1"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--genome-dir", default=None, help="Override BASE/genome_seqs")
    ap.add_argument("--taxonomy", default=None, help="Override BASE/files/f13_dataset.txt")
    ap.add_argument("--gca-source", default=None, help="Override BASE/files/master_mges.gff3")
    ap.add_argument("--out-dir", default=None, help="Override BASE/genome_seqs_div_by_speci")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    genome_dir = Path(args.genome_dir) if args.genome_dir else paths["genome_dir"]
    out_dir = Path(args.out_dir) if args.out_dir else paths["genome_seqs_div_by_speci"]
    taxonomy_file = Path(args.taxonomy) if args.taxonomy else files_dir / "f13_dataset.txt"
    gca_source = Path(args.gca_source) if args.gca_source else files_dir / "master_mges.gff3"
    section_cfg = cfg["section_2_2"][SECTION_KEY]

    log = get_logger("2.2.1_divide_genomes_by_speci", paths["logs_dir"])

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Divide genomes by specI")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_dir, log):
        n_species = len(list(out_dir.iterdir())) if out_dir.exists() else 0
        log.info(f"{out_dir} now has {n_species:,} specI subfolders (via placeholder)")
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_2_1_divide_by_speci"):
        return

    for p, label in (
        (genome_dir, "genome_dir (Section 2.1.4)"),
        (taxonomy_file, "taxonomy file (Section 2.1.3)"),
        (gca_source, "GCA source GFF3 (Section 2.1.1)"),
    ):
        if not p.exists():
            log.error(f"{p} not found -- {label} must exist first.")
            raise SystemExit(1)

    t0 = time.time()

    log.info(f"Loading taxonomy map from {taxonomy_file} ...")
    taxonomy_map = load_taxonomy_map(taxonomy_file)
    log.info(f"  {len(taxonomy_map):,} GCA -> specI entries")

    log.info(f"Building taxid.sampleid -> GCA map from {gca_source} ...")
    taxid_to_gca = build_taxid_sample_to_gca_map(gca_source)
    log.info(f"  {len(taxid_to_gca):,} taxid.sampleid -> GCA entries")

    genome_files = [
        f for f in genome_dir.iterdir()
        if f.is_file() and f.suffix.lower() in (".fasta", ".fna")
    ]
    log.info(f"Found {len(genome_files):,} genome files in {genome_dir}")

    out_dir.mkdir(parents=True, exist_ok=True)
    stats = Counter()
    species_counts: Counter = Counter()

    for f in genome_files:
        taxid_sample = f.stem  # filenames are named exactly by taxid.sampleid
        gca = taxid_to_gca.get(taxid_sample)
        if gca is None:
            stats["no_gca_mapping"] += 1
            continue

        info = taxonomy_map.get(gca)
        if info is None:
            stats["no_specI_mapping"] += 1
            continue

        specI, _label = info
        species_dir = out_dir / specI
        species_dir.mkdir(parents=True, exist_ok=True)
        dest = species_dir / f.name
        if not dest.exists() and not dest.is_symlink():
            dest.symlink_to(f.resolve())
        stats["assigned"] += 1
        species_counts[specI] += 1

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Genome files scanned      : {len(genome_files):,}")
    log.info(f"Assigned to a specI       : {stats['assigned']:,}")
    log.info(f"No taxid->GCA mapping     : {stats['no_gca_mapping']:,}")
    log.info(f"No GCA->specI mapping     : {stats['no_specI_mapping']:,}")
    log.info(f"Distinct specI folders    : {len(species_counts):,}")
    if species_counts:
        top = species_counts.most_common(3)
        log.info(f"Largest specI folders (top 3): {top}")
    log.info(f"-> {out_dir}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()