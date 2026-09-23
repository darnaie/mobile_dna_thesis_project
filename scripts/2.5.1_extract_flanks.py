#!/usr/bin/env python3
"""
2.5.1_extract_flanks.py  —  Section 2.5.1
=============================================
input : BASE/files/non_nested_is_tn.gff3 (2.1.1, MGE coordinates),
        BASE/genome_seqs_repr/ (2.2.3, representative genome FASTAs),
        BASE/files/c_mge_recombinase_table.tsv (2.4.2 -- NOT explicitly
        named when this section was requested, but structurally required:
        it's the only place each mge_id's recombinase_type/cluster_number/
        duplicate_group_number/singleton live, and that's what decides
        WHERE on disk this script files each MGE's flank pair)
output: BASE/flanks_split_acc_to_duplicates/<rec_type>/cluster_<N>/
        <duplicate_group_M | singleton>/<mge_id>_left.fasta, _right.fasta
        BASE/is_tns_seqs/<mge_id>.fasta  (only if extract_is_tn_seqs:
        true in config.yaml -- see below)

BUG FIX (see chat -- LEFT flank off-by-one): re-derived left_end_excl
from first principles and checked it against the simplest possible case
-- mge_overlap=0, which should mean the left flank touches ZERO bases of
the MGE. It didn't: left_end_excl = left_anchor + 1 always included
exactly one MGE base beyond whatever mge_overlap actually specified,
regardless of its value (reproduced directly: with mge_overlap=0, the
old formula still reached one base into the MGE). RIGHT flank extraction
was checked the same way and has no such issue -- right_anchor already
lands exactly on the MGE's own last mge_overlap bases, confirmed against
concrete coordinates. Root cause: left_anchor := (start-1)+mge_overlap
is ALREADY the correct exclusive upper bound (0-based) for a slice
capturing exactly the MGE's first mge_overlap bases -- the extra "+ 1"
on top of it was the bug. Fixed by using left_anchor directly as
left_end_excl (removing the redundant, incorrect +1); left_anchor itself
is no longer a separate variable, since nothing else used it.

Net effect of the bug: every LEFT flank ever extracted by the old
formula was composed of (mge_overlap + 1) nt of the MGE's own start and
(flank_size - mge_overlap - 1) nt of upstream genomic sequence, instead
of the intended mge_overlap / (flank_size - mge_overlap) split -- total
flank length was still always exactly flank_size (unaffected), only the
INTERNAL composition was off by one base. RIGHT flanks were never
affected. If you've already extracted flanks with the old version, this
means every LEFT flank shifted by exactly 1 nt versus what mge_overlap
actually specified -- re-running this script is the only way to pick up
the fix; nothing downstream can correct for it after the fact.

NEW: extract_is_tn_seqs (config.yaml, default false). When true, this
script ALSO writes the full MGE (is_tn) body sequence itself -- not just
its flanks -- to BASE/is_tns_seqs/<mge_id>.fasta, one flat file per MGE
(no rec_type/cluster/group subfolders, unlike the flanks tree). It's
extracted from the exact same GFF3 start/end coordinates already used to
anchor the flank windows, for every MGE that successfully gets a flank
pair written (same genome/contig/coordinate validity requirements --
nothing extra to fail on). The SAME +/- orientation correction already
applied to flanks is applied here too, using the SAME per-mge_id sign
already looked up for that MGE's flanks:
  sign == '-' (or singleton/unsigned, i.e. no flip needed either way):
      saved exactly as it reads in the genome file.
  sign == '+' : saved as its reverse complement. Unlike the flanks
      (which also SWAP left/right on '+', since reversing orientation
      turns "left of the element" into "right of the element"), the
      is_tn sequence is a single body, not a left/right pair -- there
      is nothing to swap, so revcomp is the whole transformation here.

Adapted from the original idying_flanks_13.py, with two changes:
  1. NO plotting -- that lived in the separate split_and_plot_all.py,
     which this section deliberately does not reproduce.
  2. NO grouping/merging -- split_and_plot_all.py combined every MGE in a
     duplicate group into ONE shared left.fasta/right.fasta pair. This
     script keeps files fully per-MGE (1 mge_id = 1 left.fasta + 1
     right.fasta, exactly as idying_flanks_13.py originally wrote them)
     -- only the directory LOCATION changes, now filed under
     <rec_type>/cluster_<N>/<duplicate_group | singleton>/ using the
     cluster/group assignment Sections 2.3/2.4 already computed, instead
     of landing flat in one directory.

Flank windows are offset mge_overlap nt INTO the MGE from each boundary
(anchors the flank on the immediately adjacent sequence rather than the
hard boundary itself):
  - left  flank : flank_size nt ending mge_overlap nt past the MGE start
  - right flank : flank_size nt starting mge_overlap nt before the MGE end
(the LEFT case above was mis-implemented by one base until the fix
described above; this is the corrected, intended behaviour.)

Genome matching uses each MGE's own embedded taxid.sampleid (parsed from
its mge_id via utils.taxid_sample_from_mge_id, case-preserved) rather than
re-deriving it from the GFF3 contig column -- genome_seqs_repr's files are
themselves named exactly by taxid.sampleid (see 2.2.1/2.2.3), so this
keeps matching exact-case and consistent with how the rest of the
pipeline already indexes that directory.

ORIENTATION CORRECTION (integrated from the pre-thesis_rep flank-extraction
script): duplicate-group members can be inserted on either genomic
strand, so comparing their flanks meaningfully requires normalizing them
all to the SAME relative orientation first. The sign for each mge_id
comes from BASE/vclust_deduplication/*.duplicates.txt (Section 2.4.1) --
'+' or '-', the orientation of that MGE's match relative to its group's
first/reference member (see utils.load_mge_sign_map). An mge_id that
never appears in the duplicates file at all (a true singleton) defaults
to '-', same as an unsigned token within the file (typically the
reference member itself) -- "no flip needed" in both cases, since there's
nothing to compare against.
  sign == '-' : left_flank/right_flank saved exactly as extracted (this
                was already this script's only behavior before this update).
  sign == '+' : BOTH flanks are reverse-complemented, THEN swapped --
                new_left = revcomp(right_flank), new_right = revcomp(left_flank)
                -- matching the original script's exact transformation.
                Only the post-processing changed here, not the extraction
                math itself: this still uses the same fixed flank_size/
                mge_overlap window as the '-' case, not the original
                script's separate dynamic column-identity-based extension
                logic (that's a materially different algorithm and wasn't
                asked for -- only the orientation-correction step was).

VERIFIED (see chat) this '+' transformation is correct, not just
plausible: revcomp(right_flank) -- right_flank's own shape is
[mge_overlap nt of MGE][flank_size-mge_overlap nt of downstream genomic],
so after reversal its order flips to [revcomp'd genomic][revcomp'd MGE
portion] -- i.e. genomic-first, MGE-portion-last, EXACTLY the same
positional shape as a genuine '-' LEFT flank (boundary-adjacent region
at the END). The mirror argument holds for revcomp(left_flank) landing
in the '+' RIGHT slot. This is what makes every downstream script's
"LEFT flanks have their boundary-adjacent region at the end, RIGHT
flanks at the start" assumption (2.6.1's pooling, 2.6.2/2.7.1's
side-aware scanning) hold true regardless of original strand -- confirmed
by direct derivation with concrete coordinates, not assumed.

Run standalone:
    python 2.5.1_extract_flanks.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

from utils import (
    get_logger, get_paths, load_config, parse_gff_attributes,
    taxid_sample_from_mge_id, flank_output_dir,
    load_mge_sign_map, find_duplicates_file, revcomp,
    resolve_placeholder, maybe_submit_and_exit,
)

SECTION_KEY = "2.5.1"


def load_genome(path: Path) -> dict:
    seqs: dict = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def load_placement_table(table_path: Path, log) -> dict:
    """mge_id -> (rec_type, cluster_number, duplicate_group_number, singleton),
    taking the FIRST row seen per mge_id. A multi-recombinase MGE could in
    principle have rows of more than one recombinase_type -- warns if so
    (rare in practice), files under whichever type was seen first."""
    placement: dict = {}
    seen_types: dict = defaultdict(set)
    with open(table_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mge_id = row.get("mge_id", "")
            if not mge_id:
                continue
            rt = row.get("recombinase_type") or "UNKNOWN"
            seen_types[mge_id].add(rt)
            if mge_id not in placement:
                placement[mge_id] = (
                    rt,
                    row.get("cluster_number", "NA"),
                    row.get("duplicate_group_number", "NA"),
                    row.get("singleton", "yes"),
                )
    n_conflicts = sum(1 for types in seen_types.values() if len(types) > 1)
    if n_conflicts:
        log.warning(f"{n_conflicts:,} mge_id(s) have rows with more than one "
                    f"distinct recombinase_type -- filed under the first type seen for each.")
    return placement


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--gff", default=None, help="Override BASE/files/non_nested_is_tn.gff3")
    ap.add_argument("--genome-dir", default=None, help="Override BASE/genome_seqs_repr")
    ap.add_argument("--table", default=None, help="Override BASE/files/c_mge_recombinase_table.tsv")
    ap.add_argument("--sign-map", default=None, help="Override the vclust *.duplicates.txt path (Section 2.4.1)")
    ap.add_argument("--out-dir", default=None, help="Override BASE/flanks_split_acc_to_duplicates")
    ap.add_argument("--is-tn-out-dir", default=None,
                     help="Override BASE/is_tns_seqs (only used if extract_is_tn_seqs: true)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_5"][SECTION_KEY]

    log = get_logger("2.5.1_extract_flanks", paths["logs_dir"])

    gff_path = Path(args.gff) if args.gff else files_dir / "non_nested_is_tn.gff3"
    genome_dir = Path(args.genome_dir) if args.genome_dir else paths["genome_seqs_repr"]
    table_path = Path(args.table) if args.table else files_dir / "c_mge_recombinase_table.tsv"
    sign_map_path = Path(args.sign_map) if args.sign_map else find_duplicates_file(paths["base_dir"])
    out_dir = Path(args.out_dir) if args.out_dir else paths["base_dir"] / "flanks_split_acc_to_duplicates"
    extract_is_tn_seqs = bool(section_cfg.get("extract_is_tn_seqs", False))
    is_tn_out_dir = Path(args.is_tn_out_dir) if args.is_tn_out_dir else paths["base_dir"] / "is_tns_seqs"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Flank extraction (per-MGE, filed by rec_type/cluster/duplicate group)")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_dir, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_5_1_extract_flanks"):
        return

    for p, label in (
        (gff_path, "non_nested_is_tn.gff3 (Section 2.1.1)"),
        (genome_dir, "representative genomes (Section 2.2.3)"),
        (table_path, "c_mge_recombinase_table.tsv (Section 2.4.2)"),
    ):
        if not p.exists():
            log.error(f"{p} not found -- {label} must exist first.")
            raise SystemExit(1)

    if sign_map_path is None or not sign_map_path.exists():
        log.error(
            f"No *.duplicates.txt found under {paths['base_dir'] / 'vclust_deduplication'} "
            f"-- run 2.4.1_vclust_dedup.sh first (needed for +/- orientation correction)."
        )
        raise SystemExit(1)

    flank_size = section_cfg.get("flank_size", 530)
    mge_overlap = section_cfg.get("mge_overlap", 30)

    t0 = time.time()

    log.info(f"[0] Loading +/- orientation sign map from {sign_map_path}")
    sign_map = load_mge_sign_map(sign_map_path)
    log.info(f"    {len(sign_map):,} mge_ids have an explicit sign; anything else defaults to '-'")

    log.info(f"[1] Indexing genomes in {genome_dir}")
    genome_index = {
        f.stem: f for f in genome_dir.iterdir()
        if f.is_file() and f.suffix.lower() in (".fasta", ".fna")
    }
    log.info(f"    Indexed {len(genome_index):,} genomes")

    log.info(f"[2] Loading placement info from {table_path}")
    placement = load_placement_table(table_path, log)
    log.info(f"    {len(placement):,} mge_ids have a rec_type/cluster/group assignment")

    log.info(f"[3] Parsing GFF3: {gff_path}")
    mge_records = []
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
            contig = cols[0]
            start, end = int(cols[3]), int(cols[4])
            mge_records.append((mge_id, contig, start, end))
    log.info(f"    Parsed {len(mge_records):,} MGE records")

    log.info(f"Flank size  : {flank_size} nt")
    log.info(f"MGE overlap : {mge_overlap} nt into the MGE on each side")
    log.info(f"Output root : {out_dir}")
    log.info(f"extract_is_tn_seqs : {extract_is_tn_seqs}" +
             (f"  -> {is_tn_out_dir}" if extract_is_tn_seqs else ""))
    log.info("-" * 70)

    if extract_is_tn_seqs:
        is_tn_out_dir.mkdir(parents=True, exist_ok=True)

    stats: Counter = Counter()
    genome_cache: dict = {}

    for mge_id, contig, start, end in mge_records:
        stats["total"] += 1

        if mge_id not in placement:
            stats["no_placement"] += 1
            continue

        taxid_sample = taxid_sample_from_mge_id(mge_id)
        if not taxid_sample or taxid_sample not in genome_index:
            stats["no_genome"] += 1
            continue

        if taxid_sample not in genome_cache:
            genome_cache[taxid_sample] = load_genome(genome_index[taxid_sample])
        genome = genome_cache[taxid_sample]

        if contig not in genome:
            stats["no_contig"] += 1
            continue
        contig_seq = genome[contig]

        if start < 1 or end > len(contig_seq):
            stats["bad_coords"] += 1
            continue

        # All arithmetic is 0-based Python slice notation; GFF3 coords are
        # 1-based inclusive (MGE occupies indices [start-1 : end]).
        #
        # BUG FIX (see module docstring): left_end_excl used to be
        # (start-1)+mge_overlap+1 -- one too many. (start-1)+mge_overlap is
        # ALREADY the correct exclusive upper bound for a slice capturing
        # exactly the MGE's first mge_overlap bases; verified directly
        # against the mge_overlap=0 case (should touch zero MGE bases, and
        # now does).
        left_end_excl = (start - 1) + mge_overlap
        left_start_incl = max(0, left_end_excl - flank_size)
        left_flank = contig_seq[left_start_incl:left_end_excl]

        right_anchor = end - mge_overlap
        right_end_excl = min(len(contig_seq), right_anchor + flank_size)
        right_flank = contig_seq[right_anchor:right_end_excl]

        if len(left_flank) < flank_size:
            stats["short_left"] += 1
        if len(right_flank) < flank_size:
            stats["short_right"] += 1

        rec_type, cluster_number, duplicate_group_number, singleton = placement[mge_id]
        mge_dir = flank_output_dir(out_dir, rec_type, cluster_number, duplicate_group_number, singleton)
        mge_dir.mkdir(parents=True, exist_ok=True)

        sign = sign_map.get(mge_id, "-")
        if sign == "+":
            stats["plus_strand"] += 1
            out_left = revcomp(right_flank)
            out_right = revcomp(left_flank)
        else:
            stats["minus_strand_or_singleton"] += 1
            out_left = left_flank
            out_right = right_flank

        with open(mge_dir / f"{mge_id}_left.fasta", "w") as f:
            f.write(f">{mge_id}_left\n{out_left}\n")
        with open(mge_dir / f"{mge_id}_right.fasta", "w") as f:
            f.write(f">{mge_id}_right\n{out_right}\n")

        if extract_is_tn_seqs:
            # Same [start-1:end] convention as the flank anchors above --
            # the MGE's own full annotated body, untrimmed.
            is_tn_seq = contig_seq[start - 1:end]
            out_is_tn = revcomp(is_tn_seq) if sign == "+" else is_tn_seq
            with open(is_tn_out_dir / f"{mge_id}.fasta", "w") as f:
                f.write(f">{mge_id}\n{out_is_tn}\n")
            stats["is_tn_extracted"] += 1

        stats["success"] += 1

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Total MGEs in GFF3                       : {stats['total']:,}")
    log.info(f"  no rec_type/cluster/group assignment   : {stats['no_placement']:,}")
    log.info(f"  no representative genome                : {stats['no_genome']:,}")
    log.info(f"  contig not found in genome               : {stats['no_contig']:,}")
    log.info(f"  coordinates out of range                 : {stats['bad_coords']:,}")
    log.info(f"Successfully extracted (fasta pairs written): {stats['success']:,}")
    log.info(f"  '+' strand (revcomp + swapped)           : {stats['plus_strand']:,}")
    log.info(f"  '-' strand or singleton (saved as-is)    : {stats['minus_strand_or_singleton']:,}")
    log.info(f"  truncated left flank                     : {stats['short_left']:,}")
    log.info(f"  truncated right flank                    : {stats['short_right']:,}")
    if extract_is_tn_seqs:
        log.info(f"  is_tn sequences also extracted           : {stats['is_tn_extracted']:,}")
        log.info(f"-> {is_tn_out_dir}")
    log.info(f"-> {out_dir}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()