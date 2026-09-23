#!/usr/bin/env python3
"""
2.6.2_consensus_sanity_check.py  —  Section 2.6.2 (NEW)
============================================================
input : BASE/consensus/duplicate_group_level/<rec_type>/cluster_<N>/
        duplicate_group_<M>/{left,right}.raw.fasta  (Section 2.6.1)
output: BASE/consensus/sanity_check/sanity_check_long.log
        BASE/consensus/sanity_check/sanity_check_short.log

WHAT THIS CHECKS, AND WHY
------------------------------
2.5.1_extract_flanks.py deliberately extracts flanks so that mge_overlap
nt of every member's flank come from INSIDE the MGE itself, not the
flanking host genome -- specifically so that, WITHIN one duplicate group
(members already known to be ~100% identical copies of the same
element, per Section 2.4's own vclust dedup), that mge_overlap-derived
portion should be IDENTICAL across every member. If it isn't, something
upstream is wrong -- a genome mismatch, a coordinate error, a
duplicate-group assignment that shouldn't have been made, or an
orientation/sign bug -- and it would be silently invisible to
everything downstream, which only ever looks at consensus OUTSIDE this
region (see 2.6.3_build_consensus.py's own mge_overlap exclusion). This
script exists purely to catch that category of problem directly, by
brute-force comparison, before it can hide inside a consensus result.

SIDE-AWARE ANCHOR POSITION (see chat -- read this before trusting the
output): the request this was built from said "the first 30nt" without
qualifying by side. Taken completely literally, checking the literal
first mge_overlap characters of every LEFT sequence would NOT be
checking the MGE-derived anchor at all -- re-derived directly from
2.5.1_extract_flanks.py's own coordinate math (and re-verified multiple
times earlier in this conversation): a LEFT flank's boundary-adjacent,
MGE-derived region is at the END of its extracted sequence, not the
start; only RIGHT flanks have it at the start. Checking the literal
first mge_overlap characters of a LEFT flank would be comparing
genomic, non-MGE sequence -- exactly the region that's EXPECTED to
differ between independent insertion sites, so a "mismatch" there would
be meaningless noise, not a real finding. Given the stated PURPOSE is
explicitly to catch the MGE-derived part deviating, this script anchors
side-aware by default (last mge_overlap nt for LEFT, first mge_overlap
nt for RIGHT) -- set anchor_side_aware: false in config.yaml to instead
check the literal first mge_overlap characters of every file regardless
of side, if that's what was actually meant.

LOGIC, per (rec_type, cluster, duplicate_group, side):
  1. Read {side}.raw.fasta. If it doesn't exist or has zero sequences,
     logged distinctly (long log only) and excluded from every short-log
     tally -- neither "lonely" nor "multi", since there's nothing there
     to characterize either way.
  2. Exactly 1 sequence ("2 or more '>' signs" in the request's own
     parenthetical is the operative definition of the "multi" bucket
     below, i.e. >= 2 -- 1 is its own, separate "lonely" bucket): no
     second copy to compare against, logged as lonely, tallied
     separately, not further checked.
  3. >= 2 sequences: each sequence's anchor (mge_overlap nt, position
     per anchor_side_aware above) is extracted. If a sequence is
     shorter than mge_overlap, its anchor is however much of it exists
     (never crashes; logged in the long log if this happens, since a
     short flank is itself worth knowing about). If every sequence's
     anchor is EXACTLY identical (case-sensitive, straight string
     equality): logged as matching. Otherwise: logged as a MISMATCH,
     with every distinct anchor variant found, its count, and which
     mge_ids carry it -- specifically so a real problem is
     immediately actionable (which sequences to go look at), not just
     flagged as "something's wrong somewhere in this group".

rec_types (config, default "all"): restrict to specific rec_type folder
names under duplicate_group_level/, same convention as the rest of
Section 2.6/2.7.

Run standalone:
    python 2.6.2_consensus_sanity_check.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.6.2"


def parse_rec_types(arg: Optional[str]) -> Optional[set]:
    if arg is None or str(arg).strip().lower() == "all":
        return None
    return {t.strip() for t in str(arg).split(",") if t.strip()}


def read_fasta_pairs(path: Path) -> List[Tuple[str, str]]:
    """[(header_id, sequence), ...] -- keeps the header (mge_id) alongside
    each sequence, unlike the plain read_fasta() used elsewhere, since
    a mismatch report needs to name WHICH mge_ids carry which variant."""
    pairs: List[Tuple[str, str]] = []
    name = None
    buf: List[str] = []
    if not path.exists():
        return pairs
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    pairs.append((name, "".join(buf)))
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip().upper())
        if name is not None:
            pairs.append((name, "".join(buf)))
    return pairs


def anchor_of(seq: str, side: str, mge_overlap: int, side_aware: bool) -> str:
    if not side_aware:
        return seq[:mge_overlap]
    return seq[-mge_overlap:] if side == "left" else seq[:mge_overlap]


def find_duplicate_groups(dup_group_level_dir: Path, rec_types_filter: Optional[set]):
    """Yields (rec_type, cluster_name, group_name, group_dir) for every
    duplicate_group_* folder found -- singletons (which never get their
    own duplicate_group_* folder in the first place, per Section 2.4/2.5's
    own convention) are structurally excluded just by this glob, same as
    every other Section 2.6/2.7 script that walks this tree."""
    for rt_dir in sorted(dup_group_level_dir.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue
        for cluster_dir in sorted(rt_dir.glob("cluster_*")):
            if not cluster_dir.is_dir():
                continue
            for group_dir in sorted(cluster_dir.glob("duplicate_group_*")):
                if not group_dir.is_dir():
                    continue
                yield rec_type, cluster_dir.name, group_dir.name, group_dir


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--consensus-dir", default=None, help="Override BASE/consensus")
    ap.add_argument("--out-dir", default=None, help="Override BASE/consensus/sanity_check")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    section_cfg = cfg.get("section_2_6", {}).get(SECTION_KEY, {})

    log = get_logger("2.6.2_consensus_sanity_check", paths["logs_dir"])

    consensus_dir = Path(args.consensus_dir) if args.consensus_dir else paths["base_dir"] / "consensus"
    dup_group_level_dir = consensus_dir / "duplicate_group_level"
    out_dir = Path(args.out_dir) if args.out_dir else consensus_dir / "sanity_check"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Consensus sanity check (mge_overlap anchor consistency)")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_dir, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_6_2_consensus_sanity_check"):
        return

    if not dup_group_level_dir.exists():
        log.error(f"{dup_group_level_dir} not found -- run 2.6.1_align_flanks.sh first.")
        raise SystemExit(1)

    mge_overlap = int(section_cfg.get("mge_overlap", 30))
    side_aware = bool(section_cfg.get("anchor_side_aware", True))
    rec_types_filter = parse_rec_types(section_cfg.get("rec_types", "all"))

    log.info(f"duplicate_group_level dir : {dup_group_level_dir}")
    log.info(f"mge_overlap               : {mge_overlap}")
    log.info(f"anchor_side_aware         : {side_aware}"
             + ("  (last {}nt for LEFT, first {}nt for RIGHT)".format(mge_overlap, mge_overlap) if side_aware
                else "  (literal first {}nt of every file, regardless of side)".format(mge_overlap)))
    log.info(f"rec_types filter          : {'ALL' if rec_types_filter is None else sorted(rec_types_filter)}")
    log.info("-" * 70)

    out_dir.mkdir(parents=True, exist_ok=True)
    long_log_path = out_dir / "sanity_check_long.log"
    short_log_path = out_dir / "sanity_check_short.log"

    # per rec_type tallies for the short log
    tallies: Dict[str, Counter] = defaultdict(Counter)  # rec_type -> {lonely, multi_matching, multi_mismatching, empty}

    n_groups_seen = 0
    n_units_seen = 0

    with open(long_log_path, "w") as long_log:
        for rec_type, cluster_name, group_name, group_dir in find_duplicate_groups(dup_group_level_dir, rec_types_filter):
            n_groups_seen += 1
            for side in ("left", "right"):
                n_units_seen += 1
                fasta_path = group_dir / f"{side}.raw.fasta"
                pairs = read_fasta_pairs(fasta_path)
                label = f"{rec_type} / {cluster_name} / {group_name} / {side}"

                if len(pairs) == 0:
                    tallies[rec_type]["empty"] += 1
                    long_log.write(f"{label}: NO SEQUENCES (file missing or empty: {fasta_path})\n")
                    continue

                if len(pairs) == 1:
                    tallies[rec_type]["lonely"] += 1
                    mge_id, seq = pairs[0]
                    long_log.write(f"{label}: lonely (1 sequence -- {mge_id}, {len(seq)}nt)\n")
                    continue

                # >= 2 sequences -- check anchor consistency
                short_flank_notes = []
                anchor_to_mge_ids: Dict[str, List[str]] = defaultdict(list)
                for mge_id, seq in pairs:
                    if len(seq) < mge_overlap:
                        short_flank_notes.append(f"{mge_id} (only {len(seq)}nt, shorter than mge_overlap={mge_overlap})")
                    anchor = anchor_of(seq, side, mge_overlap, side_aware)
                    anchor_to_mge_ids[anchor].append(mge_id)

                if len(anchor_to_mge_ids) == 1:
                    tallies[rec_type]["multi_matching"] += 1
                    (only_anchor,) = anchor_to_mge_ids.keys()
                    long_log.write(f"{label}: n={len(pairs)} sequences, anchor MATCHES across all -- {only_anchor!r}\n")
                else:
                    tallies[rec_type]["multi_mismatching"] += 1
                    long_log.write(f"{label}: n={len(pairs)} sequences, anchor MISMATCH -- "
                                    f"{len(anchor_to_mge_ids)} distinct variant(s):\n")
                    for variant, mge_ids in sorted(anchor_to_mge_ids.items(), key=lambda kv: -len(kv[1])):
                        long_log.write(f"    {variant!r}  (n={len(mge_ids)}): {', '.join(mge_ids)}\n")

                if short_flank_notes:
                    long_log.write(f"    NOTE -- shorter than mge_overlap: {'; '.join(short_flank_notes)}\n")

    # ── short log ────────────────────────────────────────────────────────
    with open(short_log_path, "w") as short_log:
        header = f"{'rec_type':30s} {'lonely':>8s} {'multi_matching':>15s} {'multi_mismatching':>18s} {'empty':>8s}"
        short_log.write(header + "\n")
        short_log.write("-" * len(header) + "\n")
        grand = Counter()
        for rec_type in sorted(tallies):
            t = tallies[rec_type]
            short_log.write(f"{rec_type:30s} {t['lonely']:8,d} {t['multi_matching']:15,d} "
                             f"{t['multi_mismatching']:18,d} {t['empty']:8,d}\n")
            grand.update(t)
        short_log.write("-" * len(header) + "\n")
        short_log.write(f"{'TOTAL':30s} {grand['lonely']:8,d} {grand['multi_matching']:15,d} "
                         f"{grand['multi_mismatching']:18,d} {grand['empty']:8,d}\n")

    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"duplicate groups scanned      : {n_groups_seen:,}")
    log.info(f"(rec_type,cluster,group,side) units checked : {n_units_seen:,}")
    grand = Counter()
    for t in tallies.values():
        grand.update(t)
    log.info(f"  lonely (1 sequence)          : {grand['lonely']:,}")
    log.info(f"  multi, anchor matching       : {grand['multi_matching']:,}")
    log.info(f"  multi, anchor MISMATCH       : {grand['multi_mismatching']:,}")
    log.info(f"  empty/missing file           : {grand['empty']:,}")
    if grand["multi_mismatching"] > 0:
        log.warning(f"*** {grand['multi_mismatching']:,} unit(s) have a mismatching mge_overlap anchor within "
                    f"a duplicate group that's supposed to be near-identical -- see {long_log_path} for detail.")
    log.info(f"-> {long_log_path}")
    log.info(f"-> {short_log_path}")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()