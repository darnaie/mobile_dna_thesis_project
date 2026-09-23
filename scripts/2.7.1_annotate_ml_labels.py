#!/usr/bin/env python3
"""
2.7.1_annotate_ml_labels.py  —  Section 2.7.1
=================================================
input : BASE/flanks_split_acc_to_duplicates/ (2.5.1, per-MGE flank FASTAs)
        BASE/consensus/cluster_level/ (2.6.1, cluster-level alignments)
output: BASE/files/ml_labels.tsv (one row per duplicate-group MGE per
        side -- this is the labeled dataset 2.7.2 splits and 2.7.3 trains on)

MGE-OVERLAP EXCLUSION (see chat -- applies throughout Section 2.7 now,
not just here)
------------------------------------------------------------------------
2.5.1_extract_flanks.py deliberately extracts each flank so that its
LAST mge_overlap nt (LEFT) or FIRST mge_overlap nt (RIGHT) sit INSIDE
the MGE itself, not in genuinely flanking genomic sequence (see
2.5.1's own docstring: this is what anchors the window on a shared
reference point across independent insertions). That's the right thing
for 2.5.1's own purpose, but it means the raw flank_sequence Section 2.7
used to work with was never PURELY flanking sequence -- the boundary-
adjacent mge_overlap nt of it are, literally, still part of the element.
Labeling, training, cross-validation, and visualization all now
exclude that zone, so "the flank" means only genuinely flanking
sequence throughout Section 2.7: this script strips it from BOTH the
alignment sequences used for the consensus scan AND the flank_sequence
column written to ml_labels.tsv (via strip_mge_overlap() below), so
every downstream script (2.7.2 onward) inherits the exclusion for free
just by reading flank_sequence as before -- none of them need their own
copy of this logic. New config key: mge_overlap (default 30, matching
2.5.1's own default -- keep these two in sync).

One direct consequence: sequences are now ~(flank_size - mge_overlap)
nt long (e.g. 500, not 530) instead of the full flank_size -- update
seq_length in 2.7.3/2.7.4/2.7.5's config accordingly, or the model will
zero-pad every input by mge_overlap nt for no reason.

WHAT THIS DOES, AND WHY (see chat)
------------------------------------
Builds a soft-probability training label for "is there a boundary-marking
motif (TIR-like) right at this flank's edge" for every duplicate-group
MGE in the rec_types listed below -- SIX real families plus one
deliberate negative control, all cross-checked against the empirical
family-level reference table you provided (from TnCentral/ISfinder-style
summaries of Ends/DR/IRs per family and sub-group):

  DDE_Tnp_IS1        DDE_Tnp_IS1595     DDE_Tnp_IS66
  DDE_Tnp_ISL3       DDE_Tnp_Tn3        DDE_Tnp_ISAZ013
  DEDD_Tnp_IS110     <- negative control: marked "N" for IRs in the
                         reference table, no Ends pattern given at all.
                         Goes through the EXACT SAME procedure as the
                         other six -- it just structurally can never
                         reach the top label tier, because there is no
                         literature pattern for it to match. That's a
                         cleaner negative control than hand-overriding
                         its labels: the model has to learn "this family
                         doesn't have the signal" from the data itself,
                         same as it has to learn the positive families.

FIXED: RIGHT-SIDE trim_start ASYMMETRY (see chat)
------------------------------------------------------
Re-derived directly from 2.5.1_extract_flanks.py's own coordinate math:
the LEFT flank's boundary-adjacent region (the mge_overlap nt that sit
INSIDE the MGE) is at the END of its extracted sequence; the RIGHT
flank's boundary-adjacent region is at the START. Naively trimming the
first N columns of every alignment (as consensus_building_evaluation.py's
plain scan_alignment_multi_segment_with_skip() does) is therefore only
correct for LEFT -- for RIGHT, with trim_start >= mge_overlap (true of
your production config: both are 30), it was eating directly into the
boundary-adjacent region itself, up to and including truncating or
destroying real embedded motifs before this script ever got a chance to
find them. Reproduced directly, not just reasoned about: an 8 nt motif
sitting at a right flank's own start, scanned with trim=5, returned only
its last 3 nt instead of the full 8.

Fixed via scan_alignment_multi_segment_with_skip_side_aware() (new,
purely additive function in consensus_building_evaluation.py -- the
plain, side-blind function it wraps is completely untouched, so nothing
about the LEFT side, or any other existing caller, changed). For
side="right" it reverses each sequence's character order (NOT
reverse-complement -- purely about which end trim_start skips, nothing
to do with strand) before scanning, then maps the result back into the
flank's own original coordinates. This script now uses that side-aware
version for both sides.

2.6.3_build_consensus.py's own alignment branch, and
consensus_building_evaluation.py's own Strategy C/D/E row computation,
have NOT been switched over to the side-aware version -- that's a
deliberate scope decision (see chat), since it would silently change
cluster_consensus_results.tsv / e_mge_recombinase_table.tsv and this
script's own comparison table, both of which you may already be relying
on elsewhere. Worth doing as a separate, explicit step if you want
consistency across the whole pipeline, not bundled into this fix.

SECOND FIX, FOUND WHILE VERIFYING THE FIRST ONE (see chat): even with
trim_start correctly side-aware, a genuinely conserved boundary motif
often still wasn't being recovered -- because it sits at the literal
EDGE of the extracted flank window (true for BOTH sides, by construction
of 2.5.1_extract_flanks.py), there's frequently nothing beyond it, WITHIN
the window, to trigger a low-identity column and formally CLOSE the
scan's current run. scan_alignment_multi_segment_with_skip()'s
`segments` list, by design, never includes an unclosed ("maxed") run --
correct for that function's own purpose, but it meant a real
boundary-adjacent motif reaching the edge of the window was silently
invisible to this script's label. Plausibly a large part of why so many
real clusters come back maxed_out (829 of 1,054 -- 78.7% -- in the
six-family evaluation run already produced): that maxed_out status may
often mean "conservation held all the way to the edge of what we
extracted," not "nothing here."

Fixed via scan_alignment_multi_segment_with_skip_trailing() /
scan_alignment_multi_segment_with_skip_full_side_aware() (new, in
consensus_building_evaluation.py, alongside the first fix), which also
returns that trailing run's own (start, text) when present. This script
now checks the trailing run FIRST -- it's always the single closest
possible candidate to the boundary, by construction, whenever it
qualifies (length in [4,20]) -- and only falls back to the best CLOSED
segment if there is no trailing run, or it doesn't qualify. The
label_reason column records which path a given label came from
(*_trailing_at_boundary vs *_closed_segment) so this is auditable
per-row, not just asserted.

LABELING PROCEDURE, per (rec_type, cluster, side) -- computed ONCE per
cluster/side and then applied to every duplicate-group MGE in that
cluster (mirrors 2.6.2's own cluster-level sharing exactly; singletons
are never touched at all, since this script only ever iterates
duplicate_group_* folders in the first place -- same exclusion as
2.6.1/2.6.3/find_terminal_inverted_repeats.py's multi-copy branch, for
the same reason: a singleton has no duplicate to build a shared
alignment from):
  1. Re-scan that cluster's *.aln.fasta with the EXACT SAME strategy as
     your real 2.6.2 config (Strategy D: scan_alignment_multi_segment_
     with_skip, threshold=fixed_threshold, run_len=a_run_length,
     trim=trim_start) -- imported directly from
     consensus_building_evaluation.py, not reimplemented, so this always
     matches whatever cluster_consensus_results.tsv itself would show.
  2. Unlike 2.6.2 (which keeps the LONGEST closed segment as "the"
     consensus), this keeps every segment with length in [4, 20]
     (inclusive both ends -- deliberately NOT the strict "> 4 and < 20"
     bucket the evaluation table's summary columns use, which silently
     drops exactly-4 and exactly-20 nt segments; that's fine for a
     summary statistic, not fine for something you're about to use as a
     training label) and picks whichever ONE of those sits CLOSEST to
     the actual MGE boundary -- the largest-start segment for LEFT (its
     boundary-adjacent region is at the end), the smallest-start segment
     for RIGHT (its boundary-adjacent region is at the start).
  3. No such segment (cluster fully maxed-out, or no closed segment ever
     lands in the [4,20] window near the boundary) -> label = 0.3.
  4. A segment was found -> label = 0.7 baseline. If that segment's own
     sequence fuzzy-matches ANY of the family's reference Ends patterns
     from the table (case in the pattern is meaningful: UPPERCASE
     position must match exactly, lowercase/n = wildcard; checked in
     both the segment's own orientation and its reverse complement,
     since it isn't known in advance which strand the reference pattern
     was written against) -> label = 1.0.
These three numbers (0.3 / 0.7 / 1.0) are a starting point, not a
calibrated result -- worth a sensitivity check later (retrain once with
different values, see how much the trained model's behaviour actually
moves) rather than treating them as fixed truth.

Run standalone:
    python 2.7.1_annotate_ml_labels.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _ilu
import re
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.7.1"

# ══════════════════════════════════════════════════════════════════════════
#  REFERENCE ENDS PATTERNS -- transcribed from the family/sub-group tables
#  you provided (TnCentral/ISfinder-style summaries). UPPERCASE = must
#  match exactly at that position; lowercase or 'n'/'N' = wildcard.
#  DEDD_Tnp_IS110 is deliberately absent -- no Ends pattern was given for
#  it in the reference table (it's also marked "N" for IRs there), which
#  is exactly what makes it a usable negative control: it goes through
#  the identical procedure below and can still never reach label 1.0.
# ══════════════════════════════════════════════════════════════════════════

REFERENCE_ENDS_PATTERNS: Dict[str, List[str]] = {
    "DDE_Tnp_IS1": ["GGnnnTG"],
    "DDE_Tnp_IS1595": [
        "GGCnnTG",      # ISPna2
        "CGCTCTT",      # ISH4
        "GGGgctg",      # IS1016
        "CcTGATT",      # IS1595 (the sub-group sharing the family's own name)
        "nnnGcnTATC",   # ISSod11
        "ggnnatTAT",    # ISNwi1
        "CGGnnTT",      # ISNha5
    ],
    "DDE_Tnp_IS66": ["GTAA"],
    "DDE_Tnp_ISL3": ["GG"],
    "DDE_Tnp_Tn3": ["GGGG"],
    "DDE_Tnp_ISAZ013": ["Gagg", "GaaG"],  # "Ga/g" in the table -- both alternatives at that position
    # DEDD_Tnp_IS110: intentionally no entry -- see docstring.
}

TARGET_REC_TYPES = [
    "DDE_Tnp_IS1", "DDE_Tnp_IS1595", "DDE_Tnp_IS66",
    "DDE_Tnp_ISL3", "DDE_Tnp_Tn3", "DDE_Tnp_ISAZ013",
    "DEDD_Tnp_IS110",
]

MIN_LEN, MAX_LEN = 4, 20  # inclusive both ends -- see docstring on why this differs from the eval table's ">4 and <20"




def revcomp(seq: str) -> str:
    # BUG FIX (see chat): this was missing the [::-1] reversal -- it only
    # complemented, never reversed, so pattern_matches()'s "check the
    # reverse-complement orientation too" branch was silently checking the
    # plain COMPLEMENT instead for as long as this function has existed.
    # Caught by testing a fixture where a genuine TIR pair (right-side
    # motif stored as the true reverse complement of the left-side one --
    # the biologically realistic case) needed real revcomp matching to
    # succeed; it didn't, until this fix.
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def pattern_matches(seq: str, pattern: str) -> bool:
    """Slide a window the length of `pattern` across `seq` (and across
    revcomp(seq)); a match at any offset, in either orientation, counts.
    Case in `pattern` is meaningful (see module docstring); case in
    `seq` is not (compared uppercase)."""
    L = len(pattern)
    if len(seq) < L:
        return False

    def try_one(s: str) -> bool:
        for i in range(len(s) - L + 1):
            window = s[i:i + L]
            ok = True
            for w_ch, p_ch in zip(window, pattern):
                if p_ch.isupper() and p_ch != "N" and w_ch != p_ch:
                    ok = False
                    break
            if ok:
                return True
        return False

    return try_one(seq.upper()) or try_one(revcomp(seq.upper()))


def matches_any_reference_pattern(seq: str, rec_type: str) -> bool:
    patterns = REFERENCE_ENDS_PATTERNS.get(rec_type, [])
    return any(pattern_matches(seq, p) for p in patterns)


# ══════════════════════════════════════════════════════════════════════════
#  DISCOVERY -- duplicate_group_* folders only, matching the multi-copy
#  branch's exclusion of singletons everywhere else in this pipeline
# ══════════════════════════════════════════════════════════════════════════

def discover_duplicate_group_mges(flanks_dir: Path, rec_types: List[str]) -> Dict[Tuple[str, str], List[Tuple[str, str]]]:
    """Returns {(rec_type, cluster_num): [(duplicate_group_num, mge_id), ...]}."""
    result: Dict[Tuple[str, str], List[Tuple[str, str]]] = defaultdict(list)
    for rec_type in rec_types:
        rec_dir = flanks_dir / rec_type
        if not rec_dir.is_dir():
            continue
        for cluster_dir in sorted(rec_dir.glob("cluster_*")):
            if not cluster_dir.is_dir():
                continue
            cm = re.fullmatch(r"cluster_(\d+)", cluster_dir.name)
            if not cm:
                continue
            cluster_num = cm.group(1)
            for group_dir in sorted(cluster_dir.glob("duplicate_group_*")):
                if not group_dir.is_dir():
                    continue
                gm = re.fullmatch(r"duplicate_group_(\d+)", group_dir.name)
                if not gm:
                    continue
                dup_group = gm.group(1)
                mge_ids = sorted({
                    fp.name[: -len("_left.fasta")] for fp in group_dir.glob("*_left.fasta")
                } | {
                    fp.name[: -len("_right.fasta")] for fp in group_dir.glob("*_right.fasta")
                })
                for mge_id in mge_ids:
                    result[(rec_type, cluster_num)].append((dup_group, mge_id))
    return result


def read_fasta_single(path: Path) -> str:
    seq_lines: List[str] = []
    with open(path, encoding="ascii", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line.upper())
    return "".join(seq_lines)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--flanks-dir", default=None, help="Override BASE/flanks_split_acc_to_duplicates")
    ap.add_argument("--consensus-dir", default=None, help="Override BASE/consensus")
    ap.add_argument("--out", default=None, help="Override BASE/files/ml_labels.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_7"][SECTION_KEY]

    log = get_logger("2.7.1_annotate_ml_labels", paths["logs_dir"])

    flanks_dir = Path(args.flanks_dir) if args.flanks_dir else paths["base_dir"] / "flanks_split_acc_to_duplicates"
    consensus_dir = Path(args.consensus_dir) if args.consensus_dir else paths["base_dir"] / "consensus"
    cluster_align_dir = consensus_dir / "cluster_level"
    out_path = Path(args.out) if args.out else files_dir / "ml_labels.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Annotate ML training labels")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_7_1_annotate_ml_labels"):
        return

    # Reuse the exact same, already-tested consensus-scanning logic --
    # not reimplemented, same pattern 2.6.2 already uses.
    eval_spec = _ilu.spec_from_file_location(
        "consensus_building_evaluation", Path(__file__).resolve().parent / "consensus_building_evaluation.py"
    )
    eval_mod = _ilu.module_from_spec(eval_spec)
    eval_spec.loader.exec_module(eval_mod)
    read_fasta = eval_mod.read_fasta
    scan_alignment_multi_segment_with_skip_full_side_aware = eval_mod.scan_alignment_multi_segment_with_skip_full_side_aware
    strip_mge_overlap = eval_mod.strip_mge_overlap  # moved to consensus_building_evaluation.py -- see chat; 2.6.2 imports the same copy

    threshold = float(section_cfg.get("fixed_threshold", 1.0))
    run_len = int(section_cfg.get("a_run_length", 2))
    trim_start = int(section_cfg.get("trim_start", 30))
    mge_overlap = int(section_cfg.get("mge_overlap", 30))
    rec_types_raw = section_cfg.get("rec_types", "all")
    if rec_types_raw == "all":
        rec_types = TARGET_REC_TYPES
    else:
        rec_types = [t.strip() for t in str(rec_types_raw).split(",") if t.strip()]

    log.info(f"rec_types    : {rec_types}")
    log.info(f"threshold    : {threshold}   run_len: {run_len}   trim_start: {trim_start}   mge_overlap: {mge_overlap}")
    log.info(f"length bucket: [{MIN_LEN}, {MAX_LEN}] inclusive")
    log.info("-" * 70)

    if not flanks_dir.is_dir():
        log.error(f"{flanks_dir} not found -- run 2.5.1_extract_flanks.py first.")
        raise SystemExit(1)
    if not cluster_align_dir.is_dir():
        log.error(f"{cluster_align_dir} not found -- run 2.6.1_align_flanks.sh first.")
        raise SystemExit(1)

    t0 = time.time()

    log.info("Discovering duplicate-group MGEs ...")
    cluster_members = discover_duplicate_group_mges(flanks_dir, rec_types)
    log.info(f"  {sum(len(v) for v in cluster_members.values()):,} MGEs across "
             f"{len(cluster_members):,} (rec_type, cluster) combinations")

    rows = []
    n_label_10 = n_label_07 = n_label_03 = 0

    for (rec_type, cluster_num), members in sorted(cluster_members.items()):
        prefix = f"{rec_type}_cluster_{cluster_num}"
        cluster_dir = cluster_align_dir / rec_type

        # One boundary-closest segment decision PER (cluster, side) --
        # then shared across every duplicate-group MGE in that cluster,
        # exactly mirroring 2.6.2's own cluster-level sharing.
        side_result: Dict[str, dict] = {}
        for side in ("left", "right"):
            aln_path = cluster_dir / f"{prefix}_{side}.aln.fasta"
            if not aln_path.exists():
                side_result[side] = {"label": 0.3, "reason": "no_alignment_file",
                                      "seg_start": "NA", "seg_len": "NA", "seg_seq": "NA"}
                continue
            seqs = read_fasta(aln_path)
            seqs = [strip_mge_overlap(s, side, mge_overlap) for s in seqs]
            seqs = [s for s in seqs if s]  # a member entirely inside the overlap zone contributes nothing
            if not seqs:
                side_result[side] = {"label": 0.3, "reason": "no_sequence_outside_overlap",
                                      "seg_start": "NA", "seg_len": "NA", "seg_seq": "NA"}
                continue
            segments, _maxed, trailing = scan_alignment_multi_segment_with_skip_full_side_aware(
                seqs, threshold, run_len, trim_start, side,
            )

            # The trailing (unclosed) run, when present, is BY CONSTRUCTION
            # already the boundary-closest candidate available -- for LEFT
            # it's whatever extends all the way to the sequence's own end
            # (the boundary-adjacent end); for RIGHT, after coordinate
            # correction, it's whatever extends back to position 0 (also
            # the boundary-adjacent end). No closed segment could ever be
            # closer to the boundary than that, so it's checked FIRST, and
            # only falls back to the best closed segment if it doesn't
            # qualify (wrong length, or no trailing run at all -- see
            # scan_alignment_multi_segment_with_skip_trailing()'s
            # docstring for why a genuine boundary motif very often shows
            # up here rather than as a closed segment in the first place).
            chosen = None
            chosen_source = None
            if trailing is not None and MIN_LEN <= len(trailing[1]) <= MAX_LEN:
                chosen = trailing
                chosen_source = "trailing_at_boundary"
            else:
                candidates = [(s, t) for s, t in segments if MIN_LEN <= len(t) <= MAX_LEN]
                if candidates:
                    # LEFT: boundary-adjacent region is at the END -> largest start.
                    # RIGHT: boundary-adjacent region is at the START -> smallest start.
                    chosen = max(candidates, key=lambda st: st[0]) if side == "left" else min(candidates, key=lambda st: st[0])
                    chosen_source = "closed_segment"

            if chosen is None:
                side_result[side] = {"label": 0.3, "reason": "no_segment_in_length_bucket",
                                      "seg_start": "NA", "seg_len": "NA", "seg_seq": "NA"}
                continue

            seg_start, seg_seq = chosen
            if matches_any_reference_pattern(seg_seq, rec_type):
                side_result[side] = {"label": 1.0, "reason": f"found_literature_match_{chosen_source}",
                                      "seg_start": seg_start, "seg_len": len(seg_seq), "seg_seq": seg_seq}
            else:
                side_result[side] = {"label": 0.7, "reason": f"found_no_literature_match_{chosen_source}",
                                      "seg_start": seg_start, "seg_len": len(seg_seq), "seg_seq": seg_seq}

        group_dirs_cache: Dict[str, Path] = {}
        for dup_group, mge_id in members:
            group_dir = flanks_dir / rec_type / f"cluster_{cluster_num}" / f"duplicate_group_{dup_group}"
            for side in ("left", "right"):
                flank_path = group_dir / f"{mge_id}_{side}.fasta"
                if not flank_path.exists():
                    continue
                flank_seq = read_fasta_single(flank_path)
                flank_seq = strip_mge_overlap(flank_seq, side, mge_overlap)
                r = side_result[side]
                rows.append({
                    "rec_type": rec_type, "cluster_number": cluster_num, "duplicate_group_number": dup_group,
                    "mge_id": mge_id, "side": side, "flank_sequence": flank_seq,
                    "label": r["label"], "label_reason": r["reason"],
                    "segment_start": r["seg_start"], "segment_length": r["seg_len"], "segment_sequence": r["seg_seq"],
                })
                if r["label"] == 1.0:
                    n_label_10 += 1
                elif r["label"] == 0.7:
                    n_label_07 += 1
                else:
                    n_label_03 += 1

    out_columns = [
        "rec_type", "cluster_number", "duplicate_group_number", "mge_id", "side", "flank_sequence",
        "label", "label_reason", "segment_start", "segment_length", "segment_sequence",
    ]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Rows written (MGE x side)     : {len(rows):,}")
    log.info(f"  label 1.0 (literature match): {n_label_10:,}")
    log.info(f"  label 0.7 (found, no match) : {n_label_07:,}")
    log.info(f"  label 0.3 (not found)       : {n_label_03:,}")
    for rt in rec_types:
        n = sum(1 for r in rows if r["rec_type"] == rt)
        log.info(f"  {rt:20s}: {n:,} rows")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()