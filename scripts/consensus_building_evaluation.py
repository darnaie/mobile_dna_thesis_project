#!/usr/bin/env python3
"""
consensus_building_evaluation.py
===================================
UTILITY script -- operates on the thesis_rep tree (flanks_split_acc_to_
duplicates, consensus/), NOT part of the config/section-numbered
pipeline itself, though it shares BASE/scripts with it.

Compares 15 consensus-extraction strategy/universe/threshold combinations
and writes one summary table (+ a text-rendered copy and extra detail in
a log file).

SCOPE CHANGE: duplicate-group-level STACKING is no longer evaluated
here at all -- Strategy A/B now only ever run on CLUSTER-LEVEL
stackings. Duplicate-group-level and cluster-level ALIGNMENT (Strategy C)
are both still evaluated, unchanged. Universe totals (total_mges_in /
total_clusters_in / total_duplicate_groups_in) are still computed from
the full duplicate-group-level tree, though -- that's still the most
complete "how much is really out there" reference count, even though no
ROWS are generated from it anymore.

Strategy A - "run-length tolerant cutoff", on cluster-level stacking data.
  A run of >N consecutive below-threshold columns triggers a cut
  (run length default 2). Tested at FIXED thresholds 0.9 / 0.7 / 0.5 (3
  rows), plus one DYNAMIC-threshold row (see below).

Strategy B - "single-column cutoff", same cluster-level stacking data.
  Cuts at the very first below-threshold column, no run tolerance.
  Tested at FIXED thresholds 0.9 / 0.7 / 0.5 (3 rows), plus one
  DYNAMIC-threshold row.

Strategy C - "alignment-based, multi-segment". Scans the ENTIRE
  alignment, closing off every maximal run of identity >= threshold as
  its own separate consensus SEGMENT (see scan_alignment_multi_segment's
  own docstring for the exact "maxed" semantics). Run at:
    - duplicate_group_level/<rec_type>/cluster_<N>/*.aln.fasta (1 FIXED row)
    - cluster_level/<rec_type>/*.aln.fasta (1 FIXED row + 1 DYNAMIC row)

DYNAMIC THRESHOLDING (+3 rows: one each for Strategy A, B, and C,
all on cluster-level data, matching 2.6.3_build_consensus.py's own
default behaviour): instead of one threshold fixed for the whole row,
each unit's OWN threshold is chosen from its sequence count --
threshold = dynamic_threshold_small (0.9) if n_sequences <=
dynamic_threshold_cutoff (3), else dynamic_threshold_large (0.7). This
needed no new worker functions -- _worker_stacking/_worker_alignment
already just consume whatever threshold value sits in their task tuple;
compute_stacking_row_dynamic()/compute_alignment_row_dynamic() just
resolve that threshold PER TASK before submitting to the pool, instead
of once for the whole row.

Strategy D - alignment-based, multi-segment, SAME as Strategy C (100%
  fixed threshold and dynamic thresholding both included) but with
  run-length tolerance added: a run of up to d_run (default 2)
  consecutive below-threshold columns no longer immediately closes a
  segment -- only once the low run itself reaches d_run does the
  segment close, retroactively trimmed back to exclude that whole low
  run (exactly consensus_strategy_a's "skipping" tolerance, applied to
  Strategy C's multi-segment scan instead of Strategy A's single-cutoff
  scan). +3 rows, directly mirroring Strategy C's three: group-level
  (fixed), cluster-level (fixed), cluster-level (dynamic). Implemented
  as scan_alignment_multi_segment_with_skip() / _worker_alignment_skip()
  -- _compute_alignment_row_from_tasks() takes the worker function as a
  parameter so Strategy C and D share the exact same aggregation code.

Strategy E - alignment-based, multi-segment, cluster-level ONLY (matches
  2.6.3_build_consensus.py's own scope -- no group-level row, unlike
  C/D). Threshold is FIXED at 100% (--c-threshold), same as C/D -- what's
  dynamic here is the run-length tolerance ITSELF, resolved per unit
  from its own sequence count: n_sequences <= --e-cutoff (default: same
  as --dynamic-cutoff) -> NO skipping at all (run_len=1, exactly
  Strategy C's strict behaviour for that unit); otherwise skipping IS
  allowed (run_len=--d-run, exactly Strategy D's behaviour for that
  unit). +1 row. Implemented as compute_alignment_row_dynamic_nt_skip()
  -- reuses _worker_alignment_skip() directly (run_len=1 there already
  behaves identically to the strict scan, since a single low column
  immediately satisfies "low_run >= 1"), so no new worker function was
  needed, only a new per-unit run_len resolution step before task
  construction (same pattern as the dynamic-THRESHOLD rows use for
  threshold instead of run_len).

SEGMENT START POSITIONS: scan_alignment_multi_segment() and
scan_alignment_multi_segment_with_skip() return (start, text) pairs per
segment instead of bare text -- `start` is the segment's starting
column in the ORIGINAL (untrimmed) alignment, i.e. already offset by
`trim`. _worker_alignment/_worker_alignment_skip unpack this
transparently (only the text is used for this script's own length-based
statistics); 2.6.3_build_consensus.py is what actually surfaces the
position, in cluster_consensus_results.tsv.

BUG FIX -- directory scanning, thesis_rep tree: scan_grouped_dir() was
written for a flat legacy tree shape where group-level files sit
directly inside cluster_N/. The real thesis_rep tree nests one level
deeper (cluster_N/duplicate_group_M/ or cluster_N/singleton/), so
scan_grouped_dir()'s `if not fp.is_file(): continue` skipped literally
everything when pointed at it -- every "universe total" (total_mges_in
etc.) silently came out as 0, and both Strategy C/D's group-level rows
came out completely empty, while the cluster-level rows (which use the
unaffected scan_cluster_level_dir(), since 2.6.1 writes THAT tree flat)
looked fine and hid the problem. Fixed with two tree-shape-aware
scanners used by default now: scan_duplicate_group_level_dir() (for
consensus/duplicate_group_level/) and
scan_flanks_tree_for_universe_totals() (for
flanks_split_acc_to_duplicates/, deliberately skipping singleton/
folders -- see that function's own docstring for why). scan_grouped_dir()
itself is left in place, unused by default, only in case you ever need
to point --stacking-dir/--alignments-dir back at the original legacy
tree shape it was written for.

BUG FIX -- trim_start was not side-aware (see chat): re-derived
directly from 2.5.1_extract_flanks.py's own coordinate arithmetic, a
LEFT flank's boundary-adjacent (mge_overlap) region sits at the END of
its extracted sequence, while a RIGHT flank's sits at the START. `trim`
unconditionally skipping the first N columns of whatever it's given is
therefore only skipping boundary-DISTANT sequence for LEFT -- for RIGHT
it eats into the boundary-adjacent region itself, and can truncate or
destroy the very motif being searched for (an 8 nt motif sitting at a
right flank's own start, scanned with trim=5, returns only its last 3
nt instead of the full 8 -- reproduced directly, not just reasoned
about). Fixed via a new, purely ADDITIVE wrapper,
scan_alignment_multi_segment_with_skip_side_aware() -- the existing
scan_alignment_multi_segment_with_skip() is completely unchanged, so
every already-validated LEFT-side (and any caller not yet updated to
pass a side) result is untouched. For side="right", the wrapper reverses
each sequence's character order (NOT reverse-complement -- this is
purely about which physical end `trim` skips, nothing to do with
strand/complementarity), scans the reversed sequences with the same
trim (now correctly skipping the boundary-distant end, which after
reversal sits at the front), then maps the result back to the ORIGINAL
sequence's own coordinates. Only 2.7.1_annotate_ml_labels.py has been
switched over to use it so far (see that script and the accompanying
note in chat) -- 2.6.3_build_consensus.py's own alignment branch, and
this script's own Strategy C/D/E row-computation, still call the
side-blind function directly and have the SAME asymmetry on their
RIGHT-side numbers. That's a deliberate scope decision, not an
oversight -- rewriting those changes what cluster_consensus_results.tsv
and this script's own comparison table already show, which is worth a
separate, explicit decision rather than a silent side-effect of an ML
pipeline fix.

GAP-HANDLING (unchanged from before): every strategy shares ONE identity
function (position_identities()), which only ever drops a sequence from
a column's denominator when it's too SHORT to have a character there at
all; a gap it legitimately holds always counts as its own value in the
majority vote.

Run standalone:
    python consensus_building_evaluation.py
    python consensus_building_evaluation.py --rec-types rve,c2_n1ser
    python consensus_building_evaluation.py --project-dir /other/path --workers 32
"""

from __future__ import annotations

import argparse
import random
import re
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ══════════════════════════════════════════════════════════════════════════════
#  MGE-OVERLAP EXCLUSION -- shared by 2.6.3_build_consensus.py and
#  2.7.1_annotate_ml_labels.py (moved here from 2.7.1 so both import the
#  SAME implementation instead of maintaining separate copies that could
#  silently drift apart). See either caller's own docstring for the full
#  rationale on WHY this exclusion exists in the first place.
# ══════════════════════════════════════════════════════════════════════════════

def strip_mge_overlap(seq: str, side: str, mge_overlap: int) -> str:
    """Removes the mge_overlap-sized zone that sits INSIDE the MGE (see
    2.5.1_extract_flanks.py's own docstring for the extraction geometry
    this assumes): for LEFT, that's the LAST mge_overlap characters; for
    RIGHT, the FIRST mge_overlap characters. Applied identically wherever
    a flank or alignment sequence is about to be used for CONSENSUS
    BUILDING specifically (not for the alignment step itself, which still
    uses the full, overlap-inclusive sequence as its anchor) -- so every
    caller agrees on what "the flank" means for this purpose. Safe
    against sequences shorter than or equal to mge_overlap (returns ""
    rather than raising or wrapping around) -- the whole sequence would
    be inside the MGE in that case."""
    if mge_overlap <= 0:
        return seq
    if len(seq) <= mge_overlap:
        return ""
    return seq[:-mge_overlap] if side == "left" else seq[mge_overlap:]


# ══════════════════════════════════════════════════════════════════════════════
#  FASTA I/O
# ══════════════════════════════════════════════════════════════════════════════

def read_fasta(path: Path) -> List[str]:
    """Return list of uppercase sequences from a FASTA file."""
    seqs, buf = [], []
    with open(path, encoding="ascii", errors="replace") as fh:
        for raw in fh:
            line = raw.rstrip()
            if line.startswith(">"):
                if buf:
                    seqs.append("".join(buf).upper())
                buf = []
            elif line:
                buf.append(line)
    if buf:
        seqs.append("".join(buf).upper())
    return seqs

# ══════════════════════════════════════════════════════════════════════════════
#  SHARED IDENTITY FUNCTION -- used by ALL strategies
# ══════════════════════════════════════════════════════════════════════════════

def position_identities(sequences: List[str]) -> List[Tuple[float, str]]:
    """
    Per-column identity. A sequence only drops out of a column's
    denominator if it's too SHORT to have any character there at all
    ("ran out"); every character it DOES hold -- including '-' gaps or
    ambiguous codes -- counts as its own real, distinct value in the
    majority vote. Returns [(identity, majority_base), ...] per column.
    """
    if not sequences:
        return []
    result = []
    for col in range(max(len(s) for s in sequences)):
        bases = [s[col] for s in sequences if len(s) > col]
        if not bases:
            break
        top_base, top_cnt = Counter(bases).most_common(1)[0]
        result.append((top_cnt / len(bases), top_base))
    return result

# ══════════════════════════════════════════════════════════════════════════════
#  STRATEGY A — run-length tolerant cutoff (STACKING data)
# ══════════════════════════════════════════════════════════════════════════════

def consensus_strategy_a(
    seqs: List[str], threshold: float, run_len: int, trim: int,
) -> Tuple[str, str, int]:
    """
    Trim the first `trim` nt, then scan columns left->right. Track a run
    of consecutive columns with identity <= threshold. When the run
    reaches `run_len`, cut BEFORE the start of that run. If no such run
    is found, return the full column-wise consensus.
    Returns (consensus_str, stop_reason, consensus_len);
    stop_reason in {'empty', 'pos1_low', 'cutoff', 'maxed_out'}.
    """
    trimmed = [s[trim:] for s in seqs if len(s) > trim]
    trimmed = [s for s in trimmed if s]
    if not trimmed:
        return ("", "empty", 0)

    pos_data = position_identities(trimmed)
    if not pos_data:
        return ("", "empty", 0)

    con_bases: List[str] = []
    low_run = 0
    low_start = 0
    pos1_low = False

    for i, (identity, base) in enumerate(pos_data):
        if i == 0 and identity <= threshold:
            pos1_low = True
        con_bases.append(base)
        if identity <= threshold:
            if low_run == 0:
                low_start = i
            low_run += 1
        else:
            low_run = 0
        if low_run >= run_len:
            final = "".join(con_bases[:low_start])
            reason = "pos1_low" if pos1_low else "cutoff"
            return (final, reason, len(final))

    final = "".join(con_bases)
    reason = "pos1_low" if pos1_low else "maxed_out"
    return (final, reason, len(final))

# ══════════════════════════════════════════════════════════════════════════════
#  STRATEGY B — single-column cutoff (STACKING data)
# ══════════════════════════════════════════════════════════════════════════════

def consensus_strategy_b(seqs: List[str], threshold: float, trim: int) -> Tuple[str, str, int]:
    """Trim first `trim` nt, then cut at the first column with identity <= threshold."""
    trimmed = [s[trim:] for s in seqs if len(s) > trim]
    trimmed = [s for s in trimmed if s]
    if not trimmed:
        return ("", "empty", 0)

    pos_data = position_identities(trimmed)
    if not pos_data:
        return ("", "empty", 0)

    for col, (ident, _base) in enumerate(pos_data):
        if ident <= threshold:
            consensus = "".join(b for _, b in pos_data[:col])
            flag = "low_start" if col == 0 else "normal"
            return (consensus, flag, len(consensus))

    consensus = "".join(b for _, b in pos_data)
    return (consensus, "maxed_out", len(consensus))

# ══════════════════════════════════════════════════════════════════════════════
#  STRATEGY C — alignment-based, multi-segment (ALIGNMENT data)
# ══════════════════════════════════════════════════════════════════════════════

def scan_alignment_multi_segment(
    seqs: List[str], threshold: float, trim: int = 0,
) -> Tuple[List[Tuple[int, str]], bool]:
    """
    Scans the WHOLE alignment (not just until the first break), closing
    off every maximal run of identity >= threshold as its own separate
    consensus segment, then continuing to look for more.

    Returns (segments, maxed):
      segments : (start, text) for every CLOSED run (a run that ended
                 because some later column dropped below threshold) --
                 there can be zero, one, or several. `start` is the
                 segment's starting column in the ORIGINAL (untrimmed)
                 alignment -- i.e. already offset by `trim` -- so it's
                 directly usable as "where the consensus starts relative
                 to the start of the alignment".
      maxed    : True only if the LAST column of the alignment was still
                 part of a run that was never closed (identity stayed
                 >= threshold all the way to the end) -- that trailing,
                 unclosed run is deliberately NOT included in `segments`
                 (its true boundary is unknown). If the scan ends on a
                 column that already broke a run, maxed is False even if
                 `segments` came out empty.

    trim skips the first `trim` ALIGNMENT COLUMNS (not raw nucleotides).
    """
    trimmed = [s[trim:] for s in seqs if len(s) > trim]
    trimmed = [s for s in trimmed if s]
    if not trimmed:
        return [], False

    pos_data = position_identities(trimmed)
    if not pos_data:
        return [], False

    segments: List[Tuple[int, str]] = []
    current: List[str] = []
    in_run = False
    seg_start = 0

    for i, (ident, base) in enumerate(pos_data):
        if ident >= threshold:
            if not in_run:
                seg_start = i
            current.append(base)
            in_run = True
        else:
            if in_run and current:
                segments.append((trim + seg_start, "".join(current)))
            current = []
            in_run = False

    maxed = in_run
    return segments, maxed

# ══════════════════════════════════════════════════════════════════════════════
#  STRATEGY D — alignment-based, multi-segment, WITH run-length tolerance
#  (same idea as Strategy C, but a short run of below-threshold columns --
#  up to run_len -- doesn't immediately close a segment; only once the low
#  run itself reaches run_len does the segment close, retroactively trimmed
#  back to before that low run started. This is Strategy A's "skipping"
#  tolerance applied to Strategy C's multi-segment scan, instead of
#  Strategy A's single-cutoff scan.)
# ══════════════════════════════════════════════════════════════════════════════

def scan_alignment_multi_segment_with_skip(
    seqs: List[str], threshold: float, run_len: int, trim: int = 0,
) -> Tuple[List[Tuple[int, str]], bool]:
    """
    Like scan_alignment_multi_segment, but tolerates a run of up to
    `run_len` consecutive below-threshold columns inside an open segment
    without closing it (matching consensus_strategy_a's run-length
    tolerance). Only once a low run reaches `run_len` does the segment
    close -- retroactively trimmed back to exclude that whole low run,
    same as consensus_strategy_a's low_start bookkeeping -- and scanning
    continues afterwards to look for more segments (unlike
    consensus_strategy_a, which returns immediately on its first cut).

    Returns (segments, maxed) with the exact same semantics as
    scan_alignment_multi_segment (including start positions already
    offset by `trim`, see its docstring): `maxed` is True only if a
    segment was still open (whether cleanly above threshold or
    mid-tolerated-low-run) when the alignment ran out of columns -- that
    trailing, unresolved segment is excluded from `segments` either way,
    since its true boundary is unknown.
    """
    trimmed = [s[trim:] for s in seqs if len(s) > trim]
    trimmed = [s for s in trimmed if s]
    if not trimmed:
        return [], False

    pos_data = position_identities(trimmed)
    if not pos_data:
        return [], False

    segments: List[Tuple[int, str]] = []
    current: List[str] = []   # bases in the currently-open segment, including any not-yet-resolved low run tail
    low_run = 0
    low_start = 0              # index into `current` where the current low run began
    in_run = False              # whether a segment is currently open at all
    seg_start = 0               # column (in trimmed coordinates) where the current segment began

    for i, (ident, base) in enumerate(pos_data):
        if ident >= threshold:
            if not in_run:
                seg_start = i
            current.append(base)
            in_run = True
            low_run = 0
        else:
            if not in_run:
                # no open segment to extend or break -- nothing to do
                continue
            if low_run == 0:
                low_start = len(current)
            current.append(base)
            low_run += 1
            if low_run >= run_len:
                closed = "".join(current[:low_start])
                if closed:
                    segments.append((trim + seg_start, closed))
                current = []
                in_run = False
                low_run = 0

    maxed = in_run
    return segments, maxed


def scan_alignment_multi_segment_with_skip_side_aware(
    seqs: List[str], threshold: float, run_len: int, trim: int, side: str,
) -> Tuple[List[Tuple[int, str]], bool]:
    """
    Side-aware wrapper around scan_alignment_multi_segment_with_skip --
    see the module docstring's "BUG FIX -- trim_start was not side-aware"
    note for the full reasoning. In short: a LEFT flank's
    boundary-adjacent region sits at the END of its sequence, a RIGHT
    flank's sits at the START (re-derived directly from
    2.5.1_extract_flanks.py's coordinate arithmetic), so unconditionally
    trimming the first `trim` columns is only correct for LEFT.

    side == "left"  : calls scan_alignment_multi_segment_with_skip
                       UNCHANGED.
    side == "right" : reverses every sequence's character order (NOT
                       reverse-complemented -- purely about which end
                       `trim` skips, nothing to do with strand), scans
                       the reversed sequences with the same
                       trim/threshold/run_len, then maps every returned
                       segment's start position AND text back into the
                       ORIGINAL sequence's own coordinates.

    Returns the same (segments, maxed) shape as the function it wraps,
    always in ORIGINAL sequence coordinates regardless of side. Only
    CLOSED segments are returned here -- see
    scan_alignment_multi_segment_with_skip_full_side_aware() below if you
    also need the trailing, unclosed run when maxed is True.
    """
    if side == "left":
        return scan_alignment_multi_segment_with_skip(seqs, threshold, run_len, trim)
    if side != "right":
        raise ValueError(f"side must be 'left' or 'right', got {side!r}")

    seq_length = max((len(s) for s in seqs), default=0)
    reversed_seqs = [s[::-1] for s in seqs]
    rev_segments, maxed = scan_alignment_multi_segment_with_skip(reversed_seqs, threshold, run_len, trim)

    orig_segments: List[Tuple[int, str]] = []
    for rev_start, rev_text in rev_segments:
        orig_start = seq_length - rev_start - len(rev_text)
        orig_text = rev_text[::-1]
        orig_segments.append((orig_start, orig_text))
    return orig_segments, maxed


def scan_alignment_multi_segment_with_skip_trailing(
    seqs: List[str], threshold: float, run_len: int, trim: int = 0,
) -> Tuple[List[Tuple[int, str]], bool, Optional[Tuple[int, str]]]:
    """
    Same scan as scan_alignment_multi_segment_with_skip, but ALSO
    returns the trailing, unclosed run (start, text) when maxed is True
    -- None when maxed is False (nothing was left open at the end).

    WHY THIS EXISTS (see chat): the boundary sits at the literal EDGE of
    every extracted flank window, on whichever side is boundary-adjacent
    (both, by construction of 2.5.1_extract_flanks.py -- see
    scan_alignment_multi_segment_with_skip_side_aware()'s docstring).
    A genuinely conserved boundary motif therefore very often reaches
    the end of the alignment while still "open" -- there's nothing
    beyond it in the extracted window to trigger a low-identity column
    and formally CLOSE the run -- and scan_alignment_multi_segment_with_
    skip()'s `segments` list, by design, never includes an unclosed run
    (its true full extent, beyond the window, is genuinely unknown).
    That's the right call for "segments" specifically, but it means a
    caller that only ever looks at `segments` will systematically miss
    real boundary-adjacent conservation whenever it happens to run all
    the way to the edge -- plausibly a large share of why so many real
    clusters here come back maxed_out (829 of 1,054 -- 78.7% -- in the
    six-family evaluation run this pipeline has already produced). This
    function exposes that trailing run explicitly so a caller CAN choose
    to treat it as evidence, instead of it being silently invisible.

    This is a genuinely separate concern from the LEFT/RIGHT trim
    asymmetry fix above -- it affects BOTH sides equally, and existed
    before that fix (it just happened to be masked by the RIGHT side
    already being broken in a different, more severe way). Not wired
    into any other function by default; callers that want it use this
    one explicitly.
    """
    trimmed = [s[trim:] for s in seqs if len(s) > trim]
    trimmed = [s for s in trimmed if s]
    if not trimmed:
        return [], False, None

    pos_data = position_identities(trimmed)
    if not pos_data:
        return [], False, None

    segments: List[Tuple[int, str]] = []
    current: List[str] = []
    low_run = 0
    low_start = 0
    in_run = False
    seg_start = 0

    for i, (ident, base) in enumerate(pos_data):
        if ident >= threshold:
            if not in_run:
                seg_start = i
            current.append(base)
            in_run = True
            low_run = 0
        else:
            if not in_run:
                continue
            if low_run == 0:
                low_start = len(current)
            current.append(base)
            low_run += 1
            if low_run >= run_len:
                closed = "".join(current[:low_start])
                if closed:
                    segments.append((trim + seg_start, closed))
                current = []
                in_run = False
                low_run = 0

    maxed = in_run
    trailing = (trim + seg_start, "".join(current)) if maxed and current else None
    return segments, maxed, trailing


def scan_alignment_multi_segment_with_skip_full_side_aware(
    seqs: List[str], threshold: float, run_len: int, trim: int, side: str,
) -> Tuple[List[Tuple[int, str]], bool, Optional[Tuple[int, str]]]:
    """
    Combines scan_alignment_multi_segment_with_skip_side_aware() (fixes
    the LEFT/RIGHT trim asymmetry) with
    scan_alignment_multi_segment_with_skip_trailing() (surfaces the
    trailing unclosed run instead of discarding it) -- both fixes
    together, in original sequence coordinates regardless of side.
    """
    if side == "left":
        return scan_alignment_multi_segment_with_skip_trailing(seqs, threshold, run_len, trim)
    if side != "right":
        raise ValueError(f"side must be 'left' or 'right', got {side!r}")

    seq_length = max((len(s) for s in seqs), default=0)
    reversed_seqs = [s[::-1] for s in seqs]
    rev_segments, maxed, rev_trailing = scan_alignment_multi_segment_with_skip_trailing(
        reversed_seqs, threshold, run_len, trim,
    )

    orig_segments: List[Tuple[int, str]] = []
    for rev_start, rev_text in rev_segments:
        orig_start = seq_length - rev_start - len(rev_text)
        orig_segments.append((orig_start, rev_text[::-1]))

    orig_trailing = None
    if rev_trailing is not None:
        rev_start, rev_text = rev_trailing
        orig_start = seq_length - rev_start - len(rev_text)
        orig_trailing = (orig_start, rev_text[::-1])

    return orig_segments, maxed, orig_trailing

# ══════════════════════════════════════════════════════════════════════════════
#  DIRECTORY SCANNING
# ══════════════════════════════════════════════════════════════════════════════

GroupKey = Tuple[str, str, str]   # (rec_type, cluster_num, dup_group)
ClusterKey = Tuple[str, str]      # (rec_type, cluster_num)


def _strip_known_suffix(name: str, suffixes: Tuple[str, ...]) -> Optional[str]:
    for suf in suffixes:
        if name.endswith(suf):
            return name[: -len(suf)]
    return None


def scan_grouped_dir(
    root: Path, rec_types_filter: Optional[set], suffixes: Tuple[str, ...],
    io_workers: int,
) -> Dict[GroupKey, Dict[str, List[str]]]:
    """
    LEGACY scanner (see module docstring's BUG FIX note) -- assumes
    group-level files sit directly inside cluster_N/, one level
    shallower than the real thesis_rep tree actually nests. Left in
    place, unused by default, only for pointing --stacking-dir/
    --alignments-dir back at a tree shaped like this on purpose.
        <root>/<rec_type>/cluster_<N>/<rec_type>_<N>_<group>_<side>{suffix}.
    Returns {(rec_type, cluster_num, dup_group): {"left": [...], "right": [...]}}.
    """
    file_tasks: List[Tuple[GroupKey, str, Path]] = []

    for rt_dir in sorted(root.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue

        for cl_dir in sorted(rt_dir.iterdir()):
            if not cl_dir.is_dir():
                continue
            m = re.fullmatch(r"cluster_(\d+)", cl_dir.name)
            if not m:
                continue
            cluster_num = m.group(1)

            for fp in sorted(cl_dir.iterdir()):
                if not fp.is_file():
                    continue
                stem = _strip_known_suffix(fp.name, suffixes)
                if stem is None:
                    continue
                parts = stem.rsplit("_", 2)
                if len(parts) < 3 or parts[-1] not in ("left", "right"):
                    continue
                dup_group, side = parts[-2], parts[-1]
                file_tasks.append(((rec_type, cluster_num, dup_group), side, fp))

    def _read_one(task):
        key, side, fp = task
        return key, side, read_fasta(fp)

    universe: Dict[GroupKey, Dict[str, List[str]]] = defaultdict(lambda: {"left": [], "right": []})
    with ThreadPoolExecutor(max_workers=io_workers) as pool:
        for key, side, seqs in pool.map(_read_one, file_tasks):
            universe[key][side] = seqs

    return dict(universe)


def scan_duplicate_group_level_dir(
    root: Path, rec_types_filter: Optional[set], suffix: str, io_workers: int,
) -> Dict[GroupKey, Dict[str, List[str]]]:
    """
    Scanner for 2.6.1_align_flanks.sh's REAL duplicate_group_level output
    tree:
        <root>/<rec_type>/cluster_<N>/duplicate_group_<M>/left<suffix>, right<suffix>
    2.6.1 never writes anything for singleton/ (there's no second copy
    to align there), so this never needs to look for or exclude a
    "singleton" folder -- unlike scan_flanks_tree_for_universe_totals()
    below, which scans a tree that DOES have singleton/ folders sitting
    right next to duplicate_group_*/ ones.
    Returns {(rec_type, cluster_num, dup_group): {"left": [...], "right": [...]}}.
    """
    file_tasks: List[Tuple[GroupKey, str, Path]] = []

    for rt_dir in sorted(root.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue

        for cl_dir in sorted(rt_dir.iterdir()):
            if not cl_dir.is_dir():
                continue
            m = re.fullmatch(r"cluster_(\d+)", cl_dir.name)
            if not m:
                continue
            cluster_num = m.group(1)

            for group_dir in sorted(cl_dir.glob("duplicate_group_*")):
                if not group_dir.is_dir():
                    continue
                gm = re.fullmatch(r"duplicate_group_(\d+)", group_dir.name)
                if not gm:
                    continue
                dup_group = gm.group(1)

                for side in ("left", "right"):
                    fp = group_dir / f"{side}{suffix}"
                    if fp.is_file():
                        file_tasks.append(((rec_type, cluster_num, dup_group), side, fp))

    def _read_one(task):
        key, side, fp = task
        return key, side, read_fasta(fp)

    universe: Dict[GroupKey, Dict[str, List[str]]] = defaultdict(lambda: {"left": [], "right": []})
    with ThreadPoolExecutor(max_workers=io_workers) as pool:
        for key, side, seqs in pool.map(_read_one, file_tasks):
            universe[key][side] = seqs

    return dict(universe)


def scan_flanks_tree_for_universe_totals(
    root: Path, rec_types_filter: Optional[set], io_workers: int,
) -> Dict[GroupKey, Dict[str, List[str]]]:
    """
    Scanner for the REAL thesis_rep flanks tree -- 2.5.1_extract_flanks.py's
    per-MGE raw flank output, used ONLY for this script's "universe
    totals" (how many MGEs/clusters/duplicate groups are really out
    there, as a fixed reference count every row's percentages divide
    into -- see module docstring):
        <root>/<rec_type>/cluster_<N>/duplicate_group_<M>/<mge_id>_left.fasta, _right.fasta
        <root>/<rec_type>/cluster_<N>/singleton/<mge_id>_left.fasta, _right.fasta
    ONLY duplicate_group_* folders are scanned -- singleton/ is
    deliberately excluded (matching every other place in this pipeline
    that treats a singleton as having no duplicate to compare against:
    2.6.1, find_terminal_inverted_repeats.py's multi-copy branch, etc.).
    No consensus strategy here can EVER find a consensus for a singleton
    in the first place, so counting them into "how many MGEs a strategy
    could possibly find a consensus for" would make every
    pct_mges_consensus_found number look artificially low against an
    unreachable ceiling, not more complete.

    Unlike scan_duplicate_group_level_dir() above, a single duplicate
    group here is made of MANY per-MGE files (one _left.fasta/_right.fasta
    pair per member), not one pre-concatenated stack -- so this collects
    (extends) sequences from every member file into that group's list,
    rather than assigning from a single file per key/side.
    Returns {(rec_type, cluster_num, dup_group): {"left": [...], "right": [...]}}.
    """
    file_tasks: List[Tuple[GroupKey, str, Path]] = []

    for rt_dir in sorted(root.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue

        for cl_dir in sorted(rt_dir.iterdir()):
            if not cl_dir.is_dir():
                continue
            m = re.fullmatch(r"cluster_(\d+)", cl_dir.name)
            if not m:
                continue
            cluster_num = m.group(1)

            for group_dir in sorted(cl_dir.glob("duplicate_group_*")):
                if not group_dir.is_dir():
                    continue
                gm = re.fullmatch(r"duplicate_group_(\d+)", group_dir.name)
                if not gm:
                    continue
                dup_group = gm.group(1)

                for fp in sorted(group_dir.iterdir()):
                    if not fp.is_file():
                        continue
                    if fp.name.endswith("_left.fasta"):
                        file_tasks.append(((rec_type, cluster_num, dup_group), "left", fp))
                    elif fp.name.endswith("_right.fasta"):
                        file_tasks.append(((rec_type, cluster_num, dup_group), "right", fp))

    def _read_one(task):
        key, side, fp = task
        return key, side, read_fasta(fp)

    universe: Dict[GroupKey, Dict[str, List[str]]] = defaultdict(lambda: {"left": [], "right": []})
    with ThreadPoolExecutor(max_workers=io_workers) as pool:
        for key, side, seqs in pool.map(_read_one, file_tasks):
            universe[key][side].extend(seqs)  # extend, not assign -- multiple files share one key/side here

    return dict(universe)


def scan_cluster_level_dir(
    root: Path, rec_types_filter: Optional[set], io_workers: int,
    suffix: str = ".aln.fasta",
) -> Dict[ClusterKey, Dict[str, List[str]]]:
    """
    alignments/cluster_level/<rec_type>/<rec_type>_cluster_<N>_<left|right><suffix>
    Two file kinds live side by side in the same folders: *.aln.fasta and
    *.raw.fasta. Only ONE suffix is read per call.
    Returns {(rec_type, cluster_num): {"left": [...], "right": [...]}}.
    """
    file_tasks: List[Tuple[ClusterKey, str, Path]] = []

    for rt_dir in sorted(root.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue

        prefix = f"{rec_type}_cluster_"
        for fp in sorted(rt_dir.iterdir()):
            if not fp.is_file() or not fp.name.endswith(suffix):
                continue
            if not fp.name.startswith(prefix):
                continue
            stem = fp.name[: -len(suffix)]
            remainder = stem[len(prefix):]          # "<N>_left" or "<N>_right"
            parts = remainder.rsplit("_", 1)
            if len(parts) != 2 or parts[1] not in ("left", "right"):
                continue
            cluster_num, side = parts
            file_tasks.append(((rec_type, cluster_num), side, fp))

    def _read_one(task):
        key, side, fp = task
        return key, side, read_fasta(fp)

    universe: Dict[ClusterKey, Dict[str, List[str]]] = defaultdict(lambda: {"left": [], "right": []})
    with ThreadPoolExecutor(max_workers=io_workers) as pool:
        for key, side, seqs in pool.map(_read_one, file_tasks):
            universe[key][side] = seqs

    return dict(universe)

# ══════════════════════════════════════════════════════════════════════════════
#  PER-UNIT WORKERS (run in a process pool)
# ══════════════════════════════════════════════════════════════════════════════

def _worker_stacking(task):
    """One duplicate group OR cluster, Strategy A or B.
    Returns (key, n_members, group_maxed, n_found, total_len_found, n_in_range, n_under_4, n_over_20)."""
    key, left_seqs, right_seqs, strategy, threshold, run_len, trim = task
    n_members = max(len(left_seqs), len(right_seqs))
    if strategy == "A":
        _, left_reason, left_len = consensus_strategy_a(left_seqs, threshold, run_len, trim)
        _, right_reason, right_len = consensus_strategy_a(right_seqs, threshold, run_len, trim)
    else:
        _, left_reason, left_len = consensus_strategy_b(left_seqs, threshold, trim)
        _, right_reason, right_len = consensus_strategy_b(right_seqs, threshold, trim)
    group_maxed = (left_reason == "maxed_out") or (right_reason == "maxed_out")

    n_found = 0
    total_len_found = 0
    n_in_range = n_under_4 = n_over_20 = 0
    for reason, length in ((left_reason, left_len), (right_reason, right_len)):
        if reason != "maxed_out":
            n_found += 1
            total_len_found += length
            if 4 < length < 20:
                n_in_range += 1
            elif length < 4:
                n_under_4 += 1
            elif length > 20:
                n_over_20 += 1

    return key, n_members, group_maxed, n_found, total_len_found, n_in_range, n_under_4, n_over_20


def _worker_alignment(task):
    """One duplicate group OR cluster, Strategy C.
    Returns (key, n_members, unit_maxed, multi, n_segments, total_length,
              n_in_range, n_under_4, n_over_20)."""
    key, left_seqs, right_seqs, threshold, trim = task
    n_members = max(len(left_seqs), len(right_seqs))
    left_segs, left_maxed = scan_alignment_multi_segment(left_seqs, threshold, trim)
    right_segs, right_maxed = scan_alignment_multi_segment(right_seqs, threshold, trim)
    unit_maxed = left_maxed or right_maxed
    multi = len(left_segs) > 1 or len(right_segs) > 1
    all_segs = [text for _start, text in left_segs] + [text for _start, text in right_segs]
    n_segments = len(all_segs)
    total_length = sum(len(s) for s in all_segs)
    n_in_range = sum(1 for s in all_segs if 4 < len(s) < 20)
    n_under_4 = sum(1 for s in all_segs if len(s) < 4)
    n_over_20 = sum(1 for s in all_segs if len(s) > 20)
    return key, n_members, unit_maxed, multi, n_segments, total_length, n_in_range, n_under_4, n_over_20


def _worker_alignment_skip(task):
    """One duplicate group OR cluster, Strategy D (Strategy C + run-length
    tolerance) -- OR Strategy E when the caller resolved run_len
    dynamically per unit before building the task tuple; this worker
    doesn't know or care which, it just consumes whatever run_len it's
    given. Same return shape as _worker_alignment.
    task = (key, left_seqs, right_seqs, threshold, run_len, trim)."""
    key, left_seqs, right_seqs, threshold, run_len, trim = task
    n_members = max(len(left_seqs), len(right_seqs))
    left_segs, left_maxed = scan_alignment_multi_segment_with_skip(left_seqs, threshold, run_len, trim)
    right_segs, right_maxed = scan_alignment_multi_segment_with_skip(right_seqs, threshold, run_len, trim)
    unit_maxed = left_maxed or right_maxed
    multi = len(left_segs) > 1 or len(right_segs) > 1
    all_segs = [text for _start, text in left_segs] + [text for _start, text in right_segs]
    n_segments = len(all_segs)
    total_length = sum(len(s) for s in all_segs)
    n_in_range = sum(1 for s in all_segs if 4 < len(s) < 20)
    n_under_4 = sum(1 for s in all_segs if len(s) < 4)
    n_over_20 = sum(1 for s in all_segs if len(s) > 20)
    return key, n_members, unit_maxed, multi, n_segments, total_length, n_in_range, n_under_4, n_over_20

# ══════════════════════════════════════════════════════════════════════════════
#  ROW COMPUTATION
# ══════════════════════════════════════════════════════════════════════════════

NA = "NA"


def _length_stats(raw: List[Tuple], found_idx: int, len_idx: int, inrange_idx: int,
                   under4_idx: int, over20_idx: int) -> dict:
    total_found = sum(r[found_idx] for r in raw)
    total_len = sum(r[len_idx] for r in raw)
    total_in_range = sum(r[inrange_idx] for r in raw)
    total_under_4 = sum(r[under4_idx] for r in raw)
    total_over_20 = sum(r[over20_idx] for r in raw)
    return {
        "avg_consensus_length": (total_len / total_found) if total_found else None,
        "pct_length_4_to_20": (total_in_range / total_found * 100) if total_found else None,
        "pct_length_under_4": (total_under_4 / total_found * 100) if total_found else None,
        "pct_length_over_20": (total_over_20 / total_found * 100) if total_found else None,
    }


def compute_stacking_row(
    universe: Dict, strategy: str, threshold: float,
    run_len: int, trim: int, pool: ProcessPoolExecutor, chunksize: int,
) -> dict:
    """FIXED threshold: same value for every unit in the row."""
    tasks = [
        (key, data["left"], data["right"], strategy, threshold, run_len, trim)
        for key, data in universe.items()
    ]
    raw = list(pool.map(_worker_stacking, tasks, chunksize=chunksize))
    results = [(key, n, maxed) for key, n, maxed, *_ in raw]
    stats = _aggregate(results)
    stats.update(_length_stats(raw, 3, 4, 5, 6, 7))
    return stats


def compute_stacking_row_dynamic(
    universe: Dict, strategy: str, run_len: int, trim: int,
    small_thr: float, large_thr: float, cutoff: int,
    pool: ProcessPoolExecutor, chunksize: int,
) -> dict:
    """DYNAMIC threshold: resolved PER UNIT from its own sequence count,
    before tasks are built -- _worker_stacking itself needs no changes."""
    tasks = []
    for key, data in universe.items():
        n_members = max(len(data["left"]), len(data["right"]))
        thr = small_thr if n_members <= cutoff else large_thr
        tasks.append((key, data["left"], data["right"], strategy, thr, run_len, trim))
    raw = list(pool.map(_worker_stacking, tasks, chunksize=chunksize))
    results = [(key, n, maxed) for key, n, maxed, *_ in raw]
    stats = _aggregate(results)
    stats.update(_length_stats(raw, 3, 4, 5, 6, 7))
    return stats


def compute_alignment_row(
    universe: Dict, threshold: float, trim: int, pool: ProcessPoolExecutor, chunksize: int,
) -> Tuple[dict, dict, int, Dict[Tuple[str, str], int]]:
    """FIXED threshold, Strategy C. Returns (aggregated_stats, multi_consensus_info, total_segments, segments_by_cluster)."""
    tasks = [(key, data["left"], data["right"], threshold, trim) for key, data in universe.items()]
    return _compute_alignment_row_from_tasks(tasks, pool, chunksize, _worker_alignment)


def compute_alignment_row_dynamic(
    universe: Dict, trim: int, small_thr: float, large_thr: float, cutoff: int,
    pool: ProcessPoolExecutor, chunksize: int,
) -> Tuple[dict, dict, int, Dict[Tuple[str, str], int]]:
    """DYNAMIC threshold, Strategy C, resolved per unit -- same pattern as the stacking version."""
    tasks = []
    for key, data in universe.items():
        n_members = max(len(data["left"]), len(data["right"]))
        thr = small_thr if n_members <= cutoff else large_thr
        tasks.append((key, data["left"], data["right"], thr, trim))
    return _compute_alignment_row_from_tasks(tasks, pool, chunksize, _worker_alignment)


def compute_alignment_row_skip(
    universe: Dict, threshold: float, run_len: int, trim: int, pool: ProcessPoolExecutor, chunksize: int,
) -> Tuple[dict, dict, int, Dict[Tuple[str, str], int]]:
    """FIXED threshold, Strategy D (Strategy C + run-length tolerance)."""
    tasks = [(key, data["left"], data["right"], threshold, run_len, trim) for key, data in universe.items()]
    return _compute_alignment_row_from_tasks(tasks, pool, chunksize, _worker_alignment_skip)


def compute_alignment_row_dynamic_skip(
    universe: Dict, run_len: int, trim: int, small_thr: float, large_thr: float, cutoff: int,
    pool: ProcessPoolExecutor, chunksize: int,
) -> Tuple[dict, dict, int, Dict[Tuple[str, str], int]]:
    """DYNAMIC threshold, Strategy D, resolved per unit."""
    tasks = []
    for key, data in universe.items():
        n_members = max(len(data["left"]), len(data["right"]))
        thr = small_thr if n_members <= cutoff else large_thr
        tasks.append((key, data["left"], data["right"], thr, run_len, trim))
    return _compute_alignment_row_from_tasks(tasks, pool, chunksize, _worker_alignment_skip)


def compute_alignment_row_dynamic_nt_skip(
    universe: Dict, threshold: float, trim: int, cutoff: int, skip_run_len: int,
    pool: ProcessPoolExecutor, chunksize: int,
) -> Tuple[dict, dict, int, Dict[Tuple[str, str], int]]:
    """Strategy E: threshold is FIXED (unlike the dynamic-THRESHOLD rows
    above) at 100% same as Strategy C/D, but run_len (the skip tolerance
    itself) is what's resolved PER UNIT this time -- n_members <= cutoff
    -> run_len=1 (no tolerance at all, identical behaviour to Strategy C),
    otherwise run_len=skip_run_len (Strategy D's tolerance). Reuses
    _worker_alignment_skip directly: run_len=1 there is already exactly
    equivalent to the strict scan (a single low column immediately
    reaches "low_run >= 1" and closes), so no separate strict worker is
    needed for the n_members<=cutoff case."""
    tasks = []
    for key, data in universe.items():
        n_members = max(len(data["left"]), len(data["right"]))
        run_len = 1 if n_members <= cutoff else skip_run_len
        tasks.append((key, data["left"], data["right"], threshold, run_len, trim))
    return _compute_alignment_row_from_tasks(tasks, pool, chunksize, _worker_alignment_skip)


def _compute_alignment_row_from_tasks(tasks, pool: ProcessPoolExecutor, chunksize: int, worker_fn):
    """Shared by Strategy C (worker_fn=_worker_alignment) and Strategy D
    (worker_fn=_worker_alignment_skip) -- both workers return the exact
    same tuple shape, so the aggregation logic is identical either way."""
    raw = list(pool.map(worker_fn, tasks, chunksize=chunksize))
    results = [(key, n, maxed) for key, n, maxed, *_ in raw]
    multi_keys = [key for key, _n, _maxed, multi, *_ in raw if multi]
    stats = _aggregate(results)
    stats.update(_length_stats(raw, 4, 5, 6, 7, 8))  # n_segments, total_length, n_in_range, n_under_4, n_over_20

    total_segments = sum(r[4] for r in raw)
    multi_info = {"n_units": len(multi_keys), "all_keys": multi_keys}

    segments_by_cluster: Dict[Tuple[str, str], int] = defaultdict(int)
    for key, _n, _maxed, _multi, nseg, *_ in raw:
        segments_by_cluster[(key[0], key[1])] += nseg

    return stats, multi_info, total_segments, dict(segments_by_cluster)


def _aggregate(results: List[Tuple]) -> dict:
    """
    Shared aggregation for both stacking and alignment rows.
    results: [(key, n_members, unit_maxed), ...] where key is either a
    GroupKey (rec_type, cluster, dup_group) or a ClusterKey (rec_type, cluster).
    """
    n_mges_found = n_mges_maxed = 0
    n_units_found = n_units_maxed = 0
    cluster_any_found: Dict[Tuple[str, str], bool] = defaultdict(bool)
    cluster_all_maxed: Dict[Tuple[str, str], bool] = {}

    for key, n_members, unit_maxed in results:
        rec_type, cluster_num = key[0], key[1]
        cluster_key = (rec_type, cluster_num)

        if unit_maxed:
            n_mges_maxed += n_members
            n_units_maxed += 1
        else:
            n_mges_found += n_members
            n_units_found += 1
            cluster_any_found[cluster_key] = True

        if cluster_key not in cluster_all_maxed:
            cluster_all_maxed[cluster_key] = True
        if not unit_maxed:
            cluster_all_maxed[cluster_key] = False

    n_clusters_all_maxed = sum(1 for v in cluster_all_maxed.values() if v)
    n_clusters_any_found = sum(1 for v in cluster_any_found.values() if v)

    return {
        "n_mges_found": n_mges_found,
        "n_mges_maxed": n_mges_maxed,
        "n_clusters_all_maxed": n_clusters_all_maxed,
        "n_clusters_any_found": n_clusters_any_found,
        "n_units_found": n_units_found,
        "n_units_maxed": n_units_maxed,
    }

# ══════════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════════

def parse_rec_types(arg: Optional[str]) -> Optional[set]:
    if arg is None or arg.strip().lower() == "all":
        return None
    return {t.strip() for t in arg.split(",") if t.strip()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--project-dir", default="/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep")
    ap.add_argument("--stacking-dir", default=None, help="Default: <project-dir>/flanks_split_acc_to_duplicates "
                     "(used ONLY for universe totals now -- no rows are generated from it)")
    ap.add_argument("--alignments-dir", default=None, help="Default: <project-dir>/consensus")
    ap.add_argument("--eval-dir", default=None, help="Default: <project-dir>/files/cons_evaluation")
    ap.add_argument("--rec-types", default="all",
                     help='"all" (default), or a comma-separated subset matching rec_type '
                          'folder names exactly, e.g. "rve,c2_n1ser"')
    ap.add_argument("--a-thresholds", default="0.9,0.7,0.5")
    ap.add_argument("--a-run", type=int, default=2, help='Strategy A "skipping" run length (default 2)')
    ap.add_argument("--b-thresholds", default="0.9,0.7,0.5")
    ap.add_argument("--c-threshold", type=float, default=1.0, help="Strategy C/D identity threshold (default 100%%)")
    ap.add_argument("--d-run", type=int, default=2, help='Strategy D "skipping" run length (default 2)')
    ap.add_argument("--e-cutoff", type=int, default=None,
                     help="Strategy E: n_sequences <= this -> no skipping, else skipping allowed "
                          "(default: same as --dynamic-cutoff)")
    ap.add_argument("--dynamic-cutoff", type=int, default=3, help="n_sequences <= this -> dynamic-small threshold")
    ap.add_argument("--dynamic-small", type=float, default=0.9)
    ap.add_argument("--dynamic-large", type=float, default=0.7)
    ap.add_argument("--trim-start", type=int, default=10, help="nt/columns trimmed from the 5' end of every sequence")
    ap.add_argument("--sample-size", type=int, default=10, help="how many multi-consensus examples to show per level")
    ap.add_argument("--sample-seed", type=int, default=42, help="fixed seed for the representative sample")
    ap.add_argument("--workers", type=int, default=None, help="Process pool size (default: os.cpu_count())")
    ap.add_argument("--io-workers", type=int, default=64, help="Thread pool size for reading files")
    args = ap.parse_args()

    import os
    n_workers = args.workers or os.cpu_count() or 4
    if args.e_cutoff is None:
        args.e_cutoff = args.dynamic_cutoff

    project_dir = Path(args.project_dir)
    stacking_dir = Path(args.stacking_dir) if args.stacking_dir else project_dir / "flanks_split_acc_to_duplicates"
    alignments_dir = Path(args.alignments_dir) if args.alignments_dir else project_dir / "consensus"
    dup_group_align_dir = alignments_dir / "duplicate_group_level"
    cluster_align_dir = alignments_dir / "cluster_level"
    eval_dir = Path(args.eval_dir) if args.eval_dir else project_dir / "files/cons_evaluation"
    eval_dir.mkdir(parents=True, exist_ok=True)

    rec_types_filter = parse_rec_types(args.rec_types)
    subset_prefix = "" if rec_types_filter is None else "subset_" + "_".join(sorted(rec_types_filter)) + "_"

    table_out = eval_dir / f"{subset_prefix}consensus_strategy_comparison.tsv"
    log_out = eval_dir / f"{subset_prefix}consensus_strategy_comparison.log"

    for p, label in ((stacking_dir, "stacking (universe totals)"), (dup_group_align_dir, "duplicate_group_level alignments"),
                      (cluster_align_dir, "cluster_level alignments/stackings")):
        if not p.exists():
            sys.exit(f"{label} directory not found: {p}")

    print(f"rec_types filter: {'ALL' if rec_types_filter is None else sorted(rec_types_filter)}", flush=True)

    print("Scanning stacking directory (universe totals ONLY -- no rows generated from it) ...", flush=True)
    stacking_universe = scan_flanks_tree_for_universe_totals(stacking_dir, rec_types_filter, args.io_workers)
    print(f"  {len(stacking_universe):,} duplicate groups", flush=True)

    print("Scanning duplicate_group_level alignments (Strategy C/D, group-level, fixed) ...", flush=True)
    dup_align_universe = scan_duplicate_group_level_dir(dup_group_align_dir, rec_types_filter, ".aln.fasta", args.io_workers)
    print(f"  {len(dup_align_universe):,} duplicate groups", flush=True)

    print("Scanning cluster_level alignments (Strategy C/D, cluster-level) ...", flush=True)
    cluster_align_universe = scan_cluster_level_dir(cluster_align_dir, rec_types_filter, args.io_workers,
                                                      suffix=".aln.fasta")
    print(f"  {len(cluster_align_universe):,} clusters", flush=True)

    print("Scanning cluster_level RAW stackings (Strategy A/B, cluster-level) ...", flush=True)
    cluster_raw_universe = scan_cluster_level_dir(cluster_align_dir, rec_types_filter, args.io_workers,
                                                    suffix=".raw.fasta")
    print(f"  {len(cluster_raw_universe):,} clusters", flush=True)

    # ── universe totals (fixed across every row, computed from the FULL
    # duplicate-group-level tree regardless of which rows actually use it) ──
    total_mges_in = sum(max(len(d["left"]), len(d["right"])) for d in stacking_universe.values())
    total_dup_groups_in = len(stacking_universe)
    total_clusters_in = len({(k[0], k[1]) for k in stacking_universe.keys()})
    print(f"Universe totals: {total_mges_in:,} MGEs, {total_clusters_in:,} clusters, "
          f"{total_dup_groups_in:,} duplicate groups", flush=True)

    a_thresholds = [float(x) for x in args.a_thresholds.split(",")]
    b_thresholds = [float(x) for x in args.b_thresholds.split(",")]

    rows = []
    multi_consensus_report = {}

    with ProcessPoolExecutor(max_workers=n_workers) as pool:
        chunksize_cluster_raw = max(1, len(cluster_raw_universe) // (n_workers * 4)) if cluster_raw_universe else 1

        for thr in a_thresholds:
            label = f"Strategy A {thr} + {args.a_run}nt skipping on cluster-level stackings"
            print(f"  Running: {label} ...", flush=True)
            stats = compute_stacking_row(cluster_raw_universe, "A", thr, args.a_run, args.trim_start,
                                          pool, chunksize_cluster_raw)
            rows.append({"strategy": label, "universe": "stacking_cluster_level", **stats})

        for thr in b_thresholds:
            label = f"Strategy B {thr} + no skipping on cluster-level stackings"
            print(f"  Running: {label} ...", flush=True)
            stats = compute_stacking_row(cluster_raw_universe, "B", thr, 1, args.trim_start,
                                          pool, chunksize_cluster_raw)
            rows.append({"strategy": label, "universe": "stacking_cluster_level", **stats})

        # dynamic-threshold rows (Strategy A, B -- cluster-level)
        label = (f"Strategy A DYNAMIC (<= {args.dynamic_cutoff} seqs: {args.dynamic_small}, "
                 f"else {args.dynamic_large}) + {args.a_run}nt skipping on cluster-level stackings")
        print(f"  Running: {label} ...", flush=True)
        stats = compute_stacking_row_dynamic(cluster_raw_universe, "A", args.a_run, args.trim_start,
                                              args.dynamic_small, args.dynamic_large, args.dynamic_cutoff,
                                              pool, chunksize_cluster_raw)
        rows.append({"strategy": label, "universe": "stacking_cluster_level", **stats})

        label = (f"Strategy B DYNAMIC (<= {args.dynamic_cutoff} seqs: {args.dynamic_small}, "
                 f"else {args.dynamic_large}) + no skipping on cluster-level stackings")
        print(f"  Running: {label} ...", flush=True)
        stats = compute_stacking_row_dynamic(cluster_raw_universe, "B", 1, args.trim_start,
                                              args.dynamic_small, args.dynamic_large, args.dynamic_cutoff,
                                              pool, chunksize_cluster_raw)
        rows.append({"strategy": label, "universe": "stacking_cluster_level", **stats})

        chunksize_dup = max(1, len(dup_align_universe) // (n_workers * 4)) if dup_align_universe else 1
        label = f"Strategy C {int(args.c_threshold*100)} + no skipping on alignments group-level"
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_dup, segments_by_cluster_dup = compute_alignment_row(
            dup_align_universe, args.c_threshold, args.trim_start, pool, chunksize_dup)
        stats["avg_consensuses_per_cluster"] = (
            sum(segments_by_cluster_dup.values()) / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = (
            total_segments_dup / total_dup_groups_in if total_dup_groups_in else None
        )
        rows.append({"strategy": label, "universe": "alignment_group", **stats})
        multi_consensus_report["group_level"] = multi_info

        chunksize_clu = max(1, len(cluster_align_universe) // (n_workers * 4)) if cluster_align_universe else 1
        label = f"Strategy C {int(args.c_threshold*100)} + no skipping on alignments cluster-level"
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_clu, _sbc = compute_alignment_row(
            cluster_align_universe, args.c_threshold, args.trim_start, pool, chunksize_clu)
        stats["avg_consensuses_per_cluster"] = (
            total_segments_clu / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = None
        rows.append({"strategy": label, "universe": "alignment_cluster", **stats})
        multi_consensus_report["cluster_level"] = multi_info

        # dynamic-threshold Strategy C, cluster-level alignment
        label = (f"Strategy C DYNAMIC (<= {args.dynamic_cutoff} seqs: {args.dynamic_small}, "
                 f"else {args.dynamic_large}) on alignments cluster-level")
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_clu_dyn, _sbc2 = compute_alignment_row_dynamic(
            cluster_align_universe, args.trim_start, args.dynamic_small, args.dynamic_large, args.dynamic_cutoff,
            pool, chunksize_clu)
        stats["avg_consensuses_per_cluster"] = (
            total_segments_clu_dyn / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = None
        rows.append({"strategy": label, "universe": "alignment_cluster", **stats})
        multi_consensus_report["cluster_level_dynamic"] = multi_info

        # ── Strategy D -- same 3 universes as Strategy C above (group-level
        # fixed, cluster-level fixed, cluster-level dynamic), but
        # scan_alignment_multi_segment_with_skip instead of the strict
        # single-column-cutoff scan, tolerating args.d_run nt of low
        # identity before actually closing a segment. +3 rows total.
        label = f"Strategy D {int(args.c_threshold*100)} + {args.d_run}nt skipping on alignments group-level"
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_dup_d, segments_by_cluster_dup_d = compute_alignment_row_skip(
            dup_align_universe, args.c_threshold, args.d_run, args.trim_start, pool, chunksize_dup)
        stats["avg_consensuses_per_cluster"] = (
            sum(segments_by_cluster_dup_d.values()) / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = (
            total_segments_dup_d / total_dup_groups_in if total_dup_groups_in else None
        )
        rows.append({"strategy": label, "universe": "alignment_group_d", **stats})
        multi_consensus_report["group_level_d"] = multi_info

        label = f"Strategy D {int(args.c_threshold*100)} + {args.d_run}nt skipping on alignments cluster-level"
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_clu_d, _sbc3 = compute_alignment_row_skip(
            cluster_align_universe, args.c_threshold, args.d_run, args.trim_start, pool, chunksize_clu)
        stats["avg_consensuses_per_cluster"] = (
            total_segments_clu_d / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = None
        rows.append({"strategy": label, "universe": "alignment_cluster_d", **stats})
        multi_consensus_report["cluster_level_d"] = multi_info

        label = (f"Strategy D DYNAMIC (<= {args.dynamic_cutoff} seqs: {args.dynamic_small}, "
                 f"else {args.dynamic_large}) + {args.d_run}nt skipping on alignments cluster-level")
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_clu_dyn_d, _sbc4 = compute_alignment_row_dynamic_skip(
            cluster_align_universe, args.d_run, args.trim_start,
            args.dynamic_small, args.dynamic_large, args.dynamic_cutoff, pool, chunksize_clu)
        stats["avg_consensuses_per_cluster"] = (
            total_segments_clu_dyn_d / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = None
        rows.append({"strategy": label, "universe": "alignment_cluster_d", **stats})
        multi_consensus_report["cluster_level_dynamic_d"] = multi_info

        # ── Strategy E -- fixed 100% threshold (same as C/D), but run_len
        # (skip tolerance) itself resolved PER UNIT: n_members <= e_cutoff
        # -> no skipping at all (run_len=1, Strategy C-equivalent), else
        # skipping allowed (run_len=d_run, Strategy D-equivalent).
        # Cluster-level only (matches 2.6.3_build_consensus.py's own
        # scope) -- +1 row, not +3 like C/D got.
        label = (f"Strategy E {int(args.c_threshold*100)} DYNAMIC nt-skipping "
                 f"(<= {args.e_cutoff} seqs: no skipping, else +{args.d_run}nt skipping) on alignments cluster-level")
        print(f"  Running: {label} ...", flush=True)
        stats, multi_info, total_segments_clu_e, _sbc5 = compute_alignment_row_dynamic_nt_skip(
            cluster_align_universe, args.c_threshold, args.trim_start, args.e_cutoff, args.d_run,
            pool, chunksize_clu)
        stats["avg_consensuses_per_cluster"] = (
            total_segments_clu_e / total_clusters_in if total_clusters_in else None
        )
        stats["avg_consensuses_per_duplicate_group"] = None
        rows.append({"strategy": label, "universe": "alignment_cluster_e", **stats})
        multi_consensus_report["cluster_level_e"] = multi_info

    # ── build final table rows ──────────────────────────────────────────────
    def pct(numerator, denominator):
        if not denominator:
            return NA
        return round(numerator / denominator * 100, 2)

    def num_or_na(value):
        return NA if value is None else round(value, 2)

    # Every row is cluster-level EXCEPT the two group-level Strategy C/D rows.
    GROUP_LEVEL_UNIVERSES = ("alignment_group", "alignment_group_d")

    table_rows = []
    for r in rows:
        is_cluster_level = r["universe"] not in GROUP_LEVEL_UNIVERSES
        dup_groups_found = r["n_units_found"]
        table_rows.append({
            "strategy": r["strategy"],
            "total_mges_in": total_mges_in,
            "n_mges_consensus_found": r["n_mges_found"],
            "n_mges_maxed_out": r["n_mges_maxed"],
            "total_clusters_in": total_clusters_in,
            "n_clusters_all_maxed_out": r["n_clusters_all_maxed"],
            "n_clusters_at_least_one_consensus": r["n_clusters_any_found"],
            "total_duplicate_groups_in": NA if is_cluster_level else total_dup_groups_in,
            "n_duplicate_groups_consensus_found": NA if is_cluster_level else dup_groups_found,
            "n_duplicate_groups_maxed_out": NA if is_cluster_level else r["n_units_maxed"],
            "pct_mges_consensus_found": pct(r["n_mges_found"], total_mges_in),
            "pct_clusters_at_least_one_consensus": pct(r["n_clusters_any_found"], total_clusters_in),
            "pct_duplicate_groups_consensus_found": (
                NA if is_cluster_level else pct(dup_groups_found, total_dup_groups_in)
            ),
            "avg_consensuses_per_cluster": num_or_na(r.get("avg_consensuses_per_cluster")),
            "avg_consensuses_per_duplicate_group": num_or_na(r.get("avg_consensuses_per_duplicate_group")),
            "consensus_length_nt": num_or_na(r.get("avg_consensus_length")),
            "pct_length_under_4": num_or_na(r.get("pct_length_under_4")),
            "pct_length_4_to_20": num_or_na(r.get("pct_length_4_to_20")),
            "pct_length_over_20": num_or_na(r.get("pct_length_over_20")),
        })

    columns = [
        "strategy", "total_mges_in", "n_mges_consensus_found", "n_mges_maxed_out",
        "total_clusters_in", "n_clusters_all_maxed_out", "n_clusters_at_least_one_consensus",
        "total_duplicate_groups_in", "n_duplicate_groups_consensus_found", "n_duplicate_groups_maxed_out",
        "pct_mges_consensus_found", "pct_clusters_at_least_one_consensus", "pct_duplicate_groups_consensus_found",
        "avg_consensuses_per_cluster", "avg_consensuses_per_duplicate_group",
        "consensus_length_nt", "pct_length_under_4", "pct_length_4_to_20", "pct_length_over_20",
    ]

    import csv
    with open(table_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(table_rows)
    print(f"\nSummary table -> {table_out}", flush=True)

    # ── text-rendered table + multi-consensus report, in the log ───────────
    col_widths = {c: max(len(c), max(len(str(r[c])) for r in table_rows)) for c in columns}

    def fmt_row(r) -> str:
        return "  ".join(str(r[c]).ljust(col_widths[c]) for c in columns)

    log_lines = [
        "=" * 100,
        "  Consensus Strategy Comparison",
        f"  Generated       : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"  rec_types       : {'ALL' if rec_types_filter is None else sorted(rec_types_filter)}",
        f"  stacking_dir    : {stacking_dir} (universe totals only)",
        f"  alignments_dir  : {alignments_dir}",
        f"  trim_start      : {args.trim_start}",
        "=" * 100,
        "",
        fmt_row({c: c for c in columns}),
        "-" * sum(col_widths[c] + 2 for c in columns),
    ]
    for r in table_rows:
        log_lines.append(fmt_row(r))

    log_lines += ["", "=" * 100,
                  "  MULTI-CONSENSUS REPORT (Strategy C/D/E only -- A/B never produce more than one)",
                  f"  Examples below are a REPRESENTATIVE SAMPLE (up to {args.sample_size}, "
                  f"random seed {args.sample_seed}) -- not just the first ones directory iteration happens to find.",
                  "=" * 100]

    rng = random.Random(args.sample_seed)
    report_levels = [
        ("group_level", "duplicate group", "Strategy C, group-level"),
        ("cluster_level", "cluster", "Strategy C, cluster-level (fixed)"),
        ("cluster_level_dynamic", "cluster (dynamic threshold)", "Strategy C, cluster-level (dynamic)"),
        ("group_level_d", "duplicate group", "Strategy D, group-level"),
        ("cluster_level_d", "cluster", "Strategy D, cluster-level (fixed)"),
        ("cluster_level_dynamic_d", "cluster (dynamic threshold)", "Strategy D, cluster-level (dynamic)"),
        ("cluster_level_e", "cluster (dynamic nt-skipping)", "Strategy E, cluster-level"),
    ]
    for level_key, key_label, section_title in report_levels:
        info = multi_consensus_report[level_key]
        all_keys = info["all_keys"]
        if len(all_keys) > args.sample_size:
            examples = rng.sample(all_keys, args.sample_size)
        else:
            examples = list(all_keys)

        log_lines += [
            "",
            f"-- {section_title} --",
            f"  {key_label.capitalize()}s with more than 1 consensus segment on left and/or right: {info['n_units']:,}",
        ]
        if examples:
            log_lines.append(f"  Representative sample ({len(examples)} of {info['n_units']:,}):")
            for key in examples:
                is_group_level = level_key.startswith("group_level")
                if is_group_level:
                    rec_type, cluster_num, dup_group = key
                    example_name = (f"{rec_type}_{cluster_num}_{dup_group}_{{left,right}}.aln.fasta  "
                                     f"(under duplicate_group_level/{rec_type}/cluster_{cluster_num}/)")
                else:
                    rec_type, cluster_num = key
                    example_name = (f"{rec_type}_cluster_{cluster_num}_{{left,right}}.aln.fasta  "
                                     f"(under cluster_level/{rec_type}/)")
                log_lines.append(f"    {example_name}")
        else:
            log_lines.append("  (none)")

    log_lines += ["", "=" * 100, f"Finished : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", "=" * 100]

    with open(log_out, "w") as fh:
        fh.write("\n".join(log_lines) + "\n")
    print(f"Log -> {log_out}", flush=True)


if __name__ == "__main__":
    main()