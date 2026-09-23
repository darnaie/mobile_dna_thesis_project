#!/usr/bin/env python3
"""
2.6.3_build_consensus.py  —  Section 2.6.3
==============================================
input : BASE/consensus/cluster_level/<rec_type>/<rec_type>_cluster_<N>_<left|right>.raw.fasta
        and/or .aln.fasta (Section 2.6.1) -- CLUSTER LEVEL ONLY, deliberately
        (no duplicate-group-level consensus building here).
output: BASE/files/cluster_consensus_results.tsv (one row per
        rec_type/cluster/side actually processed -- the structured result
        Section 2.6.4 reads to append onto the per-MGE table) and an
        EXTENSIVE per-cluster log (not a strategy-comparison table like
        consensus_building_evaluation.py's -- this is the real,
        once-per-pipeline-run consensus build, not a method comparison).

Adapted from consensus_building_evaluation.py's Strategy A/B (stacking)
and Strategy C/D (alignment, multi-segment) machinery -- imported
directly from that script rather than reimplemented, so this always uses
exactly the same, already-tested consensus logic.

MGE-OVERLAP EXCLUSION (NEW, see chat -- brings this in line with
2.7.1_annotate_ml_labels.py, which already had it): 2.5.1_extract_flanks.py
deliberately extracts each flank so that its LAST mge_overlap nt (LEFT) or
FIRST mge_overlap nt (RIGHT) sit INSIDE the MGE itself -- that's the right
thing for 2.5.1's own purpose (anchoring the window on a shared reference
point across independent insertions), and 2.6.1's alignment step is
SUPPOSED to use that overlap as an anchor too (it still does -- nothing
about 2.6.1 or the alignment files themselves changed). But CONSENSUS
BUILDING is a different question: a run of "perfect identity" inside the
mge_overlap zone doesn't tell you anything about a boundary-marking motif
-- every member's overlap nt come from the same MGE family's own body, so
of course they can look conserved there, independent of whether a real
TIR-like signal exists at the true edge. Previously this script scanned
the mge_overlap zone right along with everything else, meaning a
"consensus" could be reporting boundary-irrelevant, MGE-body sequence.
Fixed via strip_mge_overlap() (moved into consensus_building_evaluation.py
so this script and 2.7.1 share one implementation, not two copies that
could drift apart): every sequence -- both the *.raw.fasta stack and the
*.aln.fasta alignment -- has its mge_overlap zone removed before ANY
consensus-building function ever sees it. The alignment step itself
(2.6.1) is untouched; only what THIS script does with its output changed.
New config key: mge_overlap (default 30, matching 2.5.1's own default --
keep these two in sync).

SEGMENT SELECTION FOR consensus_source="alignment" -- CHANGED (see chat):
this used to keep the LONGEST closed segment as "the" consensus. Fixed to
match 2.7.1's own logic exactly: among segments (and the trailing,
unclosed run, checked first -- see scan_alignment_multi_segment_with_
skip_full_side_aware()'s own docstring for why the trailing run is
usually the single closest-to-boundary candidate available) with length
in [4, 20] nt inclusive, pick whichever ONE sits CLOSEST to the actual
MGE boundary -- not the longest one found anywhere in the window. A long
run of coincidental identity in the middle of the flank is not what
"the consensus" is supposed to mean here; the boundary-proximal signal
is. This also switches the underlying scan from the side-BLIND
scan_alignment_multi_segment_with_skip() to the side-aware, trailing-run-
aware scan_alignment_multi_segment_with_skip_full_side_aware() (the same
function 2.7.1 already uses) -- completing the "side-aware" fix that
2.7.1's own docstring had explicitly deferred doing here, back when this
was a separate, smaller-scope change. This DOES change
cluster_consensus_results.tsv's real output versus previous runs -- by
design, since the previous RIGHT-side / longest-segment behaviour was
the bug being fixed, not a preserved convention.

consensus_source (config, default "stacking"): "stacking" builds via
Strategy A ("run-length tolerant cutoff", default, run length
a_run_length, or dynamically 1 vs a_run_length if use_dynamic_nt_skipping)
or Strategy B ("single-column cutoff", strategy: "B", no tolerance
concept at all -- use_dynamic_nt_skipping has no effect on it) on the
mge_overlap-stripped *.raw.fasta stack; "alignment" builds via the
boundary-closest-segment selection described above on the mge_overlap-
stripped *.aln.fasta alignment; "both" runs both and records both.

DYNAMIC THRESHOLDING (config: use_dynamic_threshold, default true): before
extracting a cluster's consensus, count how many sequences are actually
in play (n_sequences = max(len(left_seqs), len(right_seqs))). If
n_sequences <= dynamic_threshold_and_nt_skipping_cutoff (default 3), use
dynamic_threshold_small (default 0.9); otherwise use
dynamic_threshold_large (default 0.7). The reasoning: a small group's
per-column identity is far noisier (one disagreeing sequence swings the
fraction a lot), so a stricter cutoff would end the consensus almost
immediately -- a looser one still means "most of very few sequences
agree". Set use_dynamic_threshold: false in config.yaml to use a single
fixed_threshold for every cluster instead.

DYNAMIC NT-SKIPPING (config: use_dynamic_nt_skipping, default false):
only meaningful when strategy: "A" (stacking) or for the run-length
tolerance used by the alignment scan. Rather than a fixed run-length
tolerance (a_run_length nt of below-threshold columns tolerated before
closing a consensus/segment), the tolerance ITSELF is resolved PER
CLUSTER from n_sequences, using the SAME cutoff as dynamic thresholding:
n_sequences <= cutoff -> NO tolerance at all (strict cutoff at the first
below-threshold column/run); n_sequences > cutoff -> a_run_length nt of
tolerance allowed. use_dynamic_threshold and use_dynamic_nt_skipping are
independent switches and combine freely.

CONSENSUS_START (in cluster_consensus_results.tsv, between n_sequences
and consensus_sequence): where the reported consensus starts, in that
side's OWN mge_overlap-stripped sequence's coordinate numbering (0-based)
-- NOT the original, overlap-inclusive flank's numbering, and NOT the
alignment file's own column numbering either, now that stripping happens
before scanning. For stacking sources this is always trim_start itself
(or "NA" if no consensus was produced). For alignment sources it's the
actual starting position of whichever segment was picked -- "NA" if none
qualified.

rec_types (config, default "all"): restrict to specific rec_type folder
names, same convention as consensus_building_evaluation.py.

Run standalone:
    python 2.6.3_build_consensus.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Optional

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

# Reuse the exact same, already-tested consensus logic -- not reimplemented.
import importlib.util as _ilu

_eval_spec = _ilu.spec_from_file_location(
    "consensus_building_evaluation", Path(__file__).resolve().parent / "consensus_building_evaluation.py"
)
_eval_mod = _ilu.module_from_spec(_eval_spec)
_eval_spec.loader.exec_module(_eval_mod)
read_fasta = _eval_mod.read_fasta
consensus_strategy_a = _eval_mod.consensus_strategy_a
consensus_strategy_b = _eval_mod.consensus_strategy_b
strip_mge_overlap = _eval_mod.strip_mge_overlap
scan_alignment_multi_segment_with_skip_full_side_aware = _eval_mod.scan_alignment_multi_segment_with_skip_full_side_aware

SECTION_KEY = "2.6.3"
MIN_LEN, MAX_LEN = 4, 20  # inclusive both ends -- matches 2.7.1's own bucket exactly


def parse_rec_types(arg: Optional[str]) -> Optional[set]:
    if arg is None or str(arg).strip().lower() == "all":
        return None
    return {t.strip() for t in str(arg).split(",") if t.strip()}


def dynamic_or_fixed_threshold(n_sequences: int, section_cfg: dict) -> float:
    if not section_cfg.get("use_dynamic_threshold", True):
        return float(section_cfg.get("fixed_threshold", 0.9))
    cutoff = int(section_cfg.get("dynamic_threshold_and_nt_skipping_cutoff", 3))
    small = float(section_cfg.get("dynamic_threshold_small", 0.9))
    large = float(section_cfg.get("dynamic_threshold_large", 0.7))
    return small if n_sequences <= cutoff else large


def resolve_run_len(n_sequences: int, section_cfg: dict, a_run_length: int) -> int:
    """Only meaningful when strategy == 'A' (see caller). If
    use_dynamic_nt_skipping is off, the run length is just a_run_length
    for every unit, unconditionally (the pre-existing behaviour). If on,
    it's resolved PER UNIT from n_sequences: <= cutoff -> 1 (no
    tolerance at all); otherwise a_run_length (tolerance allowed). Uses
    the SAME cutoff as dynamic_or_fixed_threshold, so setting both
    use_dynamic_threshold and use_dynamic_nt_skipping true applies both
    dynamic behaviours off of one shared n_sequences check, per-unit."""
    if not section_cfg.get("use_dynamic_nt_skipping", False):
        return a_run_length
    cutoff = int(section_cfg.get("dynamic_threshold_and_nt_skipping_cutoff", 3))
    return 1 if n_sequences <= cutoff else a_run_length


def find_clusters(cluster_dir: Path, rec_types_filter: Optional[set]) -> list[tuple[str, str]]:
    """Returns [(rec_type, cluster_num), ...] discovered from *_left/right.raw.fasta
    OR *.aln.fasta filenames under cluster_dir (whichever exist)."""
    found: set[tuple[str, str]] = set()
    for rt_dir in sorted(cluster_dir.iterdir()):
        if not rt_dir.is_dir():
            continue
        rec_type = rt_dir.name
        if rec_types_filter is not None and rec_type not in rec_types_filter:
            continue
        prefix = f"{rec_type}_cluster_"
        for fp in sorted(rt_dir.iterdir()):
            if not fp.is_file() or not fp.name.startswith(prefix):
                continue
            for suffix in (".raw.fasta", ".aln.fasta"):
                if fp.name.endswith(suffix):
                    remainder = fp.name[len(prefix):-len(suffix)]  # "<N>_left" or "<N>_right"
                    parts = remainder.rsplit("_", 1)
                    if len(parts) == 2 and parts[1] in ("left", "right"):
                        found.add((rec_type, parts[0]))
    return sorted(found)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--consensus-dir", default=None, help="Override BASE/consensus")
    ap.add_argument("--out", default=None, help="Override BASE/files/cluster_consensus_results.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_6"][SECTION_KEY]

    log = get_logger("2.6.3_build_consensus", paths["logs_dir"])

    consensus_dir = Path(args.consensus_dir) if args.consensus_dir else paths["base_dir"] / "consensus"
    cluster_dir = consensus_dir / "cluster_level"
    out_path = Path(args.out) if args.out else files_dir / "cluster_consensus_results.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Build cluster-level consensus")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_6_3_build_consensus"):
        return

    if not cluster_dir.exists():
        log.error(f"{cluster_dir} not found -- run 2.6.1_align_flanks.sh first.")
        raise SystemExit(1)

    consensus_source = section_cfg.get("consensus_source", "stacking")
    if consensus_source not in ("stacking", "alignment", "both"):
        log.error(f"consensus_source must be 'stacking', 'alignment', or 'both' -- got {consensus_source!r}")
        raise SystemExit(1)
    strategy = section_cfg.get("strategy", "A")
    a_run_length = int(section_cfg.get("a_run_length", 2))
    trim_start = int(section_cfg.get("trim_start", 10))
    mge_overlap = int(section_cfg.get("mge_overlap", 30))
    use_dynamic_nt_skipping = bool(section_cfg.get("use_dynamic_nt_skipping", False))
    rec_types_filter = parse_rec_types(section_cfg.get("rec_types", "all"))

    log.info(f"consensus_source          : {consensus_source}")
    log.info(f"strategy                  : {strategy}" + (f" (base run length {a_run_length}nt)" if strategy == "A" else ""))
    log.info(f"trim_start                : {trim_start}")
    log.info(f"mge_overlap (excluded from consensus, kept for alignment): {mge_overlap}")
    log.info(f"use_dynamic_threshold     : {section_cfg.get('use_dynamic_threshold', True)}")
    log.info(f"use_dynamic_nt_skipping   : {use_dynamic_nt_skipping}")
    if section_cfg.get("use_dynamic_threshold", True) or use_dynamic_nt_skipping:
        log.info(f"  cutoff (n_sequences <=) : {section_cfg.get('dynamic_threshold_and_nt_skipping_cutoff', 3)}")
    if section_cfg.get("use_dynamic_threshold", True):
        log.info(f"    -> threshold {section_cfg.get('dynamic_threshold_small', 0.9)}, else "
                 f"{section_cfg.get('dynamic_threshold_large', 0.7)}")
    else:
        log.info(f"  fixed_threshold          : {section_cfg.get('fixed_threshold', 0.9)}")
    if use_dynamic_nt_skipping:
        log.info(f"    -> no skipping (run_len=1), else +{a_run_length}nt skipping")
    log.info(f"alignment-source segment length bucket : [{MIN_LEN}, {MAX_LEN}] inclusive, boundary-closest wins (not longest)")
    log.info(f"rec_types filter          : {'ALL' if rec_types_filter is None else sorted(rec_types_filter)}")
    log.info("-" * 70)

    t0 = time.time()
    clusters = find_clusters(cluster_dir, rec_types_filter)
    log.info(f"Found {len(clusters):,} rec_type/cluster combinations to process")
    log.info("-" * 70)

    out_rows = []
    n_processed = n_skipped_no_files = 0

    for rec_type, cluster_num in clusters:
        rt_dir = cluster_dir / rec_type
        prefix = f"{rec_type}_cluster_{cluster_num}"

        log.info(f"\n=== {rec_type} / cluster_{cluster_num} ===")

        sources_to_run = []
        if consensus_source in ("stacking", "both"):
            sources_to_run.append("stacking")
        if consensus_source in ("alignment", "both"):
            sources_to_run.append("alignment")

        any_written = False
        for source in sources_to_run:
            suffix = ".raw.fasta" if source == "stacking" else ".aln.fasta"
            left_path = rt_dir / f"{prefix}_left{suffix}"
            right_path = rt_dir / f"{prefix}_right{suffix}"

            left_seqs_raw = read_fasta(left_path) if left_path.exists() else []
            right_seqs_raw = read_fasta(right_path) if right_path.exists() else []

            if not left_seqs_raw and not right_seqs_raw:
                log.info(f"  [{source}] no {suffix} files found -- skipping")
                continue

            # MGE-overlap exclusion (see module docstring): applied here,
            # AFTER reading the (overlap-inclusive) alignment/stack, so
            # 2.6.1's own output is completely untouched -- only what this
            # script does with it changes. Side-aware (LEFT strips its
            # last mge_overlap chars, RIGHT its first) -- same convention
            # 2.7.1 already uses.
            left_seqs = [s for s in (strip_mge_overlap(s, "left", mge_overlap) for s in left_seqs_raw) if s]
            right_seqs = [s for s in (strip_mge_overlap(s, "right", mge_overlap) for s in right_seqs_raw) if s]

            n_sequences = max(len(left_seqs), len(right_seqs))
            if n_sequences == 0:
                log.info(f"  [{source}] every sequence was entirely inside the mge_overlap zone -- skipping")
                continue

            threshold = dynamic_or_fixed_threshold(n_sequences, section_cfg)
            run_len = resolve_run_len(n_sequences, section_cfg, a_run_length) if strategy == "A" else 1
            log.info(f"  [{source}] n_sequences={n_sequences}  threshold={threshold}"
                     + (f"  run_len={run_len}" if strategy == "A" else "")
                     + f"  (left={len(left_seqs)} seqs, right={len(right_seqs)} seqs, mge_overlap={mge_overlap} already stripped)")

            for side, seqs in (("left", left_seqs), ("right", right_seqs)):
                if not seqs:
                    log.info(f"    {side}: no sequences -- skipping")
                    continue

                if source == "stacking":
                    if strategy == "A":
                        cons, reason, length = consensus_strategy_a(seqs, threshold, run_len, trim_start)
                        method_label = f"stacking_A_run{run_len}" + ("_dyn" if use_dynamic_nt_skipping else "")
                    else:
                        cons, reason, length = consensus_strategy_b(seqs, threshold, trim_start)
                        method_label = "stacking_B"
                    maxed = (reason == "maxed_out")
                    cons_start = trim_start if length > 0 else "NA"
                    log.info(f"    {side}: reason={reason}  length={length}  "
                             f"consensus={cons[:30]}{'...' if len(cons) > 30 else ''}")
                else:
                    # Side-aware, trailing-run-aware scan (same function
                    # 2.7.1 already uses) -- see module docstring for why
                    # this replaced the side-blind scanner.
                    segments, _maxed, trailing = scan_alignment_multi_segment_with_skip_full_side_aware(
                        seqs, threshold, run_len, trim_start, side,
                    )
                    method_label = f"alignment_boundary_closest_run{run_len}" + ("_dyn" if use_dynamic_nt_skipping else "")

                    # The trailing (unclosed) run is, by construction,
                    # already the boundary-closest candidate available
                    # whenever it qualifies -- checked first, same as
                    # 2.7.1. Only falls back to the best CLOSED segment
                    # (boundary-closest among those, NOT longest) if the
                    # trailing run doesn't qualify.
                    chosen = None
                    chosen_kind = None
                    if trailing is not None and MIN_LEN <= len(trailing[1]) <= MAX_LEN:
                        chosen = trailing
                        chosen_kind = "trailing_at_boundary"
                    else:
                        candidates = [(s, t) for s, t in segments if MIN_LEN <= len(t) <= MAX_LEN]
                        if candidates:
                            chosen = max(candidates, key=lambda st: st[0]) if side == "left" else min(candidates, key=lambda st: st[0])
                            chosen_kind = "closed_segment"

                    maxed = trailing is not None and chosen_kind != "trailing_at_boundary"
                    if chosen is not None:
                        cons_start, cons = chosen
                        length = len(cons)
                    else:
                        cons_start, cons, length = "NA", "", 0

                    all_lengths = [len(t) for _s, t in segments] + ([len(trailing[1])] if trailing else [])
                    log.info(f"    {side}: {len(segments)} closed segment(s) + "
                             f"{'1 trailing run' if trailing else 'no trailing run'}, "
                             f"chosen={chosen_kind or 'none'}  length={length}  all_lengths={all_lengths}")

                out_rows.append({
                    "rec_type": rec_type, "cluster_number": cluster_num, "side": side,
                    "consensus_source": source, "method": method_label, "threshold_used": threshold,
                    "n_sequences": n_sequences, "consensus_start": cons_start,
                    "consensus_sequence": cons if cons else "NA",
                    "consensus_length": length, "maxed": "yes" if maxed else "no",
                })
                any_written = True

        if any_written:
            n_processed += 1
        else:
            n_skipped_no_files += 1

    log.info("\n" + "-" * 70)
    log.info(f"Writing results table: {out_path}")
    out_columns = [
        "rec_type", "cluster_number", "side", "consensus_source", "method",
        "threshold_used", "n_sequences", "consensus_start", "consensus_sequence", "consensus_length", "maxed",
    ]
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"rec_type/cluster combinations found : {len(clusters):,}")
    log.info(f"  processed (>=1 result row)        : {n_processed:,}")
    log.info(f"  skipped (no matching files)        : {n_skipped_no_files:,}")
    log.info(f"Result rows written                  : {len(out_rows):,}")
    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()