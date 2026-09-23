#!/usr/bin/env python3
"""
pick_representatives.py  —  Section 2.2.3
=============================================
input : BASE/files/skani_output/*.tsv (Section 2.2.2)
output: BASE/genome_seqs_repr/  (a directory of representative genome FASTAs)

Adapted from the original parallel_ani.sh's Python driver (Steps 1-3):

  Step 1 -- for each skani TSV, accept a genome pair only when
              ANI                  >= ani_threshold
              Align_fraction_ref   >= align_threshold
              Align_fraction_query >= align_threshold
            and union them via union-find. Also records each genome's
            FASTA path from the TSV's Ref_file/Query_file columns.
  Step 2 -- merge all per-TSV groups into one global union-find (kept for
            robustness/parity with the original script, even though in
            practice groups never span TSVs here, since Section 2.2.1
            already partitioned genomes by specI and skani only ever
            compares within one specI folder).
  Step 3 -- pick one representative per group (alphabetically first
            accession) and symlink its FASTA into BASE/genome_seqs_repr/.

Run standalone:
    python 2.2.3_pick_representatives.py --config path/to/config.yaml
    python 2.2.3_pick_representatives.py --input-dir ./skani_output --out-dir ./genome_seqs_repr --ani-threshold 99.8 --align-threshold 90
"""

from __future__ import annotations

import argparse
import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.2.3"


# ── Union-Find ──────────────────────────────────────────────────────────────

def uf_find(parent: dict, x: str) -> str:
    root = x
    while parent.get(root, root) != root:
        root = parent[root]
    while parent.get(x, x) != root:
        parent[x], x = root, parent.get(x, x)
    return root


def uf_union(parent: dict, a: str, b: str) -> None:
    ra, rb = uf_find(parent, a), uf_find(parent, b)
    if ra != rb:
        if ra < rb:
            parent[rb] = ra
        else:
            parent[ra] = rb


def uf_ensure(parent: dict, x: str) -> None:
    parent.setdefault(x, x)


# ── Step 1 ───────────────────────────────────────────────────────────────────

def process_tsv(tsv_path: Path, ani_threshold: float, align_threshold: float):
    parent: dict = {}
    path_map: dict = {}
    accepted = rejected_ani = rejected_align = 0

    with open(tsv_path, newline="", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            ref_name = (row.get("Ref_file") or "").strip()
            query_name = (row.get("Query_file") or "").strip()
            if not ref_name or not query_name:
                continue

            path_map[ref_name] = ref_name
            path_map[query_name] = query_name
            uf_ensure(parent, ref_name)
            uf_ensure(parent, query_name)

            try:
                ani = float(row.get("ANI", 0))
                af_r = float(row.get("Align_fraction_ref", 0))
                af_q = float(row.get("Align_fraction_query", 0))
            except (ValueError, TypeError):
                continue

            if ani < ani_threshold:
                rejected_ani += 1
                continue
            if af_r < align_threshold or af_q < align_threshold:
                rejected_align += 1
                continue

            accepted += 1
            uf_union(parent, ref_name, query_name)

    return parent, path_map, accepted, rejected_ani, rejected_align


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--input-dir", default=None, help="Override BASE/files/skani_output")
    ap.add_argument("--out-dir", default=None, help="Override BASE/genome_seqs_repr")
    ap.add_argument("--ani-threshold", type=float, default=None)
    ap.add_argument("--align-threshold", type=float, default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    input_dir = Path(args.input_dir) if args.input_dir else paths["files_dir"] / "skani_output"
    out_dir = Path(args.out_dir) if args.out_dir else paths["genome_seqs_repr"]
    section_cfg = cfg["section_2_2"][SECTION_KEY]

    log = get_logger("2.2.3_pick_representatives", paths["logs_dir"])

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Pick representative genomes")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_dir, log):
        n_fasta = len(list(out_dir.glob("*.fasta"))) if out_dir.exists() else 0
        log.info(f"{out_dir} now has {n_fasta:,} representative genomes (via placeholder)")
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    extra_sbatch = [f"#SBATCH --cpus-per-task={section_cfg.get('slurm_cpus_per_task', 16)}"]
    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log,
                              job_name="2_2_3_pick_representatives", extra_sbatch_lines=extra_sbatch):
        return

    ani_threshold = args.ani_threshold if args.ani_threshold is not None else section_cfg["ani_threshold"]
    align_threshold = (
        args.align_threshold if args.align_threshold is not None else section_cfg["align_threshold"]
    )

    if not input_dir.exists():
        log.error(f"{input_dir} not found -- run skani_dedup.py (Section 2.2.2) first.")
        raise SystemExit(1)

    tsv_files = sorted(input_dir.glob("*.tsv"))
    log.info(f"ANI threshold      : {ani_threshold}")
    log.info(f"Align threshold    : {align_threshold}")
    log.info(f"Found {len(tsv_files):,} TSV(s) in {input_dir}")
    if not tsv_files:
        log.error("No TSV files found -- nothing to do.")
        raise SystemExit(1)

    t0 = time.time()

    # ── Step 1 ───────────────────────────────────────────────────────────
    parent_global: dict = {}
    path_map_global: dict = {}
    total_accepted = total_rej_ani = total_rej_align = 0

    for tsv_path in tsv_files:
        parent, path_map, accepted, rej_ani, rej_align = process_tsv(
            tsv_path, ani_threshold, align_threshold
        )
        total_accepted += accepted
        total_rej_ani += rej_ani
        total_rej_align += rej_align
        path_map_global.update(path_map)
        for node in parent:
            uf_ensure(parent_global, node)
        for node in parent:
            uf_union(parent_global, node, uf_find(parent, node))
        log.info(f"  {tsv_path.name}: {len(parent)} genomes | accepted={accepted} "
                 f"rejected_ani={rej_ani} rejected_align={rej_align}")

    log.info(f"[Step 1] Done -- {len(parent_global):,} unique genomes seen, "
             f"{total_accepted:,} accepted pairs (rejected_ani={total_rej_ani:,}, "
             f"rejected_align={total_rej_align:,})")

    # ── Step 2 ───────────────────────────────────────────────────────────
    global_groups: dict = defaultdict(list)
    for node in parent_global:
        global_groups[uf_find(parent_global, node)].append(node)

    representatives = [sorted(members)[0] for _root, members in sorted(global_groups.items())]
    log.info(f"[Step 2] Global groups: {len(global_groups):,}  "
             f"Representatives: {len(representatives):,}")

    # ── Step 3 ───────────────────────────────────────────────────────────
    out_dir.mkdir(parents=True, exist_ok=True)
    n_copied = n_missing = 0
    for rep in representatives:
        fasta_path = Path(path_map_global.get(rep, rep))
        if not fasta_path.exists():
            log.warning(f"  WARN: representative FASTA not found on disk: {fasta_path}")
            n_missing += 1
            continue
        dest = out_dir / fasta_path.name
        if not dest.exists() and not dest.is_symlink():
            dest.symlink_to(fasta_path.resolve())
        n_copied += 1

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY")
    log.info("-" * 70)
    log.info(f"Unique genomes seen        : {len(parent_global):,}")
    log.info(f"Identity groups            : {len(global_groups):,}")
    log.info(f"Representatives selected   : {len(representatives):,}")
    log.info(f"  linked successfully      : {n_copied:,}")
    log.info(f"  missing on disk          : {n_missing:,}")
    log.info(f"-> {out_dir}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()