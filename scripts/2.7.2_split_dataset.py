#!/usr/bin/env python3
"""
2.7.2_split_dataset.py  —  Section 2.7.2
============================================
input : BASE/files/ml_labels.tsv (2.7.1)
        BASE/files/f13_dataset.txt (2.1.3 -- taxonomy table, GCA -> specI;
        read directly by THIS script, not propagated through 2.7.1, so
        this whole fix is contained to this one file -- see SPECIES-LEVEL
        LEAKAGE PREVENTION below)
output: BASE/files/ml_labels_split.tsv (same rows, + two new columns:
        "split" -- train/val/test, for single-split training (2.7.3)
        "fold"  -- 0..n_folds-1, for k-fold cross-validation (2.7.4)
        Both are computed independently and both are always written --
        2.7.3 only ever reads "split"; 2.7.4 only ever reads "fold". You
        don't need to choose one mode over the other at split time.)

SPECIES-LEVEL LEAKAGE PREVENTION (NEW, see chat)
------------------------------------------------------
The previous version of this script split at the duplicate-group level
only: every MGE in one duplicate group always landed in the same split,
but two DIFFERENT duplicate groups from the SAME species could land in
different splits. That's a real leakage risk this fix closes: two
distinct duplicate groups from the same species can still share
species-level, non-transposon-specific sequence characteristics (GC
content, codon usage, other host-genome biases) that have nothing to do
with the boundary-motif signal the model is supposed to be learning --
if the test set contains a species the model already saw (via a
different duplicate group) during training, a good test score could
partly reflect having learned "what this species' DNA generally looks
like" rather than "what a boundary motif looks like". This version
additionally guarantees NO SPECIES (specI) appears in more than one of
{train, val, test} (and, separately, in more than one fold) at all.

Species are read directly here (Path to f13_dataset.txt, the same
taxonomy table used elsewhere in this pipeline -- see 2.2.5's own
usage), by deriving each row's GCA accession from its own mge_id (via
utils.gca_from_mge_id(), the same helper the rest of the pipeline
already uses for this) and looking it up. This keeps the whole fix
contained to this one script -- ml_labels.tsv itself is read, never
written to, and 2.7.1 doesn't need to change or re-run.

WHY THIS ISN'T AS SIMPLE AS "GROUP BY SPECIES INSTEAD OF BY GROUP"
--------------------------------------------------------------------------
A duplicate group's members can, in principle, come from MORE THAN ONE
species: duplicate-group membership (Section 2.4) is based on the
transposon sequence being ~100% identical, which horizontal transfer
can produce across species boundaries -- two different host species
carrying an identical copy of the same mobile element is exactly the
kind of event this pipeline is built to detect duplicates of. This means
species and duplicate groups don't nest cleanly one inside the other:
naively assigning "each species to one split" could still leave a
single duplicate group with some members' species assigned to train and
others assigned to test.

Resolved with a union-find over species: for every duplicate group,
every species it touches is unioned into one connected component (any
two species that ever co-occur in the same duplicate group -- in ANY
rec_type, not just the current one, since the leakage concern is
species-level, not family-level -- must end up in the same split).
Splitting then operates on these SPECIES COMPONENTS, not on individual
species or individual duplicate groups directly: a component is the
smallest unit that can be assigned to a split without risk of splitting
a duplicate group's own members across two different splits. In
practice the large majority of components are expected to be
single-species (most duplicate groups are single-species), but the
algorithm is correct regardless of how much cross-species duplication
the real data actually has.

An MGE whose GCA can't be resolved from its own ID, or whose GCA isn't
in the taxonomy table, is given a synthetic, uniquely-its-own "species"
(so it can never be silently unioned with a real species it doesn't
actually share) -- counted and reported, not silently dropped.

WHY THE RESULT IS NOT EXACTLY 80:10:10, AND WHY THAT'S EXPECTED
--------------------------------------------------------------------------
Once species components (rather than individual duplicate groups) are
the unit being placed into train/val/test, hitting the target ratios
exactly -- overall AND per rec_type -- is no longer always possible: a
single large, highly-connected species component might contribute a
big, indivisible chunk of one rec_type's groups, forcing that rec_type's
achieved split ratio away from 80:10:10 regardless of which split it's
assigned to. Assignment here is GREEDY: species components are placed
one at a time (largest first), each into whichever of {train, val,
test} currently leaves the smallest total weighted deviation from
target proportions, considering both the overall ratio and every
rec_type the component touches. This gets close to the target ratios in
the typical case (most components are small and single-species) without
ever violating the hard species-exclusivity constraint. The ACTUAL
achieved ratios -- not just the requested ones -- are always logged,
per rec_type, so a real deviation is visible rather than assumed away
(see the summary this script prints).

Fold assignment (assign_folds() below) uses the exact same species-
component machinery, generalised from a 3-way (train/val/test) to an
n_folds-way partition -- for the same reason: a fold-based cross-
validation test set has exactly the same species-leakage risk a
single held-out test set does, once you notice it.

STRATIFIED BY rec_type
--------------------------
As before, the target proportions themselves are still meant to give
proportional representation to all 7 rec_types, including the negative
control (DEDD_Tnp_IS110) -- see the greedy assignment above for how this
interacts with the new species constraint.

Reproducible: shuffling/ordering is seeded (--seed, default 42 for the
train/val/test split; --fold-seed, default 42 as well, for fold
assignment -- separate seeds so tuning one doesn't silently perturb the
other) and done on the sorted list of species-component keys, so
re-running with the same seed(s) on the same ml_labels.tsv and the same
taxonomy table always produces the same split and the same fold
assignment.

Run standalone:
    python 2.7.2_split_dataset.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit, gca_from_mge_id, load_taxonomy_map

SECTION_KEY = "2.7.2"

GroupKey = Tuple[str, str, str]  # (rec_type, cluster_number, duplicate_group_number)


# ══════════════════════════════════════════════════════════════════════════
#  SPECIES RESOLUTION
# ══════════════════════════════════════════════════════════════════════════

def resolve_species(mge_id: str, taxonomy_map: Dict[str, tuple]) -> str:
    """mge_id -> specI, via its own embedded GCA accession. Falls back to
    a synthetic, uniquely-its-own species key (never collides with a real
    specI, and never collides with another unresolved mge_id either) when
    the GCA can't be parsed or isn't in the taxonomy table -- see module
    docstring for why this is the safe default rather than dropping the
    row or lumping all unresolved rows together."""
    gca = gca_from_mge_id(mge_id)
    if gca is not None:
        info = taxonomy_map.get(gca)
        if info is not None:
            return info[0]  # specI
        return f"__unresolved_gca__:{gca}"
    return f"__unresolved_id__:{mge_id}"


# ══════════════════════════════════════════════════════════════════════════
#  UNION-FIND over species
# ══════════════════════════════════════════════════════════════════════════

class UnionFind:
    def __init__(self):
        self.parent: Dict[str, str] = {}

    def find(self, x: str) -> str:
        self.parent.setdefault(x, x)
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]  # path halving
            x = self.parent[x]
        return x

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[ra] = rb


def build_species_components(
    group_species: Dict[GroupKey, Set[str]],
) -> Dict[GroupKey, str]:
    """Every duplicate group's own set of member species is unioned
    together (they must end up in the same split -- see module
    docstring). Returns {group_key: component_id} -- component_id is
    just one representative species string per component, stable given
    the same input (union-find's own path-halving is deterministic for
    a fixed processing order, and callers always iterate group_species
    in sorted key order)."""
    uf = UnionFind()
    for key in sorted(group_species):
        species_set = group_species[key]
        species_list = sorted(species_set)
        for s in species_list:
            uf.find(s)  # ensure registered even if this group is single-species
        for a, b in zip(species_list, species_list[1:]):
            uf.union(a, b)

    group_to_component: Dict[GroupKey, str] = {}
    for key in sorted(group_species):
        # component id = the (deterministic) union-find root of any one of
        # this group's own species -- all its species share the same root
        # by construction, so any one of them identifies the component.
        any_species = sorted(group_species[key])[0]
        group_to_component[key] = uf.find(any_species)
    return group_to_component


# ══════════════════════════════════════════════════════════════════════════
#  GREEDY, SPECIES-COMPONENT-AWARE ASSIGNMENT
# ══════════════════════════════════════════════════════════════════════════

def _component_footprints(
    groups_by_component: Dict[str, List[GroupKey]], rec_type_of: Dict[GroupKey, str],
) -> Dict[str, Dict[str, int]]:
    """component_id -> {"__total__": n, rec_type_a: n_a, rec_type_b: n_b, ...}
    -- how many duplicate groups this component contributes overall and
    per rec_type. Used by the greedy assignment to score candidate
    placements against both the overall and the per-rec_type target
    proportions at once."""
    footprints: Dict[str, Dict[str, int]] = {}
    for comp_id, keys in groups_by_component.items():
        counts: Dict[str, int] = defaultdict(int)
        for k in keys:
            counts[rec_type_of[k]] += 1
            counts["__total__"] += 1
        footprints[comp_id] = dict(counts)
    return footprints


def greedy_assign_components(
    groups_by_component: Dict[str, List[GroupKey]],
    rec_type_of: Dict[GroupKey, str],
    bucket_names: List[str],
    bucket_target_fracs: List[float],
    seed: int,
) -> Dict[str, str]:
    """
    Assigns every species component to exactly one of bucket_names (e.g.
    ["train","val","test"], or ["fold0",...,"foldN-1"]), respecting the
    hard species constraint (a component is an atomic unit -- see module
    docstring) while trying to approximate bucket_target_fracs, both
    overall and per rec_type.

    Greedy, largest-component-first (classic bin-packing heuristic: the
    biggest, least-flexible items are placed while the most options are
    still open; small components fill in the remaining gaps at the end,
    where they can do the least damage to the achieved ratios). Ties
    broken by a seeded shuffle so the result is reproducible but not
    dependent on incidental dict/sort ordering.

    Scoring: for a candidate bucket, computes the ratio of (that bucket's
    count after tentatively adding this component) to (that bucket's own
    target count) -- for the OVERALL total and for EVERY rec_type this
    component touches -- and takes the WORST (highest) of these ratios as
    that bucket's score; the component is placed in whichever bucket
    minimizes this worst-case ratio (a standard greedy load-balancing /
    makespan-minimization heuristic). This is deliberately RATIO-based,
    not a squared deviation from the absolute target count: with very
    unequal target fractions (80:10:10), an absolute-deviation score
    systematically avoids the LARGEST bucket for early, large components
    -- overshooting the small 10% buckets looks numerically "cheaper" in
    absolute terms than undershooting the big 80% one, even though it is
    proportionally far worse. Caught by testing directly: an earlier,
    absolute-deviation version of this function put ZERO components in
    "train" across an entire test run, precisely because of this effect.
    Not an optimal solver (exact balanced multi-way partitioning under a
    hard grouping constraint is NP-hard in general) -- see module
    docstring for why the achieved ratios are always logged rather than
    assumed to match the targets.
    """
    footprints = _component_footprints(groups_by_component, rec_type_of)
    total_by_rec_type: Dict[str, int] = defaultdict(int)
    for k, rt in rec_type_of.items():
        total_by_rec_type[rt] += 1
    grand_total = len(rec_type_of)

    target_total = {b: f * grand_total for b, f in zip(bucket_names, bucket_target_fracs)}
    target_by_rec_type = {
        b: {rt: f * total_by_rec_type[rt] for rt in total_by_rec_type} for b, f in zip(bucket_names, bucket_target_fracs)
    }

    current_total: Dict[str, int] = {b: 0 for b in bucket_names}
    current_by_rec_type: Dict[str, Dict[str, int]] = {b: defaultdict(int) for b in bucket_names}

    rng = random.Random(seed)
    comp_ids = list(groups_by_component.keys())
    rng.shuffle(comp_ids)  # break ties randomly, reproducibly
    comp_ids.sort(key=lambda c: -footprints[c].get("__total__", 0))  # largest first (stable re: shuffle above)

    assignment: Dict[str, str] = {}
    for comp_id in comp_ids:
        fp = footprints[comp_id]
        best_bucket, best_score = None, None
        for b in bucket_names:
            tentative_total = current_total[b] + fp.get("__total__", 0)
            # target_total[b] > 0 whenever this bucket's own fraction > 0 and
            # there's at least one group overall -- both true in practice, but
            # guarded anyway rather than risking a ZeroDivisionError on some
            # future degenerate config (e.g. a bucket given a 0% target).
            ratio_overall = tentative_total / target_total[b] if target_total[b] > 0 else float("inf")
            worst_ratio = ratio_overall
            for rt, n in fp.items():
                if rt == "__total__":
                    continue
                tentative_rt = current_by_rec_type[b][rt] + n
                target_rt = target_by_rec_type[b].get(rt, 0.0)
                ratio_rt = tentative_rt / target_rt if target_rt > 0 else float("inf")
                worst_ratio = max(worst_ratio, ratio_rt)
            if best_score is None or worst_ratio < best_score:
                best_bucket, best_score = b, worst_ratio
        assignment[comp_id] = best_bucket
        current_total[best_bucket] += fp.get("__total__", 0)
        for rt, n in fp.items():
            if rt == "__total__":
                continue
            current_by_rec_type[best_bucket][rt] += n

    return assignment


# ══════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--labels", default=None, help="Override BASE/files/ml_labels.tsv")
    ap.add_argument("--taxonomy-file", default=None, help="Override BASE/files/f13_dataset.txt")
    ap.add_argument("--out", default=None, help="Override BASE/files/ml_labels_split.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_7"][SECTION_KEY]

    log = get_logger("2.7.2_split_dataset", paths["logs_dir"])

    labels_path = Path(args.labels) if args.labels else files_dir / "ml_labels.tsv"
    taxonomy_path = Path(args.taxonomy_file) if args.taxonomy_file else files_dir / "f13_dataset.txt"
    out_path = Path(args.out) if args.out else files_dir / "ml_labels_split.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Split dataset by duplicate group, species-exclusive (stratified by rec_type)")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_path, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_7_2_split_dataset"):
        return

    if not labels_path.exists():
        log.error(f"{labels_path} not found -- run 2.7.1_annotate_ml_labels.py first.")
        raise SystemExit(1)
    if not taxonomy_path.exists():
        log.error(f"{taxonomy_path} not found -- needed to resolve each MGE's species "
                  f"(specI) for leakage-safe splitting. Run 2.1.3_download_f13.py first, "
                  f"or pass --taxonomy-file explicitly.")
        raise SystemExit(1)

    train_frac = float(section_cfg.get("train_frac", 0.8))
    val_frac = float(section_cfg.get("val_frac", 0.1))
    test_frac = 1.0 - train_frac - val_frac
    seed = int(section_cfg.get("seed", 42))
    n_folds = int(section_cfg.get("n_folds", 5))
    fold_seed = int(section_cfg.get("fold_seed", 42))

    log.info(f"Split fractions: train={train_frac}  val={val_frac}  test={test_frac:.2f}  (seed={seed})")
    log.info(f"K-fold: n_folds={n_folds}  (fold_seed={fold_seed})")
    log.info(f"Taxonomy file  : {taxonomy_path}")
    log.info("-" * 70)

    t0 = time.time()

    with open(labels_path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    log.info(f"Loaded {len(rows):,} rows from {labels_path}")

    taxonomy_map = load_taxonomy_map(taxonomy_path)
    log.info(f"Loaded {len(taxonomy_map):,} GCA -> specI entries from {taxonomy_path}")

    # ── group -> set of member species, and group -> rec_type ──────────────
    group_species: Dict[GroupKey, Set[str]] = defaultdict(set)
    rec_type_of: Dict[GroupKey, str] = {}
    n_unresolved = 0
    for row in rows:
        key = (row["rec_type"], row["cluster_number"], row["duplicate_group_number"])
        rec_type_of[key] = row["rec_type"]
        species = resolve_species(row["mge_id"], taxonomy_map)
        if species.startswith("__unresolved"):
            n_unresolved += 1
        group_species[key].add(species)

    n_multi_species_groups = sum(1 for s in group_species.values() if len(s) > 1)
    log.info(f"{len(group_species):,} unique duplicate groups")
    log.info(f"  multi-species duplicate groups : {n_multi_species_groups:,} "
             f"(members span >1 species -- their species are forced into the same split component)")
    log.info(f"  rows with unresolved species    : {n_unresolved:,} "
             f"(given a synthetic, uniquely-its-own species so they're never wrongly grouped)")

    # ── species components (the atomic unit for both split and fold assignment) ─
    group_to_component = build_species_components(group_species)
    groups_by_component: Dict[str, List[GroupKey]] = defaultdict(list)
    for key, comp in group_to_component.items():
        groups_by_component[comp].append(key)
    log.info(f"  {len(groups_by_component):,} species component(s) after union-find "
             f"(a single-species duplicate group with a species seen nowhere else in the dataset "
             f"contributes its own singleton component)")
    log.info("-" * 70)

    # ── train/val/test assignment, per species component ───────────────────
    split_assignment_by_component = greedy_assign_components(
        groups_by_component, rec_type_of, ["train", "val", "test"], [train_frac, val_frac, test_frac], seed,
    )
    group_split: Dict[GroupKey, str] = {
        key: split_assignment_by_component[comp] for key, comp in group_to_component.items()
    }

    # ── fold assignment, per species component (independent seed/partition) ─
    fold_bucket_names = [f"fold{i}" for i in range(n_folds)]
    fold_assignment_by_component = greedy_assign_components(
        groups_by_component, rec_type_of, fold_bucket_names, [1.0 / n_folds] * n_folds, fold_seed,
    )
    group_fold: Dict[GroupKey, int] = {
        key: int(fold_assignment_by_component[comp][len("fold"):]) for key, comp in group_to_component.items()
    }

    # ── sanity check: hard species-exclusivity constraint actually holds ────
    species_to_splits: Dict[str, Set[str]] = defaultdict(set)
    species_to_folds: Dict[str, Set[int]] = defaultdict(set)
    for key, species_set in group_species.items():
        for s in species_set:
            species_to_splits[s].add(group_split[key])
            species_to_folds[s].add(group_fold[key])
    n_species_split_violations = sum(1 for s in species_to_splits.values() if len(s) > 1)
    n_species_fold_violations = sum(1 for s in species_to_folds.values() if len(s) > 1)
    if n_species_split_violations or n_species_fold_violations:
        log.error(f"*** INTERNAL ERROR: {n_species_split_violations} species span >1 split, "
                  f"{n_species_fold_violations} species span >1 fold -- this should be impossible "
                  f"given the union-find construction above; please report this.")
        raise SystemExit(1)
    log.info(f"Verified: no species appears in more than one split, and no species appears in more than one fold.")

    out_rows = []
    for row in rows:
        key = (row["rec_type"], row["cluster_number"], row["duplicate_group_number"])
        new_row = dict(row)
        new_row["split"] = group_split[key]
        new_row["fold"] = group_fold[key]
        out_rows.append(new_row)

    out_columns = list(rows[0].keys()) + ["split", "fold"]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=out_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info("SUMMARY -- ACHIEVED (not just requested) proportions, since the species constraint")
    log.info("can force real deviation from the nominal targets -- see module docstring.")
    log.info("-" * 70)
    groups_by_rec_type: Dict[str, List[GroupKey]] = defaultdict(list)
    for key, rt in rec_type_of.items():
        groups_by_rec_type[rt].append(key)

    log.info(f"{'rec_type':22s} {'groups':>8s} {'train':>8s} {'val':>8s} {'test':>8s}")
    any_empty_split = False
    for rec_type in sorted(groups_by_rec_type):
        keys = groups_by_rec_type[rec_type]
        n_train = sum(1 for k in keys if group_split[k] == "train")
        n_val = sum(1 for k in keys if group_split[k] == "val")
        n_test = sum(1 for k in keys if group_split[k] == "test")
        log.info(f"{rec_type:22s} {len(keys):>8,} {n_train:>8,} {n_val:>8,} {n_test:>8,}")
        if n_val == 0 or n_test == 0:
            any_empty_split = True
    if any_empty_split:
        log.info("")
        log.info("WARNING: at least one rec_type has 0 groups in val or test -- either too few")
        log.info("duplicate groups for that family, or its species are entangled (via cross-species")
        log.info("duplicate groups) with a large component that landed entirely in another split.")
        log.info("Check whether that's acceptable before training.")

    log.info("")
    log.info(f"{'rec_type':22s} {'groups':>8s} " + " ".join(f"{'fold'+str(f):>7s}" for f in range(n_folds)))
    any_empty_fold = False
    for rec_type in sorted(groups_by_rec_type):
        keys = groups_by_rec_type[rec_type]
        counts = [sum(1 for k in keys if group_fold[k] == f) for f in range(n_folds)]
        log.info(f"{rec_type:22s} {len(keys):>8,} " + " ".join(f"{c:>7,}" for c in counts))
        if any(c == 0 for c in counts):
            any_empty_fold = True
    if any_empty_fold:
        log.info("")
        log.info(f"WARNING: at least one rec_type has 0 groups in at least one of the {n_folds} folds "
                  f"-- same causes as the split warning above. That family's per-fold test metrics "
                  f"will be missing (not zero) for whichever fold(s) this affects.")

    log.info(f"-> {out_path}")
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()