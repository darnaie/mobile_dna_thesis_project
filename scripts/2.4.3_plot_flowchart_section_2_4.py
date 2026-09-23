#!/usr/bin/env python3
"""
2.4.3_plot_flowchart_section_2_4.py  —  Section 2.4.3
==========================================================
input : BASE/files/c_mge_recombinase_table.tsv (Section 2.4.2),
        BASE/vclust_deduplication/*.fasta.duplicates.txt (Section 2.4.1)
output: BASE/plots/flowchart_section_2_4.png,
        BASE/plots/histogram_duplicate_group_sizes_2_4.png,
        BASE/plots/histogram_duplicate_groups_vs_singletons_by_rectype_2_4.png

Same shape as the earlier standalone "flowchart3_duplicate_groups"
figure, now wired into the real pipeline:

  Node 1  "is_tn"  -- LIVE, unique mge_id count in c_mge_recombinase_table.tsv
                      (was a given/hardcoded figure before -- no longer
                      needed now that the real upstream table exists).
  Node 2  "# of duplicate groups" / "# of singletons" -- LIVE, straight off
          the singleton/duplicate_group_number columns Section 2.4.2 added
          (no need to re-cross-reference the raw duplicates file for this
          part -- that bookkeeping already happened in 2.4.2).
  Node 3  every recombinase type with, PER TYPE, THREE numbers now (see
          "NEW" below): # of clusters / # of duplicate groups / # of
          singletons.
  Node 4  (new) how many non-singleton (singleton == "no") mge_ids are
          present in c_mge_recombinase_table.tsv overall.

NEW (see chat): Node 3's per-rec_type lines used to show two numbers
(duplicate groups / singletons). A third, PRECEDING number -- # of
clusters, i.e. how many Section 2.3 (MMseqs2, ~85% identity) clusters
that rec_type has, upstream of and coarser than Section 2.4's own
duplicate-group deduplication (~99-100% identity, WITHIN each cluster)
-- is now shown too, since one cluster routinely contains several
duplicate groups plus possibly some singletons, and reporting only the
finer-grained numbers on their own hid that nesting. No new input file
is needed for this: c_mge_recombinase_table.tsv already carries its own
cluster_number column, inherited unchanged from Section 2.3.2's
b_mge_recombinase_table.tsv (2.4.2 appends onto b's rows rather than
replacing them -- see the table lineage this pipeline has followed
throughout). Rows with cluster_number "NA" (unclustered) are excluded
from the cluster count, same "NA sentinel" convention used everywhere
else in this pipeline.

Histograms:
  1. Group-size distribution over ALL groups in the raw duplicates file
     (not just is_tn) -- this one DOES need the raw file, since
     c_mge_recombinase_table.tsv only ever reflects is_tn-relevant groups.
  2. Per-recombinase-type duplicate groups vs. singletons, annotated with
     total is_tn count per type (from c_mge_recombinase_table.tsv alone).

STYLING (see chat, same change applied uniformly across 2.2.5/2.3.3/2.4.3):
bigger overall image (graph-level dpi=150, up from graphviz's own default
of 96), bigger font relative to box size (fontsize 16 for nodes / 13 for
edges, up from 11 / 10), and tighter box padding + less inter-node
whitespace (node margin 0.08,0.06 down from 0.15; nodesep/ranksep
0.35/0.45 down from 0.6/0.7).

Requires: pip install graphviz matplotlib  (+ system 'dot' binary)

Run standalone:
    python 2.4.3_plot_flowchart_section_2_4.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import graphviz
import matplotlib.pyplot as plt

from utils import get_logger, get_paths, load_config, parse_duplicate_groups_file

SECTION_KEY = "2.4.3"
OUT_STEM = "flowchart_section_2_4"
GROUP_SIZE_HIST_NAME = "histogram_duplicate_group_sizes_2_4.png"
TYPE_HIST_NAME = "histogram_duplicate_groups_vs_singletons_by_rectype_2_4.png"
OUT_FORMAT = "png"

# Shared graphviz styling -- see module docstring for the reasoning behind each value.
GRAPH_ATTR = {"rankdir": "TB", "splines": "polyline", "nodesep": "0.35", "ranksep": "0.45", "dpi": "150"}
NODE_ATTR = {"shape": "box", "fontname": "Courier", "fontsize": "16", "margin": "0.08,0.06"}
EDGE_ATTR = {"fontname": "Helvetica", "fontsize": "13", "arrowsize": "0.7"}


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def find_duplicates_file(base_dir: Path) -> Path | None:
    vclust_dir = base_dir / "vclust_deduplication"
    if not vclust_dir.exists():
        return None
    candidates = sorted(vclust_dir.glob("*.duplicates.txt"))
    return candidates[0] if candidates else None


def analyze_c_table(rows: list[dict]) -> dict:
    mge_ids = {r["mge_id"] for r in rows if r.get("mge_id")}

    # dedupe to one row per mge_id for group/singleton/cluster bookkeeping
    # (a multi-recombinase mge_id would otherwise be counted more than once)
    seen: dict = {}
    for r in rows:
        mid = r.get("mge_id")
        if mid and mid not in seen:
            seen[mid] = r

    n_groups_seen = set()
    n_singletons = 0
    n_non_singleton = 0
    groups_by_type: dict = defaultdict(set)
    singletons_by_type: dict = defaultdict(set)
    mge_ids_by_type: dict = defaultdict(set)
    # NEW: cluster_number is inherited unchanged from Section 2.3.2's
    # b_mge_recombinase_table.tsv (2.4.2 appends onto b's rows, doesn't
    # replace them) -- so it's already right here in every row, no
    # separate file needed. Keyed on (rec_type, cluster_number), same
    # "scope cluster identity to its own rec_type" reasoning 2.3.3 already
    # applies -- cluster "1" of one type and cluster "1" of another are
    # not the same cluster. "NA" (unclustered) is excluded, same
    # NA-sentinel convention used everywhere else in this pipeline.
    clusters_by_type: dict = defaultdict(set)

    for mid, r in seen.items():
        t = r.get("recombinase_type", "")
        mge_ids_by_type[t].add(mid)

        cluster_num = r.get("cluster_number")
        if cluster_num and cluster_num != "NA":
            clusters_by_type[t].add(cluster_num)

        if r.get("singleton") == "yes":
            n_singletons += 1
            singletons_by_type[t].add(mid)
        else:
            n_non_singleton += 1
            grp = r.get("duplicate_group_number")
            if grp and grp != "NA":
                n_groups_seen.add(grp)
                groups_by_type[t].add(grp)

    return {
        "n_mge_ids": len(mge_ids),
        "n_groups": len(n_groups_seen),
        "n_singletons": n_singletons,
        "n_non_singleton": n_non_singleton,
        "groups_by_type": groups_by_type,
        "singletons_by_type": singletons_by_type,
        "mge_ids_by_type": mge_ids_by_type,
        "clusters_by_type": clusters_by_type,
    }


def fmt(n) -> str:
    return f"{n:,}" if isinstance(n, int) else str(n)


def box_label(lines: list[str]) -> str:
    return "".join(f"{line}\\l" for line in lines)


def build_flowchart(result: dict, out_stem: str) -> None:
    dot = graphviz.Digraph(
        "flowchart_section_2_4",
        format=OUT_FORMAT,
        graph_attr=GRAPH_ATTR,
        node_attr=NODE_ATTR,
        edge_attr=EDGE_ATTR,
    )

    dot.node("node1", box_label([f"{fmt(result['n_mge_ids'])} is_tn"]),
              style="filled", fillcolor="#fbdada", color="#c1494a")

    dot.node("node2", box_label([
        f"# of duplicate groups = {fmt(result['n_groups'])}",
        f"# of singletons = {fmt(result['n_singletons'])}",
    ]), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    types_sorted = sorted(
        set(result["clusters_by_type"]) | set(result["groups_by_type"]) | set(result["singletons_by_type"]),
        key=lambda t: (
            len(result["clusters_by_type"].get(t, ()))
            + len(result["groups_by_type"].get(t, ()))
            + len(result["singletons_by_type"].get(t, ()))
        ),
        reverse=True,
    )
    # NEW: three numbers per rec_type now -- clusters / groups / singletons,
    # in that order (coarsest to finest grouping), see module docstring.
    node3_lines = ["recombinase type : clusters / groups / singletons"]
    for t in types_sorted:
        nc = len(result["clusters_by_type"].get(t, ()))
        ng = len(result["groups_by_type"].get(t, ()))
        ns = len(result["singletons_by_type"].get(t, ()))
        node3_lines.append(f"{t:<20} {nc:>6} / {ng:>6} / {ns:<6}")
    dot.node("node3", box_label(node3_lines), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    dot.node("node4", box_label([
        "non-singletons present",
        f"{fmt(result['n_non_singleton'])}",
    ]), style="filled", fillcolor="#d9f2d9", color="#4a934a")

    dot.edge("node1", "node2", label="duplicate groups")
    dot.edge("node2", "node3")
    dot.edge("node3", "node4")

    rendered = dot.render(filename=out_stem, cleanup=True)
    dot.save(filename=f"{out_stem}.gv")
    print(f"Wrote {rendered}")
    print(f"Wrote {out_stem}.gv")


def build_group_size_histogram(groups: list[set], out_path: str) -> None:
    sizes = [len(g) for g in groups]
    if not sizes:
        print("  No duplicate groups found -- skipping group-size histogram.")
        return
    size_counts = Counter(sizes)
    min_s, max_s = min(sizes), max(sizes)

    fig, ax = plt.subplots(figsize=(9, 6))
    xs = list(range(min_s, max_s + 1))
    ys = [size_counts.get(x, 0) for x in xs]
    ax.bar(xs, ys, color="#4a75c4")
    ax.set_yscale("log")
    ax.set_xlabel("duplicate group size (# of MGEs in the group)")
    ax.set_ylabel("number of duplicate groups (log scale)")
    ax.set_title(f"Duplicate group size distribution (all {len(groups):,} groups, all MGE types)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def build_type_histogram(result: dict, out_path: str) -> None:
    rec_types = sorted(
        set(result["groups_by_type"]) | set(result["singletons_by_type"]) | set(result["mge_ids_by_type"])
    )
    if not rec_types:
        print("  No recombinase types found -- skipping type histogram.")
        return

    n_groups = [len(result["groups_by_type"].get(t, ())) for t in rec_types]
    n_singletons = [len(result["singletons_by_type"].get(t, ())) for t in rec_types]
    n_total_is_tn = [len(result["mge_ids_by_type"].get(t, ())) for t in rec_types]

    x = range(len(rec_types))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(8, len(rec_types) * 1.3), 6))
    ax.bar([i - width / 2 for i in x], n_groups, width, label="# duplicate groups", color="#4a75c4")
    ax.bar([i + width / 2 for i in x], n_singletons, width, label="# singletons", color="#c99a2e")

    ax.set_xticks(list(x))
    ax.set_xticklabels(rec_types, rotation=45, ha="right")
    ax.set_ylabel("count")
    ax.set_title("Duplicate groups vs. singletons per recombinase type")
    ax.legend()

    top_scale = max(n_groups + n_singletons + [1])
    for i, n_mge in zip(x, n_total_is_tn):
        top = max(n_groups[i], n_singletons[i])
        ax.text(i, top + top_scale * 0.02, f"{n_mge:,} is_tns", ha="center", va="bottom", fontsize=9)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
    print(f"Wrote {out_path}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    plots_dir = paths["plots_dir"]
    section_cfg = cfg["section_2_4"][SECTION_KEY]

    log = get_logger("2.4.3_plot_flowchart_section_2_4", paths["logs_dir"])

    if str(section_cfg.get("create_plot", "no")).lower() != "yes":
        log.info(f"Section {SECTION_KEY}: create_plot != 'yes' -> plotting skipped.")
        return

    c_table = files_dir / "c_mge_recombinase_table.tsv"
    duplicates_file = find_duplicates_file(paths["base_dir"])

    if not c_table.exists():
        log.error(f"{c_table} not found -- run append_duplicate_info.py (Section 2.4.2) first.")
        raise SystemExit(1)
    if duplicates_file is None:
        log.error(f"No *.duplicates.txt found under {paths['base_dir'] / 'vclust_deduplication'} "
                   f"-- run vclust_dedup.py (Section 2.4.1) first.")
        raise SystemExit(1)

    log.info(f"Analyzing {c_table}")
    rows = read_tsv(c_table)
    result = analyze_c_table(rows)
    log.info(f"  is_tn={result['n_mge_ids']:,} groups={result['n_groups']:,} "
             f"singletons={result['n_singletons']:,} non_singleton={result['n_non_singleton']:,}")
    n_clusters_total = len({(t, c) for t, cs in result["clusters_by_type"].items() for c in cs})
    log.info(f"  clusters (all types, from cluster_number)={n_clusters_total:,}")

    log.info(f"Loading all duplicate groups (all MGE types) from {duplicates_file}")
    all_groups = parse_duplicate_groups_file(duplicates_file)
    log.info(f"  {len(all_groups):,} total groups (all MGE types)")

    out_stem = str(plots_dir / OUT_STEM)
    build_flowchart(result, out_stem)
    build_group_size_histogram(all_groups, str(plots_dir / GROUP_SIZE_HIST_NAME))
    build_type_histogram(result, str(plots_dir / TYPE_HIST_NAME))

    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()