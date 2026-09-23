#!/usr/bin/env python3
"""
2.3.3_plot_flowchart_section_2_3.py  —  Section 2.3.3
==========================================================
input : BASE/files/a_mge_recombinase_table.tsv (2.2.4, used as "filtered"),
        BASE/files/b_mge_recombinase_table.tsv (2.3.2, used as "rearranged"),
        BASE/files/f13_dataset.txt (2.1.3)
output: BASE/plots/flowchart_section_2_3.png,
        BASE/plots/histogram_clusters_vs_singletons_by_rectype_2_3.png

Same shape as the earlier standalone "flowchart2_and_histogram" figure,
now wired into the real pipeline:

  Node 1  "is_tn"  -- LIVE from a_mge_recombinase_table.tsv (unique
                      mge_ids -> count; genomes/species via GCA -> f13).
                      Previously a given/hardcoded figure -- now that
                      Section 2.2.4 produces this table for real, hardcoding
                      it no longer made sense.
  Node 2  "is_tn"  -- LIVE from a_mge_recombinase_table.tsv:
                        * # mge_ids with MORE THAN ONE recombinase row
                        * # mge_ids with EXACTLY ONE recombinase row
                        * # unique recombinase_type values
  Node 3  "is_tn"  -- reached via an edge labeled "<ani>% ANI clustering"
                      (the actual configured threshold, not a fixed
                      label): from b_mge_recombinase_table.tsv,
                        * # of clusters
                        * # of singletons (clusters with exactly one member)

BUG FIX vs. the original standalone version: cluster_number is assigned
*within each recombinase_type* (restarts at 1 for every type), so
"cluster 1" of one type and "cluster 1" of another type are NOT the same
cluster. The original singleton-counting logic keyed on bare
cluster_number alone, which could conflate them. This version keys
cluster identity on (recombinase_type, cluster_number) throughout.

Additional, new: after all three nodes, also reports and annotates the
number of is_tn mge_ids that were NOT assigned to any cluster at all
(every one of their recombinase rows has is_representative == "NA|NA").

STYLING (see chat, same change applied uniformly across 2.2.5/2.3.3/2.4.3):
bigger overall image (graph-level dpi=150, up from graphviz's own default
of 96), bigger font relative to box size (fontsize 16 for nodes / 13 for
edges, up from 11 / 10), and tighter box padding + less inter-node
whitespace (node margin 0.08,0.06 down from 0.15; nodesep/ranksep
0.35/0.45 down from 0.6/0.7).

Requires: pip install graphviz matplotlib  (+ system 'dot' binary)

Run standalone:
    python 2.3.3_plot_flowchart_section_2_3.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import graphviz
import matplotlib.pyplot as plt

from utils import get_logger, get_paths, load_config, load_taxonomy_map, gca_from_mge_id

SECTION_KEY = "2.3.3"
OUT_STEM = "flowchart_section_2_3"
HIST_NAME = "histogram_clusters_vs_singletons_by_rectype_2_3.png"
OUT_FORMAT = "png"

# Shared graphviz styling -- see module docstring for the reasoning behind each value.
GRAPH_ATTR = {"rankdir": "TB", "splines": "polyline", "nodesep": "0.35", "ranksep": "0.45", "dpi": "150"}
NODE_ATTR = {"shape": "box", "fontname": "Courier", "fontsize": "16", "margin": "0.08,0.06"}
EDGE_ATTR = {"fontname": "Helvetica", "fontsize": "13", "arrowsize": "0.7"}


def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def analyze_a_table(rows: list[dict], taxonomy_map: dict) -> dict:
    mge_ids = {r["mge_id"] for r in rows if r.get("mge_id")}
    rows_per_mge = Counter(r["mge_id"] for r in rows)
    n_multi = sum(1 for c in rows_per_mge.values() if c > 1)
    n_single = sum(1 for c in rows_per_mge.values() if c == 1)
    rec_types = {r["recombinase_type"] for r in rows if r.get("recombinase_type")}

    mge_ids_by_type: dict = defaultdict(set)
    for r in rows:
        t = r.get("recombinase_type")
        if t:
            mge_ids_by_type[t].add(r["mge_id"])

    gca_set = set()
    specI_set = set()
    for mid in mge_ids:
        gca = gca_from_mge_id(mid)
        if not gca:
            continue
        gca_set.add(gca)
        info = taxonomy_map.get(gca)
        if info:
            specI_set.add(info[0])

    return {
        "n_mge_ids": len(mge_ids), "n_multi": n_multi, "n_single": n_single,
        "rec_types": rec_types, "mge_ids_by_type": mge_ids_by_type,
        "gca_set": gca_set, "specI_set": specI_set,
    }


def analyze_b_table(rows: list[dict]) -> dict:
    """Cluster identity is scoped to (recombinase_type, cluster_number) --
    see module docstring for why bare cluster_number isn't safe."""
    cluster_size: Counter = Counter()
    clusters_by_type: dict = defaultdict(set)
    mge_ids_by_type: dict = defaultdict(set)
    unassigned_mge_ids: set = set()
    mge_id_has_cluster: dict = defaultdict(bool)

    for r in rows:
        t = r.get("recombinase_type", "")
        c = r.get("cluster_number", "NA")
        mge_id = r.get("mge_id", "")
        is_rep = r.get("is_representative", "NA|NA")

        if t:
            mge_ids_by_type[t].add(mge_id)

        if is_rep != "NA|NA" and c and c != "NA":
            cluster_size[(t, c)] += 1
            clusters_by_type[t].add(c)
            mge_id_has_cluster[mge_id] = True
        else:
            mge_id_has_cluster.setdefault(mge_id, False)

    for mge_id, has_cluster in mge_id_has_cluster.items():
        if not has_cluster:
            unassigned_mge_ids.add(mge_id)

    singleton_keys = {k for k, n in cluster_size.items() if n == 1}
    singleton_clusters_by_type: dict = defaultdict(set)
    for (t, c) in singleton_keys:
        singleton_clusters_by_type[t].add(c)

    return {
        "cluster_size": cluster_size,
        "n_clusters": len(cluster_size),
        "n_singletons": len(singleton_keys),
        "clusters_by_type": clusters_by_type,
        "singleton_clusters_by_type": singleton_clusters_by_type,
        "mge_ids_by_type": mge_ids_by_type,
        "unassigned_mge_ids": unassigned_mge_ids,
    }


def fmt(n) -> str:
    return f"{n:,}" if isinstance(n, int) else str(n)


def box_label(lines: list[str]) -> str:
    return "".join(f"{line}\\l" for line in lines)


def build_flowchart(node1, node2, b_result, ani_threshold, out_stem) -> None:
    dot = graphviz.Digraph(
        "flowchart_section_2_3",
        format=OUT_FORMAT,
        graph_attr=GRAPH_ATTR,
        node_attr=NODE_ATTR,
        edge_attr=EDGE_ATTR,
    )

    dot.node("node1", box_label([
        f"{fmt(node1['n_mge_ids'])} is_tn",
        f"{fmt(len(node1['gca_set']))} genomes",
        f"{fmt(len(node1['specI_set']))} species (specI)",
    ]), style="filled", fillcolor="#fbdada", color="#c1494a")

    dot.node("node2", box_label([
        f"{fmt(node2['n_multi'])} mge_ids with >1 recombinase",
        f"{fmt(node2['n_single'])} mge_ids with 1 recombinase",
        f"{fmt(len(node2['rec_types']))} unique recombinase types",
    ]), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    dot.node("node3", box_label([
        f"# of clusters = {fmt(b_result['n_clusters'])}",
        f"# of singletons = {fmt(b_result['n_singletons'])}",
    ]), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    dot.edge("node1", "node2")
    dot.edge("node2", "node3", label=f"{ani_threshold * 100:.0f}% ANI clustering")

    rendered = dot.render(filename=out_stem, cleanup=True)
    dot.save(filename=f"{out_stem}.gv")
    print(f"Wrote {rendered}")
    print(f"Wrote {out_stem}.gv")


def build_histogram(b_result: dict, a_result: dict, out_path: str) -> None:
    rec_types = sorted(b_result["clusters_by_type"].keys())
    n_clusters = [len(b_result["clusters_by_type"][t]) for t in rec_types]
    n_singletons = [len(b_result["singleton_clusters_by_type"][t]) for t in rec_types]
    n_unique_mges = [
        len(a_result["mge_ids_by_type"].get(t, set()) | b_result["mge_ids_by_type"].get(t, set()))
        for t in rec_types
    ]

    x = range(len(rec_types))
    width = 0.38

    fig, ax = plt.subplots(figsize=(max(8, len(rec_types) * 1.3), 6))
    ax.bar([i - width / 2 for i in x], n_clusters, width, label="# clusters", color="#4a75c4")
    ax.bar([i + width / 2 for i in x], n_singletons, width, label="# singletons", color="#c99a2e")

    ax.set_xticks(list(x))
    ax.set_xticklabels(rec_types, rotation=45, ha="right")
    ax.set_ylabel("count")
    ax.set_title("Clusters vs. singletons per recombinase type")
    ax.legend()

    top_scale = max(n_clusters + n_singletons + [1])
    for i, n_mge in zip(x, n_unique_mges):
        top = max(n_clusters[i], n_singletons[i])
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
    section_cfg = cfg["section_2_3"][SECTION_KEY]

    log = get_logger("2.3.3_plot_flowchart_section_2_3", paths["logs_dir"])

    if str(section_cfg.get("create_plot", "no")).lower() != "yes":
        log.info(f"Section {SECTION_KEY}: create_plot != 'yes' -> plotting skipped.")
        return

    a_table = files_dir / "a_mge_recombinase_table.tsv"
    b_table = files_dir / "b_mge_recombinase_table.tsv"
    taxonomy_file = files_dir / "f13_dataset.txt"
    ani_threshold = cfg["section_2_3"]["2.3.1"]["ani_threshold"]

    for p in (a_table, b_table, taxonomy_file):
        if not p.exists():
            log.error(f"{p} not found -- run the earlier Section 2.3 scripts first.")
            raise SystemExit(1)

    log.info("Loading taxonomy map ...")
    taxonomy_map = load_taxonomy_map(taxonomy_file)

    log.info(f"Analyzing Node 1/2 from {a_table}")
    a_rows = read_tsv(a_table)
    node1 = analyze_a_table(a_rows, taxonomy_map)
    node2 = node1  # same source table drives both nodes' underlying rows
    log.info(f"  is_tn={node1['n_mge_ids']:,} genomes={len(node1['gca_set']):,} "
             f"species={len(node1['specI_set']):,} n_multi={node1['n_multi']:,} "
             f"n_single={node1['n_single']:,} rec_types={len(node1['rec_types']):,}")

    log.info(f"Analyzing Node 3 from {b_table}")
    b_rows = read_tsv(b_table)
    b_result = analyze_b_table(b_rows)
    log.info(f"  clusters={b_result['n_clusters']:,} singletons={b_result['n_singletons']:,} "
             f"unassigned_is_tns={len(b_result['unassigned_mge_ids']):,}")

    out_stem = str(plots_dir / OUT_STEM)
    build_flowchart(node1, node2, b_result, ani_threshold, out_stem)

    hist_path = str(plots_dir / HIST_NAME)
    build_histogram(b_result, node1, hist_path)

    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()