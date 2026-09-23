#!/usr/bin/env python3
"""
2.2.5_plot_flowchart_section_2_2.py  —  Section 2.2.5
=========================================================
input : BASE/files/master_mges.gff3, non_nested_is_tn.gff3 (2.1.1),
        BASE/files/f13_dataset.txt (2.1.3),
        BASE/files/a_mge_recombinase_table.tsv (2.2.4)
output: BASE/plots/flowchart_section_2_2.png

A 4-node vertical flowchart summarising Sections 2.1-2.2:

  Node 1  "MGEs"    -- live from master_mges.gff3, ALL MGE types, summing
                       every type:count pair per row.
  Node 2  "is_tn"   -- live from master_mges.gff3, is_tn only (both
                       nested and non-nested).
  Node 3  "is_tn"   -- live from non_nested_is_tn.gff3.
  Node 4  "is_tn"   -- NOW COMPUTED FOR REAL (no more given/placeholder
                       figures): from a_mge_recombinase_table.tsv (Section
                       2.2.4, already restricted to representative
                       genomes) -- is_tn count = # unique mge_id values
                       (column 1), genomes/species resolved from those
                       mge_ids' GCA accessions via the f13 taxonomy table.
                       The arrow into Node 4 states the actual ANI
                       threshold used for representative-picking
                       (Section 2.2.3's ani_threshold).

Side box (attached to Node 3): genomes-per-species distribution -- max and
min genome count, each with the actual specI ID(s), plus median and SD.

STYLING (see chat): bigger overall image (graph-level dpi=150, up from
graphviz's own default of 96), bigger font relative to box size
(fontsize 16 for nodes / 13 for edges, up from 11 / 10), and tighter box
padding + less inter-node whitespace (node margin 0.08,0.06 down from
0.15; nodesep/ranksep 0.35/0.45 down from 0.6/0.7) so the boxes and
arrows read as smaller relative to the text they carry, and less of the
image is empty space around the actual content.

create_plot: "yes"/"no" in config.yaml's section_2_2["2.2.5"] controls
whether this script does anything at all -- "no" skips it entirely.

Requires: pip install graphviz  (+ system 'dot' binary)

Run standalone:
    python 2.2.5_plot_flowchart_section_2_2.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import statistics
from pathlib import Path

import graphviz

from utils import (
    parse_gff_attributes, parse_mge_field, gca_from_mge_id, load_taxonomy_map,
    load_config, get_paths, get_logger,
)
from is_tn_helpers import count_is_tn_from_gff

SECTION_KEY = "2.2.5"
MAX_SPECI_NAMES_SHOWN = 3
OUT_STEM = "flowchart_section_2_2"
OUT_FORMAT = "png"

# Shared graphviz styling -- see module docstring for the reasoning behind each value.
GRAPH_ATTR = {"rankdir": "TB", "splines": "polyline", "nodesep": "0.35", "ranksep": "0.45", "dpi": "150"}
NODE_ATTR = {"shape": "box", "fontname": "Courier", "fontsize": "16", "margin": "0.08,0.06"}
EDGE_ATTR = {"fontname": "Helvetica", "fontsize": "13", "arrowsize": "0.7"}


def analyze_gff_all_types(path: Path, taxonomy_map: dict) -> dict:
    n_rows = 0
    n_instances = 0
    gca_set: set = set()
    specI_to_gca: dict = {}

    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_gff_attributes(cols[8])
            mge_id = attrs.get("ID", "")
            gca = gca_from_mge_id(mge_id) if mge_id else None

            n_rows += 1
            n_instances += sum(parse_mge_field(attrs.get("mge", "")).values())

            if not gca:
                continue
            gca_set.add(gca)
            info = taxonomy_map.get(gca)
            if info:
                specI, _label = info
                specI_to_gca.setdefault(specI, set()).add(gca)

    return {"n_rows": n_rows, "n_instances": n_instances, "gca_set": gca_set,
            "specI_to_gca": specI_to_gca}


def analyze_recombinase_table(path: Path, taxonomy_map: dict) -> dict:
    """Node 4: unique mge_id count from column 1, genomes/species via GCA lookup."""
    mge_ids: set = set()
    with open(path) as fh:
        header = fh.readline()
        for line in fh:
            if not line.strip():
                continue
            mge_id = line.split("\t")[0]
            mge_ids.add(mge_id)

    gca_set: set = set()
    specI_to_gca: dict = {}
    for mge_id in mge_ids:
        gca = gca_from_mge_id(mge_id)
        if not gca:
            continue
        gca_set.add(gca)
        info = taxonomy_map.get(gca)
        if info:
            specI, _label = info
            specI_to_gca.setdefault(specI, set()).add(gca)

    return {"n_instances": len(mge_ids), "gca_set": gca_set, "specI_to_gca": specI_to_gca}


def species_extremes(specI_to_gca: dict, top_n: int) -> dict:
    counts = {specI: len(gcas) for specI, gcas in specI_to_gca.items()}
    if not counts:
        return {}
    values = list(counts.values())
    max_v, min_v = max(values), min(values)
    max_specI = sorted([s for s, c in counts.items() if c == max_v])
    min_specI = sorted([s for s, c in counts.items() if c == min_v])
    return {
        "max_count": max_v, "max_specI": max_specI, "n_at_max": len(max_specI),
        "min_count": min_v, "min_specI": min_specI, "n_at_min": len(min_specI),
        "median": statistics.median(values),
        "stdev": statistics.stdev(values) if len(values) > 1 else 0.0,
        "n_species": len(counts),
    }


def fmt(n) -> str:
    return f"{n:,}" if isinstance(n, int) else str(n)


def box_label(lines: list[str]) -> str:
    return "".join(f"{line}\\l" for line in lines)


def specI_ids_str(specI_list: list[str], top_n: int, n_total: int) -> str:
    shown = specI_list[:top_n]
    s = ", ".join(shown)
    if n_total > len(shown):
        s += f" (+{n_total - len(shown)} more)"
    return s


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    plots_dir = paths["plots_dir"]
    section_cfg = cfg["section_2_2"][SECTION_KEY]

    log = get_logger("2.2.5_plot_flowchart_section_2_2", paths["logs_dir"])

    if str(section_cfg.get("create_plot", "no")).lower() != "yes":
        log.info(f"Section {SECTION_KEY}: create_plot != 'yes' -> plotting skipped, "
                 f"nothing generated or saved.")
        return

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Sections 2.1-2.2 overview flowchart")
    log.info("=" * 70)

    master_gff = files_dir / "master_mges.gff3"
    non_nested_gff = files_dir / "non_nested_is_tn.gff3"
    taxonomy_file = files_dir / "f13_dataset.txt"
    recomb_table = files_dir / "a_mge_recombinase_table.tsv"

    for p in (master_gff, non_nested_gff, taxonomy_file, recomb_table):
        if not p.exists():
            log.error(f"{p} not found -- run the earlier Section 2.1/2.2 scripts first.")
            raise SystemExit(1)

    log.info("Loading taxonomy map ...")
    taxonomy_map = load_taxonomy_map(taxonomy_file)
    log.info(f"  {len(taxonomy_map):,} GCA -> specI entries")

    log.info(f"Analyzing Node 1 (all MGE types): {master_gff}")
    node1 = analyze_gff_all_types(master_gff, taxonomy_map)
    log.info(f"  instances={node1['n_instances']:,} genomes={len(node1['gca_set']):,} "
             f"species={len(node1['specI_to_gca']):,}")

    log.info(f"Analyzing Node 2 (is_tn, nested+non-nested): {master_gff}")
    node2 = count_is_tn_from_gff(str(master_gff), taxonomy_map)
    log.info(f"  instances={node2['n_instances']:,} genomes={len(node2['gca_set']):,} "
             f"species={len(node2['specI_to_gca']):,}")

    log.info(f"Analyzing Node 3 (non-nested is_tn): {non_nested_gff}")
    node3 = analyze_gff_all_types(non_nested_gff, taxonomy_map)
    log.info(f"  instances={node3['n_instances']:,} genomes={len(node3['gca_set']):,} "
             f"species={len(node3['specI_to_gca']):,}")

    log.info(f"Analyzing Node 4 (post-dedup, from 2.2.4's table): {recomb_table}")
    node4 = analyze_recombinase_table(recomb_table, taxonomy_map)
    log.info(f"  instances={node4['n_instances']:,} genomes={len(node4['gca_set']):,} "
             f"species={len(node4['specI_to_gca']):,}")

    ani_threshold = cfg.get("section_2_2", {}).get("2.2.3", {}).get("ani_threshold", "?")
    edge4_label = f"reduced genomic redundancy ({ani_threshold}% ANI)"

    dist = species_extremes(node3["specI_to_gca"], MAX_SPECI_NAMES_SHOWN)
    if dist:
        log.info(
            f"  Node3 genomes/species: max={dist['max_count']} specI={dist['max_specI']} "
            f"min={dist['min_count']} specI={dist['min_specI']} "
            f"median={dist['median']:.1f} stdev={dist['stdev']:.2f}"
        )

    # ── build graph ─────────────────────────────────────────────────────
    dot = graphviz.Digraph(
        "flowchart_section_2_2",
        format=OUT_FORMAT,
        graph_attr=GRAPH_ATTR,
        node_attr=NODE_ATTR,
        edge_attr=EDGE_ATTR,
    )

    dot.node("node1", box_label([
        f"{fmt(node1['n_instances'])} MGEs",
        f"{fmt(len(node1['gca_set']))} genomes",
        f"{fmt(len(node1['specI_to_gca']))} species (specI)",
    ]), style="filled", fillcolor="#d9f2d9", color="#4a934a")

    dot.node("node2", box_label([
        f"{fmt(node2['n_instances'])} is_tn",
        f"{fmt(len(node2['gca_set']))} genomes",
        f"{fmt(len(node2['specI_to_gca']))} species (specI)",
    ]), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    dot.node("node3", box_label([
        f"{fmt(node3['n_instances'])} is_tn",
        f"{fmt(len(node3['gca_set']))} genomes",
        f"{fmt(len(node3['specI_to_gca']))} species (specI)",
    ]), style="filled", fillcolor="#dbe9fb", color="#4a75c4")

    dot.node("node4", box_label([
        f"{fmt(node4['n_instances'])} is_tn",
        f"{fmt(len(node4['gca_set']))} genomes",
        f"{fmt(len(node4['specI_to_gca']))} species (specI)",
    ]), style="filled", fillcolor="#fbdada", color="#c1494a")

    if dist:
        side_lines = [
            "genomes per species",
            f"max - {dist['max_count']}  "
            f"specI: {specI_ids_str(dist['max_specI'], MAX_SPECI_NAMES_SHOWN, dist['n_at_max'])}",
            f"min - {dist['min_count']}  "
            f"specI: {specI_ids_str(dist['min_specI'], MAX_SPECI_NAMES_SHOWN, dist['n_at_min'])}",
            f"median  {dist['median']:.1f}",
            f"SD      {dist['stdev']:.2f}",
            f"n species = {dist['n_species']:,}",
        ]
    else:
        side_lines = ["genomes per species", "(no data)"]

    dot.node("sidebox", box_label(side_lines),
              shape="note", style="filled", fillcolor="#fdf1cf", color="#c99a2e")

    dot.edge("node1", "node2")
    dot.edge("node2", "node3", label="non-nested")
    dot.edge("node3", "node4", label=edge4_label)
    dot.edge("node3", "sidebox", constraint="false")
    with dot.subgraph() as s:
        s.attr(rank="same")
        s.node("node3")
        s.node("sidebox")

    out_stem = str(plots_dir / OUT_STEM)
    rendered = dot.render(filename=out_stem, cleanup=True)
    dot.save(filename=f"{out_stem}.gv")
    log.info(f"Wrote {rendered}")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()