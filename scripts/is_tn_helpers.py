#!/usr/bin/env python3
"""
is_tn_helpers.py
===================
Reusable helpers for flowchart1 (and anything else that needs the same
two operations):

  count_is_tn_from_gff(gff_path, taxonomy_map)
      Walk a master_mges.gff3-style file and total up every island whose
      'mge=' field contains an is_tn entry, REGARDLESS of mge_type (i.e.
      both nested and non-nested islands are included). Counting rule is
      the same one used everywhere else in this pipeline: an island's
      'mge=' field can bundle more than one recombinase of a type, e.g.
      'mge=is_tn:2,ser:1' is ONE row but represents 2 is_tn instances (the
      count after 'is_tn:') -- so n_instances sums that count across rows,
      it does not just count matching rows.

  species_from_mge_id_table(tsv_path, taxonomy_map, mge_id_column="mge_id")
      For any TSV that carries an mge_id column (e.g.
      filtered_mge_recombinase_table.tsv), pull the GCA accession out of
      each row's mge_id (text between the 1st and 3rd underscore, e.g.
      'MGE_GCA_001760395.1_...' -> 'GCA_001760395.1') and relate it to
      specI via the taxonomy table. Each row here is already a single
      recombinase (not an island with a count suffix), so no summing is
      needed -- just de-duplicate on mge_id and look up its genome.

Not runnable on its own -- import from a flowchart script.
"""

from __future__ import annotations

import csv

from utils import parse_gff_attributes, parse_mge_field, gca_from_mge_id


def count_is_tn_from_gff(gff_path: str, taxonomy_map: dict) -> dict:
    n_rows = 0
    n_instances = 0
    gca_set: set[str] = set()
    specI_to_gca: dict[str, set[str]] = {}
    specI_label: dict[str, str] = {}
    unmapped_gca = 0

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_gff_attributes(cols[8])
            counts = parse_mge_field(attrs.get("mge", ""))
            if "is_tn" not in counts:
                continue  # not an is_tn island -- nested status doesn't matter

            n_rows += 1
            n_instances += counts["is_tn"]

            mge_id = attrs.get("ID", "")
            gca = gca_from_mge_id(mge_id) if mge_id else None
            if not gca:
                continue
            gca_set.add(gca)
            info = taxonomy_map.get(gca)
            if info is None:
                unmapped_gca += 1
                continue
            specI, label = info
            specI_to_gca.setdefault(specI, set()).add(gca)
            specI_label.setdefault(specI, label)

    return {
        "n_rows": n_rows, "n_instances": n_instances, "gca_set": gca_set,
        "specI_to_gca": specI_to_gca, "specI_label": specI_label,
        "unmapped_gca": unmapped_gca,
    }


def species_from_mge_id_table(
    tsv_path: str, taxonomy_map: dict, mge_id_column: str = "mge_id"
) -> dict:
    gca_set: set[str] = set()
    specI_to_gca: dict[str, set[str]] = {}
    specI_label: dict[str, str] = {}
    unmapped_gca = 0
    seen_mge_ids: set[str] = set()

    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            mge_id = row.get(mge_id_column, "")
            if not mge_id or mge_id in seen_mge_ids:
                continue
            seen_mge_ids.add(mge_id)

            gca = gca_from_mge_id(mge_id)
            if not gca:
                continue
            gca_set.add(gca)
            info = taxonomy_map.get(gca)
            if info is None:
                unmapped_gca += 1
                continue
            specI, label = info
            specI_to_gca.setdefault(specI, set()).add(gca)
            specI_label.setdefault(specI, label)

    return {
        "gca_set": gca_set, "specI_to_gca": specI_to_gca,
        "specI_label": specI_label, "unmapped_gca": unmapped_gca,
        "n_unique_mge_ids": len(seen_mge_ids),
    }
