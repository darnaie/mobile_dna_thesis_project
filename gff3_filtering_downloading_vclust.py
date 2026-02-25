#!/usr/bin/env python3

import os
import subprocess
from collections import Counter

# =========================
# Paths
# =========================
INPUT_GFF = "/net/bq-storage/ag-khedkar/Sofia/project_folder/master_mges.gff3"
BASE_DIR = os.path.dirname(INPUT_GFF)

FILTERED_GFF = os.path.join(BASE_DIR, "filtered_is_tns_length.gff3")

GENOME_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/genome_seqs"
DEDUP_DIR = os.path.join(GENOME_DIR, "deduplicated")

VCLUST = "/home/bq_slipaeva/.local/bin/vclust"

os.makedirs(GENOME_DIR, exist_ok=True)
os.makedirs(DEDUP_DIR, exist_ok=True)

# =========================
# Helpers
# =========================
def parse_attributes(attr_string):
    """
    Parse GFF3 attribute column into a dict
    """
    attrs = {}
    for item in attr_string.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs


# =========================
# Step 1: Read and filter GFF3
# =========================
rows = []
lengths = []

total_rows = 0

with open(INPUT_GFF) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        total_rows += 1
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue

        length = int(cols[5])
        attrs = parse_attributes(cols[8])

        # Filter: ONLY pure is_tn
        if "mge" not in attrs:
            continue

        mge_val = attrs["mge"]

        # must contain is_tn AND no comma (no other MGE types)
        if "is_tn" not in mge_val:
            continue
        if "," in mge_val:
            continue

        rows.append((line, length))
        lengths.append(length)

after_is_tn = len(rows)

# Keep only lengths that occur more than once
length_counts = Counter(lengths)
valid_lengths = {l for l, c in length_counts.items() if c > 1}

filtered_rows = [line for line, l in rows if l in valid_lengths]

with open(FILTERED_GFF, "w") as out:
    for line in filtered_rows:
        out.write(line)

print("=== GFF3 filtering summary ===")
print(f"Total rows read: {total_rows}")
print(f"Rows with pure is_tn: {after_is_tn}")
print(f"Rows kept (duplicate lengths): {len(filtered_rows)}")
print(f"Rows filtered out: {total_rows - len(filtered_rows)}")
print(f"Filtered file written to: {FILTERED_GFF}")

# =========================
# Step 2: Download genome sequences
# =========================
unique_links = set()

for line in filtered_rows:
    first_col = line.split("\t")[0]
    parts = first_col.split(".")
    if len(parts) < 2:
        continue

    full_code = f"{parts[0]}.{parts[1]}"
    first_code = parts[0]

    url = (
        "http://progenomes3.embl.de/dumpSequence.cgi"
        f"?p={full_code}&t=c&a={first_code}"
    )
    unique_links.add(url)

print(f"Unique genome links to download: {len(unique_links)}")

for url in sorted(unique_links):
    print(f"Downloading: {url}")
    subprocess.run(
        ["wget", "-c", url],
        cwd=GENOME_DIR,
        check=False
    )

# =========================
# Step 3: vclust deduplication
# =========================
fasta_files = [
    f for f in os.listdir(GENOME_DIR)
    if f.endswith(".fa") or f.endswith(".fasta") or f.endswith(".fna")
]

if not fasta_files:
    print("No FASTA files found. Skipping vclust deduplication.")
else:
    print("Running vclust deduplication...")
    subprocess.run(
        [
            VCLUST,
            "deduplicate",
            "-i", GENOME_DIR,
            "-o", DEDUP_DIR
        ],
        check=True
    )
    print("vclust deduplication completed.")
    print(f"Deduplicated genomes written to: {DEDUP_DIR}")
