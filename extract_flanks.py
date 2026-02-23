#!/usr/bin/env python3
import os
from Bio import SeqIO

# ---------- SETTINGS ----------
DUP_FILE = "/net/bq-storage/ag-khedkar/Sofia/project_folder/mge_seqs_deduplicated/deduplicated.fasta.duplicates.txt"
FASTA_FILE = "/net/bq-storage/ag-khedkar/Sofia/project_folder/duplicates_with_flanks_200bp.fasta"
FLANK_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/flanks"
MASTER_FLANKS_FASTA = "/net/bq-storage/ag-khedkar/Sofia/project_folder/master_flanks.fasta"

LEFT_FLANK_LEN = 200
RIGHT_FLANK_LEN = 200

# Ensure output directory exists; creates directory tree if missing
os.makedirs(FLANK_DIR, exist_ok=True)

# ---------- LOAD FASTA INTO DICTIONARY ----------
# Map only the "core ID" (without coordinates) to SeqRecord
"""
Load sequences from FASTA_FILE into a lookup dictionary.

Dictionary structure
--------------------
key   : core sequence identifier (without coordinate suffix)
value : Bio.SeqRecord object

Rationale
---------
Duplicate entries contain coordinates appended to IDs.
To match duplicates file entries with sequences, coordinates are removed
so that only the stable sequence identifier is used.
"""
seq_dict = {}
for rec in SeqIO.parse(FASTA_FILE, "fasta"):
    # remove the coordinates at the end after last underscore
    core_id = "_".join(rec.id.split("_")[:-1])
    seq_dict[core_id] = rec

# ---------- HELPER FUNCTION TO EXTRACT FLANK ----------
def extract_flanks(seq_record):
    """
    Extract fixed-length left and right flanking sequences.

    Parameters
    ----------
    seq_record : Bio.SeqRecord
        Sequence record containing full sequence with flanking region already included.

    Returns
    -------
    tuple (str, str)
        left  : first LEFT_FLANK_LEN bases
        right : last RIGHT_FLANK_LEN bases

    Notes
    -----
    Assumes input sequences are long enough to contain both flanks.
    No length validation is performed.
    Flanks are returned as plain strings, not Seq objects.
    """
    seq = str(seq_record.seq)
    left = seq[:LEFT_FLANK_LEN]
    right = seq[-RIGHT_FLANK_LEN:]
    return left, right

# ---------- PROCESS DUPLICATES FILE ----------
"""
Iterate through duplicate sequence pairs and extract flanking regions.

Workflow
--------
1. Read duplicate identifiers from DUP_FILE.
2. Normalize identifier formatting to match FASTA IDs.
3. Look up corresponding sequence in seq_dict.
4. Extract left and right flanks.
5. Save each flank as individual FASTA file.
6. Store flanks in memory for master FASTA output.

Output
------
• Individual flank FASTA files in FLANK_DIR
• Combined master FASTA containing all flanks
"""
master_records = []

with open(DUP_FILE) as f:
    for line in f:
        toks = line.strip().split()
        for seq_id in toks:
            # Normalize identifier:
            # • remove strand orientation prefix (+ or -)
            # • replace coordinate separator ':' with '_'
            seq_id_clean = seq_id.lstrip("+-").replace(":", "_")

            # remove coordinates at the end after last underscore
            core_id = "_".join(seq_id_clean.split("_")[:-1])

            # Verify sequence exists in FASTA dictionary
            if core_id not in seq_dict:
                print(f"Warning: {seq_id} (core ID {core_id}) not found in FASTA")
                continue

            record = seq_dict[core_id]
            left, right = extract_flanks(record)

            # Save left flank
            left_name = f"{core_id}_left"
            left_path = os.path.join(FLANK_DIR, f"{left_name}.fasta")
            with open(left_path, "w") as out:
                out.write(f">{left_name}\n{left}\n")
            master_records.append((left_name, left))

            # Save right flank
            right_name = f"{core_id}_right"
            right_path = os.path.join(FLANK_DIR, f"{right_name}.fasta")
            with open(right_path, "w") as out:
                out.write(f">{right_name}\n{right}\n")
            master_records.append((right_name, right))


# ---------- WRITE MASTER FASTA FOR ALIGNMENT ----------
"""
Write all extracted flanking sequences into a single FASTA file.

Purpose
-------
Provides a unified dataset for downstream alignment or comparative analysis.

Structure
---------
Each record name:
    <core_id>_left
    <core_id>_right
Sequence:
    Corresponding flank sequence
"""
with open(MASTER_FLANKS_FASTA, "w") as out:
    for name, seq in master_records:
        out.write(f">{name}\n{seq}\n")

print(f"Extraction complete. Flanks saved in {FLANK_DIR}")
print(f"Master FASTA for alignment: {MASTER_FLANKS_FASTA}")
