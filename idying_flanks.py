#!/usr/bin/env python3
import os

GENOME_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/genome_seqs"
DUP_FILE = "/net/bq-storage/ag-khedkar/Sofia/project_folder/mge_seqs_deduplicated/deduplicated.fasta.duplicates.txt"
FLANK_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/flanks"
MGES_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/mges"

LEFT_FLANK_LEN = 200
RIGHT_FLANK_LEN = 200
EXTRACT_FLANK = 200
MIN_GROUP_SIZE = 100

os.makedirs(FLANK_DIR, exist_ok=True)
os.makedirs(MGES_DIR, exist_ok=True)

# ---------- helpers ----------

def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def parse_entry(entry):
    e = entry.lstrip("+-")
    parts = e.split("_", 3)
    if len(parts) < 4:
        raise ValueError(f"Unexpected entry format: {entry}")
    core = parts[3]
    locus, coords = core.split(":", 1)
    start_s, end_s = coords.split("-", 1)
    return locus, int(start_s), int(end_s)


def genome_filename(entry):
    e = entry.lstrip("+-")
    parts = e.split("_", 3)
    if len(parts) < 4:
        return None
    locus = parts[3].split(":", 1)[0]
    p = locus.split(".")
    prefix = ".".join(p[:2]) if len(p) >= 2 else locus
    return f"{prefix}.fasta"


def load_genome(path):
    seqs = {}
    name = None
    buf = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name:
            seqs[name] = "".join(buf)
    return seqs


def get_genome(filename):
    if filename is None:
        return None
    path = os.path.join(GENOME_DIR, filename)
    if not os.path.exists(path):
        return None
    return load_genome(path)


def fetch_region_from_genome(genome_dict, contig, start, end, flank):
    if genome_dict is None:
        return ""
    seq = genome_dict.get(contig)
    if seq is None:
        return ""
    left1 = max(1, start - flank)
    right1 = end + flank
    s0 = left1 - 1
    e0 = min(len(seq), right1)
    return seq[s0:e0]


# ---------- master FASTA streams ----------

master_left_path = os.path.join(FLANK_DIR, "master_left.fasta")
master_right_path = os.path.join(FLANK_DIR, "master_right.fasta")
master_mges_path = os.path.join(MGES_DIR, "master_mges.fasta")

ml = open(master_left_path, "w")
mr = open(master_right_path, "w")
mm = open(master_mges_path, "w")

# ---------- processing ----------

filtered_dup_path = os.path.join(FLANK_DIR, f"duplicates_min{MIN_GROUP_SIZE}.txt")

group_i = 0
total_written = 0
missing_genomes = 0
missing_contigs = 0
bad_entries = 0

with open(DUP_FILE) as inf, open(filtered_dup_path, "w") as outf:
    for line in inf:
        toks = line.strip().split()
        if len(toks) < MIN_GROUP_SIZE:
            continue

        outf.write(line)

        group_i += 1
        group_name = f"dup_group_{group_i:05d}"

        left_path = os.path.join(FLANK_DIR, f"{group_name}_left.fasta")
        right_path = os.path.join(FLANK_DIR, f"{group_name}_right.fasta")
        mges_path = os.path.join(MGES_DIR, f"{group_name}_mges.fasta")

        with open(left_path, "w") as left_fh, \
             open(right_path, "w") as right_fh, \
             open(mges_path, "w") as mge_fh:

            for entry in toks:
                try:
                    contig, s, e = parse_entry(entry)
                except Exception:
                    bad_entries += 1
                    continue

                genome = get_genome(genome_filename(entry))
                if genome is None:
                    missing_genomes += 1
                    continue

                region = fetch_region_from_genome(genome, contig, s, e, EXTRACT_FLANK)
                if not region:
                    missing_contigs += 1
                    continue

                seq_full = genome.get(contig, "")
                if not seq_full:
                    missing_contigs += 1
                    continue

                mge_seq = seq_full[s-1:e]

                if entry.startswith("-"):
                    region = revcomp(region)
                    mge_seq = revcomp(mge_seq)

                idbase = entry.lstrip("+-").replace(":", "_")
                left_seq = region[:LEFT_FLANK_LEN]
                right_seq = region[-RIGHT_FLANK_LEN:] if len(region) >= 1 else ""

                left_name = f"{idbase}_left"
                right_name = f"{idbase}_right"
                mge_name = f"{idbase}_mge"

                # write group files
                left_fh.write(f">{left_name}\n{left_seq}\n")
                right_fh.write(f">{right_name}\n{right_seq}\n")
                mge_fh.write(f">{mge_name}\n{mge_seq}\n")

                # write master files immediately
                ml.write(f">{left_name}\n{left_seq}\n")
                mr.write(f">{right_name}\n{right_seq}\n")
                mm.write(f">{mge_name}\n{mge_seq}\n")

                total_written += 3

ml.close()
mr.close()
mm.close()

# ---------- summary ----------

print(f"Filtered duplicates written: {filtered_dup_path}")
print(f"Groups processed: {group_i}")
print(f"Sequences written (left+right+mge): {total_written}")
print(f"Bad entries: {bad_entries}")
print(f"Missing genome files: {missing_genomes}")
print(f"Missing contigs: {missing_contigs}")
print("Flanks saved in", FLANK_DIR)
print("MGEs saved in", MGES_DIR)
print("Master left:", master_left_path)
print("Master right:", master_right_path)
print("Master mges:", master_mges_path)