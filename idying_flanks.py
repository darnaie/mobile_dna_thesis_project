import os

GENOME_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/genome_seqs"
DUP_FILE = "/net/bq-storage/ag-khedkar/Sofia/project_folder/mge_seqs_deduplicated/deduplicated.fasta.duplicates.txt"
OUT_FASTA = "duplicates_with_flanks_200bp.fasta"

FLANK = 200

# ---------- helpers ----------

def revcomp(seq):
    """
    Return the reverse complement of a DNA sequence.

    Parameters
    ----------
    seq : str
        DNA sequence containing A, C, G, T, N (case-insensitive).

    Returns
    -------
    str
        Reverse-complemented sequence with original case preserved.

    Notes
    -----
    Uses translation table for nucleotide complement and reverses string.
    """
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def parse_entry(entry):
    """
    Example entry:
    MGE_GCA_900083375.1_1280.SAMEA2445630.FKXQ01000003:261901-264243

    After removing first 3 underscore blocks:
    1280.SAMEA2445630.FKXQ01000003:261901-264243

    Extract contig identifier and genomic coordinates from entry string.

    Parameters
    ----------
    entry : str
        Duplicate entry identifier, possibly prefixed with '+' or '-'.

    Returns
    -------
    tuple (str, int, int)
        contig : full contig ID matching FASTA header
        start  : start coordinate (int)
        end    : end coordinate (int)

    Notes
    -----
    Removes orientation symbol if present.
    Assumes format where locus and coordinates are separated by ':'.
    """
    entry = entry.lstrip("+-")

    core = entry.split("_", 3)[3]
    locus, coords = core.split(":")

    contig = locus  # full contig ID exactly as in FASTA
    start, end = map(int, coords.split("-"))

    return contig, start, end


def genome_filename(entry):
    """
    From:
    1280.SAMEA2445630.FKXQ01000003
    -> 1280.SAMEA2445630.fasta

    Generate genome FASTA filename from entry identifier.

    Parameters
    ----------
    entry : str
        Duplicate entry identifier with optional orientation prefix.

    Returns
    -------
    str
        Genome FASTA filename corresponding to entry.

    Notes
    -----
    Extracts first two dot-separated components of locus ID.
    Used to locate genome file in GENOME_DIR.
    """
    entry = entry.lstrip("+-")
    core = entry.split("_", 3)[3]
    locus = core.split(":")[0]

    parts = locus.split(".")
    prefix = ".".join(parts[:2])

    return f"{prefix}.fasta"


def load_genome(path):
    """
    Load a multi-FASTA genome file into memory.

    Parameters
    ----------
    path : str
        Path to FASTA file.

    Returns
    -------
    dict
        Dictionary mapping contig name -> full sequence string.

    Notes
    -----
    FASTA headers are truncated at first whitespace.
    Entire file is stored in memory.
    """
    seqs = {}
    name = None
    seq = []

    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(seq)
                name = line[1:].split()[0]
                seq = []
            else:
                seq.append(line.strip())

        if name:
            seqs[name] = "".join(seq)

    return seqs
    

def fetch_region_from_fasta(fasta_path, target_contig, start, end, flank):
    """
    Extract a genomic region with flanking bases directly from FASTA file.

    Parameters
    ----------
    fasta_path : str
        Path to genome FASTA file.
    target_contig : str
        Contig identifier to extract from.
    start : int
        Region start coordinate.
    end : int
        Region end coordinate.
    flank : int
        Number of bases to include upstream and downstream.

    Returns
    -------
    str
        Extracted sequence region including flanks.
        Returns empty string if contig not found or no overlap.

    Notes
    -----
    Reads FASTA sequentially without loading full genome into memory.
    Coordinates are treated as 0-based positions.
    Stops reading once region end is passed.
    """
    start = max(0, start - flank)
    end = end + flank

    seq_chunks = []
    current = None
    pos = 0

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                current = line[1:].split()[0]
                pos = 0
                continue

            if current != target_contig:
                continue

            line = line.strip()
            line_len = len(line)

            line_start = pos
            line_end = pos + line_len

            if line_end >= start and line_start <= end:
                s = max(0, start - line_start)
                e = min(line_len, end - line_start)
                seq_chunks.append(line[s:e])

            pos += line_len

            if pos > end:
                break

    return "".join(seq_chunks)


# ---------- read duplicate pairs ----------

pairs = []
with open(DUP_FILE) as f:
    for line in f:
        toks = line.strip().split()
        if len(toks) != 2:
            continue
        ref, dup = toks
        orientation = "+" if dup.startswith("+") else "-"
        pairs.append((ref, dup, orientation))

print("Pairs:", len(pairs))


# ---------- genome cache ----------

genome_cache = {}

def get_genome(filename):
    """
    Retrieve genome from cache or load it if not already loaded.

    Parameters
    ----------
    filename : str
        Genome FASTA filename.

    Returns
    -------
    dict or None
        Dictionary of contig -> sequence if genome exists.
        None if genome file not found.

    Notes
    -----
    Uses global genome_cache to avoid repeated disk reads.
    """
    if filename not in genome_cache:
        path = os.path.join(GENOME_DIR, filename)
        if not os.path.exists(path):
            return None
        genome_cache[filename] = load_genome(path)
    return genome_cache[filename]


# ---------- extraction ----------

written = 0
missing_genome = 0
missing_contig = 0

with open(OUT_FASTA, "w") as out:

    for ref, dup, orientation in pairs:

        contig1, s1, e1 = parse_entry(ref)
        contig2, s2, e2 = parse_entry(dup)

        genome1 = os.path.join(GENOME_DIR, genome_filename(ref))
        genome2 = os.path.join(GENOME_DIR, genome_filename(dup))

        if not os.path.exists(genome1) or not os.path.exists(genome2):
            missing_genome += 1
            continue

        region1 = fetch_region_from_fasta(genome1, contig1, s1, e1, FLANK)
        region2 = fetch_region_from_fasta(genome2, contig2, s2, e2, FLANK)

        if not region1 or not region2:
            missing_contig += 1
            continue

        if orientation == "-":
            region2 = revcomp(region2)

        id1 = ref.replace(":", "_")
        id2 = dup.lstrip("+-").replace(":", "_")

        out.write(f">{id1}\n{region1}\n")
        out.write(f">{id2}\n{region2}\n")

        written += 2

print("Sequences written:", written)
print("Pairs missing genome:", missing_genome)
print("Pairs missing contig:", missing_contig)
print("Done.")
