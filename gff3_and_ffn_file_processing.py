#!/usr/bin/env python3

import os
import gzip
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed

BASE_DIR = "/net/bq-storage/ag-khedkar/khedkar/pg3_data/pg3_genomes/"
OUTPUT_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder"

GFF_OUT = os.path.join(OUTPUT_DIR, "master_mges.gff3")
FFN_NON_NESTED_OUT = os.path.join(OUTPUT_DIR, "master_mges_non-nested.ffn.gz")
FFN_NESTED_OUT = os.path.join(OUTPUT_DIR, "master_mges_nested.ffn.gz")
LOG_FILE = os.path.join(OUTPUT_DIR, "master_mges.log")

os.makedirs(OUTPUT_DIR, exist_ok=True)

logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s"
)

def log(msg):
    print(msg, flush=True)
    logging.info(msg)

# -------------------------
# FASTA PARSER
# -------------------------

def parse_fasta_gz(path):
    header = None
    seq = []
    with gzip.open(path, "rt") as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq)
                header = line
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, "".join(seq)

# -------------------------
# WORKER FUNCTION
# -------------------------

def process_file(path):
    gff_lines = []
    non_nested = []
    nested = []

    if path.endswith(".gff3"):
        with open(path) as f:
            for line in f:
                if "\tmobile_genetic_element\t" in line:
                    gff_lines.append(line)

    elif path.endswith(".ffn.gz"):
        for header, seq in parse_fasta_gz(path):
            h = header.lower()
            if "mge=is_tn" not in h:
                continue
            record = f"{header}\n{seq}\n"
            if "mge_type=non-nested" in h:
                non_nested.append(record)
            elif "mge_type=nested" in h:
                nested.append(record)

    return gff_lines, non_nested, nested

# -------------------------
# MAIN
# -------------------------

def main():
    log("Starting parallel MGE extraction")

    files = []
    for root, _, fs in os.walk(BASE_DIR):
        for f in fs:
            if f.endswith(".gff3") or f.endswith(".ffn.gz"):
                files.append(os.path.join(root, f))

    log(f"Discovered {len(files)} files")

    with open(GFF_OUT, "w") as gff_out, \
         gzip.open(FFN_NON_NESTED_OUT, "wt") as ffn_non, \
         gzip.open(FFN_NESTED_OUT, "wt") as ffn_nest:

        with ProcessPoolExecutor() as exe:
            futures = [exe.submit(process_file, f) for f in files]

            for i, fut in enumerate(as_completed(futures), 1):
                gff, non, nest = fut.result()

                for l in gff:
                    gff_out.write(l)
                for r in non:
                    ffn_non.write(r)
                for r in nest:
                    ffn_nest.write(r)

                if i % 50 == 0:
                    log(f"Processed {i}/{len(files)} files")

    log("Finished successfully")

if __name__ == "__main__":
    main()
