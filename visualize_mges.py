#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt

INPUT_GFF = "/net/bq-storage/ag-khedkar/Sofia/project_folder/master_mges.gff3"
OUTPUT_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder"

os.makedirs(OUTPUT_DIR, exist_ok=True)

# ======================
# STORAGE
# ======================

lengths = []
n_genes = []

nested_lengths = []
non_nested_lengths = []

nested_genes = []
non_nested_genes = []

# pie counters
nested = 0
non_nested = 0
single_once = 0
single_multiple = 0
mixed = 0
is_tn_single_gene = 0
is_tn_2_only = 0
is_tn_1_only = 0
is_tn_complex = 0
other_mge = 0

# ======================
# HELPERS
# ======================

def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            attrs[k] = v
    return attrs

# ======================
# PARSE GFF
# ======================

with open(INPUT_GFF) as f:
    for line in f:
        if line.startswith("#"):
            continue

        cols = line.rstrip().split("\t")
        if len(cols) < 9:
            continue

        length = int(cols[5])
        attrs = parse_attributes(cols[8])
        mge = attrs.get("mge", "")
        mge_type = attrs.get("mge_type", "")
        genes = int(attrs.get("n_genes", 0))

        # ------------------
        # NEW PIE LOGIC
        # ------------------

        if "is_tn" in mge:
            if genes == 1:
                is_tn_single_gene += 1
            elif mge == "is_tn:2":
                is_tn_2_only += 1
            elif "," in mge :
                is_tn_complex += 1
            elif mge == "is_tn:1":
                is_tn_1_only += 1
        else:
            other_mge += 1

        # ------------------
        # EXISTING IS/TN FILTER
        # ------------------

        if "is_tn" not in mge:
            continue

        lengths.append(length)
        n_genes.append(genes)

        if mge_type == "nested":
            nested += 1
            nested_lengths.append(length)
            nested_genes.append(genes)
        elif mge_type == "non-nested":
            non_nested += 1
            non_nested_lengths.append(length)
            non_nested_genes.append(genes)

        # mgeR logic (unchanged)
        if "mgeR" in attrs:
            mgeR = attrs["mgeR"]
            if "," in mgeR:
                mixed += 1
            else:
                try:
                    count = int(mgeR.split(":")[1])
                    if count == 1:
                        single_once += 1
                    else:
                        single_multiple += 1
                except Exception:
                    mixed += 1

# ======================
# HISTOGRAM FUNCTION
# ======================

def log_histogram(data, xlabel, title, outfile, n_bins=40):
    data = np.array(data)
    data = data[data > 0]

    min_exp = np.floor(np.log10(data.min()))
    max_exp = np.ceil(np.log10(data.max()))
    bins = np.logspace(min_exp, max_exp, n_bins)

    plt.figure(figsize=(7, 5))
    plt.hist(data, bins=bins)
    plt.xscale("log")
    plt.xlabel(xlabel)
    plt.ylabel("Count")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

# ======================
# NEW PLOT 5: MGE TYPE PIE
# ======================

labels = [
    f"is_tn & n_genes=1 ({is_tn_single_gene})",
    f"is_tn:2 only ({is_tn_2_only})",
    f"is_tn complex ({is_tn_complex})",
    f"is_tn:1 ({is_tn_1_only})",
    f"Other MGEs ({other_mge})"
]

sizes = [
    is_tn_single_gene,
    is_tn_2_only,
    is_tn_complex,
    is_tn_1_only,
    other_mge
]

plt.figure()
plt.pie(sizes, labels=labels, startangle=90)
plt.title("Distribution of mobile genetic element types")
plt.tight_layout()
plt.savefig(os.path.join(OUTPUT_DIR, "mge_type_overview_pie.png"))
plt.close()

# ======================
# NEW HISTOGRAMS
# ======================

log_histogram(
    nested_lengths,
    "MGE length",
    "Length distribution of nested IS/Tn MGEs",
    os.path.join(OUTPUT_DIR, "is_tn_nested_length_histogram.png")
)

log_histogram(
    non_nested_lengths,
    "MGE length",
    "Length distribution of non-nested IS/Tn MGEs",
    os.path.join(OUTPUT_DIR, "is_tn_non_nested_length_histogram.png")
)

log_histogram(
    nested_genes,
    "Number of genes",
    "Gene count distribution of nested IS/Tn MGEs",
    os.path.join(OUTPUT_DIR, "is_tn_nested_n_genes_histogram.png")
)

log_histogram(
    non_nested_genes,
    "Number of genes",
    "Gene count distribution of non-nested IS/Tn MGEs",
    os.path.join(OUTPUT_DIR, "is_tn_non_nested_n_genes_histogram.png")
)

print("All plots (original + 5 new ones) saved to:", OUTPUT_DIR)
