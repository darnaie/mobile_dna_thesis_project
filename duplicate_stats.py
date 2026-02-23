import matplotlib.pyplot as plt
from collections import Counter
import os
import numpy as np

DUP_FILE = "/net/bq-storage/ag-khedkar/Sofia/project_folder/mge_seqs_deduplicated/deduplicated.fasta.duplicates.txt"
OUT_DIR = "/net/bq-storage/ag-khedkar/Sofia/project_folder/"

# ----------------- helpers -----------------

def read_duplicates(file_path):
    """Read duplicates file and return list of lines (retained + duplicates)"""
    with open(file_path) as f:
        lines = [line.strip() for line in f if line.strip()]
    return lines

def count_duplicates_per_mge(lines):
    """Return list of duplicate counts per MGE"""
    counts = [len(line.split()) - 1 for line in lines]  # minus retained sequence
    return counts

def species_jump_counts(lines):
    """Return tuple: (same_species_count, cross_species_count) using GCA accession to detect species jump"""
    same_species = 0
    cross_species = 0

    for line in lines:
        toks = line.split()
        # Extract GCA code from retained sequence
        retained_gca = toks[0].split("_")[2] 
        duplicates = toks[1:]
        jump = False
        for d in duplicates:
            dup_gca = d.lstrip("+-").split("_")[2] 
            if dup_gca != retained_gca:
                jump = True
                break
        if jump:
            cross_species += 1
        else:
            same_species += 1

    return same_species, cross_species


def orientation_counts(lines):
    """Return tuple: (same_direction, opposing, mixed)"""
    same_dir = 0
    opp_dir = 0
    mixed_dir = 0

    for line in lines:
        toks = line.split()
        duplicates = toks[1:]
        if not duplicates:
            continue
        orientations = [d[0] for d in duplicates]  # '+' or '-'
        if all(o=='+' for o in orientations):
            same_dir += 1
        elif all(o=='-' for o in orientations):
            opp_dir += 1
        else:
            mixed_dir += 1
    return same_dir, opp_dir, mixed_dir

# ----------------- plotting functions -----------------

def plot_duplicates_hist(counts, out_dir):
    """Bar chart: number of duplicates per MGE with log-scaled bins and nice color"""
    plt.figure(figsize=(10,6))
    
    # Log-spaced bins from 1 to max(counts)
    bins = np.logspace(np.log10(1), np.log10(max(counts)+1), num=50)
    
    plt.hist(counts, bins=bins, color='cyan', edgecolor='black')
    plt.xscale('log')  # x-axis logarithmic
    plt.yscale('log')  # y-axis logarithmic, optional for skewed distribution
    plt.xlabel("Number of duplicates per MGE (log scale)")
    plt.ylabel("Count of MGEs (log scale)")
    plt.title("Distribution of duplicates per MGE")
    
    out_path = os.path.join(out_dir, "duplicates_per_mge_hist.png")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved histogram: {out_path}")

def plot_species_pie(same_species, cross_species, out_dir):
    """Pie chart: same species vs cross-species duplicates"""
    total = same_species + cross_species
    labels = [
        f"Same species\n{same_species} ({same_species/total*100:.1f}%)",
        f"Jump between species\n{cross_species} ({cross_species/total*100:.1f}%)"
    ]
    plt.figure(figsize=(7,6))
    plt.pie([same_species, cross_species], labels=labels, autopct='', colors=["#66c2a5","#fc8d62"])
    plt.title(f"Species-level distribution of duplicates (Total MGEs: {total})")
    out_path = os.path.join(out_dir, "species_jump_pie_GCA.png")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved species pie chart: {out_path}")

def plot_orientation_pie(same_dir, opp_dir, mixed_dir, out_dir):
    """Pie chart: orientation distribution of duplicates"""
    total = same_dir + opp_dir + mixed_dir
    labels = [
        f"Same direction\n{same_dir} ({same_dir/total*100:.1f}%)",
        f"Opposing direction\n{opp_dir} ({opp_dir/total*100:.1f}%)",
        f"Mixed direction\n{mixed_dir} ({mixed_dir/total*100:.1f}%)"
    ]
    plt.figure(figsize=(7,6))
    plt.pie([same_dir, opp_dir, mixed_dir], labels=labels, autopct='', colors=["#8da0cb","#e78ac3","#a6d854"])
    plt.title(f"Orientation distribution of duplicates (Total MGEs: {total})")
    out_path = os.path.join(out_dir, "orientation_pie.png")
    plt.savefig(out_path)
    plt.close()
    print(f"Saved orientation pie chart: {out_path}")

# ----------------- main -----------------

def main():
    lines = read_duplicates(DUP_FILE)
    print("Total lines in duplicates file:", len(lines))

    # duplicates per MGE histogram
    counts = count_duplicates_per_mge(lines)

    # species jump pie
    same_species, cross_species = species_jump_counts(lines)

    # orientation pie
    same_dir, opp_dir, mixed_dir = orientation_counts(lines)

    # ----------------- plotting -----------------
    # comment out any of the following lines if you don't need a specific graph
    plot_duplicates_hist(counts, OUT_DIR)
    plot_species_pie(same_species, cross_species, OUT_DIR)
    plot_orientation_pie(same_dir, opp_dir, mixed_dir, OUT_DIR)

if __name__ == "__main__":
    main()
