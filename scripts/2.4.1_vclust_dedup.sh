#!/usr/bin/env bash
# =============================================================================
# 2.4.1_vclust_dedup.sh  —  Section 2.4.1
# =============================================================================
# Runs `vclust deduplicate` on the whole non-nested is_tn transposon
# sequences (not just their recombinase genes -- that's what makes this
# "duplicate GROUPING" of whole elements, distinct from Section 2.3's
# recombinase-level clustering).
#
# Converted from a Python wrapper to pure bash -- see 2.2.2_skani_dedup.sh's
# header comment for why (same reasoning applies here).
#
# NAMING FIX: the duplicates file vclust produces is named
# "<prefix>.duplicates.txt", NOT "<prefix>.fasta.duplicates.txt" (an
# earlier version of this pipeline assumed the latter, based on a
# real-world example file whose "-o" prefix itself happened to already
# end in ".fasta" -- that ".fasta" was part of the prefix chosen at the
# time, not something vclust appends automatically). OUTPUT_PREFIX below
# deliberately does NOT include ".fasta", so the real duplicates file is
# OUTPUT_PREFIX + ".duplicates.txt". The summary section below also
# double-checks by globbing for *.duplicates.txt, so it reports the real
# file regardless of exactly how vclust names things on your version.
#
# THIS IS AN ACTUAL SLURM JOB -- submit it with sbatch:
#     sbatch 2.4.1_vclust_dedup.sh
#
# LOGGING: --output/--error point directly at BASE_DIR/logs/ -- see
# 2.2.2_skani_dedup.sh's header comment for why.
#
# BASE_DIR / VCLUST_BIN below should match config.yaml's paths.base_dir /
# paths.vclust_path -- kept as plain bash variables here, not read from
# the YAML; keep them in sync by hand if either changes.
# =============================================================================

#SBATCH --job-name=2_4_1_vclust_dedup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=22:30:00
#SBATCH --mem=32G
#SBATCH --output=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.4.1_vclust_dedup_%j.out
#SBATCH --error=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.4.1_vclust_dedup_%j.err
#SBATCH --mail-type=END,FAIL
# Uncomment and set if you need a specific node:
##SBATCH --nodelist=cln132

set -euo pipefail

# ── Config (keep in sync with config.yaml) ─────────────────────────────────
BASE_DIR="/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep"
VCLUST_BIN="/home/bq_slipaeva/.local/bin/vclust"

INPUT_FASTA="${BASE_DIR}/files/non_nested_is_tn.fasta"
OUTPUT_DIR="${BASE_DIR}/vclust_deduplication"
OUTPUT_PREFIX="${OUTPUT_DIR}/non_nested_is_tn_deduplicated"   # NOT ".fasta" -- see header note

# PLACEHOLDER: set to 1 to symlink an already-existing vclust_deduplication/
# instead of actually running vclust; set PLACEHOLDER_PATH accordingly.
PLACEHOLDER=0
PLACEHOLDER_PATH="BASE/vclust_deduplication/"

# Module loads / environment activation, matching the original ad-hoc script:
module purge
module load compiler/gcc/10.2.0
module load devel/python/3.10.0
source ~/uv_env/bin/activate

echo "======================================================================"
echo "Section 2.4.1 — vclust deduplication of non-nested is_tn transposons"
echo "======================================================================"
echo "Started : $(date)"
echo "Host    : $(hostname)"
echo "Job ID  : ${SLURM_JOB_ID:-(not running under slurm)}"

if [[ "${PLACEHOLDER}" == "1" ]]; then
    if [[ -z "${PLACEHOLDER_PATH}" ]]; then
        echo "ERROR: PLACEHOLDER=1 but PLACEHOLDER_PATH is empty -- set it, or set PLACEHOLDER=0." >&2
        exit 1
    fi
    if [[ ! -e "${PLACEHOLDER_PATH}" ]]; then
        echo "ERROR: PLACEHOLDER_PATH does not exist: ${PLACEHOLDER_PATH}" >&2
        exit 1
    fi
    mkdir -p "$(dirname "${OUTPUT_DIR}")"
    if [[ -e "${OUTPUT_DIR}" || -L "${OUTPUT_DIR}" ]]; then
        if [[ "$(readlink -f "${OUTPUT_DIR}")" == "$(readlink -f "${PLACEHOLDER_PATH}")" ]]; then
            echo "PLACEHOLDER=1 -- ${OUTPUT_DIR} already IS ${PLACEHOLDER_PATH}; nothing to stage."
            exit 0
        fi
        rm -rf "${OUTPUT_DIR}"
    fi
    ln -s "${PLACEHOLDER_PATH}" "${OUTPUT_DIR}"
    echo "PLACEHOLDER=1 -- symlinked ${OUTPUT_DIR} -> ${PLACEHOLDER_PATH} (real computation skipped)."
    exit 0
fi

if [[ ! -x "${VCLUST_BIN}" ]]; then
    echo "ERROR: vclust binary not found or not executable at ${VCLUST_BIN}" >&2
    exit 1
fi
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: ${INPUT_FASTA} not found -- run 2.1.1_extract_mges.py first." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

n_seqs=$(grep -c '^>' "${INPUT_FASTA}")
echo "Input FASTA   : ${INPUT_FASTA} (${n_seqs} sequences)"
echo "Output prefix : ${OUTPUT_PREFIX}"
echo "----------------------------------------------------------------------"

"${VCLUST_BIN}" deduplicate -i "${INPUT_FASTA}" -o "${OUTPUT_PREFIX}"

# Don't assume an exact filename -- glob for it, same as the downstream
# scripts (append_duplicate_info.py, plot_flowchart_section_2_4.py) do.
duplicates_file=$(find "${OUTPUT_DIR}" -maxdepth 1 -name "*.duplicates.txt" | head -n1 || true)
dedup_fasta="${OUTPUT_PREFIX}.fasta"

n_dedup=0
n_groups=0
[[ -f "${dedup_fasta}" ]] && n_dedup=$(grep -c '^>' "${dedup_fasta}" || true)
if [[ -n "${duplicates_file}" && -f "${duplicates_file}" ]]; then
    n_groups=$(wc -l < "${duplicates_file}")
fi

echo "----------------------------------------------------------------------"
echo "SUMMARY"
echo "----------------------------------------------------------------------"
echo "Input sequences        : ${n_seqs}"
echo "Deduplicated sequences  : ${n_dedup}"
echo "Duplicate groups        : ${n_groups}"
echo "-> ${dedup_fasta}"
if [[ -n "${duplicates_file}" ]]; then
    echo "-> ${duplicates_file}"
else
    echo "WARNING: no *.duplicates.txt found under ${OUTPUT_DIR} -- check vclust's actual output naming."
fi
echo "Finished : $(date)"
echo "Section 2.4.1 complete."