#!/usr/bin/env bash
# =============================================================================
# 2.3.1_cluster_recombinases.sh  —  Section 2.3.1
# =============================================================================
# Runs `mmseqs easy-cluster` on the representative-genome-restricted
# recombinase sequences (Section 2.2.4's output). ANI (--min-seq-id) and
# coverage are both user-chosen below -- no fixed default.
#
# Converted from a Python wrapper to pure bash -- see 2.2.2_skani_dedup.sh's
# header comment for why (same reasoning applies here).
#
# THIS IS AN ACTUAL SLURM JOB -- submit it with sbatch:
#     sbatch 2.3.1_cluster_recombinases.sh
#
# LOGGING: --output/--error point directly at BASE_DIR/logs/ -- see
# 2.2.2_skani_dedup.sh's header comment for why this is done via SLURM's
# own redirect rather than an in-script `tee`.
#
# BASE_DIR / MMSEQS_BIN below should match config.yaml's paths.base_dir /
# paths.mmseqs2_path -- kept as plain bash variables here, not read from
# the YAML; keep them in sync by hand if either changes.
# =============================================================================

#SBATCH --job-name=2_3_1_cluster_recombinases
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --output=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.3.1_cluster_recombinases_%j.out
#SBATCH --error=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.3.1_cluster_recombinases_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ── Config (keep in sync with config.yaml) ─────────────────────────────────
BASE_DIR="/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep"
MMSEQS_BIN="/net/bq-storage/ag-khedkar/Sofia/project_folder/mmseqs/bin/mmseqs"

# Section 2.2.4's representative-genome-only recombinase fasta (NOT
# 2.1.5's is_tn_recombinase_non_nested.fasta, which covers ALL is_tn
# recombinases -- see 2.2.4's docstring for why these are two separate files).
INPUT_FASTA="${BASE_DIR}/files/is_tn_recombinase_non_nested_repr.fasta"

ANI_THRESHOLD="0.85"     # --min-seq-id
COVERAGE="0.85"          # -c
COV_MODE="0"             # bidirectional coverage
EXACT_KMER_MATCHING="0"  # only set to "1" if ANI_THRESHOLD is exactly "1.0"

OUTPUT_DIR="${BASE_DIR}/files/recombinase_cluster_${ANI_THRESHOLD}"
OUTPUT_PREFIX="${OUTPUT_DIR}/recombinase_cluster"
TMP_DIR="${OUTPUT_DIR}/tmp_mmseqs"
KEEP_TMP=0   # set to 1 to keep tmp_mmseqs/ after a successful run (debugging only)

# PLACEHOLDER: set to 1 to symlink an already-existing cluster output dir
# instead of actually running mmseqs; set PLACEHOLDER_PATH accordingly.
PLACEHOLDER=0
PLACEHOLDER_PATH="/net/bq-storage/ag-khedkar/Sofia/project_folder/recombinase_cluster_relevant_85"

# (Optional) module loads / environment activation for your cluster:
# module purge
# module load devel/python/3.10.0

echo "======================================================================"
echo "Section 2.3.1 — Recombinase clustering (mmseqs2 easy-cluster)"
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

if [[ ! -x "${MMSEQS_BIN}" ]]; then
    echo "ERROR: mmseqs binary not found or not executable at ${MMSEQS_BIN}" >&2
    exit 1
fi
if [[ ! -f "${INPUT_FASTA}" ]]; then
    echo "ERROR: ${INPUT_FASTA} not found -- run 2.2.4_build_recombinase_table.py first." >&2
    exit 1
fi

if [[ "${EXACT_KMER_MATCHING}" == "1" && "${ANI_THRESHOLD}" != "1.0" ]]; then
    echo "WARNING: EXACT_KMER_MATCHING=1 with ANI_THRESHOLD=${ANI_THRESHOLD} (not 1.0) --" \
         "this is only appropriate at min-seq-id=1.0 and may cause real near-duplicates" \
         "to be missed (under-clustering). Consider setting it back to 0."
fi

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"

n_seqs=$(grep -c '^>' "${INPUT_FASTA}")
echo "Input FASTA        : ${INPUT_FASTA} (${n_seqs} sequences)"
echo "Output prefix      : ${OUTPUT_PREFIX}"
echo "ANI (--min-seq-id) : ${ANI_THRESHOLD}"
echo "Coverage (-c)      : ${COVERAGE} (mode ${COV_MODE})"
echo "exact-kmer-matching: ${EXACT_KMER_MATCHING}"
echo "----------------------------------------------------------------------"

"${MMSEQS_BIN}" easy-cluster \
    "${INPUT_FASTA}" \
    "${OUTPUT_PREFIX}" \
    "${TMP_DIR}" \
    --min-seq-id "${ANI_THRESHOLD}" \
    -c "${COVERAGE}" \
    --cov-mode "${COV_MODE}" \
    --exact-kmer-matching "${EXACT_KMER_MATCHING}" \
    --threads "${SLURM_CPUS_PER_TASK:-20}" \
    -v 3

if [[ "${KEEP_TMP}" == "1" ]]; then
    tmp_size=$(du -sh "${TMP_DIR}" 2>/dev/null | cut -f1)
    echo "KEEP_TMP=1 -- leaving ${TMP_DIR} in place (${tmp_size})."
else
    tmp_size=$(du -sh "${TMP_DIR}" 2>/dev/null | cut -f1)
    rm -rf "${TMP_DIR}"
    echo "Removed ${TMP_DIR} (${tmp_size}) -- set KEEP_TMP=1 above to keep it for debugging."
fi

rep_fasta="${OUTPUT_PREFIX}_rep_seq.fasta"
cluster_tsv="${OUTPUT_PREFIX}_cluster.tsv"
n_clusters=0
n_members=0
[[ -f "${rep_fasta}" ]] && n_clusters=$(grep -c '^>' "${rep_fasta}")
[[ -f "${cluster_tsv}" ]] && n_members=$(wc -l < "${cluster_tsv}")

echo "----------------------------------------------------------------------"
echo "SUMMARY"
echo "----------------------------------------------------------------------"
echo "Sequences clustered       : ${n_seqs}"
echo "Clusters (representatives): ${n_clusters}"
echo "Cluster member rows       : ${n_members}"
echo "-> ${rep_fasta}"
echo "-> ${cluster_tsv}"
echo "-> ${OUTPUT_PREFIX}_all_seqs.fasta"
echo "Finished : $(date)"
echo "Section 2.3.1 complete."