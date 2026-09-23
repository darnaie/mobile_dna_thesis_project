#!/usr/bin/env bash
# =============================================================================
# 2.6.1_align_flanks.sh  —  Section 2.6.1
# =============================================================================
# REWRITE (see chat): the previous version treated every subdirectory under
# a cluster_<N>/ folder -- "duplicate_group_<M>" AND "singleton" alike -- as
# something to align, and Level 2 pooled files from ANY subdirectory two
# levels down. That's wrong for what 2.5.1 actually produces:
#
#   BASE/flanks_split_acc_to_duplicates/<rec_type>/cluster_<N>/
#       duplicate_group_<M>/<mge_id>_left.fasta   <mge_id>_right.fasta   (>=2 members)
#       singleton/<mge_id>_left.fasta             <mge_id>_right.fasta   (every singleton
#                                                                          MGE in this
#                                                                          cluster, pooled
#                                                                          into ONE shared
#                                                                          "singleton" folder)
#
# Each .fasta is one transposon's flank -- never a pre-made stacking. A
# "singleton" folder holds MGEs that were never observed as an exact
# duplicate of anything, so there is no second copy to compare against and
# therefore no conserved-boundary signal to align in the first place --
# singletons are excluded from BOTH levels below, entirely:
#
#   LEVEL 1 -- duplicate_group_level: ONLY "duplicate_group_<M>" subfolders
#   are processed (an explicit `-type d -name "duplicate_group_*"` filter,
#   not "every subfolder") -- every *_left.fasta inside ONE such folder is
#   concatenated into a raw stack, then aligned:
#     BASE/consensus/duplicate_group_level/<rec_type>/cluster_<N>/duplicate_group_<M>/
#       left.raw.fasta   left.aln.fasta
#       right.raw.fasta  right.aln.fasta
#   "singleton" folders are never touched at this level.
#
#   LEVEL 2 -- cluster_level: for a given cluster, first find its
#   "duplicate_group_*" subfolders (if any); pool every *_left.fasta from
#   ONLY those subfolders (never from "singleton") into one stack, then
#   align. If a cluster has NO duplicate_group_* subfolders at all (i.e. it
#   contains only a "singleton" bucket, or nothing) that whole cluster is
#   SKIPPED -- no output directory is created for it at this level.
#     BASE/consensus/cluster_level/<rec_type>/
#       <rec_type>_cluster_<N>_left.raw.fasta   _left.aln.fasta
#       <rec_type>_cluster_<N>_right.raw.fasta  _right.aln.fasta
#
# THIS IS AN ACTUAL SLURM JOB -- submit it with sbatch:
#     sbatch 2.6.1_align_flanks.sh
#
# LOGGING: --output/--error point directly at BASE/logs/ -- see
# 2.2.2_skani_dedup.sh's header comment for why (SLURM's own redirect,
# not an in-script tee, to avoid the hang risk documented there).
#
# BASE_DIR / MAFFT_BIN below should match config.yaml's paths.base_dir --
# kept as plain bash variables here, not read from the YAML.
# =============================================================================

#SBATCH --job-name=2_6_1_align_flanks
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=60
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --partition=single
#SBATCH --output=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.6.1_align_flanks_%j.out
#SBATCH --error=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.6.1_align_flanks_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail
shopt -s nullglob   # unmatched globs (e.g. no cluster_* dirs at all) expand to nothing, not a literal "*"

# ── Config (keep in sync with config.yaml) ─────────────────────────────────
BASE_DIR="/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep"
MAFFT_BIN="mafft"   # override to a full path if it's not on PATH in your sbatch environment

INPUT_BASE="${BASE_DIR}/flanks_split_acc_to_duplicates"
OUT_BASE="${BASE_DIR}/consensus"
DUP_OUT="${OUT_BASE}/duplicate_group_level"
CLU_OUT="${OUT_BASE}/cluster_level"

THREADS="${SLURM_CPUS_PER_TASK:-16}"

# PLACEHOLDER: set to 1 to symlink an already-existing consensus/ dir
# instead of actually running mafft; set PLACEHOLDER_PATH accordingly.
PLACEHOLDER=0
PLACEHOLDER_PATH=""

# (Optional) module loads / environment activation for your cluster:
# module purge
# module load bio/mafft/7.490

echo "======================================================================"
echo "Section 2.6.1 — Flank alignment (mafft)"
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
    mkdir -p "$(dirname "${OUT_BASE}")"
    if [[ -e "${OUT_BASE}" || -L "${OUT_BASE}" ]]; then
        if [[ "$(readlink -f "${OUT_BASE}")" == "$(readlink -f "${PLACEHOLDER_PATH}")" ]]; then
            echo "PLACEHOLDER=1 -- ${OUT_BASE} already IS ${PLACEHOLDER_PATH}; nothing to stage."
            exit 0
        fi
        rm -rf "${OUT_BASE}"
    fi
    ln -s "${PLACEHOLDER_PATH}" "${OUT_BASE}"
    echo "PLACEHOLDER=1 -- symlinked ${OUT_BASE} -> ${PLACEHOLDER_PATH} (real computation skipped)."
    exit 0
fi

if ! command -v "${MAFFT_BIN}" &>/dev/null; then
    echo "ERROR: mafft ('${MAFFT_BIN}') not found in PATH. Load the relevant module or set MAFFT_BIN." >&2
    exit 1
fi
if [[ ! -d "${INPUT_BASE}" ]]; then
    echo "ERROR: ${INPUT_BASE} not found -- run 2.5.1_extract_flanks.py first." >&2
    exit 1
fi

mkdir -p "${DUP_OUT}" "${CLU_OUT}"

echo "Input   : ${INPUT_BASE}"
echo "DupOut  : ${DUP_OUT}"
echo "CluOut  : ${CLU_OUT}"
echo "Threads : ${THREADS}"
echo "----------------------------------------------------------------------"

align_one_side() {
    # $1 = raw output path   $2 = aln output path   $3.. = input fasta files
    local raw_out="$1" ali_out="$2"
    shift 2
    cat "$@" > "${raw_out}"
    "${MAFFT_BIN}" --auto --thread "${THREADS}" --quiet "${raw_out}" > "${ali_out}"
}

# ==============================================================================
# LEVEL 1 — duplicate_group_level  (duplicate_group_* subfolders ONLY)
# ==============================================================================
echo ""
echo "--- LEVEL 1: duplicate_group_level (stacking + alignment, per group; singletons excluded) ---"

n_groups_done=0
n_groups_skipped_empty=0

for rec_dir in "${INPUT_BASE}"/*/; do
    [ -d "${rec_dir}" ] || continue
    rec_type="$(basename "${rec_dir}")"

    for cluster_dir in "${rec_dir}"cluster_*/; do
        [ -d "${cluster_dir}" ] || continue
        cluster_name="$(basename "${cluster_dir}")"

        # ONLY duplicate_group_* subfolders -- "singleton" (and anything
        # else) is never a candidate here, by construction of this glob.
        for group_dir in "${cluster_dir}"duplicate_group_*/; do
            [ -d "${group_dir}" ] || continue
            group_name="$(basename "${group_dir}")"

            out_dir="${DUP_OUT}/${rec_type}/${cluster_name}/${group_name}"

            mapfile -t left_fastas < <(find "${group_dir}" -maxdepth 1 -name "*_left.fasta" | sort)
            mapfile -t right_fastas < <(find "${group_dir}" -maxdepth 1 -name "*_right.fasta" | sort)
            if [ "${#left_fastas[@]}" -eq 0 ] && [ "${#right_fastas[@]}" -eq 0 ]; then
                n_groups_skipped_empty=$((n_groups_skipped_empty + 1))
                continue
            fi

            mkdir -p "${out_dir}"

            if [ "${#left_fastas[@]}" -gt 0 ]; then
                echo "  [L1-left]  ${rec_type}/${cluster_name}/${group_name}  (${#left_fastas[@]} seqs)"
                align_one_side "${out_dir}/left.raw.fasta" "${out_dir}/left.aln.fasta" "${left_fastas[@]}"
            fi
            if [ "${#right_fastas[@]}" -gt 0 ]; then
                echo "  [L1-right] ${rec_type}/${cluster_name}/${group_name}  (${#right_fastas[@]} seqs)"
                align_one_side "${out_dir}/right.raw.fasta" "${out_dir}/right.aln.fasta" "${right_fastas[@]}"
            fi
            n_groups_done=$((n_groups_done + 1))
        done
    done
done

echo "Level 1 done: ${n_groups_done} duplicate_group units aligned, ${n_groups_skipped_empty} empty ones skipped -- $(date)"

# ==============================================================================
# LEVEL 2 — cluster_level  (pooled from duplicate_group_* subfolders ONLY;
# a cluster with no duplicate_group_* subfolders at all -- e.g. singletons
# only -- is skipped entirely, no output directory created for it)
# ==============================================================================
echo ""
echo "--- LEVEL 2: cluster_level (stacking + alignment, whole cluster pooled; singleton-only clusters skipped) ---"

n_clusters_done=0
n_clusters_skipped_no_dupgroups=0

for rec_dir in "${INPUT_BASE}"/*/; do
    [ -d "${rec_dir}" ] || continue
    rec_type="$(basename "${rec_dir}")"
    out_rec="${CLU_OUT}/${rec_type}"

    for cluster_dir in "${rec_dir}"cluster_*/; do
        [ -d "${cluster_dir}" ] || continue
        cluster_name="$(basename "${cluster_dir}")"   # cluster_<N>

        # This cluster's duplicate_group_* subfolders -- the ONLY thing
        # LEVEL 2 pools from. "singleton" is structurally excluded: it
        # simply never matches this glob.
        dup_group_dirs=("${cluster_dir}"duplicate_group_*/)
        if [ "${#dup_group_dirs[@]}" -eq 0 ]; then
            # cluster has no duplicate groups at all (singletons only, or
            # nothing) -- skip the WHOLE cluster at this level.
            n_clusters_skipped_no_dupgroups=$((n_clusters_skipped_no_dupgroups + 1))
            continue
        fi

        left_fastas=()
        right_fastas=()
        for group_dir in "${dup_group_dirs[@]}"; do
            [ -d "${group_dir}" ] || continue
            while IFS= read -r -d '' f; do left_fastas+=("$f"); done \
                < <(find "${group_dir}" -maxdepth 1 -name "*_left.fasta" -print0)
            while IFS= read -r -d '' f; do right_fastas+=("$f"); done \
                < <(find "${group_dir}" -maxdepth 1 -name "*_right.fasta" -print0)
        done
        IFS=$'\n' left_fastas=($(sort <<<"${left_fastas[*]-}")); unset IFS
        IFS=$'\n' right_fastas=($(sort <<<"${right_fastas[*]-}")); unset IFS

        if [ "${#left_fastas[@]}" -eq 0 ] && [ "${#right_fastas[@]}" -eq 0 ]; then
            n_clusters_skipped_no_dupgroups=$((n_clusters_skipped_no_dupgroups + 1))
            continue
        fi

        mkdir -p "${out_rec}"

        if [ "${#left_fastas[@]}" -gt 0 ]; then
            raw_left="${out_rec}/${rec_type}_${cluster_name}_left.raw.fasta"
            ali_left="${out_rec}/${rec_type}_${cluster_name}_left.aln.fasta"
            echo "  [L2-left]  ${rec_type}/${cluster_name}  (${#left_fastas[@]} seqs, pooled from ${#dup_group_dirs[@]} duplicate group(s))"
            align_one_side "${raw_left}" "${ali_left}" "${left_fastas[@]}"
        fi
        if [ "${#right_fastas[@]}" -gt 0 ]; then
            raw_right="${out_rec}/${rec_type}_${cluster_name}_right.raw.fasta"
            ali_right="${out_rec}/${rec_type}_${cluster_name}_right.aln.fasta"
            echo "  [L2-right] ${rec_type}/${cluster_name}  (${#right_fastas[@]} seqs, pooled from ${#dup_group_dirs[@]} duplicate group(s))"
            align_one_side "${raw_right}" "${ali_right}" "${right_fastas[@]}"
        fi
        n_clusters_done=$((n_clusters_done + 1))
    done
done

echo "Level 2 done: ${n_clusters_done} clusters aligned, ${n_clusters_skipped_no_dupgroups} singleton-only/empty clusters skipped -- $(date)"
echo ""
echo "======================================================================"
echo "All alignments complete -- $(date)"
echo "======================================================================"