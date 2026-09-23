#!/usr/bin/env bash
# =============================================================================
# 2.2.2_skani_dedup.sh  —  Section 2.2.2
# =============================================================================
# Runs `skani triangle` within each specI subfolder (never across species),
# producing one all-vs-all pairwise ANI TSV per folder.
#
# Converted from a Python wrapper to pure bash: this step is just
# orchestrating calls to an external binary (skani) -- Python added
# subprocess-management overhead without adding anything bash can't do
# natively, and native bash + sbatch is the more idiomatic way to run a
# tool like this on this cluster.
#
# THIS IS AN ACTUAL SLURM JOB -- submit it with sbatch, don't just run it:
#     sbatch 2.2.2_skani_dedup.sh
#
# LOGGING: --output/--error below point directly at BASE_DIR/logs/, so
# everything this job prints is captured there by SLURM itself, with the
# job ID in the filename (2.2.2_skani_dedup_<jobid>.out/.err). This is
# deliberately NOT done via `exec > >(tee ...)` inside the script body --
# that pattern hangs when combined with the background (`&`) parallelism
# below (the tee subprocess's pipe doesn't close cleanly while background
# children still hold its stdout open), which testing caught before this
# shipped. Letting SLURM's own --output/--error do the capturing sidesteps
# that failure mode entirely.
#
# BASE_DIR / SKANI_BIN below should match paths.base_dir / paths.skani_path
# in config.yaml -- kept as plain bash variables here (not read from the
# YAML) so this script has zero Python/PyYAML dependency; keep them in
# sync by hand if either changes.
#
# SCRIPT_DIR-style path detection is deliberately NOT used (no
# $BASH_SOURCE tricks) -- sbatch copies the submitted script into a spool
# directory and runs it from there, so path detection based on the
# script's own location would resolve to that spool directory, not to
# wherever this script actually lives. Every path below is hardcoded.
# =============================================================================

#SBATCH --job-name=2_2_2_skani_dedup
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --mem=40G
#SBATCH --output=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.2.2_skani_dedup_%j.out
#SBATCH --error=/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep/logs/2.2.2_skani_dedup_%j.err
#SBATCH --mail-type=END,FAIL

set -euo pipefail

# ── Config (keep in sync with config.yaml) ─────────────────────────────────
BASE_DIR="/net/bq-storage/ag-khedkar/Sofia/project_folder/thesis_rep"
SKANI_BIN="/net/bq-storage/ag-khedkar/Sofia/project_folder/skani"

INPUT_DIR="${BASE_DIR}/genome_seqs_div_by_speci"
OUTPUT_DIR="${BASE_DIR}/files/skani_output"

MAX_PARALLEL_SPECIES=4   # how many specI folders to process at once
THREADS_PER_SPECIES=8    # (MAX_PARALLEL_SPECIES * THREADS_PER_SPECIES) should <= --cpus-per-task above

# PLACEHOLDER: set to 1 to symlink an already-existing skani_output/ instead
# of actually running skani; set PLACEHOLDER_PATH accordingly.
PLACEHOLDER=0
PLACEHOLDER_PATH=""

# (Optional) module loads / environment activation for your cluster:
# module purge
# module load devel/python/3.10.0

echo "======================================================================"
echo "Section 2.2.2 — skani 100% ANI dedup within each specI folder"
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

if [[ ! -x "${SKANI_BIN}" ]]; then
    echo "ERROR: skani binary not found or not executable at ${SKANI_BIN}" >&2
    exit 1
fi
if [[ ! -d "${INPUT_DIR}" ]]; then
    echo "ERROR: ${INPUT_DIR} not found -- run 2.2.1_divide_genomes_by_speci.py first." >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"

echo "Input dir           : ${INPUT_DIR}"
echo "Output dir          : ${OUTPUT_DIR}"
echo "Max parallel species: ${MAX_PARALLEL_SPECIES}"
echo "Threads per species : ${THREADS_PER_SPECIES}"
echo "----------------------------------------------------------------------"

n_species_total=0

process_species() {
    local dir="$1"
    local name
    name=$(basename "$dir")
    local fasta_count
    fasta_count=$(ls "${dir}"/*.fasta 2>/dev/null | wc -l)

    if [[ "${fasta_count}" -gt 1 ]]; then
        "${SKANI_BIN}" triangle "${dir}"/*.fasta \
            -t "${THREADS_PER_SPECIES}" \
            -E -s 0 --min-af 0 \
            -o "${OUTPUT_DIR}/${name}_ani_results.tsv"
        echo "  [${name}] ${fasta_count} genomes -> ${name}_ani_results.tsv"
    else
        echo "  [${name}] skipped (only ${fasta_count} genome(s), nothing to dedup)"
    fi
}
export -f process_species
export SKANI_BIN THREADS_PER_SPECIES OUTPUT_DIR

for SPECIES_DIR in "${INPUT_DIR}"/*/; do
    n_species_total=$((n_species_total + 1))
    process_species "${SPECIES_DIR}" &

    if [[ $(jobs -r -p | wc -l) -ge ${MAX_PARALLEL_SPECIES} ]]; then
        wait -n
    fi
done
wait

n_tsv=$(find "${OUTPUT_DIR}" -maxdepth 1 -name "*.tsv" | wc -l)

echo "----------------------------------------------------------------------"
echo "SUMMARY"
echo "----------------------------------------------------------------------"
echo "specI folders scanned : ${n_species_total}"
echo "TSV outputs written   : ${n_tsv}"
echo "-> ${OUTPUT_DIR}"
echo "Finished : $(date)"
echo "Section 2.2.2 complete."