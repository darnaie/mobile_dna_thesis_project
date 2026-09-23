# Welcome to Transposon boundary identification master's project!

This repository contains the complete data-processing and machine-learning pipeline used to localise bacterial transposon boundaries directly from flanking DNA sequence, developed as part of a master's thesis. It combines an alignment-based consensus approach, a convolutional neural network, and a terminal-inverted-repeat (TIR) search adapted from ISEScan, all built on top of the proGenomes3 database.

For the full scientific write-up -- motivation, methods, results, and discussion -- see the thesis itself. This README covers only how to actually run the code.

## Repository structure

Scripts are named `<section>.<subsection>_<description>.py` (or `.sh` for standalone SLURM/bash steps), e.g. `2.1.1_extract_mges.py`, `2.6.1_align_flanks.sh`. They correspond to Methods Sections 2.1-2.7 of the thesis, run roughly in that order. `run_pipeline.py` runs every Section 2.1-2.7 step in sequence as a single command.

A few files are not numbered pipeline steps:
- `utils.py` -- shared helpers (config loading, GFF3 parsing, path resolution) every script imports from.
- `is_tn_helpers.py`, `consensus_building_evaluation.py`, `ml_model_config.py` -- shared logic imported by more than one numbered step; see each file's own docstring for which ones.
- `find_terminal_inverted_repeats.py` and `run_find_terminal_inverted_repeats.sh` -- the standalone TIR-identification analysis (Methods 2.3 / Results 3.8). This one is deliberately **not** part of `config.yaml` or `run_pipeline.py` -- it reads directly from two paths you set at the top of the `.sh` wrapper (see below), rather than through the shared config.

## Requirements

See `requirements.txt` for Python packages. You'll also need four external command-line tools this pipeline calls but doesn't install: **skani**, **mmseqs2**, **vclust**, and **mafft**. Install each yourself, then point `config.yaml`'s `paths.skani_path` / `paths.mmseqs2_path` / `paths.vclust_path` / `paths.mafft_path` at their binaries.

## Using `config.yaml`

This is the one file you need to edit before running anything.

**1. Set your real paths.** At the top, under `paths:`, set `base_dir` (where all pipeline output should live -- leave `""` to default to wherever `config.yaml` itself sits), `progenomes3_dir` (your local copy of the proGenomes3 database), and the four tool paths mentioned above. Every other directory the pipeline uses (`files/`, `logs/`, `plots/`, `genome_seqs/`, etc.) is derived automatically from `base_dir` -- you don't set those individually.

**2. Understand the placeholder mechanism.** Every data-producing step has two keys:
```yaml
placeholder: false
path_to_output: "NA"
```
- `placeholder: false` (the default in this repo, as delivered) means the step actually computes its output for real.
- `placeholder: true` means: instead of computing anything, copy or symlink an *already-existing* file/directory (given in `path_to_output`) into place, and skip the real computation. This is meant for re-running a later step without recomputing an expensive earlier one you already have output for -- point `path_to_output` at your own prior output, not at anything in this repo (there isn't any committed output here).

If you're running the pipeline for the first time, leave every `placeholder` as `false` and every `path_to_output` as `"NA"` -- that's the state this repo ships in.

**3. A few steps don't read this file at all.** `2.2.2_skani_dedup.sh`, `2.3.1_cluster_recombinases.sh`, `2.4.1_vclust_dedup.sh`, and `2.6.1_align_flanks.sh` are standalone SLURM/bash scripts with their **own** hardcoded paths and their **own** `PLACEHOLDER=0/1` toggle near the top of each `.sh` file -- editing `config.yaml`'s corresponding block has no effect on these four. Edit the variables inside the `.sh` file itself instead (`BASE_DIR`, the relevant `*_BIN` variable, and `PLACEHOLDER`/`PLACEHOLDER_PATH`).

**4. Adjust method parameters as needed.** Most sections also carry their own tunable parameters in the same block -- e.g. `ani_threshold` and `align_threshold` for genome deduplication (2.2.3), `fixed_threshold` and `a_run_length` for consensus building (2.6.3), or `seq_length`/`batch_size`/`learning_rate` for CNN training (2.7.3). See the comment above each section in `config.yaml` for what it reads and what it writes.

## Running the pipeline

Run everything in order:
```bash
python run_pipeline.py --config path/to/config.yaml
```
Run a specific range only:
```bash
python run_pipeline.py --config path/to/config.yaml --from 2.4.1 --to 2.5.2
```
Or run a single script directly (every script also works standalone):
```bash
python 2.1.1_extract_mges.py --config path/to/config.yaml
```

If `slurm.use_slurm: true` in `config.yaml`, a step submits itself via `sbatch` and returns immediately rather than blocking -- `run_pipeline.py` does not wait for a submitted job to finish before starting the next step. If you're using SLURM, it's safer to run the pipeline in chunks via `--from`/`--to`, confirming each SLURM job has actually completed before starting the next chunk, rather than running the whole thing in one call.

## Running TIR identification

This one is separate from the rest, since it isn't part of `config.yaml`:
```bash
sbatch run_find_terminal_inverted_repeats.sh
```
Edit `SCRIPT_DIR`, `FLANKS_DIR`, `IS_TN_DIR`, and `OUT_TSV` near the top of that file to your own paths first. See the script's own docstring, or Methods Section 2.3 of the thesis, for what it does and how its thresholds were chosen.

## A note on file numbering

The numbering of files in this repository (`2.1.1`, `2.6.3`, etc.) broadly follows the Methods sections of the thesis, but doesn't correspond to it one-to-one -- the thesis went through several rounds of editing after the code was written, and section numbers shifted in ways the filenames didn't always follow. This doesn't affect what any script does. As a rough map: files `2.1.x` through `2.7.x` correspond to Methods Sections 2.1, 2.2, and Results Sections 3.1-3.7 of the thesis; TIR identification (`find_terminal_inverted_repeats.py`) corresponds to Methods Section 2.3 and Results Section 3.8.
