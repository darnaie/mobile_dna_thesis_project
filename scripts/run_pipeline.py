#!/usr/bin/env python3
"""
run_pipeline.py
==================
Runs every Section 2.1-2.4 script in order, as separate subprocesses
(each exactly as if you'd run it by hand), stopping at the first failure.

logs/run_pipeline.log contains the SAME content that scrolls past in the
terminal -- every line each child step prints is streamed live to both
stdout and the log file, verbatim, as it happens (not buffered until the
child finishes). Orchestrator-level bookkeeping (step boundaries,
start/stop) is written the same way, tagged [run_pipeline] so it's easy
to tell apart from a step's own output.

Steps come in two kinds now:
  .py steps  -- run as `python3 <script>.py --config <config>`, exactly
                as before. Whatever's configured for that section
                (placeholder: true/false, create_plot: yes/no,
                use_slurm: true/false) still applies, since this just
                calls the script -- it doesn't bypass any of that logic.
  .sh steps  -- the skani/mmseqs2/vclust wrapper scripts (Sections 2.2.2,
                2.3.1, 2.4.1) are plain bash/sbatch scripts, not Python,
                so they're invoked differently depending on
                slurm.use_slurm in config.yaml:
                  use_slurm: true  -> `sbatch <script>.sh` (submitted,
                                       asynchronous -- see caveat below)
                  use_slurm: false -> `bash <script>.sh` (run directly,
                                       blocking, e.g. for local testing)

CAVEAT: if a step submits itself to SLURM (either a .py step with
use_slurm: true, or a .sh step submitted via sbatch), this runner does
NOT wait for that job to finish before moving on to the next step, since
it has no way to know when a submitted job actually completes. If you're
using SLURM, run the pipeline step-by-step instead (or in chunks, via
--from/--to) so you can confirm each job has actually finished before the
next one needs its output.

Run:
    python run_pipeline.py [--config path/to/config.yaml] [--from 2.1.1] [--to 2.2.5]
"""

from __future__ import annotations

import argparse
import datetime
import subprocess
import sys
from pathlib import Path

from utils import load_config, get_paths

SCRIPTS_DIR = Path(__file__).resolve().parent

# (section_key, script filename) in run order
PIPELINE = [
    ("2.1.1", "2.1.1_extract_mges.py"),
    ("2.1.2", "2.1.2_recombinase_gff.py"),
    ("2.1.3", "2.1.3_download_f13.py"),
    ("2.1.4", "2.1.4_download_genomes.py"),
    ("2.1.5", "2.1.5_recombinase_seq.py"),
    ("2.2.1", "2.2.1_divide_genomes_by_speci.py"),
    ("2.2.2", "2.2.2_skani_dedup.sh"),
    ("2.2.3", "2.2.3_pick_representatives.py"),
    ("2.2.4", "2.2.4_build_recombinase_table.py"),
    ("2.2.5", "2.2.5_plot_flowchart_section_2_2.py"),
    ("2.3.1", "2.3.1_cluster_recombinases.sh"),
    ("2.3.2", "2.3.2_append_clustering_info.py"),
    ("2.3.3", "2.3.3_plot_flowchart_section_2_3.py"),
    ("2.4.1", "2.4.1_vclust_dedup.sh"),
    ("2.4.2", "2.4.2_append_duplicate_info.py"),
    ("2.4.3", "2.4.3_plot_flowchart_section_2_4.py"),
    ("2.5.1", "2.5.1_extract_flanks.py"),
    ("2.5.2", "2.5.2_append_flank_sequences.py"),
    ("2.6.1", "2.6.1_align_flanks.sh"),
    ("2.6.2", "2.6.2_consensus_sanity_check.py"),
    ("2.6.3", "2.6.3_build_consensus.py"),
    ("2.6.4", "2.6.4_append_consensus.py"),
    ("2.7.1", "2.7.1_annotate_ml_labels.py"),
    ("2.7.2", "2.7.2_split_dataset.py"),
    ("2.7.3", "2.7.3_train_evaluate_model.py"),
    ("2.7.4", "2.7.4_kfold_cross_validation.py"),
    ("2.7.5", "2.7.5_plot_ml_results.py")
]


def _ts() -> str:
    return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def log_and_print(fh, msg: str) -> None:
    """Orchestrator's own bookkeeping messages -- timestamped and tagged
    so they're visually distinct from a step's own (already-timestamped)
    output."""
    line = f"{_ts()} [run_pipeline] {msg}"
    print(line)
    fh.write(line + "\n")
    fh.flush()


def passthrough(fh, line: str) -> None:
    """A child step's own output, written through UNCHANGED -- what you
    see in the terminal is exactly what lands in the log file."""
    print(line)
    fh.write(line + "\n")
    fh.flush()


def run_command(fh, cmd: list[str]) -> int:
    try:
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            text=True, bufsize=1,
        )
    except FileNotFoundError as exc:
        passthrough(fh, f"ERROR: could not run {cmd!r}: {exc}")
        if cmd[0] == "sbatch":
            passthrough(fh, "`sbatch` not found -- is SLURM available on this machine? "
                            "Set slurm.use_slurm to false in config.yaml to run locally instead.")
        return 1
    for line in process.stdout:
        passthrough(fh, line.rstrip("\n"))
    process.wait()
    return process.returncode


def run_step(fh, script_path: Path, config_path: str, use_slurm: bool) -> int:
    if script_path.suffix == ".sh":
        if use_slurm:
            log_and_print(fh, f"slurm.use_slurm=true -- submitting via: sbatch {script_path}")
            return run_command(fh, ["sbatch", str(script_path)])
        else:
            log_and_print(fh, f"slurm.use_slurm=false -- running directly (blocking): bash {script_path}")
            return run_command(fh, ["bash", str(script_path)])
    else:
        return run_command(fh, [sys.executable, str(script_path), "--config", config_path])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default=None)
    ap.add_argument("--from", dest="from_key", default=None,
                     help="First section to run (e.g. 2.2.1) -- skips everything before it")
    ap.add_argument("--to", dest="to_key", default=None,
                     help="Last section to run (e.g. 2.1.5) -- stops after it")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    config_path = args.config or str(Path(paths["base_dir"]) / "config.yaml")
    log_path = paths["logs_dir"] / "run_pipeline.log"
    use_slurm = bool(cfg.get("slurm", {}).get("use_slurm", False))

    keys = [k for k, _ in PIPELINE]
    start = keys.index(args.from_key) if args.from_key else 0
    end = keys.index(args.to_key) + 1 if args.to_key else len(PIPELINE)
    steps = PIPELINE[start:end]

    with open(log_path, "w") as fh:
        log_and_print(fh, "=" * 70)
        log_and_print(fh, f"Running {len(steps)} step(s): {[k for k, _ in steps]}")
        log_and_print(fh, f"slurm.use_slurm = {use_slurm}")
        log_and_print(fh, "=" * 70)

        for section_key, script_name in steps:
            log_and_print(fh, f"--- {section_key} : {script_name} " + "-" * 40)
            script_path = SCRIPTS_DIR / script_name
            if not script_path.exists():
                log_and_print(fh, f"ERROR: {script_path} not found -- stopping pipeline here.")
                raise SystemExit(1)

            returncode = run_step(fh, script_path, config_path, use_slurm)
            if returncode != 0:
                log_and_print(
                    fh,
                    f"Section {section_key} ({script_name}) exited with code "
                    f"{returncode} -- stopping pipeline here.",
                )
                raise SystemExit(returncode)

        log_and_print(fh, "=" * 70)
        log_and_print(fh, "Pipeline run complete.")
        log_and_print(fh, "=" * 70)


if __name__ == "__main__":
    main()