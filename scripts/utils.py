"""
utils.py
========
Shared helpers for the thesis_rep pipeline scripts. Import with:

    from utils import (
        load_config, get_logger, parse_gff_attributes, parse_mge_field,
        gca_from_mge_id, taxid_sample_from_mge_id, load_taxonomy_map,
    )

Every stage script lives in thesis_rep/scripts/ and, by default, loads
thesis_rep/config.yaml (one directory up). Nothing in here talks to the
filesystem except load_config / get_logger / load_taxonomy_map.
"""

from __future__ import annotations

import logging
import os
import sys
from pathlib import Path
from typing import Dict, Tuple

try:
    import yaml
except ImportError:  # pragma: no cover
    sys.exit(
        "PyYAML is required (pip install pyyaml / conda install pyyaml) "
        "to run thesis_rep scripts."
    )

REPO_ROOT = Path(__file__).resolve().parent.parent  # thesis_rep/


# ──────────────────────────────────────────────────────────────────────────
# Config
# ──────────────────────────────────────────────────────────────────────────

def default_config_path() -> Path:
    return REPO_ROOT / "config.yaml"


def _expand_base_dir(obj, base_dir: str):
    """
    Recursively replaces any string value that is exactly 'BASE' or starts
    with 'BASE/' with the real base_dir, throughout the WHOLE config
    structure. This is what lets config.yaml write things like:

        path_to_output: "BASE/files/master_mges.gff3"

    instead of always spelling out the full absolute path -- i.e. "wherever
    this output would normally land if the section ran for real". Applied
    once, here, at load time -- resolve_placeholder() /
    resolve_placeholder_multi() and every script that calls them never
    need to know this expansion happened; they just see final real paths.
    """
    if isinstance(obj, dict):
        return {k: _expand_base_dir(v, base_dir) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_expand_base_dir(v, base_dir) for v in obj]
    if isinstance(obj, str):
        if obj == "BASE":
            return base_dir
        if obj.startswith("BASE/"):
            return str(Path(base_dir) / obj[len("BASE/"):])
    return obj


def load_config(config_path: str | Path | None = None) -> dict:
    """
    Load config.yaml. Defaults to thesis_rep/config.yaml.

    paths.base_dir defaults to the directory CONTAINING this config file
    if not set explicitly. Every other directory (files/, logs/, plots/,
    slurm_logs/, genome_seqs/, ...) is DERIVED from base_dir by
    get_paths() below -- config.yaml's paths section only ever holds
    base_dir, progenomes3_dir, and the three tool paths (skani/mmseqs2/vclust).

    Any string anywhere in the config that is exactly "BASE" or starts
    with "BASE/" is expanded to the real base_dir (see _expand_base_dir) --
    e.g. path_to_output: "BASE/files/master_mges.gff3".
    """
    path = Path(config_path) if config_path else default_config_path()
    if not path.exists():
        sys.exit(f"Config file not found: {path}")
    with open(path) as fh:
        cfg = yaml.safe_load(fh)

    cfg.setdefault("paths", {})
    if not cfg["paths"].get("base_dir"):
        cfg["paths"]["base_dir"] = str(path.parent)

    cfg = _expand_base_dir(cfg, cfg["paths"]["base_dir"])

    return cfg


def get_paths(cfg: dict) -> Dict[str, Path]:
    """
    Derives every directory this pipeline reads/writes from paths.base_dir.
    config.yaml's `paths:` section itself only ever contains base_dir,
    progenomes3_dir, skani_path, mmseqs2_path, vclust_path -- nothing else.
    Creates the always-needed output directories if they don't exist yet.
    """
    base_dir = Path(cfg["paths"]["base_dir"])
    paths = {
        "base_dir": base_dir,
        "progenomes3_dir": Path(cfg["paths"]["progenomes3_dir"]),
        "skani_path": cfg["paths"].get("skani_path", "skani"),
        "mmseqs2_path": cfg["paths"].get("mmseqs2_path", "mmseqs"),
        "vclust_path": cfg["paths"].get("vclust_path", "vclust"),
        "files_dir": base_dir / "files",
        "logs_dir": base_dir / "logs",
        "plots_dir": base_dir / "plots",
        "slurm_logs_dir": base_dir / "slurm_logs",
        "genome_dir": base_dir / "genome_seqs",
        "genome_seqs_div_by_speci": base_dir / "genome_seqs_div_by_speci",
        "genome_seqs_repr": base_dir / "genome_seqs_repr",
    }
    for key in ("files_dir", "logs_dir", "plots_dir", "slurm_logs_dir"):
        paths[key].mkdir(parents=True, exist_ok=True)
    return paths


def resolve_placeholder(section_cfg: dict, output_path: str | Path, log) -> bool:
    """
    Generic PLACEHOLDER handling shared by every data-producing Section 2.x
    script (plotting subsections use the separate create_plot: yes/no
    switch instead -- see config.yaml).

    section_cfg is one subsection's config block, e.g.
    cfg["section_2_1"]["2.1.2"], which is expected to carry:
      placeholder     : bool
      path_to_output  : str  (a real path when placeholder is true; "NA" otherwise)

    If placeholder is true: stages path_to_output as output_path and
    returns True -- the caller should skip its real computation for this
    output file (just log and move on).
    If placeholder is false (or missing): returns False -- the caller
    should go ahead and actually compute this output.

    Raises clearly if placeholder is true but path_to_output isn't a real,
    existing path (rather than silently doing nothing).

    Safe to point path_to_output at the exact same path the real
    computation would have written to (e.g. you already ran this step for
    real once and just want future runs to skip recomputing it) -- that's
    detected and treated as "already staged", not an error.
    """
    if not section_cfg.get("placeholder", False):
        return False

    src = section_cfg.get("path_to_output")
    if not src or src in ("NA", "CHANGEME") or "CHANGEME" in str(src):
        raise ValueError(
            f"placeholder=true but path_to_output is not a real path (got {src!r}) "
            f"-- fill in config.yaml, or set placeholder=false to compute this for real."
        )
    src_path = Path(src)
    if not src_path.exists():
        raise FileNotFoundError(f"placeholder path_to_output does not exist: {src_path}")

    _stage_placeholder_path(src_path, Path(output_path), log)
    return True


def resolve_placeholder_multi(section_cfg: dict, output_map: Dict[str, "str | Path"], log) -> bool:
    """
    Like resolve_placeholder(), but for a subsection with MORE THAN ONE
    output file/directory. output_map maps {config_key: output_path}, e.g.:

        resolve_placeholder_multi(section_cfg, {
            "path_to_output": master_gff_out,
            "path_to_output_non_nested_gff": non_nested_gff_out,
            "path_to_output_non_nested_fasta": non_nested_fasta_out,
        }, log)

    When placeholder is true, EVERY key in output_map must have a real,
    existing path set in section_cfg -- there is no partial/"best effort"
    placeholder anymore. If you can't supply all of them, set
    placeholder to false and compute the section for real instead.
    """
    if not section_cfg.get("placeholder", False):
        return False

    missing = []
    resolved: Dict[str, Path] = {}
    for key in output_map:
        src = section_cfg.get(key)
        if not src or src in ("NA", "CHANGEME") or "CHANGEME" in str(src):
            missing.append(key)
        else:
            resolved[key] = Path(src)

    if missing:
        raise ValueError(
            f"placeholder=true but {len(missing)} of {len(output_map)} required "
            f"path(s) are not set: {missing}. This subsection has {len(output_map)} "
            f"outputs, so ALL of them need a real path when using a placeholder -- "
            f"fill them all in, or set placeholder=false to compute everything for real."
        )

    for key, src_path in resolved.items():
        if not src_path.exists():
            raise FileNotFoundError(f"{key} does not exist: {src_path}")

    for key, output_path in output_map.items():
        _stage_placeholder_path(resolved[key], Path(output_path), log)

    return True


def _stage_placeholder_path(src_path: Path, output_path: Path, log) -> None:
    """
    Copies (files) or symlinks (directories -- too large to duplicate) 
    src_path to output_path. If src_path and output_path already refer to
    the exact same file on disk (e.g. path_to_output was pointed at the
    real computation's own usual output location), does nothing but log
    that it's already staged, instead of erroring.
    """
    import shutil

    output_path.parent.mkdir(parents=True, exist_ok=True)

    try:
        already_same = (
            (output_path.exists() or output_path.is_symlink())
            and src_path.resolve() == output_path.resolve()
        )
    except OSError:
        already_same = False

    if already_same:
        log.info(f"PLACEHOLDER=true -- {output_path} already IS {src_path}; "
                 f"nothing to stage, using it as-is.")
        return

    if src_path.is_dir():
        if output_path.exists() or output_path.is_symlink():
            if output_path.is_symlink() or output_path.is_file():
                output_path.unlink()
            elif output_path.is_dir() and not any(output_path.iterdir()):
                output_path.rmdir()
        if not output_path.exists():
            output_path.symlink_to(src_path, target_is_directory=True)
        log.info(f"PLACEHOLDER=true -- symlinked {output_path} -> {src_path} "
                 f"(real computation skipped)")
    else:
        shutil.copy2(src_path, output_path)
        log.info(f"PLACEHOLDER=true -- copied {src_path} -> {output_path} "
                 f"(real computation skipped)")


# ──────────────────────────────────────────────────────────────────────────
# SLURM — optional. Every section script calls maybe_submit_and_exit()
# right after logging is set up; if slurm.use_slurm is true and this isn't
# already the submitted job, it re-invokes itself via sbatch and tells the
# caller to return immediately (the real work happens later, inside the
# submitted job).
# ──────────────────────────────────────────────────────────────────────────

SLURM_CHILD_ENV_VAR = "THESIS_REP_SLURM_CHILD"


def maybe_submit_and_exit(
    cfg: dict, script_file: str, cli_args: list[str], log,
    job_name: str, extra_sbatch_lines: list[str] | None = None,
    extra_shell_lines: list[str] | None = None,
) -> bool:
    """
    Returns True if the caller should STOP right now (a SLURM job was just
    submitted on its behalf). Returns False if the caller should proceed
    with the real computation (either because slurm is disabled, or
    because this process IS the already-submitted SLURM job).

    extra_shell_lines: raw shell lines inserted AFTER the #SBATCH block but
    BEFORE the python invocation -- e.g. `module load ...`, `source
    ~/env/bin/activate`. If not given, falls back to
    slurm.default_extra_shell_lines in config.yaml (applies to every
    section unless that section passes its own).
    """
    import subprocess
    import tempfile

    slurm_cfg = cfg.get("slurm", {})
    if not slurm_cfg.get("use_slurm", False):
        return False

    if os.environ.get(SLURM_CHILD_ENV_VAR) == "1":
        return False  # this IS the submitted job -- do the real work

    paths = get_paths(cfg)
    slurm_logs_dir = paths["slurm_logs_dir"]
    slurm_logs_dir.mkdir(parents=True, exist_ok=True)

    header_lines = [
        line for line in str(slurm_cfg.get("sbatch_header", "")).splitlines()
        if line.strip()
    ]

    shell_lines = extra_shell_lines
    if shell_lines is None:
        shell_lines = slurm_cfg.get("default_extra_shell_lines", []) or []

    script_lines = ["#!/usr/bin/env bash"]
    script_lines += header_lines
    script_lines += [
        f"#SBATCH --job-name={job_name}",
        f"#SBATCH --output={slurm_logs_dir}/{job_name}_%j.out",
        f"#SBATCH --error={slurm_logs_dir}/{job_name}_%j.err",
    ]
    if extra_sbatch_lines:
        script_lines += extra_sbatch_lines
    script_lines.append("")
    if shell_lines:
        script_lines += list(shell_lines)
        script_lines.append("")
    script_lines.append(f"export {SLURM_CHILD_ENV_VAR}=1")
    cmd = " ".join([sys.executable, str(script_file)] + [str(a) for a in cli_args])
    script_lines.append(cmd)

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".sbatch", prefix=f"{job_name}_",
        dir=str(slurm_logs_dir), delete=False,
    ) as tf:
        tf.write("\n".join(script_lines) + "\n")
        sbatch_path = tf.name

    log.info(f"use_slurm=true -- submitting via sbatch: {sbatch_path}")
    try:
        result = subprocess.run(
            ["sbatch", sbatch_path], capture_output=True, text=True, check=False
        )
        log.info(f"sbatch stdout: {result.stdout.strip()}")
        if result.returncode != 0:
            log.error(f"sbatch failed (exit {result.returncode}): {result.stderr.strip()}")
    except FileNotFoundError:
        log.error(
            "`sbatch` command not found -- is SLURM available on this machine? "
            "Set slurm.use_slurm to false in config.yaml to run locally instead."
        )
    return True



# ──────────────────────────────────────────────────────────────────────────
# Logging — every script gets its own log file under logs_dir, plus stdout.
# get_logger() ALSO tees sys.stdout itself, so plain print() calls (not
# just log.info()) end up in the log file too -- important when a script
# is run standalone rather than via run_pipeline.py (which already
# captures everything a step prints, regardless of how it prints it, at
# the subprocess level).
# ──────────────────────────────────────────────────────────────────────────

class _TeeStream:
    """Writes to both the original stream and a log file handle. Used to
    make bare print() calls show up in a script's own log file, without
    touching logging's own handlers (which write to a DIFFERENT,
    untouched reference to the original stdout -- see get_logger below --
    so nothing here causes log.info() lines to be duplicated)."""

    def __init__(self, original, log_fh):
        self._original = original
        self._log_fh = log_fh

    def write(self, data):
        self._original.write(data)
        try:
            self._log_fh.write(data)
            self._log_fh.flush()
        except ValueError:
            pass  # log file handle already closed -- don't crash on exit

    def flush(self):
        self._original.flush()
        try:
            self._log_fh.flush()
        except ValueError:
            pass

    def __getattr__(self, item):
        return getattr(self._original, item)


def get_logger(name: str, logs_dir: str | Path) -> logging.Logger:
    logs_dir = Path(logs_dir)
    logs_dir.mkdir(parents=True, exist_ok=True)
    log_path = logs_dir / f"{name}.log"

    original_stdout = sys.stdout  # captured BEFORE any tee reassignment

    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()  # avoid duplicate handlers on re-import/re-run

    fmt = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")

    # Truncate once, explicitly, for a fresh log per run. IMPORTANT: both
    # handles opened below use APPEND ('a', i.e. O_APPEND) mode, not 'w' --
    # two separate file objects writing to the same path only interleave
    # correctly if both reseek to the true end-of-file before every write
    # (which is what O_APPEND guarantees at the OS level). Mixing one 'w'
    # handle (which just advances its OWN remembered position, blind to
    # writes made through the other handle) with one 'a' handle silently
    # clobbers whichever one writes second -- this was a real bug caught
    # by testing before delivery.
    open(log_path, "w").close()

    fh = logging.FileHandler(log_path, mode="a")
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    # Bound to the ORIGINAL stdout object, not the tee -- so log.info()
    # writes go straight to the real terminal once, not through the tee
    # (which would otherwise also copy them into the log file a second time).
    sh = logging.StreamHandler(original_stdout)
    sh.setFormatter(fmt)
    logger.addHandler(sh)

    raw_fh = open(log_path, "a")
    sys.stdout = _TeeStream(original_stdout, raw_fh)

    logger.info(f"Logging to {log_path}")
    return logger


# ──────────────────────────────────────────────────────────────────────────
# GFF3 attribute parsing
# ──────────────────────────────────────────────────────────────────────────

def parse_gff_attributes(attr_col: str) -> Dict[str, str]:
    """Parse a GFF3 column-9 attribute string into a dict."""
    result: Dict[str, str] = {}
    for item in attr_col.split(";"):
        item = item.strip()
        if "=" in item:
            k, _, v = item.partition("=")
            result[k] = v
    return result


def parse_mge_field(mge_value: str) -> Dict[str, int]:
    """
    Parse a GFF3 'mge=' attribute value into {type: count}.

    Examples
    --------
    "is_tn:1"          -> {"is_tn": 1}
    "is_tn:2,ce:1"      -> {"is_tn": 2, "ce": 1}
    "is_tn"             -> {"is_tn": 1}      # no explicit count -> assume 1
    ""                  -> {}

    This is what lets us report the true number of *recombinase instances*
    of each type (summed across the ":N" suffixes), as opposed to the
    number of MGE *island rows*, which under-counts whenever a nested
    island bundles more than one recombinase together
    (e.g. "is_tn:2,ce:1" is ONE gff3 row but represents 3 recombinases:
    2 is_tn-type + 1 ce-type).
    """
    counts: Dict[str, int] = {}
    if not mge_value:
        return counts
    for token in mge_value.split(","):
        token = token.strip()
        if not token:
            continue
        if ":" in token:
            mtype, _, n = token.partition(":")
            try:
                n = int(n)
            except ValueError:
                n = 1
        else:
            mtype, n = token, 1
        counts[mtype] = counts.get(mtype, 0) + n
    return counts


def is_tn_count(mge_value: str) -> int:
    """Convenience wrapper: how many is_tn recombinase instances does this
    'mge=' field represent? (0 if none)."""
    return parse_mge_field(mge_value).get("is_tn", 0)


def total_recombinase_count(mge_value: str) -> int:
    """Sum of all type:count pairs in a 'mge=' field."""
    return sum(parse_mge_field(mge_value).values())


def parse_duplicate_groups_file(path: str | Path) -> list[set]:
    """
    vclust-style duplicates file: one line = one duplicate group; each
    token is an mge_id, optionally prefixed with '+' or '-' (strand of
    the match relative to the group's first/representative member),
    stripped here to recover the bare mge_id. Returns a list of sets, one
    per line/group (empty lines skipped).
    """
    groups = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            members = set()
            for token in line.split():
                if token and token[0] in "+-":
                    token = token[1:]
                if token:
                    members.add(token)
            if members:
                groups.append(members)
    return groups


def load_mge_sign_map(path: str | Path) -> Dict[str, str]:
    """
    Reads a vclust-style duplicates file and returns {mge_id: sign},
    sign in {'+', '-'}. Each line is one duplicate group; each token is
    an mge_id optionally prefixed with '+' or '-' -- the sign of the
    match's orientation relative to the group's first/reference member.
    An UNSIGNED token (no '+'/'-' prefix at all -- typically the
    reference member itself, compared trivially to itself) defaults to
    '-' ("no flip needed"), matching the original pipeline script's own
    parsing convention. mge_ids that never appear in the file at all
    (true singletons) simply won't be present in the returned dict --
    callers should default those to '-' too (same reasoning: nothing to
    compare against, so no orientation flip is warranted).
    """
    sign_map: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            for token in line.split():
                if not token:
                    continue
                if token[0] in "+-":
                    sign, mge_id = token[0], token[1:]
                else:
                    sign, mge_id = "-", token
                if mge_id:
                    sign_map[mge_id] = sign
    return sign_map


def find_duplicates_file(base_dir: str | Path) -> Path | None:
    """
    Finds the vclust *.duplicates.txt file under BASE/vclust_deduplication/,
    whatever its exact prefix (glob-based -- see 2.4.1_vclust_dedup.sh's
    header comment on why an exact filename shouldn't be assumed).
    """
    vclust_dir = Path(base_dir) / "vclust_deduplication"
    if not vclust_dir.exists():
        return None
    candidates = sorted(vclust_dir.glob("*.duplicates.txt"))
    return candidates[0] if candidates else None


def revcomp(seq: str) -> str:
    """Reverse-complement, tolerant of ambiguous 'N'/lowercase (matches
    the original pipeline script's own translate table exactly)."""
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def parse_mmseqs_cluster_tsv(path: str | Path) -> Tuple[set, Dict[str, str]]:
    """
    mmseqs easy-cluster's <prefix>_cluster.tsv: representative<TAB>member
    per line (a representative also appears as a member of itself).
    Returns (representatives set, {member_id: representative_id}).
    """
    representatives: set = set()
    member_to_rep: Dict[str, str] = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split("\t")
            if len(cols) < 2:
                continue
            rep_id, member_id = cols[0].strip(), cols[1].strip()
            representatives.add(rep_id)
            member_to_rep[member_id] = rep_id
    return representatives, member_to_rep


# ──────────────────────────────────────────────────────────────────────────
# ID parsing:  MGE_GCA_001031825.1_573.SAMN03280285.KQ088489:1814719-1815699
# ──────────────────────────────────────────────────────────────────────────

def _split_mge_id(mge_id: str) -> list[str]:
    return mge_id.split("_")


def gca_from_mge_id(mge_id: str) -> str | None:
    """
    'MGE_GCA_001031825.1_573.SAMN03280285.KQ088489:...' -> 'GCA_001031825.1'
    (i.e. the text between the 1st and 3rd underscore).
    """
    parts = _split_mge_id(mge_id)
    if len(parts) < 3 or parts[1] != "GCA":
        return None
    return f"GCA_{parts[2]}"


def taxid_sample_from_mge_id(mge_id: str) -> str | None:
    """
    'MGE_GCA_001031825.1_573.SAMN03280285.KQ088489:1814719-1815699'
    -> '573.SAMN03280285'   (taxid.sampleid, used as the genome/contig key)
    """
    parts = _split_mge_id(mge_id)
    if len(parts) < 4:
        return None
    tail = parts[3]                    # '573.SAMN03280285.KQ088489:...'
    dot_parts = tail.split(".")
    if len(dot_parts) >= 2 and dot_parts[0].isdigit():
        return f"{dot_parts[0]}.{dot_parts[1]}"
    return None


def taxid_sample_from_contig_id(contig_id: str) -> str:
    """
    '573.SAMN03280285.KQ088490' -> '573.samn03280285'
    Same "first two dot-separated tokens, lower-cased" key used for
    matching contig IDs (gff3 column 1) to genome FASTA files.
    """
    parts = contig_id.split(".")
    if len(parts) >= 2:
        return ".".join(parts[:2]).lower()
    return contig_id.lower()


# ──────────────────────────────────────────────────────────────────────────
# Flank output directory layout -- SHARED between Section 2.5.1 (writes the
# per-MGE flank fastas) and 2.5.2 (reads them back to append flank_sequence
# to the table). Both import this rather than reimplementing the naming, so
# the writer and reader can never drift out of sync with each other.
# ──────────────────────────────────────────────────────────────────────────

def flank_group_dir_name(duplicate_group_number: str, singleton: str) -> str:
    """
    The innermost folder name for an MGE's flank files, based on its
    Section 2.4 duplicate-group assignment. Singletons (singleton == "yes",
    or no real duplicate_group_number) get a fixed "singleton" folder
    instead of a number, so the two cases are unambiguous on disk.
    """
    dup = "" if duplicate_group_number is None else str(duplicate_group_number).strip()
    if str(singleton).strip().lower() == "yes" or not dup or dup == "NA":
        return "singleton"
    return f"duplicate_group_{dup}"


def flank_output_dir(flanks_base: str | Path, rec_type: str, cluster_number: str,
                      duplicate_group_number: str, singleton: str) -> Path:
    """
    <flanks_base>/<rec_type>/cluster_<cluster_number>/<duplicate_group_N | singleton>/
    -- the directory a given MGE's *_left.fasta / *_right.fasta live in.
    """
    rec_type = (rec_type or "UNKNOWN").strip() or "UNKNOWN"
    cluster_number = (str(cluster_number).strip() if cluster_number else "NA") or "NA"
    group_dir = flank_group_dir_name(duplicate_group_number, singleton)
    return Path(flanks_base) / rec_type / f"cluster_{cluster_number}" / group_dir


def parse_header_attributes(header: str) -> Dict[str, str]:
    """
    Parse key=value pairs out of a FASTA header, tolerant of both
    semicolon- and whitespace-delimited attribute strings, e.g.:
      '>ID mge=is_tn:1;mge_type=non-nested;size=21454'
      '>ID mge=is_tn:1 mge_type=non-nested size=21454'
    Keys are lower-cased on the way in (values are not) so callers don't
    have to worry about 'Mge=' vs 'mge=' vs 'MGE=' -- see
    parse_ffn_header_mge_fields() below for a stricter, fully
    case-insensitive variant used specifically for the mge/mge_type pair.
    """
    import re

    body = header[1:] if header.startswith(">") else header
    attrs: Dict[str, str] = {}
    for token in re.split(r"[;\s]+", body):
        if "=" in token:
            k, _, v = token.partition("=")
            attrs[k.strip().lower()] = v.strip()
    return attrs


def parse_ffn_header_mge_fields(header: str) -> Tuple[Dict[str, int], str]:
    """
    Case-insensitive extraction of the 'mge=' and 'mge_type=' fields from a
    .ffn.gz FASTA header. Mirrors the original, known-working approach of
    matching against a FULLY LOWERCASED copy of the header (the original
    pipeline script did `h = header.lower()` before substring-checking for
    "mge=is_tn" / "mge_type=non-nested" -- if that lowering is skipped,
    any case variation in the real header, e.g. 'Mge=IS_TN:1', silently
    fails every match and produces an empty output file with no error).

    Returns (mge_counts, mge_type):
      mge_counts : {type: count} dict, e.g. {'is_tn': 1}, keys lower-cased
                   (parsed via parse_mge_field on the lower-cased value)
      mge_type   : 'non-nested' / 'nested' / '' (lower-cased)
    """
    import re

    h_lower = header.lower()
    m = re.search(r"mge=([^;\s]+)", h_lower)
    mge_counts = parse_mge_field(m.group(1)) if m else {}
    mt = re.search(r"mge_type=([^;\s]+)", h_lower)
    mge_type = mt.group(1) if mt else ""
    return mge_counts, mge_type


def build_taxid_sample_to_gca_map(gff_path: str | Path) -> Dict[str, str]:
    """
    Scans a master_mges.gff3-style file's ID attributes and builds
    {taxid.sampleid: GCA_accession} by extracting both from each MGE ID
    (e.g. 'MGE_GCA_001031825.1_573.SAMN03280285.KQ088489:...' contributes
    '573.SAMN03280285' -> 'GCA_001031825.1'). Used to relate genome FASTA
    files (named by taxid.sampleid) back to their GCA accession, which is
    what the f13 taxonomy table is keyed on.
    """
    mapping: Dict[str, str] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attrs = parse_gff_attributes(cols[8])
            mge_id = attrs.get("ID", "")
            if not mge_id:
                continue
            taxid_sample = taxid_sample_from_mge_id(mge_id)
            gca = gca_from_mge_id(mge_id)
            if taxid_sample and gca and taxid_sample not in mapping:
                mapping[taxid_sample] = gca
    return mapping


def species_and_genomes_from_mge_ids(
    mge_ids, taxonomy_map: Dict[str, Tuple[str, str]]
) -> Dict[str, object]:
    """
    Given an iterable of mge_id strings, extracts each one's GCA and
    relates it to specI via taxonomy_map. Returns:
      gca_set       : {GCA, ...}
      specI_to_gca  : {specI: {GCA, ...}}
      unmapped_gca  : count of GCAs with no taxonomy entry
    """
    gca_set: set = set()
    specI_to_gca: Dict[str, set] = {}
    unmapped_gca = 0
    for mge_id in mge_ids:
        gca = gca_from_mge_id(mge_id)
        if not gca:
            continue
        gca_set.add(gca)
        info = taxonomy_map.get(gca)
        if info is None:
            unmapped_gca += 1
            continue
        specI, _label = info
        specI_to_gca.setdefault(specI, set()).add(gca)
    return {"gca_set": gca_set, "specI_to_gca": specI_to_gca, "unmapped_gca": unmapped_gca}



#   f13_staged_data_with_speci_and_tax.txt, TAB-separated:
#   GCA  specI  taxid  species(+strain)  strain
# ──────────────────────────────────────────────────────────────────────────

def load_taxonomy_map(staged_file: str | Path) -> Dict[str, Tuple[str, str]]:
    """
    Returns {GCA_accession: (specI_id, species_label)}.
    Lines that don't parse cleanly are skipped (and counted by the caller
    if it wants to log that).
    """
    mapping: Dict[str, Tuple[str, str]] = {}
    staged_file = Path(staged_file)
    if not staged_file.exists():
        return mapping

    with open(staged_file) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 4:
                continue
            gca = cols[0].strip()
            spec_i = cols[1].strip()
            species = cols[3].strip()
            strain = cols[4].strip() if len(cols) > 4 else ""
            label = (species + " " + strain).strip()
            mapping[gca] = (spec_i, label)
    return mapping


# ──────────────────────────────────────────────────────────────────────────
# Simple FASTA reader (used by the sequence-extraction scripts)
# ──────────────────────────────────────────────────────────────────────────

def read_fasta(path: str | Path) -> Dict[str, str]:
    """Parse a (non-gzipped) multi-contig FASTA into {header_id: sequence}."""
    seqs: Dict[str, str] = {}
    name = None
    buf: list[str] = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
        if name is not None:
            seqs[name] = "".join(buf)
    return seqs


def flank_file_path(
    flanks_dir: str | Path, rec_type: str, cluster_number: str,
    duplicate_group: str, mge_id: str, side: str,
) -> Path:
    """
    Single shared definition of where one MGE's flank fasta lives:

        BASE/flanks_split_acc_to_duplicates/<rec_type>/cluster_<cluster_number>/<duplicate_group>/<mge_id>_<side>.fasta

    Used identically by Section 2.5.1 (writes these files) and 2.5.2
    (reads them back to append flank_sequence to the table) -- kept here,
    in exactly one place, so the two scripts can't silently drift apart
    on the directory convention.

    duplicate_group is used as-is, including the literal string "NA" for
    singletons (2.4.2 writes "NA" there when a row has no duplicate group;
    this function does not special-case that -- singletons just get their
    own "NA" subfolder within their cluster, same as any other group value).
    """
    dup_group = duplicate_group if duplicate_group else "NA"
    return (
        Path(flanks_dir) / str(rec_type) / f"cluster_{cluster_number}"
        / str(dup_group) / f"{mge_id}_{side}.fasta"
    )


def variant_matches(mge_type: str, variant: str) -> bool:
    """Does a row with this mge_type belong in the requested output variant?"""
    if variant == "merged":
        return mge_type in ("nested", "non-nested")
    if variant == "non_nested":
        return mge_type == "non-nested"
    if variant == "nested":
        return mge_type == "nested"
    return False