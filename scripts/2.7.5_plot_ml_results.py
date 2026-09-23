#!/usr/bin/env python3
"""
2.7.5_plot_ml_results.py  —  Section 2.7.5
==============================================
input : BASE/ml_model/train_log.tsv (2.7.3, per-epoch train/val loss)
        BASE/files/ml_test_predictions.tsv (2.7.3, per-row test predictions)
        BASE/ml_model/attributions/*.tsv (2.7.3, per-base IG attributions
        -- requires captum to have been installed when 2.7.3 ran, else
        this directory won't exist and the attribution plot is skipped)
output: BASE/plots/ml_training_curve.png / .pdf
        BASE/plots/ml_per_rec_type_accuracy.png / .pdf
        BASE/plots/ml_prediction_distribution_by_rec_type.png / .pdf
        BASE/plots/ml_attributions_grid.png / .pdf
        BASE/plots/ml_model_architecture.png / .pdf

Five plots, meant to be thesis-ready as-is (300 dpi PNG + vector PDF,
readable at print size, colourblind-safe palette):

  1. TRAINING CURVE -- train_loss and val_loss per epoch (log y-axis,
     since train_loss here spans more than a 10x range), with the
     actual checkpointed epoch (lowest val_loss -- the one
     best_model.pt was saved from) marked explicitly. Log scale also
     makes the val_loss volatility across epochs easier to read than a
     linear axis would, without exaggerating it.

  2. PER-REC_TYPE TIER ACCURACY -- one bar per rec_type (each prediction
     rounded to whichever of {0.3, 0.7, 1.0} it's nearest to, compared
     against the true label's own nearest tier), sample size (n)
     annotated on each bar, sorted descending so the ordering itself
     tells the story, with the overall (pooled) accuracy drawn as a
     reference line. A rec_type with a test split of 0 (see 2.7.2's own
     warning about strata too small to split 80:10:10 -- your
     DDE_Tnp_ISAZ013 run hit exactly this) is drawn as an explicit "no
     test data" marker rather than silently omitted, so it doesn't read
     as an oversight.

  3. PREDICTION DISTRIBUTION BY TRUE LABEL TIER, per rec_type -- small
     multiples (one panel per rec_type with test data), x-axis = true
     label tier, y-axis = the model's raw (continuous) prediction, as
     strip+box plots. Plot 2 answers "how often is the model right,
     coarsely"; this one shows WHY, at the resolution the model
     actually operates at -- whether the three tiers separate cleanly
     as three distinct predicted-value clusters (good calibration) or
     overlap heavily (confusion), and whether that differs by family.

  4. ATTRIBUTIONS GRID (NEW) -- every attribution file 2.7.3 wrote (all
     of them, not a sample of them -- one small panel each), showing
     per-position attribution MAGNITUDE (L1 norm across the A/C/G/T
     channels, since the saved attribution files record all four
     channels' scores but not which base was actually present at each
     position -- see the note in plot_attributions()'s own docstring for
     why summing magnitude, rather than guessing the observed base from
     sign, is what's actually defensible here). Each panel is shaded
     over whichever region 2.7.1 originally used as that example's own
     training label, wherever that's recoverable, so you can see
     directly -- not just as the single 0/N number 2.7.3 already logs --
     whether the model's attribution actually concentrates where the
     alignment-based scan originally found its evidence, or somewhere
     else entirely.

  5. MODEL ARCHITECTURE (NEW) -- a block diagram of the actual network
     2.7.3 trains, built directly from ml_model_config() (imported from
     2.7.3, not re-typed here -- so if you tune the architecture later,
     regenerating this plot picks the change up automatically instead of
     silently going stale).

Run standalone:
    python 2.7.5_plot_ml_results.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.7.5"

TIERS = (0.3, 0.7, 1.0)

# Colourblind-safe qualitative palette (Okabe-Ito), reused across all three plots.
PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"]

plt.rcParams.update({
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "figure.dpi": 100,
    "savefig.dpi": 300,
    "axes.spines.top": False,
    "axes.spines.right": False,
})


def nearest_tier(x: float, tiers=TIERS) -> float:
    return min(tiers, key=lambda t: abs(t - x))


def tier_accuracy(preds: List[float], labels: List[float]) -> float:
    correct = sum(1 for p, y in zip(preds, labels) if nearest_tier(p) == nearest_tier(y))
    return correct / len(preds) if preds else float("nan")


def save_both(fig, out_dir: Path, name: str) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / f"{name}.png", bbox_inches="tight")
    fig.savefig(out_dir / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 1 -- training curve
# ══════════════════════════════════════════════════════════════════════════

def plot_training_curve(train_log_rows: List[dict], out_dir: Path, log) -> None:
    epochs = [int(r["epoch"]) for r in train_log_rows]
    train_loss = [float(r["train_loss"]) for r in train_log_rows]
    val_loss = [float(r["val_loss"]) for r in train_log_rows]

    best_idx = min(range(len(val_loss)), key=lambda i: val_loss[i])
    best_epoch = epochs[best_idx]
    best_val = val_loss[best_idx]

    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(epochs, train_loss, color=PALETTE[0], marker="o", markersize=3, linewidth=1.5, label="Train loss")
    ax.plot(epochs, val_loss, color=PALETTE[1], marker="o", markersize=3, linewidth=1.5, label="Validation loss")
    ax.axvline(best_epoch, color=PALETTE[2], linestyle="--", linewidth=1.2, alpha=0.8)
    ax.annotate(
        f"checkpoint saved\n(epoch {best_epoch}, val loss {best_val:.4f})",
        xy=(best_epoch, best_val), xytext=(best_epoch + 0.06 * (max(epochs) - min(epochs) + 1), best_val * 2.2),
        fontsize=9, color=PALETTE[2],
        arrowprops=dict(arrowstyle="->", color=PALETTE[2], lw=1),
    )
    ax.set_yscale("log")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("MSE loss (log scale)")
    ax.set_title("Training and validation loss")
    ax.legend(frameon=False, loc="upper right")
    ax.grid(axis="y", alpha=0.25, which="both")
    fig.tight_layout()
    save_both(fig, out_dir, "ml_training_curve")
    log.info(f"  Best checkpoint: epoch {best_epoch}, val_loss={best_val:.5f} "
             f"(training ran {len(epochs)} epochs before early stopping)")
    log.info(f"-> {out_dir / 'ml_training_curve.png'} / .pdf")


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 2 -- per-rec_type tier accuracy
# ══════════════════════════════════════════════════════════════════════════

def plot_per_rec_type_accuracy(
    pred_rows: List[dict], all_rec_types: List[str], out_dir: Path, log,
) -> None:
    by_rt: Dict[str, List[dict]] = {}
    for row in pred_rows:
        by_rt.setdefault(row["rec_type"], []).append(row)

    overall_preds = [float(r["prediction"]) for r in pred_rows]
    overall_labels = [float(r["label"]) for r in pred_rows]
    overall_acc = tier_accuracy(overall_preds, overall_labels)

    stats = []
    for rt in all_rec_types:
        rows = by_rt.get(rt, [])
        if not rows:
            stats.append((rt, None, 0))
            continue
        preds = [float(r["prediction"]) for r in rows]
        labels = [float(r["label"]) for r in rows]
        stats.append((rt, tier_accuracy(preds, labels), len(rows)))

    # sort: rec_types WITH data first (descending accuracy), then any with no data
    with_data = sorted([s for s in stats if s[1] is not None], key=lambda s: -s[1])
    without_data = [s for s in stats if s[1] is None]
    ordered = with_data + without_data

    fig, ax = plt.subplots(figsize=(8, 4.5))
    x = range(len(ordered))
    heights = [s[1] if s[1] is not None else 0.0 for s in ordered]
    colors = [PALETTE[i % len(PALETTE)] for i in range(len(ordered))]
    bars = ax.bar(x, heights, color=colors, width=0.6)

    for i, (rt, acc, n) in enumerate(ordered):
        if acc is None:
            ax.text(i, 0.02, "no test\ndata", ha="center", va="bottom", fontsize=8.5, style="italic", color="#666666")
        else:
            ax.text(i, acc + 0.015, f"n={n}", ha="center", va="bottom", fontsize=8.5, color="#333333")

    ax.axhline(overall_acc, color="#333333", linestyle=":", linewidth=1.2)
    ax.text(len(ordered) - 0.4, overall_acc + 0.015, f"overall = {overall_acc:.3f}",
            ha="right", fontsize=9, color="#333333")

    ax.set_xticks(list(x))
    ax.set_xticklabels([s[0] for s in ordered], rotation=30, ha="right")
    ax.set_ylabel("Tier accuracy")
    ax.set_ylim(0, 1.12)
    ax.set_title("Test-set tier accuracy by recombinase type")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    save_both(fig, out_dir, "ml_per_rec_type_accuracy")

    log.info(f"  Overall test tier accuracy: {overall_acc:.3f}")
    for rt, acc, n in ordered:
        log.info(f"    {rt:20s} n={n:5d}  tier_accuracy={'NA' if acc is None else f'{acc:.3f}'}")
    log.info(f"-> {out_dir / 'ml_per_rec_type_accuracy.png'} / .pdf")


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 3 -- prediction distribution by true label tier, faceted by rec_type
# ══════════════════════════════════════════════════════════════════════════

def plot_prediction_distribution(pred_rows: List[dict], all_rec_types: List[str], out_dir: Path, log) -> None:
    """
    Box-plot summary (median/quartiles/whiskers) with jittered scatter
    points overlaid -- predicted value (y) vs. true label tier (x),
    one panel per rec_type with test data.

    NOTE on this version: an earlier iteration of this plot dropped the
    box overlay entirely, keeping only the scatter -- brought back here
    (interpreting "where did the boxes go" as wanting both together, not
    a return to box-only), since the box gives an at-a-glance summary
    the raw scatter alone doesn't, while the scatter still shows the
    real per-point spread and sample density the box alone would hide.
    If pure scatter without boxes was actually what was wanted, this is
    the one line to change back (drop the ax.boxplot call below).

    Fixed 3-column layout, not min(4, n_panels) -- with the current 6
    rec_types that have test data (IS1, IS1595, IS66, ISL3, Tn3, IS110),
    this gives exactly two full rows of 3 panels each, rather than the
    4-then-2 layout a 4-column grid produces for 6 panels. If the number
    of rec_types with test data ever changes, this still degrades
    gracefully (fewer than 3 -> that many columns; not a multiple of 3 ->
    the last row just has empty slots turned off, same as before).
    """
    by_rt: Dict[str, List[dict]] = {}
    for row in pred_rows:
        by_rt.setdefault(row["rec_type"], []).append(row)

    rec_types_with_data = [rt for rt in all_rec_types if by_rt.get(rt)]
    n_panels = len(rec_types_with_data)
    if n_panels == 0:
        log.info("  No rec_types with test data -- skipping prediction-distribution plot.")
        return

    n_cols = min(3, n_panels)
    n_rows = -(-n_panels // n_cols)  # ceil division
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4.2 * n_cols, 3.8 * n_rows), sharey=True)
    axes_flat = axes.flatten() if n_panels > 1 else [axes]

    for i, rt in enumerate(rec_types_with_data):
        ax = axes_flat[i]
        rows = by_rt[rt]
        by_tier: Dict[float, List[float]] = {t: [] for t in TIERS}
        for row in rows:
            true_tier = nearest_tier(float(row["label"]))
            by_tier[true_tier].append(float(row["prediction"]))

        present_tiers = [t for t in TIERS if by_tier[t]]
        box_data = [by_tier[t] for t in present_tiers]
        positions = list(range(len(present_tiers)))

        bp = ax.boxplot(box_data, positions=positions, widths=0.5, patch_artist=True,
                         showfliers=False, zorder=2)
        for patch, t in zip(bp["boxes"], present_tiers):
            patch.set_facecolor(PALETTE[TIERS.index(t)])
            patch.set_alpha(0.35)
        for median in bp["medians"]:
            median.set_color("#333333")
            median.set_linewidth(1.3)

        for pos, t in zip(positions, present_tiers):
            ys = by_tier[t]
            jitter = [pos + (hash((rt, t, j)) % 1000 / 1000 - 0.5) * 0.5 for j in range(len(ys))]
            ax.scatter(jitter, ys, s=9, color=PALETTE[TIERS.index(t)], alpha=0.4, linewidths=0, zorder=3)

        ax.set_xticks(positions)
        ax.set_xticklabels([str(t) for t in present_tiers])
        ax.set_xlim(-0.5, len(present_tiers) - 0.5)
        ax.set_title(f"{rt}\n(n={len(rows)})", fontsize=10)
        ax.set_ylim(-0.05, 1.05)
        ax.grid(axis="y", alpha=0.2)
        if i % n_cols == 0:
            ax.set_ylabel("Predicted value")
        if i >= n_panels - n_cols:
            ax.set_xlabel("True label tier")

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.suptitle("Predicted value vs. true label tier, by recombinase type", y=1.02, fontsize=13)
    fig.tight_layout()
    save_both(fig, out_dir, "ml_prediction_distribution_by_rec_type")
    log.info(f"-> {out_dir / 'ml_prediction_distribution_by_rec_type.png'} / .pdf")


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 4 -- attribution grid (all examples)
# ══════════════════════════════════════════════════════════════════════════

def parse_attribution_filename(stem: str, known_rec_types: List[str]) -> Optional[tuple]:
    """Attribution filenames are '<rec_type>_<mge_id>_<side>' (see
    2.7.3's compute_attributions()) -- but both rec_type (e.g.
    'DDE_Tnp_IS1') and mge_id (e.g. 'MGE_GCA_...') routinely contain
    underscores themselves, so splitting on '_' alone is ambiguous.
    Resolved by matching the KNOWN rec_type as an exact prefix (longest
    match first, in case one rec_type name is itself a prefix of
    another) and the KNOWN 'left'/'right' as the exact suffix -- whatever
    remains in between is the mge_id, whole, underscores and all."""
    for side in ("left", "right"):
        suffix = f"_{side}"
        if not stem.endswith(suffix):
            continue
        body = stem[: -len(suffix)]
        for rt in sorted(known_rec_types, key=len, reverse=True):
            prefix = f"{rt}_"
            if body.startswith(prefix):
                mge_id = body[len(prefix):]
                if mge_id:
                    return rt, mge_id, side
    return None


def load_segment_lookup(labels_split_path: Path) -> Dict[tuple, Optional[tuple]]:
    """(rec_type, mge_id, side) -> (segment_start, segment_length) in that
    row's OWN original flank coordinates, or None if 2.7.1 didn't find
    one for it. Coordinate conversion for the right side (see below) is
    deliberately NOT done here -- it depends on seq_length, which this
    function doesn't know about; done by the caller instead."""
    lookup: Dict[tuple, Optional[tuple]] = {}
    with open(labels_split_path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            key = (row["rec_type"], row["mge_id"], row["side"])
            s, l = row.get("segment_start"), row.get("segment_length")
            if s not in (None, "NA") and l not in (None, "NA"):
                try:
                    lookup[key] = (int(s), int(l))
                    continue
                except ValueError:
                    pass
            lookup[key] = None
    return lookup


def plot_attributions(
    attr_dir: Path, labels_split_path: Optional[Path], seq_length: int,
    out_dir: Path, log,
) -> None:
    """
    Every *.tsv file under attr_dir gets its own small panel -- ALL of
    them, not a sample, since 2.7.3 already did the sampling (via
    --n-attribution-examples-per-rec-type) when it decided which test
    examples to compute attributions for in the first place.

    Per position, the four saved channels (A/C/G/T) are reduced to a
    single "attribution magnitude" (L1 norm: |A|+|C|+|G|+|T|) for
    plotting. NOT reduced to "the attribution of whichever base was
    actually there" -- that would require knowing the actual input
    sequence at each position, which 2.7.3's attribution output doesn't
    currently record (only the four raw scores per position), and
    guessing the base from which channel's attribution happens to be
    largest isn't reliable enough to build a plot on. Magnitude is the
    honest reduction given what's actually saved: it still shows WHERE
    along the flank the model's prediction is sensitive to the input,
    which is the main thing a small-multiples overview like this is for.

    Each panel is shaded over the boundary segment 2.7.1 originally used
    to build that example's OWN training label, when labels_split_path
    is given and a segment is recoverable for it. IMPORTANT coordinate
    note: RIGHT-side flanks are reverse-complemented before being fed to
    the model (see 2.7.3's FlankDataset), so the attribution's own
    position axis for a right-side example is in REVCOMP'D coordinates,
    while segment_start/segment_length from ml_labels_split.tsv are in
    the flank's ORIGINAL coordinates (2.7.1 never reverses anything).
    Converted here exactly the same way 2.7.3's own compute_attributions()
    already does for its "peak_in_found_segment" check: for a segment at
    [s, s+l) in a sequence of length seq_length, its position after a
    full-sequence reversal is [seq_length-(s+l), seq_length-s). Getting
    this wrong (i.e. overlaying RIGHT-side shading using unconverted
    coordinates) would silently mis-place the shaded region for every
    right-side panel -- exactly the class of bug this pipeline has
    already hit once, on the trim_start fix.

    X-AXIS RANGE (per panel, not shared) -- BUG FIX, see chat: this used
    to restrict each panel's x-axis to that example's CLUSTER-LEVEL
    ALIGNMENT length (read from its .aln.fasta under consensus/
    cluster_level/). That was wrong: a mafft alignment's length reflects
    gap columns inserted to accommodate indels across ALL of that
    cluster's members, and is routinely LONGER than any individual raw
    flank -- some panels were ending up with x-axes stretching out to
    800, 1300, even 2500, while the actual attribution data (always
    exactly seq_length positions, since that's the model's fixed one-hot
    input width) stopped at ~530, leaving a long, flat, empty-looking
    tail that looked like missing or gap-padded data. It wasn't gap
    padding -- it was an axis limit taken from the wrong coordinate
    space entirely (the aligned "alignment" view, not the raw "stacking"
    view the model itself actually sees).

    Fixed by using each row's own flank_sequence length from
    ml_labels_split.tsv directly instead -- exactly the raw, unaligned
    sequence 2.7.3's FlankDataset one-hot-encodes (padded to seq_length
    with all-zero columns if naturally shorter, which happens only for
    flanks truncated at a contig edge in 2.5.1). No separate file read
    needed at all now; this is simply len(flank_sequence) for that
    example's own row. Falls back to the full seq_length if
    labels_split_path wasn't given or the row isn't found.
    """
    attr_files = sorted(attr_dir.glob("*.tsv"))
    if not attr_files:
        log.info(f"  No attribution files found in {attr_dir} -- skipping (captum may not have been "
                 f"installed when 2.7.3 ran; see its log for 'captum not installed').")
        return

    # known_rec_types has to come from somewhere real -- derive from the
    # labels-split file if available, else fall back to a permissive
    # heuristic (strip the last two '_'-joined segments as mge_id+side)
    # so this still degrades gracefully without it.
    labels_rows: List[dict] = []
    if labels_split_path and labels_split_path.exists():
        with open(labels_split_path, newline="") as fh:
            labels_rows = list(csv.DictReader(fh, delimiter="\t"))
    if labels_rows:
        known_rec_types = sorted({row["rec_type"] for row in labels_rows})
    else:
        known_rec_types = sorted({fp.stem.rsplit("_", 2)[0] for fp in attr_files})

    segment_lookup = load_segment_lookup(labels_split_path) if labels_split_path and labels_split_path.exists() else {}
    raw_length_lookup: Dict[tuple, int] = {
        (row["rec_type"], row["mge_id"], row["side"]): len(row.get("flank_sequence", "")) for row in labels_rows
    }

    n_panels = len(attr_files)
    n_cols = min(6, n_panels)
    n_rows = -(-n_panels // n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(2.6 * n_cols, 2.3 * n_rows), sharex=False)
    axes_flat = axes.flatten() if n_panels > 1 else [axes]

    n_shaded = 0
    n_xlim_restricted = 0
    n_raw_length_missing = 0
    for i, fp in enumerate(attr_files):
        ax = axes_flat[i]
        parsed = parse_attribution_filename(fp.stem, known_rec_types)
        rec_type, mge_id, side = parsed if parsed else ("?", fp.stem, "?")

        with open(fp, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        positions = [int(r["position"]) for r in rows]
        magnitude = [abs(float(r["A"])) + abs(float(r["C"])) + abs(float(r["G"])) + abs(float(r["T"])) for r in rows]

        ax.fill_between(positions, magnitude, color=PALETTE[0], alpha=0.6, linewidth=0)
        ax.plot(positions, magnitude, color=PALETTE[0], linewidth=0.6)

        seg = segment_lookup.get((rec_type, mge_id, side)) if parsed else None
        if seg is not None:
            s, l = seg
            if side == "left":
                s = seq_length - (s + l)  # see docstring: same conversion as 2.7.3's own check
            ax.axvspan(s, s + l, color=PALETTE[2], alpha=0.25, linewidth=0)
            n_shaded += 1

        raw_len = raw_length_lookup.get((rec_type, mge_id, side)) if parsed else None
        if raw_len:
            ax.set_xlim(0, raw_len)
            n_xlim_restricted += 1
        else:
            n_raw_length_missing += 1
            ax.set_xlim(0, seq_length)

        short_id = mge_id if len(mge_id) <= 16 else mge_id[:7] + "…" + mge_id[-6:]
        ax.set_title(f"{rec_type}\n{short_id} ({side})", fontsize=7.5)
        ax.tick_params(labelsize=7)
        ax.set_yticks([])
        ax.set_xlabel("Position", fontsize=8)  # every panel now (sharex=False, so labels no longer shared)

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.suptitle(
        "Per-base attribution magnitude (Integrated Gradients), all sampled test examples\n"
        "(shaded band = boundary segment 2.7.1 used for that example's own training label; "
        "x-axis limited to that example's own raw flank length, not the aligned length)",
        y=1.01, fontsize=11,
    )
    fig.tight_layout()
    save_both(fig, out_dir, "ml_attributions_grid")
    log.info(f"  x-axis restricted to raw flank length for {n_xlim_restricted}/{n_panels} panel(s)"
             + (f" ({n_raw_length_missing} row(s) not found in ml_labels_split.tsv, fell back to full seq_length)"
                if n_raw_length_missing else ""))
    log.info(f"  {n_panels} attribution panel(s) plotted, {n_shaded} with a shaded label-segment overlay")
    log.info(f"-> {out_dir / 'ml_attributions_grid.png'} / .pdf")


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 6 -- per-rec_type PWM / sequence logo, built from each example's
#  own highest-attribution window
# ══════════════════════════════════════════════════════════════════════════

def build_rec_type_pwm_windows(
    attr_dir: Path, labels_split_path: Path, all_rec_types: List[str],
    window_len: int, log,
) -> Dict[str, List[str]]:
    """
    For every attribution file belonging to a given rec_type, slides a
    window of `window_len` positions across that example's own
    attribution-magnitude curve, finds the window with the HIGHEST total
    magnitude (i.e. the window the model's own attribution says matters
    most for that specific example), and extracts the ACTUAL SEQUENCE at
    that window's position -- not from the attribution file itself
    (2.7.3 never saves the input sequence there, only the four per-
    channel attribution scores -- see plot_attributions()'s own
    docstring), but reconstructed from ml_labels_split.tsv's own
    flank_sequence column, reverse-complemented for right-side examples
    exactly the way 2.7.3's FlankDataset does before encoding, so the
    extracted window is drawn from the SAME sequence, in the SAME
    orientation and coordinate space, that actually produced the
    attribution scores being searched.

    Returns {rec_type: [window_seq, window_seq, ...]} -- one window per
    attribution file for that rec_type (only for files where a window of
    the full requested length was actually available and a matching
    flank_sequence was found; shorter alignments or unmatched examples
    are skipped and counted, not silently included as short/padded
    windows, since a PWM requires every row to be the same length).
    """
    attr_files = sorted(attr_dir.glob("*.tsv"))
    if not attr_files:
        return {}

    with open(labels_split_path, newline="") as fh:
        label_rows = list(csv.DictReader(fh, delimiter="\t"))
    seq_lookup = {
        (row["rec_type"], row["mge_id"], row["side"]): row["flank_sequence"] for row in label_rows
    }
    known_rec_types = sorted({row["rec_type"] for row in label_rows}) or all_rec_types

    windows_by_rec_type: Dict[str, List[str]] = {rt: [] for rt in all_rec_types}
    n_skipped_no_sequence = 0
    n_skipped_too_short = 0

    for fp in attr_files:
        parsed = parse_attribution_filename(fp.stem, known_rec_types)
        if not parsed:
            continue
        rec_type, mge_id, side = parsed

        flank_seq = seq_lookup.get((rec_type, mge_id, side))
        if flank_seq is None:
            n_skipped_no_sequence += 1
            continue
        encoded_seq = revcomp(flank_seq.upper()) if side == "left" else flank_seq.upper()

        with open(fp, newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        magnitude = [abs(float(r["A"])) + abs(float(r["C"])) + abs(float(r["G"])) + abs(float(r["T"])) for r in rows]

        if len(magnitude) < window_len or len(encoded_seq) < window_len:
            n_skipped_too_short += 1
            continue

        # sliding-window sum via cumulative sums -- O(n), not O(n * window_len)
        cumsum = [0.0]
        for m in magnitude:
            cumsum.append(cumsum[-1] + m)
        best_start, best_sum = 0, -1.0
        for start in range(0, len(magnitude) - window_len + 1):
            window_sum = cumsum[start + window_len] - cumsum[start]
            if window_sum > best_sum:
                best_sum = window_sum
                best_start = start

        if best_start + window_len > len(encoded_seq):
            n_skipped_too_short += 1
            continue
        window_seq = encoded_seq[best_start:best_start + window_len]
        if len(window_seq) == window_len and set(window_seq) <= set("ACGT"):
            windows_by_rec_type.setdefault(rec_type, []).append(window_seq)

    if n_skipped_no_sequence or n_skipped_too_short:
        log.info(f"  (PWM windows) skipped {n_skipped_no_sequence} example(s) with no matching "
                 f"flank_sequence, {n_skipped_too_short} with fewer than {window_len} positions available")

    return windows_by_rec_type


def plot_rec_type_pwms(
    attr_dir: Path, labels_split_path: Optional[Path], all_rec_types: List[str],
    window_len: int, out_dir: Path, log,
) -> None:
    """
    One sequence logo per recombinase type (small multiples), each built
    by stacking that family's own per-example highest-attribution
    windows (see build_rec_type_pwm_windows()) into a position weight
    matrix and rendering it with logomaker -- the same information-
    content-scaled letter-height convention as a standard WebLogo. This
    is the per-family, integrated-across-flanks counterpart to Figure 5
    (which shows one raw attribution curve per individual example): each
    panel here answers "if I only look at where each of this family's
    own flanks says its attribution is highest, and pool the actual
    sequence at those positions across every flank, what motif emerges" --
    directly comparable in spirit to a family's literature-reported
    terminal-repeat consensus.

    Requires labels_split_path (to recover each example's actual
    sequence, not saved in the attribution files themselves -- see
    build_rec_type_pwm_windows()) and captum-produced attribution files
    to already exist; skips (with a clear log message, not a silent
    empty plot) if either is unavailable, or if a given rec_type ends up
    with zero usable windows (e.g. every one of its examples' alignments
    was shorter than window_len).
    """
    if not labels_split_path or not labels_split_path.exists():
        log.info(f"  {labels_split_path} not found -- cannot recover flank sequences for PWM windows, skipping.")
        return
    if not attr_dir.is_dir() or not any(attr_dir.glob("*.tsv")):
        log.info(f"  No attribution files in {attr_dir} -- skipping PWM plot.")
        return

    try:
        import logomaker
    except ImportError:
        log.warning("  *** logomaker not installed -- ml_pwm_by_rec_type.png/.pdf will NOT be written. "
                    "Run `pip install logomaker` (or your cluster's equivalent) and re-run this script "
                    "to get it -- nothing else in this run is affected.")
        return

    windows_by_rec_type = build_rec_type_pwm_windows(attr_dir, labels_split_path, all_rec_types, window_len, log)
    rec_types_with_windows = [rt for rt in all_rec_types if len(windows_by_rec_type.get(rt, [])) >= 2]
    if not rec_types_with_windows:
        log.info("  No rec_type had at least 2 usable high-attribution windows -- skipping PWM plot.")
        return

    n_panels = len(rec_types_with_windows)
    n_cols = min(3, n_panels)
    n_rows = -(-n_panels // n_cols)
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(3.6 * n_cols, 2.4 * n_rows))
    axes_flat = axes.flatten() if n_panels > 1 else [axes]

    for i, rt in enumerate(rec_types_with_windows):
        ax = axes_flat[i]
        windows = windows_by_rec_type[rt]
        counts_matrix = logomaker.alignment_to_matrix(windows)
        logomaker.Logo(counts_matrix, ax=ax, color_scheme="classic")
        ax.set_title(f"{rt}\n(n={len(windows)} windows)", fontsize=9.5)
        ax.set_xlabel("Position in window", fontsize=8)
        ax.set_ylabel("bits", fontsize=8)
        ax.tick_params(labelsize=7)

    for j in range(n_panels, len(axes_flat)):
        axes_flat[j].axis("off")

    fig.suptitle(
        f"Per-recombinase-type sequence logo, built from each example's own {window_len} nt "
        f"highest-attribution window",
        y=1.02, fontsize=12,
    )
    fig.tight_layout()
    save_both(fig, out_dir, "ml_pwm_by_rec_type")
    log.info(f"  {n_panels} rec_type logo(s) plotted (window_len={window_len})")
    log.info(f"-> {out_dir / 'ml_pwm_by_rec_type.png'} / .pdf")


# ══════════════════════════════════════════════════════════════════════════
#  PLOT 5 -- model architecture diagram
# ══════════════════════════════════════════════════════════════════════════

def _architecture_boxes(n_rec_types: int, seq_length: int) -> list:
    """Shared box content/colors for both orientations -- single source so
    the vertical and horizontal diagrams can never drift apart in content,
    only in layout."""
    from ml_model_config import ml_model_config as _ml_model_config
    cfg = _ml_model_config()

    channels = cfg["conv_channels"]
    dilations = cfg["dilations"]
    kernel = cfg["kernel_size"]
    rt_dim = cfg["rec_type_embed_dim"]
    side_dim = cfg["side_embed_dim"]
    head_dim = cfg["head_hidden_dim"]

    boxes = [(f"Input flank\n(one-hot, 4 × {seq_length})", "#EDEDED")]
    for ch, dil in zip(channels, dilations):
        boxes.append((f"Conv1d {ch}ch, k={kernel}\ndilation={dil}\n+ BatchNorm + ReLU", PALETTE[0]))
    boxes.append(("Global average pool\n(over sequence axis)", PALETTE[1]))
    boxes.append((f"Concatenate with\nrec_type embed ({rt_dim}d, {n_rec_types} types)\n"
                   f"+ side embed ({side_dim}d)", PALETTE[3]))
    boxes.append((f"Linear → {head_dim}\nReLU + Dropout({cfg['dropout']})", PALETTE[4]))
    boxes.append(("Linear → 1\nSigmoid", PALETTE[4]))
    boxes.append(("Predicted confidence\n[0, 1]", "#EDEDED"))
    return boxes, channels, dilations, head_dim


def plot_model_architecture(scripts_dir: Path, n_rec_types: int, seq_length: int, out_dir: Path, log) -> None:
    """
    Block diagram of the actual network 2.7.3 trains, top-to-bottom.
    Architecture hyperparameters come from ml_model_config.py (a
    separate, dependency-light module -- see chat / that file's own
    docstring for why: this used to importlib the WHOLE 2.7.3 module to
    reach ml_model_config(), which also executed 2.7.3's own `import
    torch` in the process and broke this plot outright wherever torch
    isn't installed, even though drawing a box diagram never needed it.
    Direct import now, no torch involved -- still the same single-
    source-of-truth property (regenerating the plot always reflects
    whatever 2.7.3 would actually build), just without the accidental
    dependency.

    Sized up again per a follow-up request -- bigger text relative to
    box size, tighter vertical spacing between boxes (less unused
    whitespace), larger figure overall. See also
    plot_model_architecture_horizontal() below, the left-to-right
    counterpart with identical box content (from the same
    _architecture_boxes(), so the two can't drift apart), for when a
    wide layout suits the thesis page better than a tall one.

    Input box label now reads the real seq_length instead of a
    hardcoded 530 -- flanks are shorter now that 2.7.1 excludes the
    mge_overlap zone (see that script's docstring), so a fixed 530
    would have been wrong as soon as that landed.
    """
    boxes, channels, dilations, head_dim = _architecture_boxes(n_rec_types, seq_length)

    n = len(boxes)
    fig, ax = plt.subplots(figsize=(7.2, 2.0 * n))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, n)
    ax.axis("off")

    box_h = 0.82  # up from 0.72 -- less gap between boxes, less empty space overall
    for i, (text, color) in enumerate(boxes):
        y = n - i - 1
        is_light = color == "#EDEDED"
        rect = plt.Rectangle((0.02, y + (1 - box_h) / 2), 0.96, box_h,
                              facecolor=color, edgecolor="#333333", linewidth=1.6,
                              alpha=0.9 if not is_light else 1.0, zorder=2)
        ax.add_patch(rect)
        text_color = "#111111" if is_light else "white"
        ax.text(0.5, y + 0.5, text, ha="center", va="center", fontsize=17, fontweight="medium",
                 color=text_color, zorder=3)
        if i < n - 1:
            ax.annotate("", xy=(0.5, y - 0.015), xytext=(0.5, y + (1 - box_h) / 2 - 0.015),
                        arrowprops=dict(arrowstyle="-|>", color="#333333", lw=2.0, mutation_scale=22), zorder=1)

    ax.set_title("Boundary-motif CNN architecture", fontsize=22, pad=14)
    fig.tight_layout()
    save_both(fig, out_dir, "ml_model_architecture")

    log.info(f"  Diagram built from 2.7.3's actual ml_model_config(): "
             f"{len(channels)} conv blocks, dilations {dilations}, head_hidden_dim={head_dim}, seq_length={seq_length}")
    log.info(f"-> {out_dir / 'ml_model_architecture.png'} / .pdf")


def plot_model_architecture_horizontal(scripts_dir: Path, n_rec_types: int, seq_length: int, out_dir: Path, log) -> None:
    """
    Same content as plot_model_architecture() (see _architecture_boxes(),
    the single shared source for both), laid out left-to-right instead of
    top-to-bottom -- same information, different orientation for whatever
    fits the thesis page layout better.

    Sized up again per a follow-up request, same as the vertical version,
    PLUS the connecting arrows are deliberately more pronounced here than
    in the vertical version -- thicker shaft, bigger arrowhead
    (mutation_scale raised well above the vertical diagram's) -- since
    "more pronounced arrows" was specifically called out for this
    orientation. Tighter horizontal gap between boxes too, so less of
    the figure is empty space between them.
    """
    boxes, channels, dilations, head_dim = _architecture_boxes(n_rec_types, seq_length)

    n = len(boxes)
    fig, ax = plt.subplots(figsize=(3.1 * n, 5.2))
    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box_w = 0.94  # up from 0.90 -- tighter horizontal gap between boxes
    for i, (text, color) in enumerate(boxes):
        x = i
        is_light = color == "#EDEDED"
        rect = plt.Rectangle((x + (1 - box_w) / 2, 0.04), box_w, 0.92,
                              facecolor=color, edgecolor="#333333", linewidth=1.6,
                              alpha=0.9 if not is_light else 1.0, zorder=2)
        ax.add_patch(rect)
        text_color = "#111111" if is_light else "white"
        ax.text(x + 0.5, 0.5, text, ha="center", va="center", fontsize=15.5, fontweight="medium",
                 color=text_color, zorder=3)
        if i < n - 1:
            gap_start = x + 1 - (1 - box_w) / 2
            gap_end = x + 1 + (1 - box_w) / 2
            ax.annotate("", xy=(gap_end + 0.01, 0.5), xytext=(gap_start - 0.01, 0.5),
                        arrowprops=dict(arrowstyle="-|>", color="#222222", lw=3.2, mutation_scale=34), zorder=1)

    ax.set_title("Boundary-motif CNN architecture (left to right)", fontsize=22, pad=14)
    fig.tight_layout()
    save_both(fig, out_dir, "ml_model_architecture_horizontal")

    log.info(f"  Diagram built from 2.7.3's actual ml_model_config(): "
             f"{len(channels)} conv blocks, dilations {dilations}, head_hidden_dim={head_dim}, seq_length={seq_length}")
    log.info(f"-> {out_dir / 'ml_model_architecture_horizontal.png'} / .pdf")




def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--train-log", default=None, help="Override BASE/ml_model/train_log.tsv")
    ap.add_argument("--predictions", default=None, help="Override BASE/files/ml_test_predictions.tsv")
    ap.add_argument("--out-dir", default=None, help="Override BASE/plots")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    section_cfg = cfg["section_2_7"][SECTION_KEY]

    log = get_logger("2.7.5_plot_ml_results", paths["logs_dir"])

    train_log_path = Path(args.train_log) if args.train_log else paths["base_dir"] / "ml_model" / "train_log.tsv"
    predictions_path = Path(args.predictions) if args.predictions else paths["files_dir"] / "ml_test_predictions.tsv"
    out_dir = Path(args.out_dir) if args.out_dir else paths.get("plots_dir", paths["base_dir"] / "plots")

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Plot ML training/evaluation results")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, out_dir / "ml_training_curve.png", log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_7_4_plot_ml_results"):
        return

    if not train_log_path.exists():
        log.error(f"{train_log_path} not found -- run 2.7.3_train_evaluate_model.py first.")
        raise SystemExit(1)
    if not predictions_path.exists():
        log.error(f"{predictions_path} not found -- run 2.7.3_train_evaluate_model.py first.")
        raise SystemExit(1)

    all_rec_types_raw = section_cfg.get("rec_types", "all")
    seq_length = int(section_cfg.get("seq_length", 500))  # was 530 -- must match 2.7.3's own seq_length for the segment coordinate conversion to be correct (2.7.1 now excludes the mge_overlap zone)

    log.info(f"train_log   : {train_log_path}")
    log.info(f"predictions : {predictions_path}")
    log.info(f"out_dir     : {out_dir}")
    log.info("-" * 70)

    with open(train_log_path, newline="") as fh:
        train_log_rows = list(csv.DictReader(fh, delimiter="\t"))
    with open(predictions_path, newline="") as fh:
        pred_rows = list(csv.DictReader(fh, delimiter="\t"))

    log.info(f"Loaded {len(train_log_rows):,} training-log rows, {len(pred_rows):,} test predictions")

    if all_rec_types_raw == "all":
        all_rec_types = sorted({r["rec_type"] for r in pred_rows})
    else:
        all_rec_types = [t.strip() for t in str(all_rec_types_raw).split(",") if t.strip()]

    log.info("-" * 70)
    log.info("Plot 1/7: training curve ...")
    plot_training_curve(train_log_rows, out_dir, log)

    log.info("-" * 70)
    log.info("Plot 2/7: per-rec_type tier accuracy ...")
    plot_per_rec_type_accuracy(pred_rows, all_rec_types, out_dir, log)

    log.info("-" * 70)
    log.info("Plot 3/7: prediction distribution by true label tier ...")
    plot_prediction_distribution(pred_rows, all_rec_types, out_dir, log)

    log.info("-" * 70)
    log.info("Plot 4/7: attribution grid (all sampled test examples) ...")
    attr_dir = paths["base_dir"] / "ml_model" / "attributions"
    labels_split_path = predictions_path.parent / "ml_labels_split.tsv"
    if attr_dir.is_dir():
        plot_attributions(attr_dir, labels_split_path, seq_length, out_dir, log)
    else:
        log.info(f"  {attr_dir} not found -- 2.7.3 either hasn't run its attribution step, or captum "
                 f"wasn't installed when it did (check its own log). Skipping.")

    log.info("-" * 70)
    log.info("Plot 5/7: per-rec_type sequence logo from highest-attribution windows ...")
    pwm_window_len = int(section_cfg.get("pwm_window_len", 15))
    if attr_dir.is_dir():
        plot_rec_type_pwms(attr_dir, labels_split_path, all_rec_types, pwm_window_len, out_dir, log)
    else:
        log.info(f"  {attr_dir} not found -- skipping.")

    log.info("-" * 70)
    log.info("Plot 6/7: model architecture diagram (vertical) ...")
    plot_model_architecture(Path(__file__).resolve().parent, len(all_rec_types), seq_length, out_dir, log)

    log.info("-" * 70)
    log.info("Plot 7/7: model architecture diagram (horizontal) ...")
    plot_model_architecture_horizontal(Path(__file__).resolve().parent, len(all_rec_types), seq_length, out_dir, log)

    # ── final summary: check what's ACTUALLY on disk, independent of what
    # each plot function above logged -- deliberately not trusting each
    # function's own success message, since a plot that silently failed
    # to save (missing optional dependency, empty data, etc.) is exactly
    # the kind of thing that's easy to miss scrolling through a long log.
    # This is the last thing printed, so it's the thing you actually see.
    log.info("-" * 70)
    log.info("FINAL CHECK -- what's actually in " + str(out_dir))
    expected = [
        "ml_training_curve", "ml_per_rec_type_accuracy", "ml_prediction_distribution_by_rec_type",
        "ml_attributions_grid", "ml_pwm_by_rec_type", "ml_model_architecture", "ml_model_architecture_horizontal",
    ]
    any_missing = False
    for stem in expected:
        png_ok = (out_dir / f"{stem}.png").exists()
        pdf_ok = (out_dir / f"{stem}.pdf").exists()
        status = "OK" if (png_ok and pdf_ok) else ("PARTIAL" if (png_ok or pdf_ok) else "MISSING")
        if status != "OK":
            any_missing = True
        log.info(f"  [{status:7s}] {stem}.png / .pdf")
    if any_missing:
        log.warning("  *** at least one expected plot is missing or only partially written -- "
                    "scroll up to the matching 'Plot n/7' section above for why.")

    log.info("-" * 70)
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()