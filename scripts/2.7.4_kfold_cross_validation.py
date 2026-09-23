#!/usr/bin/env python3
"""
2.7.4_kfold_cross_validation.py  —  Section 2.7.4
=====================================================
input : BASE/files/ml_labels_split.tsv (2.7.2 -- specifically its "fold"
        column; the "split" column, used by single-split 2.7.3, is
        ignored here)
output: BASE/ml_model/cv_fold_<k>/best_model.pt   (one checkpoint per fold)
        BASE/ml_model/cv_fold_<k>/train_log.tsv   (per-epoch loss, per fold)
        BASE/files/ml_cv_results.tsv              (per-fold x per-rec_type
                                                     test metrics, long format)
        BASE/files/ml_cv_summary.tsv              (mean +/- std across
                                                     folds, overall and per
                                                     rec_type)
        BASE/plots/ml_cv_results.png / .pdf        (cross-validated
                                                     tier-accuracy spread)

WHAT THIS IS, AND WHY IT'S SEPARATE FROM 2.7.3
---------------------------------------------------
2.7.3 trains ONE model on ONE 80:10:10 split. That gives a single test-
set accuracy number per rec_type -- useful, but it conflates two
different kinds of uncertainty: how good the model really is, and how
lucky or unlucky that ONE particular split happened to be, especially
for smaller rec_types (150-250 test examples in the single-split run
already showed real fold-to-fold-sized differences even from ordinary
sampling noise). K-fold cross-validation addresses this directly: the
SAME architecture and training procedure is trained n_folds separate
times, each time with a DIFFERENT fold held out as the test set, and
the resulting spread of test-set metrics across folds gives an actual
estimate of how much the reported accuracy would vary if you'd happened
to split the data differently -- not just one more accuracy number, but
a genuine variance estimate around it.

This is why it's a SEPARATE script rather than a flag on 2.7.3: the
control flow is fundamentally different (one train/eval call per fold
in a loop, not one call total), and keeping the two separate means an
ordinary single-split run (2.7.3, for quick iteration) is unaffected by
this script's own choices, and vice versa.

FOLD STRUCTURE PER ITERATION
---------------------------------
For fold k (0..n_folds-1): fold k itself is the TEST set; fold
(k+1) % n_folds is held out as the VALIDATION set (for early stopping,
same as 2.7.3's own use of a validation split); every OTHER fold is
TRAINING data. This means every fold serves as the test set exactly
once, and no row is EVER in more than one of {train, val, test} within
a single fold's own iteration -- 2.7.2's own fold assignment already
guarantees no duplicate group crosses a fold boundary, so nothing
further to enforce here at the row level.

REUSES 2.7.3's OWN MODEL/DATASET/TRAINING CODE
---------------------------------------------------
BoundaryMotifCNN, FlankDataset, run_epoch, tier_accuracy, revcomp,
one_hot_encode are all imported directly from
2.7.3_train_evaluate_model.py (via importlib, since its filename isn't
a valid Python module name for a plain `import` statement -- same
mechanism 2.6.2 already uses for consensus_building_evaluation.py) --
NOT reimplemented here. This script needs torch regardless (it trains
real models, unlike the plotting scripts), so there's no reason to
avoid the import the way 2.7.5's architecture-diagram plot specifically
had to.

COST WARNING -- READ THIS BEFORE RUNNING ON THE REAL DATASET
------------------------------------------------------------------
This trains n_folds FULL models, one after another, each essentially a
complete 2.7.3 run. On your real ~14,000-row dataset, a single 2.7.3
run has taken on the order of an hour (60-75 minutes across your two
real logs) on CPU. n_folds=5 here means ROUGHLY FIVE TIMES THAT --
several hours of sequential wall-clock time, not minutes. Concretely
worth doing before committing a cluster job to this: reduce n_folds
(3 is a defensible, much cheaper alternative to 5, and still gives a
real variance estimate rather than a single point), and/or reduce
max_epochs / lower early_stopping_patience for the CV run specifically
(a rougher per-fold estimate is a reasonable trade for getting a
variance estimate across folds at all, rather than one very precise
single-split number). These are separate config keys from 2.7.3's own
(see below) specifically so tightening them for CV doesn't change your
existing single-split results.

Run:
    python 2.7.4_kfold_cross_validation.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import importlib.util as _ilu
import statistics
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from torch.utils.data import DataLoader

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit

SECTION_KEY = "2.7.4"

PALETTE = ["#0072B2", "#E69F00", "#009E73", "#D55E00", "#CC79A7", "#56B4E9", "#F0E442", "#999999"]


def _load_2_7_3_module(scripts_dir: Path):
    spec = _ilu.spec_from_file_location("m273_cv", scripts_dir / "2.7.3_train_evaluate_model.py")
    m = _ilu.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def train_one_fold(
    m273, fold: int, train_rows: List[dict], val_rows: List[dict], test_rows: List[dict],
    rec_type_to_idx: Dict[str, int], seq_length: int, batch_size: int, max_epochs: int,
    patience: int, learning_rate: float, seed: int, model_dir: Path, device: str, log,
) -> Optional[dict]:
    """Trains and evaluates ONE fold, exactly mirroring 2.7.3's own
    training loop (same model class, same optimizer/scheduler choice, same
    early-stopping-on-val-loss checkpointing) but scoped to this fold's
    own train/val/test rows. Returns per-rec_type test metrics for this
    fold, or None if this fold has no training or validation rows at all
    (possible only in a pathological case -- extremely few duplicate
    groups for every rec_type combined -- logged clearly if it happens)."""
    if not train_rows or not val_rows:
        log.warning(f"  Fold {fold}: empty train or val set -- skipping this fold entirely.")
        return None

    torch.manual_seed(seed + fold)  # different but reproducible per fold

    train_loader = DataLoader(m273.FlankDataset(train_rows, rec_type_to_idx, seq_length),
                               batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(m273.FlankDataset(val_rows, rec_type_to_idx, seq_length),
                             batch_size=batch_size, shuffle=False)

    model = m273.BoundaryMotifCNN(n_rec_types=len(rec_type_to_idx)).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=3)

    fold_dir = model_dir / f"cv_fold_{fold}"
    fold_dir.mkdir(parents=True, exist_ok=True)
    best_val_loss = float("inf")
    epochs_without_improvement = 0
    train_log_rows = []

    for epoch in range(1, max_epochs + 1):
        train_loss = m273.run_epoch(model, train_loader, optimizer, device, train=True)
        val_loss = m273.run_epoch(model, val_loader, optimizer, device, train=False)
        scheduler.step(val_loss)
        train_log_rows.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss,
                                "lr": optimizer.param_groups[0]["lr"]})
        if val_loss < best_val_loss - 1e-6:
            best_val_loss = val_loss
            epochs_without_improvement = 0
            torch.save(model.state_dict(), fold_dir / "best_model.pt")
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= patience:
                break
        log.info(f"  Fold {fold} epoch {epoch:3d}: train_loss={train_loss:.5f}  val_loss={val_loss:.5f}")

    with open(fold_dir / "train_log.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["epoch", "train_loss", "val_loss", "lr"],
                                 delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(train_log_rows)
    log.info(f"  Fold {fold}: best val_loss={best_val_loss:.5f} ({len(train_log_rows)} epochs run)")

    if not test_rows:
        log.warning(f"  Fold {fold}: no test rows -- no test metrics for this fold.")
        return {"fold": fold, "best_val_loss": best_val_loss, "per_rec_type": {}}

    model.load_state_dict(torch.load(fold_dir / "best_model.pt", map_location=device))
    model.eval()
    preds = m273.evaluate_with_breakdown(model, test_rows, rec_type_to_idx, seq_length, device, batch_size)
    labels = [float(r["label"]) for r in test_rows]

    per_rec_type = {}
    all_rec_types_in_test = sorted({r["rec_type"] for r in test_rows})
    for rt in all_rec_types_in_test:
        idxs = [i for i, r in enumerate(test_rows) if r["rec_type"] == rt]
        rt_preds = [preds[i] for i in idxs]
        rt_labels = [labels[i] for i in idxs]
        mse = float(np.mean([(p - y) ** 2 for p, y in zip(rt_preds, rt_labels)]))
        mae = float(np.mean([abs(p - y) for p, y in zip(rt_preds, rt_labels)]))
        acc = m273.tier_accuracy(rt_preds, rt_labels)
        per_rec_type[rt] = {"n": len(idxs), "mse": mse, "mae": mae, "tier_accuracy": acc}

    overall_mse = float(np.mean([(p - y) ** 2 for p, y in zip(preds, labels)]))
    overall_mae = float(np.mean([abs(p - y) for p, y in zip(preds, labels)]))
    overall_acc = m273.tier_accuracy(preds, labels)
    per_rec_type["__overall__"] = {"n": len(test_rows), "mse": overall_mse, "mae": overall_mae,
                                     "tier_accuracy": overall_acc}

    return {"fold": fold, "best_val_loss": best_val_loss, "per_rec_type": per_rec_type}


def plot_cv_results(fold_results: List[dict], all_rec_types: List[str], out_dir: Path, log) -> None:
    """One panel per rec_type (plus one for the pooled overall), each
    showing that rec_type's test tier_accuracy across every fold it had
    test data in, as a strip plot -- the spread of points IS the result:
    how much accuracy would have varied under a different split."""
    labels_to_plot = ["__overall__"] + [rt for rt in all_rec_types]
    data_by_label: Dict[str, List[float]] = {lbl: [] for lbl in labels_to_plot}
    for fr in fold_results:
        if fr is None:
            continue
        for lbl, stats in fr["per_rec_type"].items():
            if lbl in data_by_label:
                data_by_label[lbl].append(stats["tier_accuracy"])

    labels_with_data = [lbl for lbl in labels_to_plot if data_by_label[lbl]]
    if not labels_with_data:
        log.info("  No fold results with test metrics -- skipping CV results plot.")
        return

    fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(labels_with_data)), 5))
    for i, lbl in enumerate(labels_with_data):
        ys = data_by_label[lbl]
        color = "#333333" if lbl == "__overall__" else PALETTE[i % len(PALETTE)]
        jitter = [i + (hash((lbl, j)) % 1000 / 1000 - 0.5) * 0.3 for j in range(len(ys))]
        ax.scatter(jitter, ys, s=40, color=color, alpha=0.75, zorder=3,
                   label=f"n={len(ys)} fold(s)" if i == 0 else None)
        if len(ys) > 1:
            mean_y = statistics.mean(ys)
            ax.hlines(mean_y, i - 0.25, i + 0.25, color=color, linewidth=2, zorder=2)

    ax.set_xticks(range(len(labels_with_data)))
    ax.set_xticklabels(["overall" if lbl == "__overall__" else lbl for lbl in labels_with_data],
                        rotation=30, ha="right")
    ax.set_ylabel("Test tier accuracy")
    ax.set_ylim(0, 1.05)
    ax.set_title(f"Cross-validated test tier accuracy ({len(fold_results)} fold(s))\n"
                 f"(each point = one fold's held-out test accuracy; bar = mean across folds)")
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "ml_cv_results.png", dpi=300, bbox_inches="tight")
    fig.savefig(out_dir / "ml_cv_results.pdf", bbox_inches="tight")
    plt.close(fig)
    log.info(f"-> {out_dir / 'ml_cv_results.png'} / .pdf")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--labels-split", default=None, help="Override BASE/files/ml_labels_split.tsv")
    ap.add_argument("--model-dir", default=None, help="Override BASE/ml_model")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_7"][SECTION_KEY]

    log = get_logger("2.7.4_kfold_cross_validation", paths["logs_dir"])

    labels_split_path = Path(args.labels_split) if args.labels_split else files_dir / "ml_labels_split.tsv"
    model_dir = Path(args.model_dir) if args.model_dir else paths["base_dir"] / "ml_model"
    results_out = files_dir / "ml_cv_results.tsv"
    summary_out = files_dir / "ml_cv_summary.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — K-fold cross-validation")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, results_out, log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_7_4_kfold_cross_validation"):
        return

    if not labels_split_path.exists():
        log.error(f"{labels_split_path} not found -- run 2.7.2_split_dataset.py first "
                  f"(needs its 'fold' column).")
        raise SystemExit(1)

    seq_length = int(section_cfg.get("seq_length", 500))  # was 530 -- 2.7.1 now excludes the 30nt mge_overlap zone
    batch_size = int(section_cfg.get("batch_size", 64))
    max_epochs = int(section_cfg.get("max_epochs", 100))
    patience = int(section_cfg.get("early_stopping_patience", 10))
    learning_rate = float(section_cfg.get("learning_rate", 1e-3))
    seed = int(section_cfg.get("seed", 42))
    n_folds = int(section_cfg.get("n_folds", 5))

    device = "cuda" if torch.cuda.is_available() else "cpu"
    log.info(f"Device: {device}")
    log.info(f"n_folds={n_folds}  seq_length={seq_length}  batch_size={batch_size}  "
             f"max_epochs={max_epochs}  patience={patience}  lr={learning_rate}  seed={seed}")
    log.info("COST NOTE: this runs a full training loop per fold -- see module docstring "
             "for a realistic time estimate before running this on the full dataset.")
    log.info("-" * 70)

    m273 = _load_2_7_3_module(Path(__file__).resolve().parent)

    t0 = time.time()
    with open(labels_split_path, newline="") as fh:
        all_rows = list(csv.DictReader(fh, delimiter="\t"))
    if not all_rows or "fold" not in all_rows[0]:
        log.error(f"{labels_split_path} has no 'fold' column -- re-run the updated "
                  f"2.7.2_split_dataset.py first.")
        raise SystemExit(1)

    rec_types = sorted({r["rec_type"] for r in all_rows})
    rec_type_to_idx = {rt: i for i, rt in enumerate(rec_types)}
    log.info(f"Loaded {len(all_rows):,} rows, {n_folds} folds, rec_types: {rec_types}")

    fold_results: List[Optional[dict]] = []
    for k in range(n_folds):
        val_fold = (k + 1) % n_folds
        test_rows = [r for r in all_rows if int(r["fold"]) == k]
        val_rows = [r for r in all_rows if int(r["fold"]) == val_fold]
        train_rows = [r for r in all_rows if int(r["fold"]) not in (k, val_fold)]
        log.info("-" * 70)
        log.info(f"FOLD {k}/{n_folds - 1}: train={len(train_rows):,}  val={len(val_rows):,} "
                 f"(fold {val_fold})  test={len(test_rows):,} (fold {k})")
        result = train_one_fold(
            m273, k, train_rows, val_rows, test_rows, rec_type_to_idx, seq_length,
            batch_size, max_epochs, patience, learning_rate, seed, model_dir, device, log,
        )
        fold_results.append(result)

    # ── write long-format per-fold, per-rec_type results ────────────────────
    long_rows = []
    for fr in fold_results:
        if fr is None:
            continue
        for rt, stats in fr["per_rec_type"].items():
            long_rows.append({
                "fold": fr["fold"], "rec_type": rt, "n": stats["n"],
                "mse": round(stats["mse"], 5), "mae": round(stats["mae"], 5),
                "tier_accuracy": round(stats["tier_accuracy"], 4),
            })
    with open(results_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["fold", "rec_type", "n", "mse", "mae", "tier_accuracy"],
                                 delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(long_rows)
    log.info(f"-> {results_out}")

    # ── aggregate mean +/- std across folds, overall and per rec_type ───────
    by_label: Dict[str, List[dict]] = {}
    for row in long_rows:
        by_label.setdefault(row["rec_type"], []).append(row)

    summary_rows = []
    for label in ["__overall__"] + rec_types:
        entries = by_label.get(label, [])
        if not entries:
            continue
        accs = [e["tier_accuracy"] for e in entries]
        mses = [e["mse"] for e in entries]
        summary_rows.append({
            "rec_type": "overall" if label == "__overall__" else label,
            "n_folds_with_data": len(entries),
            "mean_tier_accuracy": round(statistics.mean(accs), 4),
            "std_tier_accuracy": round(statistics.stdev(accs), 4) if len(accs) > 1 else "NA",
            "mean_mse": round(statistics.mean(mses), 5),
            "std_mse": round(statistics.stdev(mses), 5) if len(mses) > 1 else "NA",
        })
    with open(summary_out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["rec_type", "n_folds_with_data", "mean_tier_accuracy",
                                                  "std_tier_accuracy", "mean_mse", "std_mse"],
                                 delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)

    log.info("-" * 70)
    log.info("CROSS-VALIDATION SUMMARY")
    log.info("-" * 70)
    for row in summary_rows:
        log.info(f"  {row['rec_type']:20s} folds={row['n_folds_with_data']}  "
                 f"tier_accuracy={row['mean_tier_accuracy']:.3f} +/- {row['std_tier_accuracy']}  "
                 f"mse={row['mean_mse']:.5f} +/- {row['std_mse']}")
    log.info(f"-> {summary_out}")

    plots_dir = paths.get("plots_dir", paths["base_dir"] / "plots")
    plot_cv_results(fold_results, rec_types, plots_dir, log)

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()