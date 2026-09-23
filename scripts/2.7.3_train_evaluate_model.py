#!/usr/bin/env python3
"""
2.7.3_train_evaluate_model.py  —  Section 2.7.3
===================================================
input : BASE/files/ml_labels_split.tsv (2.7.2)
output: BASE/ml_model/best_model.pt              (trained weights)
        BASE/ml_model/train_log.tsv              (per-epoch train/val loss)
        BASE/files/ml_test_predictions.tsv       (per-row test-set predictions)
        BASE/ml_model/attributions/<rec_type>_<mge_id>_<side>.tsv
                                                  (per-base Integrated-Gradients
                                                  attribution for a handful of
                                                  sampled test examples per rec_type)

WHAT THIS MODEL ACTUALLY PREDICTS, AND WHY (see chat)
----------------------------------------------------------
This is NOT a base-pair-resolution signal-regression model like the
BPNet notebook it's adapted from -- deliberately. That notebook's target
(ChIP-seq read counts) is a real per-base MEASURED quantity; the label
2.7.1 built is a single confidence SCORE per (MGE, side) -- 0.3 / 0.7 /
1.0 -- describing the whole flank, not a curve across it. Forcing that
into a per-base target would mean inventing structure the label doesn't
actually contain. So the model predicts ONE scalar per flank (global
average pooling before the regression head, not a per-position "deconv"
head), trained with MSE toward that label. Per-base importance is still
fully available afterward -- see ATTRIBUTION below -- it just comes from
interpreting the trained scalar model, exactly the way the BPNet
notebook's own attribution step does (it also collapses its per-base
output to a scalar, via a summed position mask, before running
Integrated Gradients).

INPUT ENCODING
------------------
Each flank is one-hot encoded (4 channels: A/C/G/T; any other character,
e.g. N, encodes as all-zero). flank_sequence itself, as written by
2.7.1, already excludes the mge_overlap zone (see that script's own
docstring) -- what's encoded here is genuinely flanking sequence only,
with the true MGE boundary sitting exactly at one edge of the sequence:
the END for LEFT flanks, the START for RIGHT flanks (re-derived directly
from 2.5.1_extract_flanks.py's own coordinate arithmetic).

LEFT-side flanks are reverse-complemented before encoding (SWITCHED from
RIGHT, see chat): a genuine TIR pair is a pair of reverse-complements of
each other, so presenting both flanks in a shared orientation means
picking ONE consistent "which end is the boundary" convention and
revcomp'ing whichever side doesn't naturally have it there. Either
choice (revcomp LEFT so both boundaries land at the START, or revcomp
RIGHT so both land at the END -- the previous version of this script)
is internally self-consistent; this is genuinely a convention choice,
not a correctness fix. Revcomp'ing LEFT means both flanks now present
their boundary-adjacent end at position 0 (the START of the encoded
sequence) rather than the end.

CONDITIONING ON rec_type AND side
--------------------------------------
Both are embedded (small learned lookup tables) and concatenated onto
the pooled sequence features before the final regression head, rather
than training either fully separate per-family models or one
family-blind model. Reasoning (see chat): rec_types range from ~50
duplicate groups to a handful in this dataset, too few for reliable
separate models across the board, but genuinely different in what a
"boundary motif" looks like (see the literature Ends patterns 2.7.1
matches against) -- conditioning lets one shared trunk pool statistical
strength across families while the embedding gives it a way to know
which family's expectation applies, without hand-building 7 separate
architectures.

MODEL
---------
A small 1D CNN: a stack of Conv1d -> BatchNorm -> ReLU blocks with
increasing dilation (1, 2, 4, 8) to build up receptive field over the
flank without needing pooling-induced resolution loss, global average
pooling over the sequence axis, concatenated with the rec_type/side
embeddings, then a 2-layer MLP head to one scalar output with a sigmoid
(labels live in [0.3, 1.0], sigmoid's [0,1] range comfortably covers
that). Exact layer sizes are in ml_model_config() below, in one place,
on purpose -- these are reasonable defaults for a first pass, not a
tuned result, and are the first thing worth adjusting if training
doesn't behave well on the real dataset. See DEFAULT_CONFIG for exact
values.

TRAINING
------------
MSE loss, Adam optimizer, ReduceLROnPlateau on validation loss, early
stopping (patience configurable). Batches are built directly from
ml_labels_split.tsv's own train/val/test split column -- 2.7.2 already
guaranteed no duplicate group crosses a split boundary, so nothing
further to enforce here.

EVALUATION
--------------
On the held-out test split: overall MSE/MAE, MSE/MAE broken down PER
rec_type (this is where you'd actually see if a family like
DDE_Tnp_IS1595 -- which your own reference table shows has several
visibly different sub-group Ends patterns -- is being under-served by
the shared trunk), and a coarser "tier accuracy" (each prediction
rounded to whichever of {0.3, 0.7, 1.0} it's closest to, compared
against the true label) as an easier number to explain to someone who
isn't going to read an MSE value.

ATTRIBUTION -- "which nucleotides matter"
----------------------------------------------
Integrated Gradients (Captum, same tool the BPNet notebook uses) against
a uniform (0.25, 0.25, 0.25, 0.25) baseline, run on a handful of sampled
test-set examples per rec_type (--n-attribution-examples). Output is one
TSV per example: one row per input position, with the per-base
attribution score for each of the 4 channels -- directly usable for the
"which nucleotides are important for which rec_type" comparison. Also
reports, for each sampled example, whether the position of PEAK
attribution actually falls inside the segment 2.7.1 originally found for
that MGE (segment_start/segment_length columns already in
ml_labels_split.tsv) -- a cheap, immediate sanity check on whether the
model's learned notion of "important" agrees with the hand-built
consensus-scan that generated its own training labels, before you trust
it on cases where that scan found nothing at all.

Run standalone:
    python 2.7.3_train_evaluate_model.py --config path/to/config.yaml
"""

from __future__ import annotations

import argparse
import csv
import random
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

from utils import get_logger, get_paths, load_config, resolve_placeholder, maybe_submit_and_exit
from ml_model_config import ml_model_config

SECTION_KEY = "2.7.3"

BASES = "ACGT"
BASE_TO_IDX = {b: i for i, b in enumerate(BASES)}


def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTN", "TGCAN")
    return seq.translate(comp)[::-1]


def one_hot_encode(seq: str, length: int) -> np.ndarray:
    """(4, length) float32 array. Longer sequences truncated, shorter
    ones right-padded with all-zero (N-like) columns. Unknown characters
    (N, ambiguity codes) also encode as all-zero."""
    arr = np.zeros((4, length), dtype=np.float32)
    seq = seq[:length]
    for i, ch in enumerate(seq):
        idx = BASE_TO_IDX.get(ch)
        if idx is not None:
            arr[idx, i] = 1.0
    return arr


# ══════════════════════════════════════════════════════════════════════════
#  DATASET
# ══════════════════════════════════════════════════════════════════════════

class FlankDataset(Dataset):
    def __init__(self, rows: List[dict], rec_type_to_idx: Dict[str, int], seq_length: int):
        self.rows = rows
        self.rec_type_to_idx = rec_type_to_idx
        self.seq_length = seq_length

    def __len__(self) -> int:
        return len(self.rows)

    def __getitem__(self, i: int):
        row = self.rows[i]
        seq = row["flank_sequence"].upper()
        if row["side"] == "left":
            seq = revcomp(seq)  # see module docstring: consistent boundary-region orientation both sides
        x_seq = one_hot_encode(seq, self.seq_length)
        rec_type_idx = self.rec_type_to_idx[row["rec_type"]]
        side_idx = 0 if row["side"] == "left" else 1
        y = float(row["label"])
        return (
            torch.tensor(x_seq, dtype=torch.float32),
            torch.tensor(rec_type_idx, dtype=torch.long),
            torch.tensor(side_idx, dtype=torch.long),
            torch.tensor(y, dtype=torch.float32),
        )


# ══════════════════════════════════════════════════════════════════════════
#  MODEL
# ══════════════════════════════════════════════════════════════════════════

# ml_model_config() itself now lives in ml_model_config.py, not here (see
# chat / that file's own docstring) -- so 2.7.5_plot_ml_results.py can
# import just the hyperparameters, for its architecture diagram, without
# needing torch installed at all.


class BoundaryMotifCNN(nn.Module):
    def __init__(self, n_rec_types: int, cfg: Optional[dict] = None):
        super().__init__()
        cfg = cfg or ml_model_config()
        channels = cfg["conv_channels"]
        k = cfg["kernel_size"]
        dilations = cfg["dilations"]
        assert len(channels) == len(dilations), "conv_channels and dilations must be the same length"

        blocks = []
        in_ch = 4
        for out_ch, dil in zip(channels, dilations):
            pad = (k - 1) * dil // 2  # 'same'-style padding so length is preserved through the stack
            blocks.append(nn.Conv1d(in_ch, out_ch, kernel_size=k, dilation=dil, padding=pad))
            blocks.append(nn.BatchNorm1d(out_ch))
            blocks.append(nn.ReLU())
            in_ch = out_ch
        self.conv_stack = nn.Sequential(*blocks)
        self.global_pool = nn.AdaptiveAvgPool1d(1)

        self.rec_type_embed = nn.Embedding(n_rec_types, cfg["rec_type_embed_dim"])
        self.side_embed = nn.Embedding(2, cfg["side_embed_dim"])

        head_in = channels[-1] + cfg["rec_type_embed_dim"] + cfg["side_embed_dim"]
        self.head = nn.Sequential(
            nn.Linear(head_in, cfg["head_hidden_dim"]),
            nn.ReLU(),
            nn.Dropout(cfg["dropout"]),
            nn.Linear(cfg["head_hidden_dim"], 1),
        )

    def forward(self, x_seq: torch.Tensor, rec_type_idx: torch.Tensor, side_idx: torch.Tensor) -> torch.Tensor:
        h = self.conv_stack(x_seq)                 # (B, C, L)
        h = self.global_pool(h).squeeze(-1)         # (B, C)
        rt = self.rec_type_embed(rec_type_idx)      # (B, rec_type_embed_dim)
        sd = self.side_embed(side_idx)              # (B, side_embed_dim)
        combined = torch.cat([h, rt, sd], dim=1)
        out = self.head(combined).squeeze(-1)       # (B,)
        return torch.sigmoid(out)


# ══════════════════════════════════════════════════════════════════════════
#  TRAIN / EVAL LOOPS
# ══════════════════════════════════════════════════════════════════════════

def run_epoch(model, loader, optimizer, device, train: bool) -> float:
    model.train(mode=train)
    total_loss, n = 0.0, 0
    loss_fn = nn.MSELoss(reduction="sum")
    for x_seq, rt_idx, side_idx, y in loader:
        x_seq, rt_idx, side_idx, y = x_seq.to(device), rt_idx.to(device), side_idx.to(device), y.to(device)
        if train:
            optimizer.zero_grad()
        pred = model(x_seq, rt_idx, side_idx)
        loss = loss_fn(pred, y)
        if train:
            loss.backward()
            optimizer.step()
        total_loss += loss.item()
        n += y.shape[0]
    return total_loss / max(n, 1)


def evaluate_with_breakdown(model, rows, rec_type_to_idx, seq_length, device, batch_size):
    ds = FlankDataset(rows, rec_type_to_idx, seq_length)
    loader = DataLoader(ds, batch_size=batch_size, shuffle=False)
    model.eval()
    preds = []
    with torch.no_grad():
        for x_seq, rt_idx, side_idx, y in loader:
            x_seq, rt_idx, side_idx = x_seq.to(device), rt_idx.to(device), side_idx.to(device)
            pred = model(x_seq, rt_idx, side_idx).cpu().numpy()
            preds.extend(pred.tolist())
    return preds


def tier_accuracy(preds: List[float], labels: List[float], tiers=(0.3, 0.7, 1.0)) -> float:
    correct = 0
    for p, y in zip(preds, labels):
        nearest_pred = min(tiers, key=lambda t: abs(t - p))
        nearest_true = min(tiers, key=lambda t: abs(t - y))
        if nearest_pred == nearest_true:
            correct += 1
    return correct / max(len(preds), 1)


# ══════════════════════════════════════════════════════════════════════════
#  ATTRIBUTION
# ══════════════════════════════════════════════════════════════════════════

def compute_attributions(model, rows, rec_type_to_idx, seq_length, device, n_per_rec_type, seed, log):
    try:
        from captum.attr import IntegratedGradients
    except ImportError:
        log.info("captum not installed -- skipping attribution step (pip install captum to enable it).")
        return []

    rng = random.Random(seed)
    by_rec_type: Dict[str, List[dict]] = {}
    for row in rows:
        by_rec_type.setdefault(row["rec_type"], []).append(row)

    sampled = []
    for rec_type, group_rows in sorted(by_rec_type.items()):
        k = min(n_per_rec_type, len(group_rows))
        sampled.extend(rng.sample(group_rows, k))

    model.eval()

    class ScalarWrapper(nn.Module):
        def __init__(self, model, rec_type_idx, side_idx):
            super().__init__()
            self.model = model
            self.rec_type_idx = rec_type_idx
            self.side_idx = side_idx

        def forward(self, x_seq):
            batch = x_seq.shape[0]
            rt = self.rec_type_idx.expand(batch)
            sd = self.side_idx.expand(batch)
            return self.model(x_seq, rt, sd).unsqueeze(-1)

    results = []
    for row in sampled:
        seq = row["flank_sequence"].upper()
        if row["side"] == "left":
            seq = revcomp(seq)
        x = torch.tensor(one_hot_encode(seq, seq_length), dtype=torch.float32).unsqueeze(0).to(device)
        rt_idx = torch.tensor([rec_type_to_idx[row["rec_type"]]], dtype=torch.long).to(device)
        side_idx = torch.tensor([0 if row["side"] == "left" else 1], dtype=torch.long).to(device)

        wrapped = ScalarWrapper(model, rt_idx, side_idx).to(device)
        ig = IntegratedGradients(wrapped)
        baseline = torch.full_like(x, 0.25)
        attributions = ig.attribute(x, baselines=baseline, n_steps=50)
        attributions = attributions.squeeze(0).detach().cpu().numpy()  # (4, seq_length)

        peak_pos = int(np.argmax(np.abs(attributions).sum(axis=0)))
        seg_start, seg_len = row.get("segment_start"), row.get("segment_length")
        peak_in_found_segment = "NA"
        if seg_start not in (None, "NA") and seg_len not in (None, "NA"):
            try:
                s, l = int(seg_start), int(seg_len)
                if row["side"] == "left":
                    # segment_start/segment_length (from 2.7.1) are in the ORIGINAL
                    # left-flank's own coordinates; the attribution was computed on
                    # the REVCOMP'd sequence actually fed to the model (see module
                    # docstring), so the segment's position under a full-sequence
                    # reversal has to be flipped the same way before comparing:
                    # a segment at [s, s+l) in a sequence of length seq_length lands
                    # at [seq_length-(s+l), seq_length-s) after revcomp.
                    s = seq_length - (s + l)
                peak_in_found_segment = "yes" if s <= peak_pos < s + l else "no"
            except ValueError:
                pass

        results.append({
            "rec_type": row["rec_type"], "mge_id": row["mge_id"], "side": row["side"],
            "attributions": attributions, "peak_pos": peak_pos, "peak_in_found_segment": peak_in_found_segment,
        })
    return results


# ══════════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════════

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=None)
    ap.add_argument("--labels-split", default=None, help="Override BASE/files/ml_labels_split.tsv")
    ap.add_argument("--model-dir", default=None, help="Override BASE/ml_model")
    ap.add_argument("--predictions-out", default=None, help="Override BASE/files/ml_test_predictions.tsv")
    args = ap.parse_args()

    cfg = load_config(args.config)
    paths = get_paths(cfg)
    files_dir = paths["files_dir"]
    section_cfg = cfg["section_2_7"][SECTION_KEY]

    log = get_logger("2.7.3_train_evaluate_model", paths["logs_dir"])

    labels_split_path = Path(args.labels_split) if args.labels_split else files_dir / "ml_labels_split.tsv"
    model_dir = Path(args.model_dir) if args.model_dir else paths["base_dir"] / "ml_model"
    predictions_out = Path(args.predictions_out) if args.predictions_out else files_dir / "ml_test_predictions.tsv"

    log.info("=" * 70)
    log.info(f"Section {SECTION_KEY} — Train and evaluate the boundary-motif CNN")
    log.info("=" * 70)

    if resolve_placeholder(section_cfg, model_dir / "best_model.pt", log):
        log.info(f"Section {SECTION_KEY} complete (via placeholder).")
        return

    if maybe_submit_and_exit(cfg, __file__, sys.argv[1:], log, job_name="2_7_3_train_evaluate_model"):
        return

    if not labels_split_path.exists():
        log.error(f"{labels_split_path} not found -- run 2.7.2_split_dataset.py first.")
        raise SystemExit(1)

    seq_length = int(section_cfg.get("seq_length", 500))  # was 530 -- 2.7.1 now excludes the 30nt mge_overlap zone; set to match flank_size - mge_overlap in your own config
    batch_size = int(section_cfg.get("batch_size", 64))
    max_epochs = int(section_cfg.get("max_epochs", 100))
    patience = int(section_cfg.get("early_stopping_patience", 10))
    learning_rate = float(section_cfg.get("learning_rate", 1e-3))
    seed = int(section_cfg.get("seed", 42))
    n_attribution_examples = int(section_cfg.get("n_attribution_examples_per_rec_type", 5))

    torch.manual_seed(seed)
    random.seed(seed)
    np.random.seed(seed)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    log.info(f"Device: {device}")
    log.info(f"seq_length={seq_length}  batch_size={batch_size}  max_epochs={max_epochs}  "
             f"patience={patience}  lr={learning_rate}  seed={seed}")
    log.info("-" * 70)

    t0 = time.time()

    with open(labels_split_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        all_rows = list(reader)
    log.info(f"Loaded {len(all_rows):,} rows from {labels_split_path}")

    rec_types = sorted({r["rec_type"] for r in all_rows})
    rec_type_to_idx = {rt: i for i, rt in enumerate(rec_types)}
    log.info(f"rec_types ({len(rec_types)}): {rec_types}")

    train_rows = [r for r in all_rows if r["split"] == "train"]
    val_rows = [r for r in all_rows if r["split"] == "val"]
    test_rows = [r for r in all_rows if r["split"] == "test"]
    log.info(f"train={len(train_rows):,}  val={len(val_rows):,}  test={len(test_rows):,}")
    if not train_rows or not val_rows:
        log.error("Empty train or val split -- check 2.7.2's output / rec_types filter.")
        raise SystemExit(1)

    train_loader = DataLoader(FlankDataset(train_rows, rec_type_to_idx, seq_length),
                               batch_size=batch_size, shuffle=True)
    val_loader = DataLoader(FlankDataset(val_rows, rec_type_to_idx, seq_length),
                             batch_size=batch_size, shuffle=False)

    model = BoundaryMotifCNN(n_rec_types=len(rec_types)).to(device)
    n_params = sum(p.numel() for p in model.parameters())
    log.info(f"Model parameters: {n_params:,}")

    optimizer = torch.optim.Adam(model.parameters(), lr=learning_rate)
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(optimizer, mode="min", factor=0.5, patience=3)

    model_dir.mkdir(parents=True, exist_ok=True)
    best_val_loss = float("inf")
    epochs_without_improvement = 0
    train_log_rows = []

    log.info("-" * 70)
    log.info("Training ...")
    for epoch in range(1, max_epochs + 1):
        train_loss = run_epoch(model, train_loader, optimizer, device, train=True)
        val_loss = run_epoch(model, val_loader, optimizer, device, train=False)
        scheduler.step(val_loss)
        current_lr = optimizer.param_groups[0]["lr"]
        log.info(f"  epoch {epoch:3d}: train_loss={train_loss:.5f}  val_loss={val_loss:.5f}  lr={current_lr:.2e}")
        train_log_rows.append({"epoch": epoch, "train_loss": train_loss, "val_loss": val_loss, "lr": current_lr})

        if val_loss < best_val_loss - 1e-6:
            best_val_loss = val_loss
            epochs_without_improvement = 0
            torch.save(model.state_dict(), model_dir / "best_model.pt")
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= patience:
                log.info(f"  Early stopping at epoch {epoch} (no val improvement for {patience} epochs).")
                break

    with open(model_dir / "train_log.tsv", "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["epoch", "train_loss", "val_loss", "lr"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(train_log_rows)

    log.info(f"Best val_loss: {best_val_loss:.5f}  -> {model_dir / 'best_model.pt'}")
    log.info("-" * 70)

    # ── evaluate on held-out test set with the BEST checkpoint ──────────────
    model.load_state_dict(torch.load(model_dir / "best_model.pt", map_location=device))
    model.eval()

    if test_rows:
        preds = evaluate_with_breakdown(model, test_rows, rec_type_to_idx, seq_length, device, batch_size)
        labels = [float(r["label"]) for r in test_rows]

        mse = float(np.mean([(p - y) ** 2 for p, y in zip(preds, labels)]))
        mae = float(np.mean([abs(p - y) for p, y in zip(preds, labels)]))
        tier_acc = tier_accuracy(preds, labels)

        log.info("TEST SET EVALUATION")
        log.info(f"  overall: n={len(test_rows):,}  MSE={mse:.5f}  MAE={mae:.5f}  tier_accuracy={tier_acc:.3f}")
        for rt in rec_types:
            idxs = [i for i, r in enumerate(test_rows) if r["rec_type"] == rt]
            if not idxs:
                log.info(f"  {rt:20s}: (no test examples for this rec_type)")
                continue
            rt_preds = [preds[i] for i in idxs]
            rt_labels = [labels[i] for i in idxs]
            rt_mse = float(np.mean([(p - y) ** 2 for p, y in zip(rt_preds, rt_labels)]))
            rt_mae = float(np.mean([abs(p - y) for p, y in zip(rt_preds, rt_labels)]))
            rt_acc = tier_accuracy(rt_preds, rt_labels)
            log.info(f"  {rt:20s}: n={len(idxs):,}  MSE={rt_mse:.5f}  MAE={rt_mae:.5f}  tier_accuracy={rt_acc:.3f}")

        pred_columns = ["rec_type", "cluster_number", "duplicate_group_number", "mge_id", "side",
                         "label", "prediction"]
        with open(predictions_out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=pred_columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row, pred in zip(test_rows, preds):
                writer.writerow({
                    "rec_type": row["rec_type"], "cluster_number": row["cluster_number"],
                    "duplicate_group_number": row["duplicate_group_number"], "mge_id": row["mge_id"],
                    "side": row["side"], "label": row["label"], "prediction": round(float(pred), 4),
                })
        log.info(f"-> {predictions_out}")
    else:
        log.info("No test rows -- skipping test-set evaluation.")

    # ── attribution ("which nucleotides matter"), sampled from the test set ─
    log.info("-" * 70)
    log.info(f"Computing Integrated Gradients attributions "
             f"({n_attribution_examples} example(s) per rec_type, from the test set) ...")
    attribution_rows = test_rows if test_rows else val_rows
    attr_results = compute_attributions(
        model, attribution_rows, rec_type_to_idx, seq_length, device, n_attribution_examples, seed, log,
    )
    if attr_results:
        attr_dir = model_dir / "attributions"
        attr_dir.mkdir(parents=True, exist_ok=True)
        n_peak_confirmed = 0
        n_peak_checkable = 0
        for r in attr_results:
            out_path = attr_dir / f"{r['rec_type']}_{r['mge_id']}_{r['side']}.tsv"
            with open(out_path, "w", newline="") as fh:
                writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
                writer.writerow(["position", "A", "C", "G", "T"])
                for pos in range(r["attributions"].shape[1]):
                    writer.writerow([pos] + [round(float(v), 6) for v in r["attributions"][:, pos]])
            if r["peak_in_found_segment"] != "NA":
                n_peak_checkable += 1
                if r["peak_in_found_segment"] == "yes":
                    n_peak_confirmed += 1
        log.info(f"  Wrote {len(attr_results):,} attribution file(s) -> {attr_dir}")
        if n_peak_checkable:
            log.info(f"  Peak attribution inside 2.7.1's own found segment: "
                     f"{n_peak_confirmed}/{n_peak_checkable} sampled examples "
                     f"({n_peak_confirmed/n_peak_checkable*100:.0f}%) -- a quick sanity check, "
                     f"not a formal validation.")

    elapsed = time.time() - t0
    log.info("-" * 70)
    log.info(f"Elapsed: {elapsed:.1f}s")
    log.info(f"Section {SECTION_KEY} complete.")


if __name__ == "__main__":
    main()