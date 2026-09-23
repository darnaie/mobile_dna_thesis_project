"""
ml_model_config.py
======================
Shared, dependency-light model hyperparameters for Section 2.7's
boundary-motif CNN. Deliberately has ZERO imports beyond the standard
library -- no torch, no captum -- so both 2.7.3_train_evaluate_model.py
(which needs it to actually build the model) and
2.7.5_plot_ml_results.py (which only needs the numbers, to draw the
architecture diagram) can import it without needing torch installed in
whatever environment happens to be running the plotting-only script.

BUG FIX (see chat): ml_model_config() used to live directly inside
2.7.3_train_evaluate_model.py, and 2.7.5 reached it via importlib on the
WHOLE 2.7.3 module -- which meant executing 2.7.3's own top-level
`import torch` too, even though drawing a diagram of these numbers never
actually needed torch. That broke 2.7.5 outright on any machine/
environment where torch isn't installed -- a perfectly reasonable thing
to want, since a plotting-only environment shouldn't need training's
heavy dependencies. Moved here instead; 2.7.3 now imports FROM this
file rather than defining ml_model_config() locally, and 2.7.5 imports
it directly too, with no importlib/torch involved at all.
"""

from __future__ import annotations


def ml_model_config() -> dict:
    """All architecture hyperparameters in one place -- see
    2.7.3_train_evaluate_model.py's module docstring for what each one
    means in context. First thing to adjust if training doesn't behave
    well; also the single source of truth 2.7.5_plot_ml_results.py's
    architecture diagram is generated from, so it can't silently go
    stale if these are tuned later."""
    return {
        "conv_channels": [64, 64, 128, 128],
        "kernel_size": 11,
        "dilations": [1, 2, 4, 8],
        "rec_type_embed_dim": 8,
        "side_embed_dim": 2,
        "head_hidden_dim": 64,
        "dropout": 0.2,
    }