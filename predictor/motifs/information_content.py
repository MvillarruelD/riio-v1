"""
information_content.py -- the single Shannon information-content helper (Phase 0 dedup).

The identical per-column IC formula `2 + sum_b p_b log2 p_b` (bits, 0..2 for a 4-row PWM) was copied
into `signals/operator_logo._ic` and inline in
`signals/motif_finder.predict`. One implementation now lives here so the IC numbers the Test-1
pre-check and the operator logos report are computed identically everywhere.

Pure / numpy-only. Run `python information_content.py` for a self-test.
"""
from __future__ import annotations

import numpy as np


def per_column(pwm) -> np.ndarray:
    """Per-column information content in bits for a 4xL probability matrix (base axis = 0)."""
    p = np.clip(np.asarray(pwm, dtype=float), 1e-9, 1.0)
    return 2.0 + (p * np.log2(p)).sum(axis=0)


def total(pwm) -> float:
    """Summed information content over all columns (bits)."""
    return float(per_column(pwm).sum())


def mean_per_column(pwm) -> float:
    """Mean per-column IC (the motif-emergence quantity used by operator_logo.select_by_motif)."""
    pc = per_column(pwm)
    return float(pc.mean()) if pc.size else 0.0


# --------------------------------------------------------------------------- self-test
def _demo() -> None:
    uniform = np.full((4, 5), 0.25)
    assert abs(total(uniform)) < 1e-6, "uniform PWM must carry ~0 bits"

    onehot = np.zeros((4, 3)); onehot[0] = 1.0          # all-A columns
    pc = per_column(onehot)
    assert np.allclose(pc, 2.0), f"monomorphic columns must be 2 bits, got {pc}"
    assert abs(total(onehot) - 6.0) < 1e-6 and abs(mean_per_column(onehot) - 2.0) < 1e-6

    # a half-informative column (50/50 between two bases) = 1 bit
    half = np.zeros((4, 1)); half[0] = half[1] = 0.5
    assert abs(per_column(half)[0] - 1.0) < 1e-6, "50/50 column must be 1 bit"
    print("OK: information content (uniform=0, one-hot=2, 50/50=1 bit) verified.")


if __name__ == "__main__":
    _demo()
