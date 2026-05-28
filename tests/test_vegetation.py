"""Unit tests for vegetation science functions (no engine, no grid plumbing)."""

import numpy as np
from gem.vegetation import logistic_growth_delta


def test_logistic_growth_simple():
    """Biomass should grow when below carrying capacity, shrink when above."""
    B = np.array([50.0, 150.0])     # below and above K
    r = np.array([0.1, 0.1])
    K = np.array([100.0, 100.0])
    
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    
    assert delta[0] > 0   # Growing (B < K)
    assert delta[1] < 0   # Shrinking (B > K)


def test_logistic_growth_on_full_grid():
    """End-to-end: broadcast spatial grid with per-species parameters."""
    # Full spatial grid: (X, Y, S)
    B = np.full((180, 360, 15), 50000.0)
    r = np.full((180, 360, 15), 0.1)
    K = np.full((180, 360, 1), 100000.0)  # broadcast across species
    
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    
    assert delta.shape == (180, 360, 15)
    assert np.all(delta > 0)  # All species growing (B < K)

