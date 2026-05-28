"""Vegetation science functions. Pure numpy, no engine imports.

These functions implement the biological mechanisms of vegetation growth.
No grid plumbing, no side effects. Inputs and outputs are explicit arrays.

Easy to test with hand-built data, easy to reason about for ecology contributors,
easy to prototype in a notebook with synthetic data.

All functions work on any shape (S,), (Y, S), or (X, Y, S) through numpy broadcasting.
See docs/processes_implementation_specification.md for the shape contract.
"""

import numpy as np
from numpy.typing import NDArray


def logistic_growth_delta(
    biomass: NDArray[np.float64],
    growth_rate: NDArray[np.float64],
    carrying_capacity: NDArray[np.float64],
    dt: float,
) -> NDArray[np.float64]:
    """Logistic (density-dependent) growth of vegetation biomass.
    
    Implements: dB/dt = r·B·(1 - B/K)
    
    where:
        B = biomass per species (or per location and species)
        r = intrinsic growth rate (species-dependent)
        K = carrying capacity (environment and time-dependent)
    
    Typical ranges:
        r ≈ 0.02–0.3 per day (depends on growth strategy)
        K ≈ 10^4–10^5 kg/ha (depends on environment)
        dt = 1 day (but parameterizable)
    
    Args:
        biomass: Current vegetation biomass. Shape: (S,), (Y, S), or (X, Y, S).
        growth_rate: Intrinsic growth rate r. Must broadcast with biomass.
        carrying_capacity: Environmental capacity K. Must broadcast with biomass.
        dt: Time step in days (default 1.0).
    
    Returns:
        Biomass change dB over dt. Shape broadcasts to match all inputs.
    
    Shape contract (enforced by assertion):
        All three array inputs must be broadcastable together.
        Examples:
          - (X, Y, 15) with (X, Y, 15) with (X, Y, 1)  ✓ broadcasts to (X, Y, 15)
          - (X, Y, 15) with (15,) with (X, Y, 1)       ✓ broadcasts to (X, Y, 15)
          - (X, Y, 15) with (X, Y, 10) with (X, Y, 1)  ✗ 15 != 10, cannot broadcast
    
    Example:
        >>> B = np.array([100.0, 200.0, 50.0])       # 3 species
        >>> r = np.array([0.1, 0.15, 0.2])           # per-species rates
        >>> K = np.array([500.0, 400.0, 300.0])      # per-species carrying capacity
        >>> dt = 1.0
        >>> dB = logistic_growth_delta(B, r, K, dt)
        >>> dB.shape
        (3,)
        >>> dB[0]  # positive because B < K for species 0
        array(7.2)
    """
    # Shape contract: all inputs must broadcast together.
    # This allows (X, Y, 15) + (X, Y, 1) + (15,) etc., but catches mismatches.
    try:
        biomass_bc, rate_bc, capacity_bc = np.broadcast_arrays(
            biomass, growth_rate, carrying_capacity
        )
    except ValueError as e:
        raise ValueError(
            f"Inputs do not broadcast: biomass {biomass.shape}, "
            f"growth_rate {growth_rate.shape}, carrying_capacity {carrying_capacity.shape}. "
            f"Details: {e}"
        )

    # Logistic equation: dB/dt = r * B * (1 - B/K)
    delta_biomass = dt * growth_rate * biomass * (1.0 - biomass / carrying_capacity)

    return delta_biomass
