"""
Dispersal simulation run to test the dispersal process

Adapted from the R prototypes in experiments/dispersal_experiments_design/functions/.

Dispersal mechanism:
- Density-dependent (diffuse_density_dependent): per-capita emigration rate is a
  sigmoid of (metabolism - net_growth), following Ryser et al. 2021 Eq. 10.

Assumes reflective boundary conditions (no biomass leaves)

Depends on functions in src/gem/dispersal.py
"""

import numpy as np
from numpy.typing import NDArray
from ....src.gem.dispersal import enforce_boundary_conditions
from ....src.gem.dispersal import compute_disperse_delta

def run_diffusion_simulation(
    biomass: NDArray[np.float64],
    max_disp_rate: NDArray[np.float64],
    n_time: int,
    net_growth: NDArray[np.float64],
    metabolism: NDArray[np.float64],
    b: float = 10.0,
) -> list[NDArray[np.float64]]:
    """
    Density-dependent diffusion simulation over n_time steps.

    Parameters
    ----------
    biomass : (n_rows, n_cols, S) array
        Initial biomass densities for S species across all cells.
    max_disp_rate : (S,) array
        Maximum per-capita dispersal rate, one value per species. Broadcast to
        the full grid shape once before the loop.
    n_time : int
        Number of time steps to simulate.
    net_growth : (n_rows, n_cols, S) array
        Per-capita net growth rates; held constant across steps when ATN is
        not coupled. Same shape as biomass.
    metabolism : (n_rows, n_cols, S) array
        Per-cell, per-species metabolic rate; sigmoid inflection point. Same
        shape as biomass.
    b : float
        Sigmoid steepness (default 10).

    Returns
    -------
    list of (n_time + 1) arrays
        Biomass snapshots from t = 0 to t = n_time (index 0 is the initial state).
    """
    boundary_number = enforce_boundary_conditions(biomass.shape[0], biomass.shape[1])

    # max_disp_rate is per-species (S,); broadcast it to the full (X, Y, S) grid
    # once so disperse_delta receives equal-shaped arrays (the caller's job, not
    # the science function's).
    max_disp_grid = np.broadcast_to(max_disp_rate, biomass.shape)

    results = [biomass.copy()]
    for _ in range(n_time):
        biomass = biomass + compute_disperse_delta(
            biomass, net_growth, metabolism, max_disp_grid, boundary_number, b
        )
        results.append(biomass.copy())

    return results
