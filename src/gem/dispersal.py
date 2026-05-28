"""
Dispersal functions for the GEM working group.

Adapted from the R prototypes in experiments/dispersal_experiments_design/functions/.

Dispersal mechanism:
- Density-dependent (diffuse_density_dependent): per-capita emigration rate is a
  sigmoid of (metabolism - net_growth), following Ryser et al. 2021 Eq. 10.

Assumes reflective boundary conditions (no biomass leaves the grid)
"""

import numpy as np
from numpy.typing import NDArray


def enforce_boundary_conditions(n_rows: int, n_cols: int) -> NDArray[np.int64]:
    """
    Return an (n_rows, n_cols, 1) integer array of neighbour counts.

    Interior cells = 4, edge cells = 3, corner cells = 2. Used to scale
    emigration losses at the grid boundary (reflective boundary conditions:
    no biomass flux leaves the domain). The trailing size-1 axis lets this
    array broadcast against a (n_rows, n_cols, S) biomass array.
    """
    m = np.full((n_rows, n_cols, 1), 4, dtype=np.int64)
    m[0, :, :]  -= 1    # top row: no northern neighbour
    m[-1, :, :] -= 1    # bottom row: no southern neighbour
    m[:, 0, :]  -= 1    # left col: no western neighbour
    m[:, -1, :] -= 1    # right col: no eastern neighbour
    return m


def disperse_delta(
    biomass: NDArray[np.float64],
    net_growth: NDArray[np.float64],
    metabolism: NDArray[np.float64],
    max_disp_rate: NDArray[np.float64],
    boundary_number: NDArray[np.int64],
    b: float = 10.0,
    dt: float = 1.0,
) -> NDArray[np.float64]:
    """
    One density-dependent diffusion step (Ryser et al. 2021, Eq. 10).

    Per-capita emigration rate is a sigmoid of (metabolism - net_growth):
        d = max_disp_rate / (1 + exp(-b * (metabolism - net_growth)))
    High when conditions are bad (net_growth < metabolism), low when good.
    Survival during dispersal is set to 1 (no matrix mortality).

    Parameters
    ----------
    biomass : (n_rows, n_cols, S) array
        Current biomass densities for S species across all cells.
    net_growth : (n_rows, n_cols, S) array
        Per-capita net growth rate nu_{i,z}: (feeding - losses - metabolism)
        / biomass. Computed from ATN state before calling this function.
        Same shape as biomass.
    metabolism : (n_rows, n_cols, S) array
        Per-cell, per-species metabolic rate x_i; the sigmoid inflection point
        (dispersal switches around this value). Same shape as biomass.
    max_disp_rate : (n_rows, n_cols, S) array
        Maximum per-capita dispersal rate (parameter a in Ryser et al. 2021).
        Species-specific in storage (shape (S,)); the caller broadcasts it up to
        the full grid shape before calling. Same shape as biomass.
    boundary_number : (n_rows, n_cols, 1) int array
        Neighbour counts from enforce_boundary_conditions(). Broadcasts over S.
    b : float
        Sigmoid steepness (default 10, as in Ryser et al. 2021).

    Returns
    -------
    (n_rows, n_cols, S) array
        Net biomass change (immigration - emigration) for this step.
        Add to the current biomass to advance the state.
    """

    # Ensure biomass, net_growth, metabolism and max_disp_rate share the same shape
    assert biomass.shape == net_growth.shape == metabolism.shape == max_disp_rate.shape, (
        f"Shape mismatch: biomass {biomass.shape}, net_growth {net_growth.shape}, "
        f"metabolism {metabolism.shape}, max_disp_rate {max_disp_rate.shape}"
    )

    # Ensure that boundary conditions enforcer has the same spatial dimensions as biomass
    assert biomass.shape[:2] == boundary_number.shape[:2], (
        f"Spatial shape mismatch: biomass {biomass.shape[:2]}, "
        f"boundary_number {boundary_number.shape[:2]}"
    )

    d_grid = max_disp_rate / (1 + np.exp(-b * (metabolism - net_growth)))

    flux_per_edge = biomass * (d_grid / 4)
    flux_out = flux_per_edge * boundary_number

    n_rows, n_cols, n_species = biomass.shape
    zero_row = np.zeros((1, n_cols, n_species))
    zero_col = np.zeros((n_rows, 1, n_species))

    from_north = np.concatenate([zero_row, flux_per_edge[:-1, :, :]], axis=0)
    from_south = np.concatenate([flux_per_edge[1:, :, :], zero_row],  axis=0)
    from_west  = np.concatenate([zero_col, flux_per_edge[:, :-1, :]], axis=1)
    from_east  = np.concatenate([flux_per_edge[:, 1:, :], zero_col],  axis=1)

    return (-flux_out + from_north + from_south + from_west + from_east) * dt


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
        biomass = biomass + disperse_delta(
            biomass, net_growth, metabolism, max_disp_grid, boundary_number, b
        )
        results.append(biomass.copy())

    return results
