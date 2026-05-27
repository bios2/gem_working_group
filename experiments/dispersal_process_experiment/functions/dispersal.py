"""
Dispersal functions for the GEM working group.

Adapted from the R prototypes in experiments/dispersal_experiments_design/functions/.

Two dispersal mechanisms are implemented:
- Density-independent (diffuse): each cell loses a fixed fraction of its
  biomass to neighbours every step, regardless of local conditions.
- Density-dependent (diffuse_density_dependent): per-capita emigration rate is a
  sigmoid of (metabolism - net_growth), following Ryser et al. 2021 Eq. 10.

Both use a 4-connected (rook) neighbourhood with reflective boundary conditions.
"""

import numpy as np
from numpy.typing import NDArray


def enforce_boundary_conditions(n_rows: int, n_cols: int) -> NDArray[np.int64]:
    """
    Return an (n_rows, n_cols) integer array of neighbour counts.

    Interior cells = 4, edge cells = 3, corner cells = 2. Used to scale
    emigration losses at the grid boundary (reflective boundary conditions:
    no biomass flux leaves the domain).
    """
    m = np.full((n_rows, n_cols), 4, dtype=np.int64)
    m[0, :]  -= 1    # top row: no northern neighbour
    m[-1, :] -= 1    # bottom row: no southern neighbour
    m[:, 0]  -= 1    # left col: no western neighbour
    m[:, -1] -= 1    # right col: no eastern neighbour
    return m


def diffuse(
    biomass: NDArray[np.float64],
    disp_rate: float,
    boundary_number: NDArray[np.int64],
) -> NDArray[np.float64]:
    """
    One density-independent diffusion step.

    Each cell loses (disp_rate / 4) * n_nbrs of its biomass, where n_nbrs is
    the number of neighbours. Biomass disperses equally to all neighbours.
    Boundary cells lose proportionally less (reflective boundary conditions).

    Parameters
    ----------
    biomass : (n_rows, n_cols) array
        Current biomass densities.
    disp_rate : float
        Per-capita dispersal coefficient.
    boundary_number : (n_rows, n_cols) int array
        Neighbour counts from enforce_boundary_conditions(). Same shape as biomass.

    Returns
    -------
    (n_rows, n_cols) array
        Net biomass change (immigration - emigration) for this step.
        Add to the current biomass to advance the state.
    """
    assert biomass.shape == boundary_number.shape, (
        f"biomass shape {biomass.shape} != boundary_number shape {boundary_number.shape}"
    )

    flux_per_edge = biomass * (disp_rate / 4)
    flux_out = flux_per_edge * boundary_number

    n_rows, n_cols = biomass.shape
    zero_row = np.zeros((1, n_cols))
    zero_col = np.zeros((n_rows, 1))

    from_north = np.vstack([zero_row, flux_per_edge[:-1, :]])   # row above sends down
    from_south = np.vstack([flux_per_edge[1:, :], zero_row])    # row below sends up
    from_west  = np.hstack([zero_col, flux_per_edge[:, :-1]])   # col left sends right
    from_east  = np.hstack([flux_per_edge[:, 1:], zero_col])    # col right sends left

    return -flux_out + from_north + from_south + from_west + from_east


def diffuse_density_dependent(
    biomass: NDArray[np.float64],
    net_growth: NDArray[np.float64],
    metabolism: float,
    max_disp_rate: float,
    boundary_number: NDArray[np.int64],
    b: float = 10.0,
) -> NDArray[np.float64]:
    """
    One density-dependent diffusion step (Ryser et al. 2021, Eq. 10).

    Per-capita emigration rate is a sigmoid of (metabolism - net_growth):
        d = max_disp_rate / (1 + exp(-b * (metabolism - net_growth)))
    High when conditions are bad (net_growth < metabolism), low when good.
    Survival during dispersal is set to 1 (no matrix mortality).

    Parameters
    ----------
    biomass : (n_rows, n_cols) array
        Current biomass densities.
    net_growth : (n_rows, n_cols) array
        Per-capita net growth rate nu_{i,z}: (feeding - losses - metabolism)
        / biomass. Computed from ATN state before calling this function.
        Same shape as biomass.
    metabolism : float
        Species metabolic rate x_i; sigmoid inflection point (dispersal switches
        around this value).
    max_disp_rate : float
        Maximum per-capita dispersal rate (parameter a in Ryser et al. 2021).
    boundary_number : (n_rows, n_cols) int array
        Neighbour counts from enforce_boundary_conditions(). Same shape as biomass.
    b : float
        Sigmoid steepness (default 10, as in Ryser et al. 2021).

    Returns
    -------
    (n_rows, n_cols) array
        Net biomass change (immigration - emigration) for this step.
        Add to the current biomass to advance the state.
    """
    assert biomass.shape == net_growth.shape == boundary_number.shape, (
        f"Shape mismatch: biomass {biomass.shape}, net_growth {net_growth.shape}, "
        f"boundary_number {boundary_number.shape}"
    )

    d_grid = max_disp_rate / (1 + np.exp(-b * (metabolism - net_growth)))

    flux_per_edge = biomass * (d_grid / 4)
    flux_out = flux_per_edge * boundary_number

    n_rows, n_cols = biomass.shape
    zero_row = np.zeros((1, n_cols))
    zero_col = np.zeros((n_rows, 1))

    from_north = np.vstack([zero_row, flux_per_edge[:-1, :]])
    from_south = np.vstack([flux_per_edge[1:, :], zero_row])
    from_west  = np.hstack([zero_col, flux_per_edge[:, :-1]])
    from_east  = np.hstack([flux_per_edge[:, 1:], zero_col])

    return -flux_out + from_north + from_south + from_west + from_east


def run_diffusion_simulation(
    biomass: NDArray[np.float64],
    max_disp_rate: float,
    n_time: int,
    net_growth: NDArray[np.float64],
    metabolism: float,
    b: float = 10.0,
) -> list[NDArray[np.float64]]:
    """
    Density-dependent diffusion simulation over n_time steps.

    Parameters
    ----------
    biomass : (n_rows, n_cols) array
        Initial biomass densities.
    max_disp_rate : float
        Maximum per-capita dispersal rate.
    n_time : int
        Number of time steps to simulate.
    net_growth : (n_rows, n_cols) array
        Per-capita net growth rates; held constant across steps when ATN is
        not coupled. Same shape as biomass.
    metabolism : float
        Species metabolic rate; sigmoid inflection point.
    b : float
        Sigmoid steepness (default 10).

    Returns
    -------
    list of (n_time + 1) arrays
        Biomass snapshots from t = 0 to t = n_time (index 0 is the initial state).
    """
    boundary_number = enforce_boundary_conditions(*biomass.shape)

    results = [biomass.copy()]
    for _ in range(n_time):
        biomass = biomass + diffuse_density_dependent(
            biomass, net_growth, metabolism, max_disp_rate, boundary_number, b
        )
        results.append(biomass.copy())

    return results
