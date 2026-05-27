"""Vegetation science functions. Pure numpy, no engine imports."""

import numpy as np
from numpy.typing import NDArray

def logistic_growth_delta(
    biomass: NDArray[np.float64],
    growth_rate: NDArray[np.float64],
    carrying_capacity: NDArray[np.float64],
    dt: float,
) -> NDArray[np.float64]:
    """Per-cell logistic growth applied to a biomass array.

    The function is shape-agnostic — the same code runs on a 1-D species
    vector or a full spatial grid. The team convention for which shapes
    to actually use:

      - (S,)       : one location, S species  (smallest unit, handy for tests)
      - (Y, S)     : a row of Y locations, S species each
      - (X, Y, S)  : the full X-by-Y grid, S species each

    Convention on axes: the trailing axis is always species; leading axes
    (if any) are spatial.

    All three array inputs must already have the **same shape**. The
    adapter is responsible for reshaping per-species or per-cell
    parameters (e.g. via ``np.broadcast_to``) before calling this
    function.

    Returns the biomass delta over one time step ``dt``.
    """
    # All inputs must already share the same shape — catches broadcast
    # mistakes in the adapter early, with a clear failure point.
    assert biomass.shape == growth_rate.shape == carrying_capacity.shape

    delta_biomass = dt * growth_rate * biomass * (1.0 - biomass / carrying_capacity)

    return delta_biomass
