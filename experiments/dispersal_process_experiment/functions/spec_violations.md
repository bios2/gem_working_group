# Spec violations in dispersal.py

This file documents the three ways `dispersal.py` deviates from the process contract
defined in [docs/processes_implementation_specification.md](../../../docs/processes_implementation_specification.md).
These violations are acceptable inside the experiment but must be fixed before the
code graduates to `src/gem/dispersal.py`.

---

## 1. Function names don't follow the `<process>_delta` pattern (§1.1)

The spec requires biomass-modifying science functions to be named `<process>_delta`.

| Current name | Required name |
|---|---|
| `diffuse` | `diffuse_delta` |
| `diffuse_density_dependent` | `diffuse_density_dependent_delta` |

`enforce_boundary_conditions` and `run_diffusion_simulation` are helpers, not
biomass-modifying science functions, so the `_delta` rule does not apply to them.

---

## 2. Missing `dt` as last argument (§1.1, §3)

The spec requires every biomass-modifying science function to take `dt: float` as its
last positional argument and to scale its output by it. The dispersal rate is currently
treated as a per-step rate, but in the engine it will be a rate per unit time (e.g. per
day). The corrected functions should multiply their return value by `dt`:

```python
def diffuse_delta(
    biomass: NDArray[np.float64],
    disp_rate: float,
    boundary_number: NDArray[np.int64],
    dt: float,                            # add this
) -> NDArray[np.float64]:
    ...
    return (-flux_out + from_north + from_south + from_west + from_east) * dt
```

---

## 3. Shape does not follow the `(X, Y, S)` convention (§2)

The spec requires the trailing axis to be species and the leading axes to be spatial,
so functions must accept `(X, Y, S)` arrays and the same code must run on `(S,)`,
`(Y, S)`, and `(X, Y, S)`.

The current functions take `(n_rows, n_cols)` — a single-species 2D grid with no
species axis. This means:
- The function can only handle one species at a time.
- It won't broadcast correctly against per-species parameters in the engine.

Fixing this requires restructuring the spatial shifts so they operate on the first two
axes while leaving the trailing species axis untouched. For example, the north shift
becomes:

```python
# current — works on (n_rows, n_cols) only
from_north = np.vstack([zero_row, flux_per_edge[:-1, :]])

# corrected — works on (..., n_rows, n_cols, S) via ellipsis indexing
from_north = np.concatenate([np.zeros_like(flux_per_edge[:1, ...]),
                              flux_per_edge[:-1, ...]], axis=0)
```

The shape assert would also need to check that `biomass.ndim >= 2` (at least one
spatial axis and one species axis) rather than asserting a specific 2D shape.
