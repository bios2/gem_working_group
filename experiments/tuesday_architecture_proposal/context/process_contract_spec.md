# Process contract specification

This document specifies the contract that every ecosystem process (vegetation, ATN, dispersal, fire, ...) must follow to be pluggable into the simulation engine described in the architecture proposal. It complements `process_synthesis_atn.md`, `process_synthesis_vegetation.md`, and `process_synthesis_distribution.md` by describing *how* a process is exposed to the engine, independently of *what* it computes.

The spec has four parts: the output contract (discrete deltas), the shape convention (broadcast-friendly numpy at the process's natural dimensionality), the typing requirements (dtype-checked signatures and a runtime shape guard), and the module layout (science modules independent of the engine, adapters consolidated next to it).

This file is intended to be reused as the canonical reference for `src/` contribution rules in future `README.md` and `AGENTS.md` / `CLAUDE.md` artifacts.

## 1. Output contract: a discrete biomass delta per time step

Every science function returns a **finite biomass change for one engine time step**, not an instantaneous rate of change.

```text
delta_B = process_science(state_arrays..., scalar_params..., dt)
```

- `state_arrays...`: the inputs the process reads (biomass slices, environmental layers, per-species parameters reshaped to broadcast). All are NumPy arrays. No engine objects, no `grid`, no `env`, no `registry`.
- `scalar_params...`: plain Python scalars (rates, exponents, ratios). No dicts, no objects.
- `dt`: the engine's time step, in whatever units the project agrees on (typically days). Always the last positional argument. Always a `float`.
- `delta_B`: the biomass change produced by this process over `dt`, same shape as the primary biomass input. The engine accumulates `delta_B` from every process into a shared delta layer and applies the sum at the end of the step (see criticism #1 in `../README.md`).

Processes whose scientific formulation is naturally a differential equation `dB/dt = f(B, ...)` (e.g. ATN) integrate that rate to a delta **inside their adapter**, using a shared numerical helper from `numerics.py` (e.g. a vectorised RK4 step). The engine never sees rates; it only sees deltas. Rationale: see [`../discretization.md`](../discretization.md).

## 2. Shape convention: each process at its natural dimensionality, written to broadcast

Each process is written at the dimensionality that is natural to its science:

- **Cell-local processes** (e.g. logistic growth) are written as elementwise math on whatever array shape they receive. The same code runs on a single cell, a row, or the whole `(X, Y, S)` grid.
- **Within-cell species coupling** (e.g. ATN feeding, predation) is written as numpy reductions on the species axis (`axis=-1` is the recommended convention but not required).
- **Between-cell coupling** (e.g. dispersal, fire spread) is written as numpy operations on the spatial axes — shifts, neighbour stencils, convolutions.

The discipline is: **write the science function so that the same numpy code works at any leading shape**. A function intended to apply per-cell should not assume `B.shape == (n_species,)`; it should use elementwise numpy operations (`*`, `+`, `np.where`, etc.) that broadcast naturally. A function that does within-cell reductions should reduce on a named axis (`axis=-1`), not on `axis=0`, so leading spatial axes pass through.

This discipline is what makes a single unit test on a `(3,)` array representative of a production run on `(180, 360, 3)`. It also avoids per-cell Python loops in adapters, which is the main performance trap of mixing scientific formulations with a global-scale engine.

## 3. Typing requirements

Type annotations harden the contract enough to catch real mistakes (wrong dtype, missing argument, mismatched scalar type) without forcing the team to learn advanced typing.

```python
import numpy as np
from numpy.typing import NDArray

def <process>_delta(
    biomass: NDArray[np.float64],
    *other_state_arrays: NDArray[np.float64],
    *scalar_params: float,
    dt: float,
) -> NDArray[np.float64]:
    ...
```

Rules:

1. **All array inputs and outputs are `NDArray[np.float64]`.** This catches `int` arrays slipping in, which silently change division semantics.
2. **Scalar parameters are `float` (or `int` when truly integer).** No dicts of parameters; pass values explicitly.
3. **`dt` is always the last positional argument and always `float`.** This is the only argument the engine itself supplies; everything else is fetched by the adapter from the shared state.
4. **Return type matches `biomass`'s shape.** Same dtype, same shape as the primary biomass input. Enforced by broadcasting at execution; documented in the docstring.
5. **One-line runtime shape guard at the top of the function** is recommended for any function that operates on more than one array:
   ```python
   assert biomass.shape == growth_rate.shape == carrying_capacity.shape
   ```
   This catches shape mismatches early with a clear message, instead of silently producing wrong broadcasts deeper in the math.

Static shape typing (e.g. PEP 646 variadic generics) is not required. Shapes are documented in the docstring and validated at runtime by the assert.

## 4. Module layout: flat package, science and adapters separated

The runtime package is flat — no subpackages, no nested directories:

```text
gem_working_group/
├── src/
│   └── gem_working_group/
│       ├── __init__.py
│       ├── EcosystemEngine.py        # the pipeline
│       ├── EcosystemGridState.py     # shared state (the "cart")
│       ├── EnvironmentState.py
│       ├── SpeciesRegistry.py
│       ├── vegetation.py             # science: pure functions, no engine imports
│       ├── atn.py                    # science: pure functions, no engine imports
│       ├── dispersal.py              # science: pure functions, no engine imports
│       ├── processes.py              # all engine adapters live here
│       └── numerics.py               # shared helpers (rk4_step, kernels, ...)
└── tests/
    ├── test_vegetation.py
    ├── test_atn.py
    └── test_dispersal.py
```

### 4.1 Science modules (`vegetation.py`, `atn.py`, `dispersal.py`, ...)

- Pure functions matching §3's typed signature.
- Import nothing from the engine. No `EcosystemGridState`, no `EnvironmentState`, no `SpeciesRegistry`. Only `numpy` and the standard library (and `numerics` if the function is an ODE rate that will be integrated by RK4).
- Follow the shape convention from §2.
- Can be loaded and tested in a notebook without instantiating the engine.

### 4.2 Adapters (`processes.py`)

- Single file holding all `apply_<process>(grid, env, dt)` functions registered with the engine.
- The only place that imports `EcosystemGridState`, `EnvironmentState`, and `SpeciesRegistry`.
- Each adapter is responsible for:
  - Slicing the right state arrays out of `grid` / `env` / `grid.registry` (biomass for the relevant group, environmental layers, per-species parameters reshaped to broadcast against the biomass shape).
  - Calling the matching science function.
  - For ODE-style processes: calling `numerics.rk4_step` (or another helper from `numerics.py`) to convert the returned rate into a delta over `dt`.
  - Writing the delta into the engine's shared delta layer via `grid.add_delta(...)`.
- A reader opening `processes.py` sees the full engine pipeline in one place: every adapter, every glue line, every layer it touches.

### 4.3 Shared helpers (`numerics.py`)

- Cross-process numerical utilities. Currently anticipated:
  - `rk4_step(derivative_fn, state, dt, n_substeps)`: vectorised fixed-step Runge-Kutta integrator, used by ODE-style adapters.
  - Reusable spatial kernels and neighbour stencils for dispersal-like processes.
- Depends only on NumPy. Does not import the engine.

## 5. Worked example: per-cell logistic growth

The simplest case: a process applied independently to each cell, with no spatial or cross-species coupling.

### 5.1 Science (`vegetation.py`)

```python
# src/gem_working_group/vegetation.py
import numpy as np
from numpy.typing import NDArray

def logistic_growth_delta(
    biomass: NDArray[np.float64],
    growth_rate: NDArray[np.float64],
    carrying_capacity: NDArray[np.float64],
    dt: float,
) -> NDArray[np.float64]:
    """
    Per-cell logistic growth applied to a biomass array.

    All array inputs broadcast against a common shape, typically (..., S_plants)
    where leading axes are spatial. The same function runs on a single cell
    (shape (S,)), a row (shape (Y, S)), or the whole grid (shape (X, Y, S)).

    Returns the biomass delta over one time step `dt`.
    """
    assert biomass.shape == growth_rate.shape == carrying_capacity.shape
    return dt * growth_rate * biomass * (1.0 - biomass / carrying_capacity)
```

Two lines of math, fully typed, runs at any leading dimensionality.

### 5.2 Adapter (`processes.py`)

```python
# src/gem_working_group/processes.py
import numpy as np
from . import vegetation
from .EcosystemGridState import EcosystemGridState
from .EnvironmentState import EnvironmentState

def apply_vegetation_growth(
    grid: EcosystemGridState,
    env: EnvironmentState,
    dt: float,
) -> None:
    plant_idx = grid.registry.get_group_indices("plants")
    B = grid.layers["biomass"][..., plant_idx]                            # (X, Y, S_plants)
    r = grid.registry.get_group_parameter("plants", "base_growth_rate")   # (S_plants,)
    K = env.get_layer("carrying_capacity")[..., np.newaxis]               # (X, Y, 1)
    r_b = np.broadcast_to(r, B.shape)
    K_b = np.broadcast_to(K, B.shape)
    delta = vegetation.logistic_growth_delta(B, r_b, K_b, dt)
    grid.add_delta("biomass", plant_idx, delta)
```

The adapter is the only piece that knows about `grid`, `env`, and the species registry. The science file is untouched by these concerns and can be unit-tested against hand-built arrays.

### 5.3 Unit test (`tests/test_vegetation.py`)

```python
# tests/test_vegetation.py
import numpy as np
from gem_working_group.vegetation import logistic_growth_delta

def test_logistic_growth_zero_at_carrying_capacity():
    B = np.array([100.0, 100.0, 100.0])
    r = np.array([0.1, 0.2, 0.3])
    K = np.array([100.0, 100.0, 100.0])
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    np.testing.assert_allclose(delta, 0.0)

def test_logistic_growth_runs_on_grid():
    shape = (4, 5, 3)
    B = np.full(shape, 50.0)
    r = np.full(shape, 0.1)
    K = np.full(shape, 100.0)
    delta = logistic_growth_delta(B, r, K, dt=1.0)
    assert delta.shape == shape
    assert np.all(delta > 0)
```

Both tests exercise the same function with no engine, no fixtures, no mocking.

## 6. Worked example sketch: ATN under this contract

ATN's current implementation on the `patn` branch already separates the science (`derivatives(y, t, cell_idx)`) from the integration (`run_cell` / `run_all_cells`). Porting it to this contract is a refactor of structure, not a rewrite of the equations.

### 6.1 Science (`atn.py`)

```python
# src/gem_working_group/atn.py
import numpy as np
from numpy.typing import NDArray

def atn_derivative(
    biomass: NDArray[np.float64],
    body_mass: NDArray[np.float64],
    adjacency: NDArray[np.float64],
    temperature_K: NDArray[np.float64],
    carrying_capacity: NDArray[np.float64],
    is_basal: NDArray[np.float64],
    dt: float,                # unused: this function returns a rate
    # ... additional scalar params for allometric and FR constants ...
) -> NDArray[np.float64]:
    """
    Unscaled spatial ATN model (Binzer/Schneider). Returns dB/dt at the
    current state; integration to a delta is handled by the adapter.

    `biomass` has shape (..., S). Within-cell species coupling is expressed
    via reductions on axis=-1. Leading axes pass through unchanged.
    """
    # ... feeding (Holling II on axis=-1), metabolism, growth, predation ...
    return dB_dt
```

`dt` appears in the signature for uniformity; the function ignores it because it returns a rate. The adapter is what converts that rate into a delta over `dt`.

### 6.2 Adapter (`processes.py`)

```python
# src/gem_working_group/processes.py (continued)
from . import atn, numerics

def apply_atn(grid: EcosystemGridState, env: EnvironmentState, dt: float) -> None:
    B = grid.layers["biomass"]
    body_mass = grid.registry.body_mass
    adjacency = grid.registry.adjacency
    T_K       = env.get_layer("temperature_K")[..., np.newaxis]
    K         = env.get_layer("carrying_capacity")[..., np.newaxis]
    is_basal  = grid.registry.is_basal

    def rate(B_state):
        return atn.atn_derivative(B_state, body_mass, adjacency, T_K, K, is_basal, dt)

    delta = numerics.rk4_step(rate, B, dt=dt, n_substeps=4)
    grid.add_delta("biomass", slice(None), delta)
```

The science function stays a pure NumPy expression of the published ATN equations. The adapter is where the integration scheme (RK4 with 4 substeps) is declared explicitly, in the same file as every other adapter.

## 7. Summary checklist for a new process contribution

- [ ] Add a new science module `src/gem_working_group/<process>.py`.
- [ ] Write the science function with the typed signature from §3.
- [ ] Follow the shape convention from §2 (broadcast-friendly numpy at the process's natural dimensionality).
- [ ] Add the runtime shape `assert` at the top of the function.
- [ ] Add `apply_<process>(grid, env, dt)` in `processes.py`. This is where slicing, parameter fetching, and (if applicable) RK4 integration live.
- [ ] Add a unit test in `tests/test_<process>.py` using hand-built arrays — no engine instance.
- [ ] If the science is an ODE, validate the adapter's RK4 substep count against a reference implementation (`scipy.odeint` on a small grid) and document the chosen substep count.
- [ ] Register the adapter with the engine in the pipeline definition.
