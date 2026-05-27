# Process contract specification

This document specifies the contract that every ecosystem process (vegetation, ATN, dispersal, metabolism, fire, ...) must follow to be pluggable into the simulation engine described in the architecture proposal. It complements `process_synthesis_atn.md`, `process_synthesis_vegetation.md`, and `process_synthesis_distribution.md` by describing *how* a process is exposed to the engine, independently of *what* it computes.

The spec has six parts: the output contract (biomass deltas for processes that change biomass, named intermediate outputs for dependency processes that don't), the shape convention (broadcast-friendly numpy at the process's natural dimensionality), the typing requirements (dtype-checked signatures and a runtime shape guard), the module layout (science modules independent of the engine, adapters consolidated next to it), worked examples for both process categories, and the handling of ODE-style processes (an explicit section, because their integration scheme is not part of the engine).

This file is intended to be reused as the canonical reference for `src/` contribution rules in future `README.md` and `AGENTS.md` / `CLAUDE.md` artifacts.

## 1. Output contract: what a process 

Processes fall into two categories. The category determines what the science function returns and how the adapter wires it into the shared state.

### 1.1 Biomass-modifying processes — return a biomass delta

Vegetation, ATN, dispersal, fire, and any other process whose *purpose* is to change biomass return a **biomass delta**: the finite change in the biomass array over one engine time step, not an instantaneous rate of change.

```text
biomass_delta = <process>_delta(state_arrays..., scalar_params..., dt)
```

- `state_arrays...`: the inputs the process reads (biomass slices, environmental layers, per-species parameters reshaped to broadcast, any dependency outputs it consumes — see §1.2). All are NumPy arrays. No engine objects, no `grid`, no `env`, no `registry`.
- `scalar_params...`: plain Python scalars (rates, exponents, ratios). No dicts, no objects.
- `dt`: the engine's time step, in whatever units the project agrees on (typically days). Always the last positional argument. Always a `float`.
- `biomass_delta`: the biomass change produced by this process over `dt`, same shape as the primary biomass input. The engine accumulates `biomass_delta` from every process into a shared delta layer and applies the sum at the end of the step (see criticism #1 in `../README.md`).

Processes whose scientific formulation is naturally a differential equation (`dB/dt = f(B, ...)`, e.g. ATN) are addressed in §7.

### 1.2 Dependency processes — return a shared intermediate quantity

Some processes do not themselves change biomass but compute a quantity — a rate, flux, or factor — that **multiple downstream biomass-modifying processes consume**. Metabolism is the canonical example: it appears as a biomass loss term inside ATN *and* as an input to dispersal's emigration trigger (see `process_synthesis_distribution.md`). Net primary productivity is another candidate, reused by basal growth and by herbivory budgets.

These processes follow the same shape and typing conventions (§2, §3), but their signature returns a named scientific quantity instead of a biomass delta and takes no `dt`:

```text
<quantity> = <dependency>(state_arrays..., scalar_params...)
```

- The return value is in the quantity's natural units (e.g. `metabolic_rate` in `biomass / time`, `npp` in `biomass / (area * time)`), with no engine time step baked in. Downstream biomass-modifying processes apply `dt` themselves when they convert the rate or flux into a biomass delta.
- The adapter writes the result into a **shared broadcasting-friendly layer** on the engine state (e.g. `grid.layers["metabolic_rate"]`) so subsequent adapters can fetch it the same way they fetch biomass or environmental layers.
- The engine pipeline must schedule every dependency process before any process that consumes its output. Ordering is the adapter author's responsibility, not the science function's.

A dependency output is *only* worth promoting to a shared layer once it is consumed by more than one biomass-modifying process. A quantity used by a single process should remain a private helper inside that process's science module — adding a shared layer for it is premature.

## 2. Shape convention: each process at its natural dimensionality, written to broadcast

Each process is written at the dimensionality that is natural to its science:

- **Cell-local processes** (e.g. logistic growth) are written as elementwise math on whatever array shape they receive. The same code runs on a single cell, a row, or the whole `(X, Y, S)` grid.
- **Within-cell species coupling** (e.g. ATN feeding, predation) is written as numpy reductions on the species axis (`axis=-1` is the recommended convention but not required).
- **Between-cell coupling** (e.g. dispersal, fire spread) is written as numpy operations on the spatial axes — shifts, neighbour stencils, convolutions.

The discipline is: **write the science function so that the same numpy code works at any leading shape**. A function intended to apply per-cell should not assume `B.shape == (n_species,)`; it should use elementwise numpy operations (`*`, `+`, `np.where`, etc.) that broadcast naturally. A function that does within-cell reductions should reduce on a named axis (`axis=-1`), not on `axis=0`, so leading spatial axes pass through.

This discipline is what makes a single unit test on a `(3,)` array representative of a production run on `(180, 360, 3)`. It also avoids per-cell Python loops in adapters, which is the main performance trap of mixing scientific formulations with a global-scale engine.

## 3. Typing requirements

Type annotations harden the contract enough to catch real mistakes (wrong dtype, missing argument, mismatched scalar type) without forcing the team to learn advanced typing.

Biomass-modifying processes (§1.1):

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

Dependency processes (§1.2):

```python
def <dependency>(
    *state_arrays: NDArray[np.float64],
    *scalar_params: float,
) -> NDArray[np.float64]:
    ...
```

Rules:

1. **All array inputs and outputs are `NDArray[np.float64]`.** This catches `int` arrays slipping in, which silently change division semantics.
2. **Scalar parameters are `float` (or `int` when truly integer).** No dicts of parameters; pass values explicitly.
3. **`dt` is always the last positional argument and always `float` — for biomass-modifying processes only.** Dependency processes return a rate or state value and never take `dt`; downstream consumers apply `dt` themselves. `dt` is the only argument the engine itself supplies; everything else is fetched by the adapter from the shared state.
4. **Return shape.** For biomass-modifying processes, same dtype and shape as the primary `biomass` input. For dependency processes, document the returned shape in the docstring (typically the broadcast shape of the inputs). Enforced at execution; not statically typed.
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
│       ├── metabolism.py             # dependency science (§1.2), pure functions
│       └── processes.py              # all engine adapters live here
└── tests/
    ├── test_vegetation.py
    ├── test_atn.py
    ├── test_dispersal.py
    └── test_metabolism.py
```

### 4.1 Science modules (`vegetation.py`, `atn.py`, `dispersal.py`, `metabolism.py`, ...)

- Pure functions matching §3's typed signature (biomass-modifying *or* dependency).
- Import nothing from the engine. No `EcosystemGridState`, no `EnvironmentState`, no `SpeciesRegistry`. Only `numpy` and the standard library.
- Follow the shape convention from §2.
- Can be loaded and tested in a notebook without instantiating the engine.
- Dependency science modules (e.g. `metabolism.py`) may be imported by other science modules as helpers when the quantity is also used inline (e.g. `atn.py` calls `metabolism.metabolic_rate(...)` internally). The same function is exposed to a dedicated adapter only once a second process consumes the result.

### 4.2 Adapters (`processes.py`)

- Single file holding all `apply_<process>(grid, env, dt)` functions registered with the engine.
- The only place that imports `EcosystemGridState`, `EnvironmentState`, and `SpeciesRegistry`.
- Each adapter is responsible for:
  - Slicing the right state arrays out of `grid` / `env` / `grid.registry` (biomass for the relevant group, environmental layers, per-species parameters reshaped to broadcast against the biomass shape, dependency outputs from shared layers).
  - Calling the matching science function.
  - **For biomass-modifying processes:** writing the delta into the engine's shared delta layer via `grid.add_delta(...)`.
  - **For dependency processes:** writing the returned quantity into a shared broadcasting-friendly layer (e.g. `grid.layers["metabolic_rate"]`) so downstream adapters in the same step can fetch it.
- Adapter scheduling enforces ordering: dependency adapters run before any biomass-modifying adapter that reads their output.
- A reader opening `processes.py` sees the full engine pipeline in one place: every adapter, every glue line, every layer it touches.

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

## 6. Worked example: metabolism as a dependency process

Metabolism returns a per-species mass-specific rate that is reused by ATN (as a loss term) and by dispersal (as part of the emigration trigger). It is a textbook dependency process (§1.2): no `dt`, no biomass delta, just a rate written to a shared layer for downstream adapters to read.

### 6.1 Science (`metabolism.py`)

```python
# src/gem_working_group/metabolism.py
import numpy as np
from numpy.typing import NDArray

def metabolic_rate(
    body_mass: NDArray[np.float64],
    temperature: NDArray[np.float64],
    X0: float,
    b_X: float,
    E_a: float,
    k_B: float,
    T0_K: float,
) -> NDArray[np.float64]:
    """
    Mass- and temperature-scaled metabolic rate, in units of biomass / time.

    Returns an array broadcast over the shapes of `body_mass` and `temperature`
    (typically (S,) and (X, Y) respectively, broadcasting to (X, Y, S)).
    No `dt`: this is a rate. Consumers multiply by `dt` themselves.
    """
    assert body_mass.ndim >= 1 and temperature.ndim >= 1
    boltzmann = np.exp(-E_a * (1.0 / (k_B * temperature) - 1.0 / (k_B * T0_K)))
    return X0 * body_mass**b_X * boltzmann
```

### 6.2 Adapter (`processes.py`)

```python
def compute_metabolic_rate(
    grid: EcosystemGridState,
    env: EnvironmentState,
    dt: float,  # accepted for signature uniformity; unused
) -> None:
    M = grid.registry.get_parameter("body_mass")[np.newaxis, np.newaxis, :]   # (1, 1, S)
    T = env.get_layer("temperature_K")[..., np.newaxis]                       # (X, Y, 1)
    rate = metabolism.metabolic_rate(
        np.broadcast_to(M, grid.shape),
        np.broadcast_to(T, grid.shape),
        X0=grid.params["X0"], b_X=grid.params["b_X"],
        E_a=grid.params["E_a"], k_B=grid.params["k_B"], T0_K=grid.params["T0_K"],
    )
    grid.set_layer("metabolic_rate", rate)
```

The adapter writes to a shared layer (`metabolic_rate`) instead of accumulating into the biomass-delta layer. The engine schedules this adapter before `apply_atn` and `apply_dispersal`, both of which fetch `grid.layers["metabolic_rate"]` as one of their inputs.

### 6.3 Consumption inside a biomass-modifying adapter

```python
def apply_atn(grid, env, dt):
    B = grid.layers["biomass"]
    x = grid.layers["metabolic_rate"]            # produced by compute_metabolic_rate
    # ... feeding, assimilation, etc. ...
    delta = atn.atn_delta(B, x, ..., dt)
    grid.add_delta("biomass", slice(None), delta)
```

ATN does not recompute metabolism; it reads the shared layer. Dispersal does the same, which is the whole point of making metabolism a dependency process rather than a private helper inside `atn.py`.

## 7. ODE-style processes (ATN and others)

Some processes are naturally written as ordinary differential equations: their published formulation is a rate `dB/dt = f(B, ...)`, not a per-step delta. ATN on the `patn` branch is the current example; future hydrodynamics, soil dynamics, or temperature-driven processes may also arrive in this form. The engine still only accepts a delta per time step (§1), so the contributor must choose how to bridge from the rate to a delta. There are two acceptable options, and the choice is made by the process's contributor:

### Option A — Integrate the rate inside the adapter

The science module exposes a function returning `dB/dt`, written under the same shape and typing conventions as any other science function (§2, §3). The adapter is the only place that calls an integrator to advance the state by `dt`. Two integrators are reasonable to use:

- `scipy.integrate.solve_ivp` (or the legacy `scipy.integrate.odeint`) with a Runge-Kutta method (`RK45`, `RK23`) when the contributor wants adaptive substep control and tolerance-based error handling. This was ATN's original choice on the `patn` branch.
- A small project-local vectorised RK4 step, when the contributor wants fixed-substep integration that vectorises cleanly over the full `(X, Y, S)` grid (typically faster than per-cell `scipy` calls at global scale).

If multiple ODE-style processes end up sharing the same fixed-step integrator, a small `numerics.py` module added next to the engine is the natural home for it (e.g. an `rk4_step(derivative_fn, state, dt, n_substeps)` helper). This module is **not part of the standard package layout** in §4; it is added only when a second ODE-style process arrives and reuses the integrator. The first ODE process can implement its integrator privately inside its adapter and only extract it once a second user appears.

Whichever integrator is used, the integration scheme and any substep count are written explicitly in the adapter, in the same file as every other adapter (`processes.py`). A reader sees both the science (a rate function in `<process>.py`) and the numerical choice (the integrator call in `processes.py`) without having to look into framework code.

### Option B — Rewrite the science to return a delta directly

The contributor restates the equations in discrete-time form: the science function returns `biomass_delta` directly, with `dt` baked into the math (e.g. a forward Euler step `delta = dt * f(state)`, or a higher-order discrete scheme). The adapter then has no integration step at all, and the process looks structurally identical to a non-ODE process (§5).

This option is appropriate when the discrete approximation is accurate enough at the engine's `dt`, when the contributor wants the simplest possible code path, or when the underlying dynamics are not stiff enough to need adaptive control.

### Choosing between A and B

The trade-offs are documented in [`../discretization.md`](../discretization.md). In short: Option A preserves the published continuous-time fidelity at a small implementation cost; Option B is simpler but requires the contributor to verify that the chosen `dt` (and any internal subcycling) keeps the dynamics accurate. Either way, an ODE-style process should validate its new path against a reference implementation (e.g. `scipy.integrate.solve_ivp` on a small grid) at least once before being trusted in production runs.

## 8. Summary checklist for a new process contribution

- [ ] Decide whether the process is biomass-modifying (§1.1) or a dependency process (§1.2).
- [ ] Add a new science module `src/gem_working_group/<process>.py`.
- [ ] Write the science function with the typed signature from §3 (biomass-modifying takes `dt`; dependency does not).
- [ ] Follow the shape convention from §2 (broadcast-friendly numpy at the process's natural dimensionality).
- [ ] Add the runtime shape `assert` at the top of the function.
- [ ] Add `apply_<process>(grid, env, dt)` in `processes.py`. This is where slicing, parameter fetching, and (if applicable) ODE integration live. Biomass-modifying adapters write to the shared delta layer; dependency adapters write to a named shared layer.
- [ ] Add a unit test in `tests/test_<process>.py` using hand-built arrays — no engine instance.
- [ ] If the science is ODE-style, choose Option A or B from §7, document the choice in the adapter, and validate against a reference integration on a small grid.
- [ ] Register the adapter with the engine in the pipeline definition, with dependency adapters scheduled before any consumer.
