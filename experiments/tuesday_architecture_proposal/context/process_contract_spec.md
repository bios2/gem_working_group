# Process contract specification

This document specifies the contract that every ecosystem process (vegetation, ATN, dispersal, fire, ...) must follow to be pluggable into the simulation engine described in the architecture proposal. It complements `process_synthesis_atn.md`, `process_synthesis_vegetation.md`, and `process_synthesis_distribution.md` by describing *how* a process is exposed to the engine, independently of *what* it computes.

The spec has three parts: the output contract (discrete deltas), the shape contract (broadcast-friendly dimensionality), and the module layout (pure science separated from engine adapter).

## 1. Output contract: discrete biomass deltas

Every science function returns a **finite biomass delta for one engine time step**, not an instantaneous rate.

```text
delta_B = process_science(state_arrays..., params)
```

- `state_arrays...`: the inputs the process reads (biomass slices, environmental layers, parameters, ...). All inputs are NumPy arrays or scalar parameters. No engine objects, no `grid`, no `env`.
- `delta_B`: the biomass change produced by this process over one time step `dt`, on the same shape as the biomass slice the process operates on. The engine accumulates `delta_B` from every process into a shared delta layer and applies the sum at the end of the step.

Processes whose scientific formulation is naturally a differential equation `dB/dt = f(B, ...)` (e.g. ATN) integrate that rate to a delta **inside their own adapter** before returning, using a shared numerical helper (see §3). The engine never sees rates; it only sees deltas.

Rationale and the full discussion behind this choice are in [`../discretization.md`](../discretization.md).

## 2. Shape contract: `(..., S)` with species last

Science functions are written to be **broadcast-friendly across leading spatial axes**.

- State arrays have shape `(..., S)` with species on the last axis.
- The leading `...` axes are spatial and arbitrary: `()` for a unit test on a single cell, `(X, Y)` for a production run on the full grid, `(N_cells,)` for a flattened representation.
- Coupling **along the species axis** (ATN feeding, predation, competition between species) is expressed as numpy reductions on `axis=-1` or `einsum` along that axis.
- Coupling **along spatial axes** (dispersal, fire spread) is expressed as numpy operations on the leading axes (shifts, convolutions, neighbour stencils).

This is a stylistic convention, not engine-enforced. Following it means the same code unit-tests on a `(S,)` array and runs vectorised on `(X, Y, S)` in production without per-cell Python loops. It is the convention that lets ATN's per-cell food-web equations run at global scale without rewriting the math.

## 3. Module layout: science, adapters, shared dependencies

Each process lives in its own module, with a clear split between pure science and engine glue.

```text
engine_v2/
├── processes/
│   ├── vegetation/
│   │   ├── science.py        # pure functions, no engine imports
│   │   └── adapter.py        # thin engine integration
│   ├── atn/
│   │   ├── science.py        # pure dB/dt or dB function
│   │   └── adapter.py        # owns the integration scheme if science is dB/dt
│   ├── dispersal/
│   │   ├── science.py
│   │   └── adapter.py
│   └── shared/
│       ├── numerics.py       # rk4_step, etc.
│       └── kernels.py        # reusable spatial kernels, neighbour stencils
└── EcosystemEngine.py
```

A single-file layout (`vegetation.py` exposing both `vegetation_growth` and `apply_vegetation_growth`) is acceptable for processes small enough not to warrant a directory, but the science / adapter split should still be visible inside the file.

### 3.1 Science module (`science.py`)

- Pure functions. No imports from the engine, no references to `grid`, `env`, `registry`, or any state object.
- Takes NumPy arrays and parameters; returns NumPy arrays.
- Follows the shape contract from §2.
- Output is either a delta `delta_B` (preferred, satisfies the contract directly) or a rate `dB_dt` (acceptable when the science is fundamentally an ODE; integration happens in the adapter).
- Trivial to unit-test with hand-built arrays. Easy to plot from a notebook without instantiating the engine.

### 3.2 Adapter module (`adapter.py`)

- The only piece that knows about the engine. Imports `EcosystemGridState`, `EnvironmentState`, `SpeciesRegistry`.
- Exposes one function with the engine signature `apply_<process>(grid, env)`. This is what `EcosystemEngine.add_process` registers.
- Responsibilities:
  - Slice the right inputs out of `grid` / `env` / `grid.registry` (biomass for relevant groups, environmental layers, per-species parameters).
  - Call the science function.
  - If the science function returns a rate, call the shared integrator (e.g. `numerics.rk4_step`) to convert it to a delta over `dt`.
  - Write the delta into the shared delta layer registered on `grid` (the engine sums and applies it at end of step; see criticism #1 in `../README.md`).
- The adapter is the place where any non-uniform routing lives: which groups the process applies to, which layers it reads and writes, what `dt` substep count it uses. The science function stays uniform; the adapter encodes the process-specific plumbing.

### 3.3 Shared dependencies (`shared/`)

- Cross-process helpers that several adapters reuse. Currently anticipated:
  - `numerics.rk4_step(derivative_fn, state, dt, n_substeps)`: vectorised fixed-step RK4 integrator. The adapter for any ODE-style process uses this to convert `dB/dt` to a delta over one engine time step.
  - `kernels`: reusable spatial neighbour stencils and dispersal kernels.
- These are optional. A process that does not need them ignores `shared/` entirely.
- `shared/` does **not** import the engine. It depends only on NumPy.

## 4. Example: ATN under this contract

ATN's existing implementation on the `patn` branch already separates `derivatives(y, t, cell_idx)` (the science) from `run_cell` / `run_all_cells` (the integration loop). Porting it to this contract is a refactor, not a rewrite of the science:

1. Move `derivatives` into `processes/atn/science.py` as a vectorised function `atn_derivatives(B, traits, adj_mat, env, params) -> dB_dt`, with leading axes broadcast over space (replace the per-cell Python loop in `run_all_cells` with broadcasting on `(X, Y, S)`).
2. Create `processes/atn/adapter.py` with `apply_atn(grid, env)`:
   - Pull `B`, traits, adjacency, environmental layers out of the shared state.
   - Build a closure `dBdt_fn = lambda B: atn_derivatives(B, ...)`.
   - Compute `delta = shared.numerics.rk4_step(dBdt_fn, B, dt=engine.dt, n_substeps=...)`.
   - Write `delta` into the engine's delta layer for biomass.
3. Validate once: run the new path against the existing per-cell `scipy.odeint` implementation on a small grid; confirm trajectories match to a documented tolerance; lock in the substep count.

The same shape applies to any future ODE-style process. The same shape, with `n_substeps=1` and `delta = dt * f(B)` (forward Euler) — or with the science function returning a delta directly — applies to discrete-time processes like vegetation logistic growth.
