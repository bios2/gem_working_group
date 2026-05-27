# Discretization: ODE vs discrete-time, and what the engine contract should be

This document synthesizes the findings, analysis, and recommendation from a working session on whether the simulation engine should ask processes to return **discrete biomass deltas** (`delta_B` per engine time step) or **continuous-time derivatives** (`dB/dt`). It supports criticism #2 and criticism #4 in [`README.md`](README.md) and the contract specification in [`context/process_contract_spec.md`](context/process_contract_spec.md).

## 1. Background: two ways of writing a dynamical process

Both styles describe how a quantity changes over time; they differ in how the rule is written and how time advances.

- **Discrete-time.** The rule is written as `B(t + dt) = B(t) + f(B(t), env, params)`. Time advances in fixed chunks of size `dt`. Implementation is straightforward NumPy arithmetic in a loop. This is what `engine_v2` currently does.
- **ODE (Ordinary Differential Equation).** The rule is written as `dB/dt = g(B, env, params)` — an instantaneous rate. To get values at concrete time points, a numerical solver (e.g. `scipy.integrate.odeint`, an adaptive Runge-Kutta scheme) integrates the rate continuously, internally choosing substep sizes to stay within an error tolerance. This is what the ATN implementation on the `patn` branch does.

Logistic growth in `engine_v2` and logistic growth as a term in the ATN equations are the *same equation*, written in the two styles.

## 2. Findings from the current code

### 2.1 ATN (`patn` branch) is written as an ODE

`experiments/atn_model/atn_model.py` separates the science from the integration:

- `derivatives(y, t, cell_idx)` returns `dB/dt` for one cell, expressed directly as the Binzer/Schneider unscaled ATN equations (logistic growth for basals, Holling-II feeding for consumers, metabolic loss, predation loss).
- `run_cell(...)` calls `scipy.integrate.odeint` on `derivatives` to integrate the cell's dynamics over a time grid.
- `run_all_cells(...)` loops over cells in Python, calling `run_cell` once per cell.

This structure already embodies the split between *pure science* and *engine integration* that criticism #2 in `README.md` recommends. It is also where the global-scale performance issue lives: at 180×360 cells, `run_all_cells` makes 64,800 Python-level `scipy.odeint` calls per simulated period.

### 2.2 `engine_v2` is written as discrete-time

`experiments/tuesday_architecture_proposal/engine_v2/processes.py::apply_vegetation_growth` writes:

```python
plants += (r * plants * (1 - (plants / cc)))
```

This is one Euler step of the logistic equation, applied in-place at each engine time step. The science and the integration are fused into a single mutating expression. There is no explicit `dt`; the time step is whatever a call to `model.step()` represents.

### 2.3 The two styles are not currently compatible

ATN cannot be dropped into `engine_v2` as-is: it produces a rate and expects an external solver, while `engine_v2` expects each process to advance biomass by one step. Some convention has to be chosen.

## 3. Analysis: what each style buys and costs

The choice was evaluated against four priorities stated for the project: performance at global scale and on long time series, testability of individual processes, reusability across processes, and scientific transparency for non-engineer contributors.

| Dimension | Discrete-time (delta) | ODE (derivative) |
|---|---|---|
| Performance at global scale | One NumPy operation updates `(X, Y, S)` in vectorised form. Naturally fast. | Adaptive solvers require either a per-cell Python loop (slow, as in current ATN) or a global ODE system (entangles processes). |
| Composability of process types | Rates, stochastic events, kernels, lookups all fit under "produce a delta over `dt`". Processes compose by summation of deltas. | Only continuous-rate processes fit cleanly. Stochastic (fire), kernel-based (dispersal hops), and discrete-event processes require operator splitting — which collapses back to discrete-time. |
| Testability of one process | Pure function `f(state, params) -> delta`. One assertion: hand it arrays, check the returned delta. | Pure function `f(state, params) -> dB/dt` is meaningful only paired with a numerical scheme; behaviour observed only after integration. |
| Reusability | One uniform `(state, env, params) -> delta` signature. | Either a second contract beside the delta contract (engine complexity), or every process expressed as a rate (constrains process types). |
| Scientific transparency | Numerical scheme is visible in the process's adapter (`delta = rk4_step(...)` or `delta = dt * f(...)`). | Numerical scheme is hidden inside the framework; a reader has to ask "what solver, what tolerance, what step size". |
| Coupling within one step | Needs the delta-layer trick (criticism #1) to be order-independent. | Handled natively by the solver. |
| Numerical accuracy on smooth dynamics | Depends on the chosen `dt`. | Adaptive error control out of the box. |

The only column ODE wins outright is **within-step coupling and accuracy on smooth, stiff dynamics**. Criticism #1's delta-layer trick is precisely the discrete-time fix for within-step coupling, which closes most of that gap. The remaining accuracy gap is addressed by either choosing a small enough engine `dt` or by letting a process integrate internally with a finer substep.

## 4. The performance argument in detail

The strongest single argument against an ODE engine contract is performance at global spatial scale.

- Vectorised NumPy on `(X, Y, S)` updates the entire grid in one operation. At 180×360 cells with N species, this is a fraction of a second per step on a laptop, regardless of N (within memory limits).
- Per-cell `scipy.odeint` is bound by Python call overhead and solver setup costs that scale with the number of cells, not the math. Tens of thousands of cells means tens of thousands of solver invocations per simulated period.
- An ODE *engine* (one global solver integrating a flattened state vector across all cells and species) avoids the per-cell loop, but then every process must be expressed as a contribution to a single `dy/dt` vector. Heterogeneous process types (stochastic ignition, kernel-based dispersal, discrete climate forcing) break this framing; they are not differentiable in the sense an ODE solver requires.

Discrete-time as the engine contract preserves vectorisation across the whole grid and remains agnostic to process type.

## 5. Does ODE-to-delta lose scientific quality?

The natural concern is: if ATN moves off `scipy.odeint`, do we lose accuracy or fidelity?

What is lost:

- Adaptive substep sizing (the solver picks substeps automatically based on local error).
- Tolerance-based error control (`rtol`, `atol` knobs).

What replaces it:

- A fixed-scheme integrator inside ATN's adapter, typically RK4 with a small fixed substep (engine `dt` divided into N internal substeps).
- A one-time validation: run the new vectorised RK4 implementation against the existing per-cell `scipy.odeint` implementation on a representative scenario, tune the substep count until trajectories match to a documented tolerance, lock it in.

For non-stiff food-web dynamics — the regime ATN typically operates in — RK4 with a sensible substep matches `odeint` to several decimals and runs orders of magnitude faster on a global grid. Madingley and published spatial ATN implementations both rely on fixed-step integration at scale; the trade-off is well understood in the literature.

The caveat: if ATN parameter sweeps push into genuinely stiff regimes (extreme rate ratios, near-extinction transients), fixed-step RK4 can become inaccurate or unstable where adaptive stepping would silently rescue it. The mitigation is that the substep count is a knob owned by the ATN adapter; it can be increased locally without touching the engine.

Net effect on quality: preserved, with one upfront validation exercise.

## 6. One real cost of a uniform discrete-time engine

If every process inherits the same engine `dt`, the choice of `dt` is constrained by the *fastest* process. If vegetation succession evolves over decades and metabolism over hours, the engine's `dt` has to be short enough for the fast process, which wastes work on the slow ones. ODE solvers handle multi-rate problems automatically; here, the management is manual: a fast process can subcycle internally inside its adapter (run several substeps for one engine step), or a slow process can be skipped on most steps (apply once every N calls). Both are simple to implement and both are visible at the adapter level rather than hidden in framework code.

This cost is acknowledged but considered acceptable: explicit per-process pacing is easier for an ecology team to reason about than tuning solver tolerances.

## 7. Proposal

**The engine contract is a discrete biomass delta per engine time step.** One signature, no exceptions:

```text
delta_B = process_science(state_arrays..., params)   # pure, NumPy in / NumPy out
apply_<process>(grid, env)                           # adapter, calls science and writes delta
```

ODE-style processes integrate the rate to a delta **inside their own adapter** using a shared numerical helper:

```python
# adapter for an ODE-style process
def apply_atn(grid, env):
    B = grid.get_layer_view("biomass")
    dBdt_fn = lambda B: atn_derivatives(B, traits, adj, env, params)
    delta = shared.numerics.rk4_step(dBdt_fn, B, dt=engine.dt, n_substeps=4)
    grid.add_delta("biomass", delta)
```

The integration scheme and substep count are written explicitly in the adapter — a reader sees in one file both the science (the rate) and the numerical choice (RK4 with N substeps). This satisfies the four priorities:

- **Performance**: vectorised across the whole grid; no per-cell Python loops.
- **Testability**: science function is a pure NumPy-in, NumPy-out function; trivially unit-tested.
- **Reusability**: one engine signature; processes compose by summing deltas; shared `numerics` helpers reused across adapters.
- **Scientific transparency**: numerical scheme visible in the same file as the science, not hidden in framework code.

The full module layout and signature conventions are specified in [`context/process_contract_spec.md`](context/process_contract_spec.md).

## 8. Validation work the proposal implies

If the proposal is adopted, the ATN team carries out one validation task before relying on the new path:

1. Port `derivatives` into a vectorised `atn_derivatives(B, traits, adj, env, params) -> dB_dt` operating on `(..., S)` arrays.
2. Build the adapter `apply_atn(grid, env)` using `shared.numerics.rk4_step`.
3. Run a representative scenario (e.g. one of the existing `reproduce_atnr_figures.py` cases) under both the legacy per-cell `scipy.odeint` path and the new vectorised RK4 path.
4. Tune `n_substeps` until trajectories match to a documented tolerance; commit the comparison plot.
5. Document the chosen `n_substeps` and its validation in the ATN adapter.

After this, the team can trust the new path as the production path; the legacy `scipy.odeint` path is retained as a reference implementation for future revalidation when ATN parameter regimes change.
