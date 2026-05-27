# Model engine proposal from Alexis

## The efficient janitor and his efficient cart

Alexis introduced this proposal with a janitor analogy. The old janitor calls the store every time he enters a new room to ask for the specific supplies he needs; the workflow only works because he has years of implicit knowledge of the building, and it is nearly impossible to transfer to an intern. The efficient janitor instead carries a **cart** stocked with the tools and products commonly needed across rooms. The cart can be extended when a new need appears, and the same tools are reused everywhere.

Mapped to the model:

- the **store** is the raw data sources scattered across the project (a Madingley-style architecture where each process reaches for whatever it needs);
- the **cart** is a small set of shared state objects (`EcosystemGridState`, `EnvironmentState`, `SpeciesRegistry`) that travel together through every process;
- each **room** is a process (vegetation growth, ATN, dispersal, fire, ...) that operates on the cart's contents and hands it to the next process.

The goal is to avoid the Madingley pain point of intertwined processes that are hard to isolate, swap, or test.

## Objectives

- Provide a modular ecosystem-model engine in which vegetation, ATN, and dispersal/distribution can be developed independently and plugged into the same simulation.
- Define a small set of shared state objects (the "cart") that act as the data contracts between processes, so adding or removing a process is a one-line change in the pipeline.
- Avoid duplicating the same biological information in several places by having every process read from and write to a single referenced state.
- Keep the spatial representation explicit (grid of cells) and the species representation indexed (`Species_ID` axis), so that group-level operations (plants, herbivores, ...) reduce to NumPy slicing.
- Accommodate processes whose science is naturally written at different dimensionalities — per-cell food-web math (ATN), per-species spatial dynamics (dispersal), fully `(X, Y, S)` elementwise dynamics (vegetation logistic growth) — under a single shape contract, so contributors do not have to reshape their equations to match the engine.
- Stay simple enough for the team to maintain (plain Python + NumPy, no specialist framework), while leaving room for tests on individual processes and on full pipelines.

## Modules description

- `EcosystemGridState`: Holds the dynamic 3D `(X, Y, Species_ID)` biological state. Stores a `biomass` layer by default and lets processes register additional layers they need to share (e.g. `net_growth_rate`, `metabolic_loss`). Exposes group-restricted views (`get_layer_view`) and an `edit_group_data` context manager that hands a process a writable working copy for one functional group and writes it back automatically.
- `EnvironmentState`: Holds the 2D `(X, Y)` environmental layers (temperature, carrying capacity, ...). Layers are added by name and shape-checked against the world grid.
- `SpeciesRegistry`: Owns the species list and the functional/trophic group definitions. Stores per-species parameters as 1D arrays for fast vectorised math, exposes group indices and group-sliced parameter arrays, and holds the feeding adjacency matrix (`add_feeding_link(resource_group, consumer_group)`).
- `processes`: A library of step functions with a common signature `process(grid, env)`. Each process reads from `grid` / `env` / `grid.registry`, performs its update, and writes back through `grid`. Current stubs: `apply_vegetation_growth`, `apply_atn_step`, `apply_dispersal`.
- `EcosystemEngine`: The pipeline. Holds references to `EcosystemGridState` and `EnvironmentState`, keeps an ordered list of processes registered through `add_process`, and runs them in order on each `step()`.

## Important principles

- **Single source of truth for state.** No duplication of information in memory across modules: every process reads from and writes to the same initialized `EcosystemGridState` / `EnvironmentState` / `SpeciesRegistry` objects passed by reference. Intermediate quantities a process needs to expose to other processes (e.g. ATN's `net_growth_rate` for dispersal) are added as named layers on the shared grid, not stored locally inside a process.
- **Common interface for processes.** Every process is a plain function with the same `(grid, env)` signature. Registering one is `engine.add_process(fn)`; removing one is deleting that line. Pipeline order is defined by registration order.
- **Vectorised, group-based operations.** Species belong to one or more named groups (`plants`, `herbivores`, `deciduous_trees`, ...). Processes operate on group slices of the state matrices, so a process never needs to loop over individual species.
- **Each process at its own dimensionality, written to broadcast.** Each science function works at the array shape that is natural to its biology — a per-cell process operates elementwise, a within-cell food-web process reduces over species, a between-cell process uses spatial stencils. The discipline is that the *same NumPy code* must run whether it is called on a single cell (e.g. `(S,)`) or on the full grid (e.g. `(X, Y, S)`). Concretely: rely on NumPy broadcasting for elementwise math, reduce on named axes (`axis=-1` for species when applicable) rather than positional ones, and avoid Python loops over cells. This keeps unit tests on hand-built arrays representative of production behaviour and avoids the per-cell loop that currently makes the ATN implementation slow at global scale.
- **Dependency direction points inward.** State objects (the cart) do not know about processes; processes depend on the state interface, not the other way around. Swapping the implementation of a process — or replacing one process with another — should not require touching the state classes.
- **Separation of concerns between the layers.** Processes describe *what happens here and now* in a cell; the engine handles *when* (the time loop); the state objects handle *what is where*. A process should not be in charge of advancing time or orchestrating other processes.
- **Extensibility through layers, not new modules.** When a process needs to publish a new quantity for downstream processes, it registers a new layer on `EcosystemGridState` rather than introducing a new shared object.

## Vincent's notes — concerns and suggested adjustments

### 1. Processes should produce a delta, not mutate biomass directly

In the current `engine_v2`, processes write back into `grid.layers["biomass"]` in-place (see `apply_vegetation_growth` editing `plants += ...` through `edit_group_data`). This means the order in which processes are registered changes the result, even when the underlying biological effects are intended to be commutative within a single time step.

Alexis suggested in response that the state could carry a **delta layer paired with each biomass-like state variable**, which processes would write into instead of mutating biomass directly. The engine would then integrate those deltas into the biomass at the end of the time step:

```
biomass_{t+1} = biomass_t + Σ_p dB_p(state_t, env_t)
```

Concretely, every process accumulates its contribution into the shared `dB` layer (one per group, or one global delta layer indexed by `Species_ID`) rather than into `biomass`. The engine applies and zeroes the delta at the end of each step. This makes within-step computation order-independent (every process sees the same `state_t`), matches how the underlying equations are usually written (a sum of contributions), and makes it trivial to inspect or plot the contribution of each process separately. It also fits naturally with the "extensibility through layers" principle above — the delta is just another registered layer on `EcosystemGridState`.

### 2. Testability — decouple the process function from how the engine integrates it

Right now a process *is* its integration into the engine: it takes `(grid, env)`, fetches what it needs, mutates the shared state, and returns nothing. That hides the scientific signature of the process. A vegetation contributor reading `apply_vegetation_growth` does not see "vegetation depends on NPP, structural fraction, competition, herbivory" — they see grid plumbing.

I'd prefer two layers:

1. **A pure scientific function** with an explicit signature, e.g.
   `vegetation_growth(biomass_plants, npp, competition, params) -> delta_biomass_plants`.
   No `grid`, no `env`, no side effects. Easy to unit-test with hand-built arrays, easy to prototype in a notebook from source, easy to reason about for a non-engineer contributor. It's implementation must be broadcast-friendly, but it does not need to know about the full grid or how to slice it.
2. **A thin adapter** that knows how to pull the right slices out of `grid` / `env` / `registry`, call the pure function, and write the delta back into the correct layer.

The adapter is what the engine registers. The science lives in code that does not depend on the engine at all.

The full proposed process contract specification in [`context/process_contract_spec.md`](tuesday_architecture_proposal/context/process_contract_spec.md), meant to be reusable as the canonical reference for `src/` contributions in future `README.md` and `AGENTS.md` / `CLAUDE.md` files.

### 3. Group-restricted processes

Some processes will only apply to a subset of groups (e.g. herbivory only on herbivores, fire only on plants, possibly different dispersal kernels for birds vs walkers). The current registration is uniform: every registered process runs on every step on the full grid, and it is the process body's job to slice. That works, but it pushes routing logic into each process.

Open question: should the engine support **scoped registration**, e.g. `engine.add_process(fn, groups=["herbivores"])` or `engine.add_process(fn, layer="biomass", groups=["plants"])`, so that the slicing is described at registration time and the pure scientific function only ever sees what it operates on? This would compose well with point 2 (the adapter becomes the scope declaration) and clarify which processes touch which groups when reading the pipeline definition.

This should not be designed in the abstract. Before committing to an engine-level mechanism, we need to clarify with the **vegetation, ATN, and dispersal teams** which processes are actually group-restricted in their formulations, and cross-check against the **literature** (Madingley, ATN papers, dispersal references) to see how they handle group-specific dynamics. The risk otherwise is to bake in a routing abstraction that does not match how the scientific processes are written, and end up paying the abstraction cost without the benefit. The engine should follow the processes' needs here, not the other way around.

### 4. ODE vs discrete-time — what should processes return?

The ATN implementation on the `patn` branch is written as a system of **Ordinary Differential Equations (ODEs)**: it exposes `dB/dt` and relies on `scipy.integrate.odeint` to advance time. The current `engine_v2` is **discrete-time**: each process writes a finite update to biomass for one engine time step. The two styles are not currently compatible, and the engine needs to take a position on which one is the contract.

I am not personally familiar enough with the trade-offs to decide this on my own, especially given the constraints I do care about: global-scale performance, long time series, testability of individual processes, reusability, and scientific transparency for non-engineer contributors.

The proposal, after a Claude discussion on this question, is to make the engine contract a **discrete biomass delta per engine time step**. ODE-style processes (ATN today, possibly others later) integrate their rate to a delta **inside their own adapter, or are rewritten to directly return the delta**. The engine never sees rates; it only sees deltas, which it sums and applies at the end of the step.

The full findings, comparison, and rationale are written up in [`discretization.md`](tuesday_architecture_proposal/discretization.md). The output and shape conventions this implies for every process are specified in [`context/process_contract_spec.md`](tuesday_architecture_proposal/context/process_contract_spec.md).