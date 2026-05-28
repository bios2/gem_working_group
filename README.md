# GEM Working Group

This repository supports the GEM working group focused on designing and prototyping a new General Ecosystem Model.

The project is inspired by Madingley, but the goal is not to reproduce it directly. The working group is exploring a new model structure that can better represent biodiversity change, species-level dynamics, trophic interactions, spatial structure, dispersal, vegetation, coexistence, and ecosystem processes.

A second goal is to document how AI coding agents can support collaborative ecological modelling, while keeping scientific assumptions, model structure, tests, and validation explicit.

## Objectives

- Design a new process-based ecosystem model.
- Compare the effects of biodiversity change and climate change on ecosystem functioning.
- Explore species-level or hybrid species/guild representations.
- Develop modular model components for trophic dynamics, vegetation, dispersal, mortality, reproduction, and spatial processes.
- Define ecological tests, diagnostics, and validation targets.
- Document lessons from AI-assisted collaborative model development.

## Resources and documentation

[📂 SharePoint Documents](https://usherbrooke.sharepoint.com/sites/ielabworkinggroup)

[💬 Teams Chat](https://teams.microsoft.com/l/team/19%3A5UwBzgYI52ESK9znTZVksNqlEEprrzf8AeMWIcRYWpU1%40thread.tacv2/conversations?groupId=f8341ee7-d260-4f7a-9890-6f13bbc5ec80&tenantId=3a5a8744-5935-45f9-9423-b32c3a5de082)

[📅 Logistics](https://usherbrooke.sharepoint.com/:f:/r/sites/ielabworkinggroup/Documents%20partages/logistics?csf=1&web=1&e=nsdYNf)

## General Collaboration guidelines

- This repo is where all notes, papers, and documents for the duration of the row group will live
- The `main` branch is protected, so please create a new branch for any changes and submit a pull request.

## AI coding assistants

[AGENTS.md](AGENTS.md) (mirrored as [CLAUDE.md](CLAUDE.md)) tells AI coding agents — Claude Code, Codex, Cursor, … — how to help on this project: how to communicate with you, what conventions to enforce, when to push back proactively, and what they must not do (no inventing process mechanics, no unfamiliar abstractions, no commits or pushes without asking). Read it before you let an agent touch the repo, so you know the behaviour to expect and what to correct. If you disagree with anything in there, open a PR against `AGENTS.md` — it is the source of truth, and `CLAUDE.md` is kept in sync from it.

## Project structure

```
gem-working-group/
├── data/         # Input data files (large rasters gitignored); see Input data files
├── docs/         # Reference documents — contracts, specifications, design notes
├── experiments/  # Prototypes and experiments; one subfolder per experiment
├── papers/       # Reference papers and bibliography
├── src/          # The model code package (processes + engine), added when implementation starts
└── README.md
```

Everyday work happens in `experiments/` while the model is being prototyped. Once a process or engine component is stable, it migrates into `src/` so it can be imported by every experiment.

## Python packaging and environment

We use Python. The dependencies actually used across the prototype branches are:

- `numpy` — numerical arrays, the backbone of every process.
- `scipy` — numerical integration for ODE-style processes (`scipy.integrate.solve_ivp`, `odeint`).
- `pandas` — tabular data handling (species traits, parameter tables, diagnostics).
- `matplotlib` — plotting for experiments and figure reproduction.
- `rasterio` — reading and writing GeoTIFF spatial inputs.
- `pytest` — the recommended test runner for unit tests (see *Model development*).

Additional libraries (e.g. `xarray`, `geopandas`) are added only when a process actually needs them.

Dependencies and packaging metadata live in a single `pyproject.toml` at the repo root. No `setup.py`, no `requirements.txt` — `pyproject.toml` is the only source of truth, so there is nothing to keep in sync by hand.

The package lives at [src/gem/](src/gem/) and imports as `import gem` (e.g. `from gem.vegetation import logistic_growth_delta`). The repo directory `gem-working-group/` is the *workspace*, not the package — keep the distinction in mind when writing imports.

**Getting started.** From the repo root:

```bash
python -m venv .venv             # create an isolated environment
source .venv/bin/activate        # macOS / Linux
.venv\Scripts\activate           # Windows
pip install -e .                 # install the project and its dependencies
pytest                           # confirm everything works
```

`pip install -e .[dev]` installs the project in **editable mode**: changes you make to the source code take effect immediately without reinstalling. The `.venv/` folder is local to your machine and is gitignored — never commit it. Each team member recreates it from `pyproject.toml`. [dev] is an optional extra that includes development dependencies like `pytest`.

VS Code detects automatically a virtual environment in the project folder and recommends you to activate it when you start working with a virtual environment in your workspace. 

## Style guide and naming conventions

We follow standard Python conventions (PEP 8) so the code is recognisable to anyone who has read a Python tutorial:

- **Modules and packages** (folders and `.py` files): lowercase with underscores — `vegetation.py`, `species_registry.py`. No CamelCase, no hyphens.
- **Functions and variables**: lowercase with underscores — `logistic_growth_delta`, `body_mass`.
- **Constants**: uppercase with underscores — `EXTINCTION_THRESHOLD`, `T0_K`.
- **Classes**: CamelCase — `EcosystemGridState`, `EnvironmentState`.

Names should describe **what** the thing is, not how it is implemented: prefer `metabolic_rate` over `m`, `carrying_capacity` over `K_arr`. Single-letter names are fine inside short math expressions where the meaning is local (`B`, `r`, `K` matching the equation you are encoding) but not in module-level APIs.

Branch names: lowercase with hyphens — `vegetation-logistic-growth`, `dispersal-density-dependent`. One feature per branch.


## Experiments and prototyping

The `experiments` folder is where we will develop and share code for model prototyping, testing, and experiments. Each experiment should have its own subfolder with a README describing the purpose, methods, and results.

Naming convention for experiment folders: `DAY_GROUPNAME_experimentNAME`. Example - `sunday_atn_bylot_experiment1`.

We recommend using Jupyter notebooks for prototyping and documentation, but feel free to use other formats as needed. The key is to keep everything organized and well-documented for future reference.

## Model development

The model evolves the ecosystem by composing modular **processes** — vegetation growth, ATN trophic dynamics, dispersal, metabolism, fire, and so on. Each process is implemented as a pure numpy function with a typed signature.

The full contract is specified in [docs/processes_implementation_specification.md](docs/processes_implementation_specification.md). What every contributor needs to know:

- **Two process categories.** Biomass-modifying processes (vegetation, ATN, dispersal) return a `biomass_delta` array — the finite change in biomass over one time step `dt`. Dependency processes (metabolism, NPP, ...) return a shared intermediate quantity (rate, flux, factor) consumed by *multiple* biomass-modifying processes; they take no `dt`.
- **Numpy at the process's natural dimensionality.** The same science code runs on a single cell, a row, or the full `(X, Y, S)` grid via standard broadcasting. Per-cell python loops are the main performance trap and are not acceptable.
- **Typed signatures, runtime shape asserts.** All arrays are `NDArray[np.float64]`; scalars are `float`. A one-line shape `assert` at the top of each science function catches broadcast mismatches early.
- **ODE-style processes** (ATN and any other rate-based formulation) pick an integration strategy explicitly — either integrate `dB/dt` to a delta (`scipy.solve_ivp` or a vectorised RK4 step) or rewrite the science to return `biomass_delta` directly. A continuous-time rate is never the public output of a process.

Adding a new process is the same recipe every time: write a typed science function in its own module (`vegetation.py`, `atn.py`, `dispersal.py`, `metabolism.py`, ...) that imports nothing but `numpy`, add a runtime shape assert, and add a unit test against hand-built arrays. How the function is wired into a running simulation is the engine's concern, covered below.

## Input data files

Input data lives in `data/` with two subfolders:

- `data/raw/`: untouched inputs as downloaded from their source (climate reanalyses, species traits, occurrences, ...). Never edit these by hand.
- `data/processed/`: inputs reprojected onto the engine grid (see *Geographic grid*), cleaned, or otherwise prepared for the engine to consume.

**Formats.** GeoTIFF (`.tif`) for spatial raster data (temperature maps, carrying capacity, ...). CSV or JSON for tabular data (species traits, parameter sets, ...). Avoid project-specific binary formats — they make data hard to inspect outside the engine.

**File naming.** Include the date and a short descriptor: `data/raw/era5_temperature_2024-06-01.tif`, `data/processed/species_traits_bylot_v2.csv`. This keeps reruns reproducible because the filename itself records which version was used.

**Large files.** Large rasters are gitignored. Each `data/raw/` subfolder should include a small script (`download.py` or `download.sh`) that reproduces the download, so the team can rebuild the dataset without committing gigabytes to git.

**Reprojection** of raw inputs onto the engine's projected grid happens in `data/processed/`, not inside the engine or in processes. The engine consumes already-projected data.

## Geographic grid

For all experiments and simulations, we will use a common spatial grid to ensure comparability of results. The grid will cover North America with the following specifications:

CRS: `ESRI:102008 North America Albers Equal Area Conic`  
PROJ: `+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs`  
Cell size: `100000 m × 100000 m`  
Grid origin: `x = -7000000 m`, `y = -2000000 m`  
Extent: `x = -7000000..5000000`, `y = -2000000..5500000`  
Cell ID: `NA100_R{row}_C{col}` from upper-left or lower-left origin, documented explicitly  

## Simulation engine

The engine is the runtime glue around the modular processes described in *Model development*. It is intentionally small — most of the action is in the science modules — and is built around three shared state objects ("the cart") that travel together through every process, plus a pipeline that runs the processes in order on each time step.

- **State management.** Three objects hold all of the model's state, and they are the **single source of truth** — processes read from them and write back to them rather than keeping their own copies:
  - `EcosystemGridState`: the dynamic `(X, Y, Species)` biological state. Holds a `biomass` layer by default, plus any named layers processes register — the shared `biomass_delta` layer that biomass-modifying processes accumulate into, dependency outputs like `metabolic_rate`, and so on.
  - `EnvironmentState`: the `(X, Y)` environmental layers (temperature, carrying capacity, ...). Layers are added by name and shape-checked against the grid.
  - `SpeciesRegistry`: the species list, the functional or trophic groups they belong to (`plants`, `herbivores`, ...), per-species traits stored as 1D arrays for fast vectorised math, and the feeding adjacency matrix.
- **Adapters (`processes.py`).** Science modules import nothing from the engine. All engine glue lives in a single file, `processes.py`, holding one `apply_<process>(grid, env, dt)` per process. An adapter slices the right state arrays, fetches parameters, calls the science function, and writes the result back — biomass-modifying adapters accumulate into the shared delta layer; dependency adapters write to a named shared layer.
- **Broadcasting.** Every process operates on numpy arrays at its natural shape — `(S,)` for a single cell, `(Y, S)` for a row, `(X, Y, S)` for the full grid. Adapters reshape per-species parameters to `(1, 1, S)` and environmental layers to `(X, Y, 1)` so they broadcast cleanly against the `(X, Y, S)` biomass array. The same science code then runs unchanged from a unit test on a `(3,)` array to a global simulation.
- **Initialization.** Setting up a run means building the three state objects: define the grid extent and projection (see *Geographic grid*), load environmental layers from `data/processed/`, load the species list and traits, set initial biomass. Initialization is a regular Python function — not a config file — so it can be parametrised and reused across experiments. Each experiment owns its own initialization script.
- **Pipeline registration.** Each adapter is registered with the engine in pipeline order via `engine.add_process(apply_metabolism)`. Dependency-process adapters must be registered before any adapter that consumes their output. Removing a process is deleting one line.

The end-of-step integration applies the accumulated `biomass_delta` to the biomass layer and zeroes the delta. This makes within-step computation order-independent (every process sees the same `state_t`) and matches how the underlying equations are usually written: a sum of contributions.

## Running simulations

Simulations live in `experiments/`, one subfolder per experiment (see *Experiments and prototyping*). Each experiment runs the engine for a defined set of conditions and stores its outputs locally.

**Minimum reproducibility checklist.** Every run should record:

- The **random seed** used (if any stochastic process is involved).
- The **configuration**: initial state, environmental data files, species list, parameter values, number of time steps, `dt`. A small JSON or YAML file alongside the outputs is sufficient.
- The **engine version**: the git commit hash the run was produced from.

Without these three, a result cannot be reproduced.

**Output storage.** Save outputs inside the experiment folder under a timestamped run name, e.g. `experiments/sunday_atn_bylot_experiment1/runs/2026-05-27_baseline/`. The folder should contain the biomass trajectory (NetCDF or `.npz` for arrays, CSV for tabular diagnostics), the configuration file, and any plots. Large outputs are gitignored; commit only what is small and informative (configuration, diagnostics, key plots).

**Sharing runs.** Notebooks are useful for exploring outputs but should not be the canonical record. Keep the data files self-describing (column names, units in metadata) so a teammate can reopen them without needing your notebook.