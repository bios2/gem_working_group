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

## Project structure



## Python packaging and environment

<!-- Dependancies -->

<!-- .toml -->

<!-- .venv and and getting started with venv and pip install -e . -->


## Style guide and naming conventions

<!-- Style guide, casing for modules, functions and naming -->


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


## Geographic grid

For all experiments and simulations, we will use a common spatial grid to ensure comparability of results. The grid will cover North America with the following specifications:

CRS: `ESRI:102008 North America Albers Equal Area Conic`  
PROJ: `+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs`  
Cell size: `100000 m × 100000 m`  
Grid origin: `x = -7000000 m`, `y = -2000000 m`  
Extent: `x = -7000000..5000000`, `y = -2000000..5500000`  
Cell ID: `NA100_R{row}_C{col}` from upper-left or lower-left origin, documented explicitly  

## Simulation engine

<!-- Describing Alex's engine state management species registry. There should be a doc describing that  -->

The engine drives a simulation forward by composing process functions (see *Model development*) into a pipeline that updates a shared `(X, Y, S)` biomass array over cells and species each time step.

- **State management.** <!-- TODO: EcosystemGridState, EnvironmentState, SpeciesRegistry, shared delta layer -->
- **Adapters (`processes.py`).** Science modules import nothing from the engine. All engine glue lives in a single file, `processes.py`, holding one `apply_<process>(grid, env, dt)` per process. Adapters slice the right state arrays, fetch parameters, call the science function, and write the result back — biomass-modifying adapters accumulate into the shared delta layer, dependency adapters write to a named shared layer.
- **Broadcasting.** <!-- TODO: how grid layers, environment layers, and per-species parameters are reshaped to broadcast against the (X, Y, S) biomass array -->
- **Initialization.** <!-- TODO: building the grid, loading initial conditions, species traits, environmental layers -->
- **Pipeline registration.** Each adapter is registered with the engine in pipeline order. Dependency-process adapters must run before any adapter that consumes their output.


## Running simulations

Inside experiments folder. Store relevant data.

---
<!-- 
- How to handle and store simulation runs (notebooks ?,  saved outputs ?) and how to make them accessible to the team. Make minimal requirements for reproducibility of runs (e.g. saving the random seed, saving the configuration file, etc.). Naming convention for runs and outputs with date and time and description.
- Input data. Script to download and preprocess input data (e.g. environmental data, species traits, etc.) and store them in a standardized format that can be easily accessed by the engine and processes. We recommend using geotiff to store spatial data and csv or json for tabular data. We also recommend using a standardized directory structure for input data (e.g. data/raw, data/processed, etc.) and a naming convention for files (e.g. data/raw/environmental_data_2024-06-01.tif). Local gitignores for large .tiff datasets that are not stored in the repository but can be downloaded and processed by the script.
- Initialization scripts for the engine (e.g. to set up the grid, load initial conditions, species list and traits, load the input data). To be described in engine section. Should be modular and reusable for different runs and configurations. Should also include error handling and logging to facilitate debugging and tracking of runs.

 -->