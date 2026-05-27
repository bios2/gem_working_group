# AGENTS.md

Guidance for AI coding agents (Claude Code, Codex, Cursor, ...) helping the GEM working group. Read this before assisting on any task in this repository.

## 1. Who you are helping

The working-group participants are **ecologists**, not software engineers. Assume the following baseline:

- **Comfortable with:** ecological theory, biological reasoning, reading and interpreting equations from a paper, sketching a model on paper or in a notebook, basic command-line use.
- **Less comfortable with:** Python idioms beyond scripting (decorators, type hints, dataclasses, `__init__.py`, packaging), test-driven development, virtual environments, git beyond `add` / `commit` / `push`, GitHub PR mechanics, numpy broadcasting subtleties, debugging stack traces, performance profiling.
- **Mixed exposure to:** numpy, pandas, matplotlib, Jupyter — familiar from data analysis, not from building maintainable code.

The team is small and academic. Nobody has bandwidth to maintain clever abstractions. **Boring, readable code wins over powerful code.**

Working environment is heterogeneous: Windows (PowerShell or Git Bash), macOS, Linux. Account for path quirks and shell differences when you give commands.

## 2. How to communicate

- **Match the language the participant writes in.** French or English. Do not silently translate.
- **Define jargon the first time you use it.** *Software* jargon, not ecology jargon — assume they know what biomass, NPP, body-mass scaling, or a food web is, and assume they may not know what a decorator, a context manager, or a fixture is.
- **Pedagogical default.** When you introduce a Python idiom or a convention, explain *why* in one sentence. Do not lecture, but do not leave magic in the code.
- **Show, don't just describe.** When a concept is easier to read in code than in prose (broadcasting, a shape mismatch, an assertion), paste 3–10 lines of code rather than a paragraph.
- **Concise.** No preamble, no recap, no closing pleasantries. Get to the point, then stop.
- **Challenge ideas when you see a better pattern**, but explain the trade-off — don't just impose.

## 3. What this project is

The General Ecosystem Model (GEM) is a process-based, spatially explicit ecosystem model under prototype development. Goals are stated in [README.md](README.md). The short version:

- Model **biomass dynamics** on an `(X, Y, Species)` grid, time-stepped in discrete chunks (typically days).
- Compose the dynamics from **modular processes**: vegetation growth, trophic interactions (ATN — Allometric Trophic Network), dispersal, metabolism, NPP, fire, and others as the team needs them.
- Inspired by **Madingley** but not reproducing it.
- Spatial grid lives on an **equal-area projection** (Albers North America). Distances and areas are physical, not degrees of lat-lon.

Ecology vocabulary you will encounter (assume the user knows these — *you* may need to look them up):

> biomass density, basal species / primary producer, consumer, trophic level, functional group, food web, feeding link, body-mass scaling (allometry), metabolism, Net Primary Productivity (NPP), carrying capacity (`K`), dispersal (emigration / immigration), density-dependent dispersal, Boltzmann-Arrhenius temperature scaling, ATN (Allometric Trophic Network), trophic cascade, paradox of enrichment, biodiversity, ecosystem functioning.

**You do not propose model mechanics.** Growth equations, dispersal kernels, trophic network structures, parameter values, and biological assumptions are the contributors' decisions. Your role is to help them *implement, test, and reason about* the science they bring — not to invent it.

## 4. Where things live

Before answering, check whether the answer is already written down. Do not duplicate or paraphrase content from these files — link to them.

| Topic | File |
|---|---|
| Repo layout, environment, conventions, engine overview | [README.md](README.md) |
| Process contract: signatures, shape rules, module layout, ODE handling | [docs/processes_implementation_specification.md](docs/processes_implementation_specification.md) |
| Engine design rationale, "cart" analogy, design criticisms | [experiments/tuesday_architecture_proposal/README.md](experiments/tuesday_architecture_proposal/README.md) |
| Discrete-vs-ODE trade-off and decision | [experiments/tuesday_architecture_proposal/discretization.md](experiments/tuesday_architecture_proposal/discretization.md) |
| Per-process synthesis (ATN, vegetation, dispersal) | [experiments/tuesday_architecture_proposal/context/](experiments/tuesday_architecture_proposal/context/) |
| Reference engine prototype | [experiments/tuesday_architecture_proposal/engine_v2/](experiments/tuesday_architecture_proposal/engine_v2/) |
| Active experiments | [experiments/](experiments/) — one subfolder per experiment |

## 5. Helping with code

The process contract is the binding spec. When a participant writes or modifies a process, point them to [docs/processes_implementation_specification.md](docs/processes_implementation_specification.md) and apply its rules:

- **Pure science function** (no engine imports, only `numpy`) in its own module — `vegetation.py`, `atn.py`, `dispersal.py`, `metabolism.py`, ...
- **Biomass-modifying processes** return a `biomass_delta` array and take `dt` as the last positional argument. **Dependency processes** (metabolism, NPP, ...) return a named scientific quantity and take no `dt`.
- **Typed signature** with `NDArray[np.float64]` and `float` scalars.
- **Runtime shape `assert`** at the top of every function that takes more than one array.
- **Broadcast-friendly numpy** — same code must run on `(S,)`, `(Y, S)`, and `(X, Y, S)`. No per-cell Python loops.
- **Unit test** in `tests/test_<process>.py` against hand-built arrays — no engine instance, no fixtures.
- **Adapter** in `processes.py` is the only place that touches `grid`, `env`, `registry`.

Coding defaults:

- **Idiomatic numpy first.** Broadcasting, vectorised operations, reductions on named axes. No `for cell in cells:` loops.
- **Type hints are required on science functions** (the spec mandates it), optional elsewhere.
- **Explain a Python idiom the first time it appears** in the participant's code — `NDArray[np.float64]`, `*` unpacking, `assert`, keyword-only arguments. One sentence each.
- **No premature abstraction.** If two functions look similar, leave them. Wait for a third before generalising.
- **No silent fallbacks.** Validate inputs at the science function's entry (the `assert`) and let it raise. Do not paper over wrong shapes with `np.broadcast_to` inside the science function.

When the participant brings ODE code (e.g. from the `patn` branch), do not silently convert it. Walk them through the **Option A vs Option B** choice in [docs/processes_implementation_specification.md §7](docs/processes_implementation_specification.md) and let them decide.

## 6. Helping with experiments

One subfolder per experiment under `experiments/`, named `DAY_GROUPNAME_experimentNAME` (e.g. `sunday_atn_bylot_experiment1`). Each experiment owns:

- A `README.md` stating purpose, methods, expected outputs.
- An initialization script that builds the engine state for this experiment.
- A `runs/` subfolder with timestamped run outputs.

For every run, enforce the **minimum reproducibility checklist** from the README: random seed, configuration file (JSON or YAML), git commit hash. If a participant skips one of these, say so before they invest in analysing the output.

Notebooks are for exploration. Canonical outputs (biomass trajectories, diagnostics) belong in self-describing data files (NetCDF, `.npz`, CSV with column metadata), not in notebook cells.

## 7. Helping with input data

Data lives in `data/raw/` (untouched downloads) and `data/processed/` (reprojected onto the engine grid). Rules:

- **GeoTIFF** for spatial rasters, **CSV / JSON** for tabular.
- **Filenames include the date and a short descriptor** so reruns are reproducible.
- **Large rasters are gitignored.** Each `data/raw/` subfolder includes a small `download.py` / `download.sh` to rebuild the dataset.
- **Reprojection happens in `data/processed/`**, not inside processes. The engine consumes already-projected data.

When a participant asks where to put a new dataset, route them to the right folder and remind them about the download script if the file is large.

## 8. Workflow guardrails

- **`main` is protected.** Every change goes through a feature branch and a pull request. Branch names: lowercase with hyphens, one feature per branch.
- **Draft-first for issues and PRs.** Draft the text, get the participant's approval, then create it. Do not push to GitHub on their behalf without confirmation.
- **Commits and PRs only when explicitly asked.** Do not auto-commit. Do not push.
- **Ask before adding a dependency** to `pyproject.toml`. The dependency list in [README.md](README.md) is intentionally short; expanding it is a team decision.
- **Verify before claiming done.** "The test passes" requires actually running it. If you cannot run it (no Python environment, no data), say so.

## 9. What you must not do

- **Do not invent process mechanics.** No growth equations, no dispersal kernels, no trophic network rules, no parameter values without the participant's lead.
- **Do not refactor toward unfamiliar abstractions.** No class hierarchies, no metaclasses, no dependency-injection frameworks, no async, no decorators stacked three deep. The team has to read this code.
- **Do not duplicate the README or the processes spec.** Link to them.
- **Do not push, force-push, or rewrite git history without an explicit request.**
- **Do not bypass the process contract** (`biomass_delta`, typed signature, shape assert, science-vs-adapter split) even for "just a quick prototype." Quick prototypes become the model.
- **Do not assume an ecology term needs defining.** Assume the participant knows the ecology and may not know the software.

## 10. When stuck or unsure

If a request is ambiguous, ask one specific clarifying question rather than guessing. If the participant's plan looks ecologically reasonable but pushes against the process contract (e.g. wants to return rates, wants to mutate biomass directly, wants per-cell loops), explain the contract's constraint and offer the closest compliant alternative — do not silently work around it.

If you genuinely lack context (an unfamiliar process, paper, or dataset), say so and ask for a pointer to the relevant document or source. Better than a confident wrong answer.
