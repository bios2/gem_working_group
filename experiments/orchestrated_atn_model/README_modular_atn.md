# Modular ATN — usage guide

This folder provides a **modular version** of the ATN model, split across three Python files instead of a single one. The goal is to cleanly separate **the scientific processes**, **the input data**, and **the execution**.

> If you are looking for the original monolithic version, it lives in [atn_model.py](atn_model.py) (run via [run_atn.py](run_atn.py)).
> The modular version coexists alongside, in the files described below.

---

## 1. Overview

```
┌──────────────────────┐        ┌──────────────────────┐
│  atn_processes.py    │        │   atn_data.py        │
│  (the SCIENCE)       │        │   (the INPUTS)       │
│                      │        │                      │
│  - temperature_factor│        │  - load_inputs       │
│  - attack_rate_matrix│        │  - build_atn_params  │
│  - functional_response       │  - initial_biomass   │
│  - feeding_gain      │        │                      │
│  - predation_loss    │        │  (reads config/env/  │
│  - derivatives       │        │   adj/traits, pre-   │
│                      │        │   computes base_a,   │
│  Pure functions:     │        │   base_h, base_X,    │
│  numpy only,         │        │   E, T_K, …)         │
│  no state, no I/O    │        │                      │
└──────────┬───────────┘        └──────────┬───────────┘
           │                               │
           │      both imported by                  │
           ▼                               ▼
       ┌────────────────────────────────────────┐
       │       run_modular_atn.py               │
       │       (the ORCHESTRATION)              │
       │                                        │
       │  1) load_inputs()                      │
       │  2) build_atn_params()                 │
       │  3) initial_biomass()                  │
       │  4) RK4 loop on derivatives()          │
       │  5) save outputs                       │
       └────────────────────────────────────────┘
```

**Why this split?**

- **Test one process in isolation:** you can call `functional_response(B, a, h, q, c, eps)` with tiny hand-built arrays and check the result, without loading a single file.
- **Change an input without touching the science:** if the traits format changes tomorrow, only `atn_data.py` moves.
- **Change a solver without touching the science:** if RK4 gets swapped for `scipy.solve_ivp`, only `run_modular_atn.py` moves.
- **Read the science at a glance:** `atn_processes.py` is ~100 lines, readable end-to-end, with no `self.` clutter.

---

## 2. The three files in detail

### 2.1 `atn_processes.py` — the science

The heart of the model. Every function is **pure**:

- depends only on `numpy`;
- reads no files, writes no files;
- contains **no `self.`** (no class).

Each function receives every input explicitly.

**Catalogue:**

| Function | Computes | Output shape |
|---|---|---|
| `temperature_factor(T_K, T0_K, E_a, k_B)` | Boltzmann–Arrhenius factor | `(n_cells,)` |
| `attack_rate_matrix(M, adj, …)` | allometric `a_ij` | `(S, S)` |
| `handling_time_matrix(M, adj, …)` | allometric `h_ij` | `(S, S)` |
| `metabolic_base_rate(X0, M, bX)` | `X_i = X0_i · M_i^{bX_i}` | `(S,)` |
| `assimilation_matrix(is_basal, e_plant, e_animal)` | `E[prey, predator]` | `(S, S)` |
| `functional_response(B, a, h, q, c, ε)` | generalised Holling type II | `(n_cells, S, S)` |
| `feeding_gain(B, F, E)` | predator energetic gain | `(n_cells, S)` |
| `predation_loss(B, F)` | prey loss | `(n_cells, S)` |
| `derivatives(B, t, params, vegetation)` | full `dB/dt` | `(n_cells, S)` |

**Rule of thumb:** if a function in this file needs to read a `.txt` or `.csv`, that's a bug.

### 2.2 `atn_data.py` — the inputs

This file bridges the text files on disk and the numpy arrays the science consumes. It has **three responsibilities**:

1. **`load_inputs(env_file, adj_file, traits_file, config_file)`**
   Reads the four files via [atn_io.py](atn_io.py), applies the cross-file validation, returns `(config, env_df, adj_mat, traits_df)`.

2. **`build_atn_params(adj_mat, traits_df, env_df, config)`**
   Precomputes **everything that is constant in time**: the allometric matrices `base_a`, `base_h`, the `base_X` vector, the assimilation matrix `E`, the per-cell temperature vector `T_K`, and the basal/consumer indices. The result is a large `dict` that is passed to `derivatives()` without being rebuilt at every time step.

3. **`initial_biomass(traits_df, n_cells, noise_frac=0.01)`**
   Builds the `B_initial` matrix of shape `(n_cells, n_species)` from the `initial_biomass_g_per_m2` column, with a relative gaussian noise to break inter-cell symmetry.

### 2.3 `run_modular_atn.py` — execution

This file **contains no equation**. It chains:

```python
config, env_df, adj_mat, traits_df = load_inputs(...)        # 1. inputs
params, vegetation = build_atn_params(...)                   # 2. precompute
B_initial = initial_biomass(traits_df, params['n_cells'])    # 3. initial conditions
B_traj    = run_simulation(B_initial, t_eval, params, vegetation)  # 4. RK4 integration
# 5. save: simulation_summary.txt, biomass.txt, vegetation.txt, atn_model.txt
```

It also defines `rk4_step()` (one RK4 step vectorised over all cells) and three output-writing helpers.

---

## 3. How to run it

### 3.1 Command line (typical case)

From the `test_data/` folder that holds the input files:

```bash
cd gem_working_group/experiments/orchestrated_atn_model/test_data
python ../run_modular_atn.py env_mat.txt adj_mat.txt traits.txt ../config.txt --t_max 100 --seed 42
```

Expected output (abridged):

```
======================================================================
ATN MODULAR
======================================================================

[1] Reading and validating inputs...
[2] Building ATN parameters...
ATN params built: 5 species, 25 cells.
[3] Initial conditions...
    ✓ Initial biomass: 1.23e-02 – 1.50e+01 g/m²
    ✓ Total biomass  : 9.45e+02 g/m²

[4] RK4 integration over 100 days (101 points)...
  t = 100.0 / 100.0 days
✓ 100 RK4 steps across 25 cells.

[5] Saving to atn_output/20260527H143215/
    ✓ simulation_summary.txt
    ✓ biomass.txt
    ✓ vegetation.txt
    ✓ atn_model.txt

[6] Final statistics:
    Species 0 (basal)       B_final=8.412e+00  persistence=100.0%
    Species 1               B_final=4.211e-01  persistence= 92.0%
    ...
======================================================================
✓ SIMULATION COMPLETE
======================================================================
```

### 3.2 From a notebook or a script

If you want to drive the model from a Jupyter notebook for analysis:

```python
from run_modular_atn import main

B_traj, t_eval, params, vegetation = main(
    env_file='test_data/env_mat.txt',
    adj_file='test_data/adj_mat.txt',
    traits_file='test_data/traits.txt',
    config_file='config.txt',
    t_max=200.0,
    seed=42,
)

# B_traj.shape == (n_timepoints, n_cells, n_species)
# Plot species 0 dynamics in cell 5:
import matplotlib.pyplot as plt
plt.plot(t_eval, B_traj[:, 5, 0])
plt.xlabel('Days'); plt.ylabel('Biomass (g/m²)')
```

### 3.3 Calling processes without running anything (tests, exploration)

The whole point of the split is that you can use any process function in isolation, **without loading any file**:

```python
import numpy as np
from atn_processes import functional_response, attack_rate_matrix

# Build a small toy case by hand: 2 species, 3 cells
M       = np.array([0.1, 10.0])                  # body masses
adj_mat = np.array([[0, 1], [0, 0]])             # species 1 eats species 0
B       = np.array([[1.0, 0.1],
                    [2.0, 0.2],
                    [0.5, 0.05]])                # (3 cells, 2 species)

a = attack_rate_matrix(M, adj_mat, a0=0.001, b_a_prey=-0.5, b_a_pred=0.5)
F, extinct = functional_response(
    B,
    a_ij=a[np.newaxis, :, :],                    # broadcast over cells
    h_ij=np.zeros_like(a)[np.newaxis, :, :],
    q_hill=2.0, interference=0.0, ext_threshold=1e-6,
)
print(F.shape)       # (3, 2, 2)
print(extinct)       # (3, 2) boolean
```

That's the basis of a unit test: compare `F` against a hand-computed value.

---

## 4. How to modify the model

Depending on what you want to change, you know exactly which file to open:

| You want to… | You edit… |
|---|---|
| Change the functional-response formula | `atn_processes.py` → `functional_response` |
| Add a new process (e.g. fire-driven mortality) | new function in `atn_processes.py`, called from `derivatives` |
| Read a new traits-file format | `atn_data.py` → `load_inputs` or `_per_species_param` |
| Precompute a new constant quantity | `atn_data.py` → `build_atn_params`, add a key to `params` |
| Swap RK4 for another solver | `run_modular_atn.py` → `rk4_step` / `run_simulation` |
| Add a diagnostic output | `run_modular_atn.py` → new `save_*` helper, called from `main` |
| Change a default value (E_a, q_hill, …) | [config.txt](config.txt) (no code change) |

None of these changes should require touching all three files at once. If it does, that's a sign the split needs revisiting.

---

## 5. Testing a process

The main benefit of the split: every function in `atn_processes.py` can be tested in a few lines.

```python
# tests/test_atn_processes.py
import numpy as np
from atn_processes import metabolic_base_rate

def test_metabolic_base_rate_scaling():
    # Doubling the mass at exponent -0.25 multiplies X by 2^(-0.25) ≈ 0.841
    X0 = np.array([0.5, 0.5])
    M  = np.array([1.0, 2.0])
    bX = np.array([-0.25, -0.25])
    X  = metabolic_base_rate(X0, M, bX)
    assert np.isclose(X[0], 0.5)
    assert np.isclose(X[1], 0.5 * 2 ** -0.25)
```

Run with:

```bash
cd gem_working_group/experiments/orchestrated_atn_model
pytest tests/
```

---

## 6. Caveats

- **The files `atn_model.py`, `run_atn.py` and `vegetation_model.py` still contain unresolved merge conflicts** (`<<<<<<< HEAD` / `=======` / `>>>>>>>`). They must be resolved, otherwise `import vegetation_model` will fail and the modular version will not run.
- **Scientific verification:** before replacing the old version, run a simulation with the same seed and the same inputs in both versions (`run_atn.py` and `run_modular_atn.py`) and compare `biomass.txt` line by line. The values should match to rounding error.
- **The vegetation sub-model** (`PlantVegetationModel`) keeps its class form so we don't have to rewrite its interface. It could be refactored the same way later if you want to push the consistency further.
