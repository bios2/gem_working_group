# Vegetation model — `vegetation_model.py`

NPP-driven growth of basal (plant) biomass for the spatial ATN.

This module owns one piece of the science: how much **leaf biomass** each plant
species gains per timestep, given the cell's Net Primary Productivity (NPP) and
the local competition between **herbs** and **trees**. It does *not* handle
consumers, metabolism, or feeding — those live in [atn_model.py](atn_model.py).

The biology and the original Madingley-derived equation are explained in
[context/vegetation.md](context/vegetation.md). This README documents what the
code actually computes and how to drive it.

---

## 1. The equation

For every basal species $i$ in a cell, the growth **rate** (g wet m⁻² day⁻¹) is

$$
G_i = NPP \cdot \psi \cdot (1 - f_{\text{struct},i}) \cdot C_i
$$

| Symbol | Meaning | Source |
|---|---|---|
| $NPP$ | cell Net Primary Productivity | `env_df['NPP']`, one value per cell |
| $\psi$ | carbon → wet-matter conversion factor | `config['psi']` (scalar) |
| $f_{\text{struct},i}$ | fraction of NPP allocated to structural tissue (wood, roots) | `traits['f_struct']`, else `config['f_struct_default']` |
| $C_i$ | herb/tree competitive partition coefficient | computed from biomass (below) |

The factor $(1 - f_{\text{struct},i})$ is the share of production that goes to
**leaves** rather than to structural tissue.

> **Note on units.** The full Harfoot et al. equation in
> [context/vegetation.md](context/vegetation.md) also carries a timestep factor
> $\delta t_{NPP}$ and a cell-area factor $A_{\text{cell}}$. Here $G_i$ is a
> **rate** (per day): the timestep is applied by the integrator's `dt` in
> [atn_model.py](atn_model.py), and $A_{\text{cell}}$ is dropped because NPP is
> supplied as an areal density. So this module returns $dB/dt$ contributions,
> never a finished biomass increment — consistent with the
> [process contract](../../docs/processes_implementation_specification.md).

### The competition partition $C_i$

$B_{\text{trees}}$ is the **total** tree biomass summed over all tree species in
the cell. Herbs and trees split production differently:

**Herbs** — all herb species in a cell share the same coefficient:

$$
C_{\text{herb}} = \frac{\alpha_{\text{herbs}}}{\alpha_{\text{herbs}} + B_{\text{trees}}}
$$

As trees accumulate, herbs get shaded out and $C_{\text{herb}} \to 0$.

**Trees** — each tree species gets a share proportional to its own biomass:

$$
C_{\text{tree},i} = \frac{B_i}{\alpha_{\text{herbs}} + B_{\text{trees}}}
$$

Summing $C_{\text{tree},i}$ over all tree species gives the total tree partition
$C_{\text{trees}} = B_{\text{trees}} / (\alpha_{\text{herbs}} + B_{\text{trees}}) = 1 - C_{\text{herb}}$,
matching the curve in [context/vegetation.md](context/vegetation.md). The per-species
form just divides that total among individual tree species by relative biomass.

| Symbol | Meaning | Source |
|---|---|---|
| $\alpha_{\text{herbs}}$ | half-saturation constant for the herb/tree split (g m⁻²) | `config['alpha_herbs_default']` |
| $B_{\text{trees}}$ | total tree biomass in the cell | summed at runtime |

$\alpha_{\text{herbs}}$ sets where the crossover happens: when
$B_{\text{trees}} = \alpha_{\text{herbs}}$, herbs and trees each get half the production.

---

## 2. Inputs

`PlantVegetationModel(traits_df, env_df, config)` expects:

**`traits_df`** — one row per species:

| Column | Required | Used for |
|---|---|---|
| `is_basal` | yes | `1` = plant (gets growth), `0` = consumer (growth is 0) |
| `vegetation_type` | yes, for basal | `'herb'` or `'tree'` — selects which $C_i$ formula applies |
| `f_struct` | optional | per-species $f_{\text{struct},i}$; NaN rows fall back to the config default |

**`env_df`** — one row per grid cell, must contain an `NPP` column (plus `x`, `y`
for output).

**`config`** — must contain `psi`, `alpha_herbs_default`, and `f_struct_default`
(see [atn_io.py](atn_io.py) `check_parameter_completeness`).

Consumers (`is_basal == 0`) always get $G_i = 0$ — that is what makes the same
`dB/dt` formula in [atn_model.py](atn_model.py) work for plants and animals alike.

---

## 3. The two interfaces

The module exposes the same calculation at two shapes, following the
[process contract](../../docs/processes_implementation_specification.md):

```python
growth(B, cell_idx)        # B: (S,)          -> G: (S,)          one cell
growth_all_cells(B)        # B: (n_cells, S)  -> G: (n_cells, S)  every cell at once
```

`growth_all_cells` is the vectorised path the integrator uses — no Python loop
over cells. `growth` is the single-cell convenience version. Both return a
**rate per species**, zero for consumers.

---

## 4. How it plugs into the ATN

The whole point of this module is the term $G$ in the biomass balance solved by
[atn_model.py](atn_model.py). Inside `derivatives()` the full per-cell, per-species
rate is

$$
\underbrace{\frac{dB}{dt}}_{\text{net change}}
= \underbrace{\text{feeding\_gain}}_{\text{assimilated prey}}
+ \underbrace{\,G\,}_{\textbf{vegetation growth}}
- \underbrace{X \cdot B}_{\text{metabolism}}
- \underbrace{\text{predation\_loss}}_{\text{eaten by others}}
$$

**$G$ is the only production term, and it comes entirely from this module.**
Everything else in the equation is loss or trophic transfer supplied by the ATN.

The coupling is two lines in `derivatives()` ([atn_model.py:197-201](atn_model.py#L197-L201)):

```python
# ---- Vegetation growth (basal species only, zero for consumers) ----
G = self.vegetation.growth_all_cells(B_safe)           # (n_cells, S)

dBdt = feeding_gain + G - X * B_safe - predation_loss  # G = production term
```

Why the same equation works for plants **and** animals:

| Species type | `feeding_gain` | `G` (this module) | net effect |
|---|---|---|---|
| **basal** (plant) | 0 — no prey rows in `adj_mat` | $NPP \cdot \psi \cdot (1-f_{\text{struct}}) \cdot C_i$ | grows from NPP, dies from metabolism + herbivory |
| **consumer** (animal) | assimilated prey | 0 — `is_basal == 0` | grows from feeding, dies from metabolism + predation |

Because $G_i = 0$ for consumers and $\text{feeding\_gain} = 0$ for basal species,
the two roles never overlap and a single vectorised expression covers both — no
`if basal / else consumer` branching in `derivatives()`.

---

## 5. Output

`save_output(B_traj, t_eval, output_dir)` writes `vegetation.txt` (long format,
one row per cell × basal species × timepoint):

```
pixel_id  x  y  time  species  delta_biomass
```

`delta_biomass` is $G_i$ in g m⁻² day⁻¹ — the instantaneous NPP-driven growth
contribution, evaluated along the saved trajectory. It is the production term
only, not the net change in plant biomass.

---

## 6. References

- Equation derivation and biological terms: [context/vegetation.md](context/vegetation.md)
- Original source: Harfoot et al. (2014) *PLoS Biology* Text S1, Terrestrial plant model
- Process contract (shapes, rates vs. deltas, science/adapter split): [docs/processes_implementation_specification.md](../../docs/processes_implementation_specification.md)
- ATN coupling and integration: [readme_atn.md](readme_atn.md), [README_modular_atn.md](README_modular_atn.md)
