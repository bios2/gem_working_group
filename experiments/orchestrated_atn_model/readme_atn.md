# Spatially Explicit ATN Model

Efficient Python implementation of the Allometric Trophic Network (ATN) model 
with spatial heterogeneity (Section 8, unscaled model from `ATN_model_spatiotemporal_formulas_parameters.Rmd`).

Basal species growth is driven by externally supplied NPP following the vegetation equation
from Harfoot et al. (2014) Text S1, with a Michaelis-Menten competitive partition between
herbs and trees (see `context/vegetation.md`).

## Scripts

| Script | Role | Run directly |
|---|---|---|
| `run_atn.py` | Main entry point: reads inputs, validates, integrates ODEs, saves results | Yes |
| `atn_model.py` | `ATNModel` class: all ODE dynamics | No (imported by other scripts) |
| `vegetation_model.py` | `PlantVegetationModel` class: NPP-driven basal growth and herb/tree bookkeeping | No (imported by `atn_model.py`) |
| `atn_io.py` | File reading and 20+ validation checks | No (imported by other scripts) |
| `config.txt` | Editable text file of all model parameters (read at runtime by `atn_io.read_config`) | No (read by `run_atn.py`) |

## Running with your own data

### Step 1 — Prepare three input files

| File | Format | What it contains |
|---|---|---|
| `env_mat.txt` | CSV | One row per spatial cell: `pixel_id`, `x`, `y`, `temperature_K`, `NPP` |
| `adj_mat.txt` | Space- or comma-separated | Square binary matrix (rows = resources, columns = consumers) |
| `traits.txt` | CSV | One row per species: `species_id`, `body_mass_g`, `is_basal`, `initial_biomass_g_per_m2`, `vegetation_type` |

See the **Input Files** section below for column details and examples.

### Step 2 — Run the simulation

```bash
python run_atn.py env_mat.txt adj_mat.txt traits.txt
```

Optional arguments (pass via the Python API):

```python
from run_atn import main

B_traj, t_eval, model = main(
    'env_mat.txt',
    'adj_mat.txt',
    'traits.txt',
    t_max=500,         # simulation length in days (default: 100)
)
```

### Step 3 — Check outputs

Results are saved to `atn_output/yyyymmddhhmmss/` (folder named with the run timestamp):

```
atn_output/
└── 20260525143012/
    ├── simulation_summary.txt   ← species traits, grid info, all model constants
    └── biomass.txt              ← long-format table (one row per pixel × time × species)
```

`biomass.txt` columns: `pixel_id`, `x`, `y`, `time_step`, `species_id`, `biomass`

The console prints per-species final biomass and the fraction of cells where each species persists.

## Model Overview

```
┌─────────────────────────────────────────────────────────────────┐
│              SPATIALLY EXPLICIT ATN MODEL ARCHITECTURE          │
└─────────────────────────────────────────────────────────────────┘

INPUT FILES:
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│ env_mat.txt  │  │ adj_mat.txt  │  │ traits.txt   │
│              │  │              │  │              │
│ Temp, NPP    │  │ Food-web     │  │ Body mass,   │
│ per cell     │  │ links        │  │ veg. type    │
└──────┬───────┘  └──────┬───────┘  └──────┬───────┘
       │                 │                 │
       └─────────────────┼─────────────────┘
                         │
                    ┌────▼────┐
                    │VALIDATE │ (20+ sanity checks)
                    └────┬────┘
                         │
        ┌────────────────┼────────────────┐
        │                │                │
   ┌────▼──────┐  ┌─────▼──────┐  ┌─────▼──────┐
   │atn_io.py  │  │atn_model.py│  │config.txt  │
   │           │  │            │  │            │
   │Read & validate           Initialize model
   └────┬──────┘  │            │  │parameters  │
        │         │            │  └──────────┬─┘
        │         │            │             │
        └─────────┼────────────┼─────────────┘
                  │            │
            ┌─────▼────────────▼────┐
            │   ATN Model Instance   │
            │   (ATNModel class)     │
            │                        │
            │  Stores:               │
            │  - Species traits      │
            │  - Food-web adjacency  │
            │  - Herb/tree indices   │
            │  - Allometric params   │
            └─────────┬──────────────┘
                      │
        ┌─────────────┼─────────────┐
        │             │             │
    ┌───▼───┐     ┌───▼───┐     ┌──▼────┐
    │ Cell  │     │ Cell  │     │ Cell  │  ... (independent cells)
    │   0   │     │   1   │     │  N    │
    └───┬───┘     └───┬───┘     └──┬────┘
        │             │            │
        ├─────────────┼────────────┤  For each cell:
        │             │            │  1. Get temperature T_K and NPP
        │  ODE        │  ODE       │  2. Compute allometric rates
        │ SOLVER      │ SOLVER     │  3. Integrate dB/dt for each species
        │ (scipy.     │ (scipy.    │  4. Apply NPP-driven growth (basal)
        │ odeint)     │ odeint)    │  5. Apply Holling Type II (consumers)
        │             │            │
    ┌───▼───────────────────────────▼──┐
    │  Biomass Trajectory B(t, cell, sp) │
    │  Output shape: (time, cells, spp) │
    └───┬────────────────────────────┬──┘
        │                            │
    ┌───▼────────────┐  ┌───────────▼────┐
    │ Save output    │  │ Print summary   │
    │ biomass.txt +  │  │ (persistence,   │
    │ summary.txt    │  │  final biomass) │
    └────────────────┘  └─────────────────┘


EQUATIONS (per cell, g/m²/day):

Basal species i (herb or tree):
  dB_i/dt = G_i                                    [NPP-driven growth]
            - X_i(M_i,T) * B_i                     [metabolism]
            - Σ_j B_j F_ji(B,M,T)                  [herbivory]

  where G_i = NPP * ψ * (1 - f_struct_i) * C_i

  Herb:  C_i = α / (α + B_trees)
  Tree:  C_i = B_i / (α + B_trees)

Consumer species j:
  dB_j/dt = B_j * Σ_i e_i F_ij(B,M,T)             [feeding gain]
            - X_j(M_j,T) * B_j                     [metabolism]
            - Σ_k B_k F_jk(B,M,T)                  [predation]

Functional response (Holling Type II):
  F_ij = a_ij(M_i,M_j,T) * B_i^q
         ────────────────────────────────────────
         1 + c_j B_j + Σ_k h_kj(M_k,M_j,T) a_kj B_k^q

Allometric rates:
  X_i = X0 * M_i^(-0.25) * exp[-E(T₀-T)/(k_B T T₀)]
  a_ij = a0 * M_i^(-0.5) * M_j^(0.5) * exp[-E(T₀-T)/(k_B T T₀)]
  h_ij = h0 * M_i^(0.5) * M_j^(-0.5)


PARAMETER FLOW:

config.txt (read by atn_io.read_config())
   ↓
ATNModel.__init__() extracts rates:
   psi, f_struct, alpha_herbs, X0, b_X, a0, b_a_prey, b_a_pred, h0, etc.
   ↓
ATNModel.run_all_cells()
   ↓
For each cell: ATNModel.derivatives(B, t, cell_idx)
   ├─ Compute T_K and NPP from env_df[cell_idx]
   ├─ Compute a_ij(M, T) and h_ij(M, T) matrices
   ├─ Compute X(M, T) vector
   ├─ For each basal species i:
   │  └─ dB_i/dt = vegetation_growth(B, NPP, C_i) - loss
   ├─ For each consumer species j:
   │  └─ dB_j/dt = feeding_gain(Σ e_i F_ij) - loss
   └─ Return dydt vector
       ↓
ODE solver (scipy.odeint) integrates forward in time
   ↓
Biomass trajectory saved to atn_output/yyyymmddhhmmss/biomass.txt (long format)
```

## Function Call Graph

```
ATNModel.__init__(adj_mat, traits_df, env_df, config)
│  Reads adj_mat, traits_df, env_df; extracts allometric constants from config
│  Instantiates PlantVegetationModel (owns herb_idx, tree_idx, f_struct, alpha_herbs)
│
└── run_all_cells(B_initial, t_eval)
        │  iterates over spatial cells
        │
        └── run_cell(B_initial[cell], cell_idx, t_eval)
                │  wraps scipy.odeint for one cell
                │
                └── derivatives(y, t, cell_idx)     ← ODE RHS, called each solver step
                        │
                        ├── _allometric_rate('attack')    → a_ij  (n_spp × n_spp)
                        ├── _allometric_rate('handling')  → h_ij  (n_spp × n_spp)
                        ├── _metabolic_rate()             → X     (n_spp,)
                        ├── vegetation.growth(B, cell)   → G     (n_spp,)   [PlantVegetationModel]
                        │       reads NPP from env_df[cell_idx]
                        │       herb: C_i = α/(α + B_trees)
                        │       tree: C_i = B[i]/(α + B_trees)
                        │
                        └── _functional_response(B, j)   → F_ij  (n_spp,)  per consumer j
                                └─ uses a_ij, h_ij set above
```

## Requirements

```bash
pip install numpy scipy pandas
```

## Quick Start

```bash
python run_atn.py env_mat.txt adj_mat.txt traits.txt
```

### With custom simulation parameters:

```python
from run_atn import main

B_traj, t_eval, model = main(
    'env_mat.txt',
    'adj_mat.txt', 
    'traits.txt',
    t_max=500,           # simulate for 500 days
)
```

## Input Files

### 1. **env_mat.txt** (CSV)

Environmental matrix with one row per spatial cell.

**Required columns:**
- `pixel_id` (index): unique pixel identifier
- `x`: x coordinate of the pixel (non-negative integer, no duplicate `x, y` pairs)
- `y`: y coordinate of the pixel (non-negative integer)
- `temperature_K`: temperature in Kelvin (e.g., 293.15 = 20°C)
- `NPP`: net primary productivity in g C m⁻² day⁻¹ (must be positive)

**Example:**
```
pixel_id,x,y,temperature_K,NPP
0,0,0,293.15,3.5
1,0,1,293.15,3.2
2,0,2,288.15,2.8
```

### 2. **adj_mat.txt** (space or comma-separated)

Adjacency/food-web matrix: binary, resource×consumer orientation.

**Format:**
- Rows = resource species (rows eaten)
- Columns = consumer species (columns that eat)
- Binary: 0 (no link) or 1 (feeding link exists)
- Must be square: `n_species × n_species`
- No self-loops (diagonal must be 0)

**Example:**
```
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
...
1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
```

### 3. **traits.txt** (CSV)

Species trait table: one row per species.

**Required columns:**
- `species_id` (index): unique species identifier
- `body_mass_g`: body mass in grams (positive)
- `is_basal`: 1 if basal (plant), 0 if consumer (animal)
- `initial_biomass_g_per_m2`: starting biomass density in g/m²
- `vegetation_type`: `'herb'` or `'tree'` for all basal species (ignored for consumers)

**Optional columns:**

| Column | Applies to | Units | Default if absent |
|---|---|---|---|
| `f_struct` | basal species | dimensionless (0–1) | `f_struct_default` in `config.txt` |
| `metabolic_rate_base` | all species | day⁻¹ | `X0` in `config.txt` |
| `metabolic_rate_exponent` | all species | dimensionless | `b_X` in `config.txt` |
| `assimilation_plant` | consumer species | dimensionless (0–1) | `e_plant` in `config.txt` |
| `assimilation_animal` | consumer species | dimensionless (0–1) | `e_animal` in `config.txt` |

Any column may have `NaN` for individual rows; those rows fall back to the global config default.
Per-species metabolic parameters are most useful when the community mixes ectotherms and endotherms,
or when individual species have well-constrained empirical metabolic rates.

**Example:**
```
species_id,body_mass_g,is_basal,initial_biomass_g_per_m2,vegetation_type,metabolic_rate_base,metabolic_rate_exponent,assimilation_plant,assimilation_animal
0,1.7,1,10.0,herb,,,,
1,3.4,1,8.0,herb,,,,
2,6.0,1,5.0,tree,,,,
3,16.7,1,3.0,tree,,,,
4,30.4,0,2.0,,0.5,-0.25,0.45,0.85
5,32.8,0,2.0,,0.88,-0.25,0.20,0.72
...
39,617906,0,0.1,,,,
```

> **Note:** `vegetation_type` must be `'herb'` or `'tree'` for every row where `is_basal == 1`. The column must be present; an absent column raises a `ValidationError`.
>
> **Trees at zero initial biomass:** If all trees start at zero, `C_tree_i = 0` and tree species receive no NPP input. Trees require a positive `initial_biomass_g_per_m2` to grow.

## Output Files

Saved to `atn_output/yyyymmddhhmmss/` (one timestamped folder per run).

### simulation_summary.txt

Human-readable record of the run:
- Run timestamp
- Number of species, time steps, pixels, and grid dimensions
- Full species trait table (`species_id`, `body_mass_g`, `is_basal`, `initial_biomass_g_per_m2`)
- All model constants from `config.txt` with descriptions

### biomass.txt

Long-format table with one row per pixel × time step × species combination.

**Columns:** `pixel_id`, `x`, `y`, `time_step`, `species_id`, `biomass`

**Example rows:**
```
pixel_id x y time_step species_id biomass
0 0 0 0.0000 0 1.002345e+01
0 0 0 0.0000 1 8.012456e+00
...
2 0 2 100.0000 39 5.123400e-03
```

### Load and analyze results:

```python
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('atn_output/20260525143012/biomass.txt', sep=' ')

# Plot species 0 at pixel (0, 0) over time
s0 = df[(df['species_id'] == 0) & (df['x'] == 0) & (df['y'] == 0)]
plt.plot(s0['time_step'], s0['biomass'])
plt.xlabel('Time (days)')
plt.ylabel('Biomass (g/m²)')
plt.show()

# Final mean biomass per species across all pixels
final = df[df['time_step'] == df['time_step'].max()]
print(final.groupby('species_id')['biomass'].mean())
```

## Configuration

Edit `config.txt` to modify parameters. Format: `key = value  # comment`. Lines starting with `#` are ignored. Booleans are `True` or `False`; scientific notation (`1e-6`) is supported.

```
# ===== ALLOMETRIC RATE CONSTANTS =====
X0             = 0.5       # metabolic rate normalization (day^-1)
b_X            = -0.25     # metabolic exponent
a0             = 0.001     # attack rate normalization (day^-1)
b_a_prey       = -0.5      # attack rate prey exponent
b_a_pred       = 0.5       # attack rate predator exponent
h0             = 0.01      # handling time normalization (days)
b_h_prey       = 0.5       # handling time prey exponent
b_h_pred       = -0.5      # handling time predator exponent

# ===== FUNCTIONAL RESPONSE =====
q_hill         = 2.0       # Hill exponent (Type II ~= 2)
interference   = 0.0       # consumer interference coefficient

# ===== ASSIMILATION EFFICIENCIES =====
e_plant        = 0.45      # plant assimilation (overridden by assimilation_plant in traits.txt)
e_animal       = 0.85      # animal assimilation (overridden by assimilation_animal in traits.txt)

# ===== VEGETATION GROWTH =====
psi                 = 9.813  # C-to-wet-matter conversion (g wet / g C)
f_struct_default    = 0.3    # default fraction of NPP to structural tissue
alpha_herbs_default = 1.0    # half-saturation constant for herb/tree competition (g/m2)

# ===== TEMPERATURE DEPENDENCE =====
use_temperature = True     # apply temperature scaling to rates (True or False)
T0_K            = 293.15   # reference temperature (K, 20 C)
k_B             = 8.617e-5 # Boltzmann constant (eV/K)
E_a             = 0.65     # activation energy (eV)

# ===== EXTINCTION =====
ext_threshold        = 1e-6  # biomass below this is treated as extinct (g/m2)
extinction_timescale = 0.1   # decay timescale for extinct species (days)
```

To use a non-default config path:

```bash
python run_atn.py env_mat.txt adj_mat.txt traits.txt my_config.txt
```

```python
B_traj, t_eval, model = main('env_mat.txt', 'adj_mat.txt', 'traits.txt', config_file='my_config.txt')
```

## Model Details

### Equations

**Basal species (plants):**

$$\frac{dB_i}{dt} = G_i - X_i B_i - \sum_j B_j F_{ij}$$

where $G_i$ is the NPP-driven leaf biomass growth rate (Harfoot et al. 2014 Text S1):

$$G_i = NPP \cdot \psi \cdot (1 - f_{\text{struct},i}) \cdot C_i$$

with the Michaelis-Menten competitive partition between vegetation types:

$$C_{\text{herb}} = \frac{\alpha}{\alpha + B_{\text{trees}}} \qquad C_{\text{tree}} = \frac{B_i}{\alpha + B_{\text{trees}}}$$

**Consumers (animals):**
$$\frac{dB_j}{dt} = B_j \sum_i e_i F_{ij} - X_j B_j - \sum_k B_k F_{jk}$$

where $e_i$ = assimilation efficiency on resource $i$.

### Allometric scaling

Metabolic and attack rates scale with body mass $M$ and temperature $T$:

$$p(M, T) = p_0 M^b \exp\left(\frac{-E_a (T_0 - T)}{k_B T T_0}\right)$$

### Functional response (Type II, Holling):

$$F_{ij} = \frac{a_{ij} B_i^q}{1 + c_j B_j + \sum_k h_{kj} a_{kj} B_k^q}$$

where $a_{ij}$ = attack rate, $q$ = Hill exponent, $h_{ij}$ = handling time, $c_j$ = consumer interference.

## Performance Tips

- **Shorter simulation:** Lower `t_max` for faster runs
- **Temperature off:** Set `use_temperature: False` to skip temperature calculations
- **Parallel cells:** Future: run cells in parallel with `concurrent.futures`

## Validation Features

The system includes comprehensive checks (all in `atn_io.py`):

✓ File existence & CSV parsing  
✓ Required columns in each file (`NPP` in env_mat, `vegetation_type` in traits)  
✓ `vegetation_type` values are `'herb'` or `'tree'` for all basal species  
✓ NPP values are positive (range printed so unit errors are caught early)  
✓ Temperature range printed (catches Celsius-vs-Kelvin errors before they corrupt allometric rates)  
✓ Species count consistency across files  
✓ No self-loops in adjacency matrix  
✓ Basal species don't consume  
✓ Body mass, biomass, and efficiency ranges  
✓ Parameter completeness and realism  

Run with bad inputs to see detailed error messages.

## References

- Binzer, A., Guill, C., Rall, B. C., & Brose, U. (2016). Interactive effects of warming, eutrophication and size structure: impacts on biodiversity and food-web structure. *Global Change Biology*, 22, 220–227.
- Brose, U., Williams, R. J., & Martinez, N. D. (2006). Allometric scaling enhances stability in complex food webs. *Ecology Letters*, 9, 1228–1236.
- Brown, J. H., Gillooly, J. F., Allen, A. P., Savage, V. M., & West, G. B. (2004). Toward a metabolic theory of ecology. *Ecology*, 85, 1771–1789.
- Delmas, E., Brose, U., Gravel, D., Stouffer, D. B., & Poisot, T. (2017). Simulations of biomass dynamics in community food webs. *Methods in Ecology and Evolution*, 8, 881–886.
- Gillooly, J. F., Brown, J. H., West, G. B., Savage, V. M., & Charnov, E. L. (2001). Effects of size and temperature on metabolic rate. *Science*, 293, 2248–2251.
- Harfoot, M. B. J., et al. (2014). Emergent global patterns of ecosystem structure and function from a mechanistic general ecosystem model. *PLOS Biology*, 12, e1001841. *(Text S1: terrestrial plant model equations)*
- Kattge, J., et al. (2011). TRY – a global database of plant traits. *Global Change Biology*, 17, 2905–2935. *(source of ψ = 9.813 g wet / g C)*
- Kraus, D., et al. (2022). Coupling the Madingley general ecosystem model to LPJ-GUESS vegetation model for more realistic vegetation dynamics. *Ecological Modelling*.
- Rall, B. C., et al. (2012). Universal temperature and body-mass scaling of feeding rates. *Philosophical Transactions of the Royal Society B*, 367, 2923–2934.
- Williams, R. J., & Martinez, N. D. (2004). Stabilization of chaotic and non-permanent food-web dynamics. *European Physical Journal B*, 38, 297–303.
- Yodzis, P., & Innes, S. (1992). Body size and consumer-resource dynamics. *The American Naturalist*, 139, 1151–1175.

### Parameter sources

| Parameter(s) | Value(s) | Source |
|---|---|---|
| `psi` | 9.813 g wet / g C | Kattge et al. (2011) via Harfoot et al. (2014) Text S1 |
| `f_struct_default` | 0.3 | De Kauwe et al. (2014); Harfoot et al. (2014) Text S1 |
| `alpha_herbs_default` | 1.0 g/m² | Michaelis-Menten partition (see `context/vegetation.md`) |
| `X0`, `b_X` | 0.5, −0.25 | Yodzis & Innes (1992); Brown et al. (2004) — per-species override via `metabolic_rate_base`, `metabolic_rate_exponent` |
| `a0`, `b_a_prey`, `b_a_pred` | 0.001, −0.5, 0.5 | Rall et al. (2012); Brose et al. (2006) |
| `h0`, `b_h_prey`, `b_h_pred` | 0.01, 0.5, −0.5 | Rall et al. (2012) |
| `e_plant`, `e_animal` | 0.45, 0.85 | Yodzis & Innes (1992) — per-species override via `assimilation_plant`, `assimilation_animal` |
| `q_hill` | 2.0 | Williams & Martinez (2004) |
| `E_a` | 0.65 eV | Brown et al. (2004); Gillooly et al. (2001) |
| `k_B` | 8.617 × 10⁻⁵ eV/K | Fundamental constant (NIST CODATA) |

## License

MIT
