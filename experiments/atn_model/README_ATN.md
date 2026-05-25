# Spatially Explicit ATN Model

Efficient Python implementation of the Allometric Trophic Network (ATN) model 
with spatial heterogeneity (Section 8, unscaled model from `ATN_model_spatiotemporal_formulas_parameters.Rmd`).

## Scripts

| Script | Role | Run directly |
|---|---|---|
| `run_atn.py` | Main entry point: reads inputs, validates, integrates ODEs, saves results | Yes |
| `test_atn.py` | Quick sanity check: runs a 10-day simulation on example files to verify setup | Yes |
| `atn_model.py` | `ATNModel` class: all ODE dynamics | No (imported by other scripts) |
| `atn_io.py` | File reading and 20+ validation checks | No (imported by other scripts) |
| `config.py` | Default allometric and numerical parameters | No (imported by other scripts) |

## Running with your own data

### Step 1 — Prepare three input files

| File | Format | What it contains |
|---|---|---|
| `env_mat.txt` | CSV | One row per spatial cell: `cell_id`, `temperature_K`, optional `K_plant_i` columns |
| `adj_mat.txt` | Space- or comma-separated | Square binary matrix (rows = resources, columns = consumers) |
| `traits.txt` | CSV | One row per species: `species_id`, `body_mass_g`, `is_basal`, `initial_biomass_g_per_m2` |

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
    n_timesteps=500,   # number of saved time points (default: 1000)
    output_dir='./results'
)
```

### Step 3 — Check outputs

Results are saved to `./atn_output/` (or the directory you specified):

- `biomass_trajectory.npy` — array of shape `(n_timesteps, n_cells, n_species)`
- `time_points.npy` — time vector in days

The console prints per-species final biomass and the fraction of cells where each species persists.

### Validate your setup first

Before a long run, use `test_atn.py` to confirm your files are readable and the model initialises correctly:

```bash
python test_atn.py
```

This runs a 10-day simulation and checks for NaN values and shape consistency.

## Model Overview

```
┌─────────────────────────────────────────────────────────────────┐
│              SPATIALLY EXPLICIT ATN MODEL ARCHITECTURE          │
└─────────────────────────────────────────────────────────────────┘

INPUT FILES:
┌──────────────┐  ┌──────────────┐  ┌──────────────┐
│ env_mat.txt  │  │ adj_mat.txt  │  │ traits.txt   │
│              │  │              │  │              │
│ Temp, K per  │  │ Food-web     │  │ Body mass,   │
│ cell         │  │ links        │  │ metabolism   │
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
   │atn_io.py  │  │atn_model.py│  │config.py   │
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
        │             │            │  1. Get temperature T_K
        │  ODE        │  ODE       │  2. Compute allometric rates
        │ SOLVER      │ SOLVER     │  3. Integrate dB/dt for each species
        │ (scipy.     │ (scipy.    │  4. Apply logistic growth (basal)
        │ odeint)     │ odeint)    │  5. Apply Holling Type II (consumers)
        │             │            │
    ┌───▼───────────────────────────▼──┐
    │  Biomass Trajectory B(t, cell, sp) │
    │  Output shape: (time, cells, spp) │
    └───┬────────────────────────────┬──┘
        │                            │
    ┌───▼────────────┐  ┌───────────▼────┐
    │ Save output    │  │ Print summary   │
    │ .npy files     │  │ (persistence,   │
    │                │  │  final biomass) │
    └────────────────┘  └─────────────────┘


EQUATIONS (per cell g):

Basal species i:
  dB_i/dt = r_i(M_i,T) * B_i * (1 - B_i/K_i)      [logistic growth]
            - X_i(M_i,T) * B_i                     [metabolism]
            - Σ_j B_j F_ji(B,M,T)                  [herbivory]

Consumer species j:
  dB_j/dt = B_j * Σ_i e_i F_ij(B,M,T)             [feeding gain]
            - X_j(M_j,T) * B_j                     [metabolism]
            - Σ_k B_k F_jk(B,M,T)                  [predation]

Functional response (Holling Type II):
  F_ij = a_ij(M_i,M_j,T) * B_i^q
         ────────────────────────────────────────
         1 + c_j B_j + Σ_k h_kj(M_k,M_j,T) a_kj B_k^q

Allometric rates:
  r_i = r0 * M_i^(-0.25) * exp[-E(T₀-T)/(k_B T T₀)]
  X_i = X0 * M_i^(-0.25) * exp[-E(T₀-T)/(k_B T T₀)]
  a_ij = a0 * M_i^(-0.5) * M_j^(0.5) * exp[-E(T₀-T)/(k_B T T₀)]
  h_ij = h0 * M_i^(0.5) * M_j^(-0.5)


PARAMETER FLOW:

config.py (default parameters)
   ↓
ATNModel.__init__() extracts rates:
   r0, b_r, X0, b_X, a0, b_a_prey, b_a_pred, h0, etc.
   ↓
ATNModel.run_all_cells()
   ↓
For each cell: ATNModel.derivatives(B, t, cell_idx)
   ├─ Compute T_K from env_df[cell_idx]
   ├─ Compute a_ij(M, T) and h_ij(M, T) matrices
   ├─ Compute X(M, T) and r(M, T) vectors
   ├─ For each basal species i:
   │  └─ dB_i/dt = logistic_growth(B_i, r_i, K_i) - loss
   ├─ For each consumer species j:
   │  └─ dB_j/dt = feeding_gain(Σ e_i F_ij) - loss
   └─ Return dydt vector
       ↓
ODE solver (scipy.odeint) integrates forward in time
   ↓
Biomass trajectory saved to .npy files
```

## Function Call Graph

```
ATNModel.__init__(adj_mat, traits_df, env_df, config)
│  Reads adj_mat, traits_df, env_df; extracts allometric constants from config
│
├── _L_matrix()                    body-size feeding kernel L_ij
│       └─ returns L * adj_mat     (food-web topology with size-matching)
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
                        ├── _basal_growth_rate()          → r     (n_spp,)
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
    n_timesteps=500,     # save 500 output points
    output_dir='./results'
)
```

## Input Files

### 1. **env_mat.txt** (CSV)

Environmental matrix with one row per spatial cell.

**Columns:**
- `cell_id` (index): unique cell identifier
- `temperature_K`: temperature in Kelvin (e.g., 293.15 = 20°C)
- `K_plant_0`, `K_plant_1`, ...: carrying capacity for each basal species (optional)

**Example:**
```
cell_id,temperature_K,K_plant_0,K_plant_1,K_plant_2,K_plant_3,K_plant_4,K_plant_5,K_plant_6,K_plant_7,K_plant_8,K_plant_9
0,293.15,100,100,100,100,100,100,100,100,100,100
1,293.15,100,100,100,100,100,100,100,100,100,100
2,288.15,90,90,90,90,90,90,90,90,90,90
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

**Optional but recommended columns:**
- `metabolic_rate_base`: if provided, overrides model default
- `metabolic_rate_exponent`: if provided, overrides model default
- `assimilation_plant`: efficiency when eating plants (0–1)
- `assimilation_animal`: efficiency when eating animals (0–1)
- `hill_exponent`: Hill coefficient for functional response (>0)

**Example:**
```
species_id,body_mass_g,is_basal,initial_biomass_g_per_m2
0,1.7,1,10.0
1,3.4,1,8.0
2,6.0,1,5.0
3,16.7,1,3.0
4,30.4,1,2.0
5,32.8,1,2.0
6,43.5,1,1.5
7,346.5,1,0.5
8,364.7,1,0.5
9,810.8,1,0.3
10,126.97,0,1.0
11,301.86,0,0.8
12,345.40,0,0.6
...
39,617906,0,0.1
```

## Output Files

Saved to `./atn_output/` by default:

- **biomass_trajectory.npy**: NumPy array of shape `(n_timesteps, n_cells, n_species)`
  - Axis 0: time points
  - Axis 1: spatial cells
  - Axis 2: species indices

- **time_points.npy**: NumPy array of time values (days)

### Load and analyze results:

```python
import numpy as np
import matplotlib.pyplot as plt

# Load data
B_traj = np.load('atn_output/biomass_trajectory.npy')
t = np.load('atn_output/time_points.npy')

print(f"Shape: {B_traj.shape}")  # (n_timesteps, n_cells, n_species)

# Plot species 0 in cell 0 over time
plt.plot(t, B_traj[:, 0, 0])
plt.xlabel('Time (days)')
plt.ylabel('Biomass (g/m²)')
plt.title('Species 0, Cell 0')
plt.show()

# Final average biomass per species
final_avg = B_traj[-1, :, :].mean(axis=0)
print(f"Final biomass per species: {final_avg}")
```

## Configuration

Edit `config.py` to modify parameters:

```python
CONFIG = {
    # Allometric parameters
    'r0': 0.5,            # basal growth normalization
    'b_r': -0.25,         # basal growth exponent
    'X0': 0.5,            # metabolic rate normalization
    'b_X': -0.25,         # metabolic exponent
    'a0': 0.001,          # attack rate normalization
    'b_a_prey': -0.5,     # attack rate prey exponent
    'b_a_pred': 0.5,      # attack rate predator exponent
    'h0': 0.01,           # handling time normalization
    
    # Functional response
    'q_hill': 2.0,        # Hill exponent (Type II ≈ 2)
    'R_opt': 100.0,       # optimal predator/prey mass ratio
    'gamma': 2.0,         # L-matrix sharpness
    
    # Efficiency
    'e_plant': 0.45,      # plant assimilation
    'e_animal': 0.85,     # animal assimilation
    
    # Temperature dependence
    'use_temperature': True,
    'T0_K': 293.15,       # reference temp (20°C)
    'E_a': 0.65,          # activation energy (eV)
    
    # Extinction
    'ext_threshold': 1e-6,  # biomass below this → extinct
}
```

## Model Details

### Equations

**Basal species (plants):**
$$\frac{dB_i}{dt} = B_i r_i G_i - X_i B_i - \sum_j B_j F_{ij}$$

where:
- $G_i = 1 - B_i/K_i$ (logistic growth regulation)
- $r_i$ = basal growth rate
- $X_i$ = metabolic loss rate
- $F_{ij}$ = feeding rate of consumer $j$ on resource $i$

**Consumers (animals):**
$$\frac{dB_j}{dt} = B_j \sum_i e_i F_{ij} - X_j B_j - \sum_k B_k F_{jk}$$

where:
- $e_i$ = assimilation efficiency on resource $i$

### Allometric scaling

All biological rates scale with body mass $M$ and temperature $T$:

$$p(M, T) = p_0 M^b \exp\left(\frac{-E(T_0 - T)}{k_B T T_0}\right)$$

### Functional response (Type II, Holling):

$$F_{ij} = \frac{a_{ij} B_i^q}{1 + c_j B_j + \sum_k h_{kj} a_{kj} B_k^q}$$

where:
- $a_{ij}$ = attack rate
- $q$ = Hill exponent
- $h_{ij}$ = handling time
- $c_j$ = consumer interference

## Performance Tips

- **Coarser output:** Lower `n_timesteps` for faster runs
- **Temperature off:** Set `use_temperature: False` to skip temperature calculations
- **High connectivity:** Increase `link_threshold` to skip weak interactions
- **Parallel cells:** Future: run cells in parallel with `concurrent.futures`

## Validation Features

The system includes comprehensive checks:

✓ File existence & CSV parsing  
✓ Required columns in each file  
✓ Species count consistency across files  
✓ No self-loops in adjacency matrix  
✓ Basal species don't consume  
✓ Body mass, biomass, and efficiency ranges  
✓ Parameter completeness and realism  
✓ Temperature range checking  

Run with bad inputs to see detailed error messages.

## References

- Binzer, A., Guill, C., Rall, B. C., & Brose, U. (2016). Interactive effects of warming, eutrophication and size structure: impacts on biodiversity and food-web structure. *Global Change Biology*, 22, 220–227.
- Brose, U., Williams, R. J., & Martinez, N. D. (2006). Allometric scaling enhances stability in complex food webs. *Ecology Letters*, 9, 1228–1236.
- Brown, J. H., Gillooly, J. F., Allen, A. P., Savage, V. M., & West, G. B. (2004). Toward a metabolic theory of ecology. *Ecology*, 85, 1771–1789.
- Delmas, E., Brose, U., Gravel, D., Stouffer, D. B., & Poisot, T. (2017). Simulations of biomass dynamics in community food webs. *Methods in Ecology and Evolution*, 8, 881–886.
- Gillooly, J. F., Brown, J. H., West, G. B., Savage, V. M., & Charnov, E. L. (2001). Effects of size and temperature on metabolic rate. *Science*, 293, 2248–2251.
- Rall, B. C., Brose, U., Hartvig, M., Kalinkat, G., Schwarzmüller, F., Vucic-Pestic, O., & Petchey, O. L. (2012). Universal temperature and body-mass scaling of feeding rates. *Philosophical Transactions of the Royal Society B*, 367, 2923–2934.
- Williams, R. J., & Martinez, N. D. (2004). Stabilization of chaotic and non-permanent food-web dynamics. *European Physical Journal B*, 38, 297–303.
- Yodzis, P., & Innes, S. (1992). Body size and consumer-resource dynamics. *The American Naturalist*, 139, 1151–1175.

### Parameter sources

| Parameter(s) | Value(s) | Source |
|---|---|---|
| `r0`, `X0` | 0.5 | Yodzis & Innes (1992) |
| `b_r`, `b_X` | −0.25 | Brown et al. (2004); Yodzis & Innes (1992) |
| `a0`, `b_a_prey`, `b_a_pred` | 0.001, −0.5, 0.5 | Rall et al. (2012); Brose et al. (2006) |
| `h0`, `b_h_prey`, `b_h_pred` | 0.01, 0.5, −0.5 | Rall et al. (2012) |
| `e_plant`, `e_animal` | 0.45, 0.85 | Yodzis & Innes (1992) |
| `R_opt` | 100 | Brose et al. (2006) |
| `q_hill` | 2.0 | Williams & Martinez (2004) |
| `E_a` | 0.65 eV | Brown et al. (2004); Gillooly et al. (2001) |
| `k_B` | 8.617 × 10⁻⁵ eV/K | Fundamental constant (NIST CODATA) |

## License

MIT
