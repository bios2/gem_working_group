"""
Read and validate input files for the spatially explicit ATN model.
Includes detailed sanity checks and diagnostic messages.
"""
# Import required libraries
import numpy as np  # numerical arrays and operations
import pandas as pd  # tabular data handling (CSV reading)
from typing import Dict, Tuple, Optional  # type hints
import sys  # system utilities for exit codes

# Custom exception class for validation errors
class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass

class ValidationWarning:
    """Warning message (non-fatal)."""
    def __init__(self, msg: str):
        # Store the warning message
        self.msg = msg
    def __str__(self):
        # Return formatted warning string with symbol
        return f"⚠ WARNING: {self.msg}"

def read_config(filepath: str) -> Dict:
    """
    Read model configuration from a plain-text file.

    Format: one 'key = value' pair per non-blank line; lines starting with '#'
    and inline '#' comments are ignored.  Values are parsed as bool (True/False)
    or float (handles scientific notation such as 1e-6 and 8.617e-5).
    """
    config: Dict = {}
    try:
        with open(filepath) as f:
            for lineno, raw in enumerate(f, 1):
                line = raw.strip()
                if not line or line.startswith('#'):
                    continue
                if '#' in line:
                    line = line[:line.index('#')].strip()
                if '=' not in line:
                    continue
                key, _, val_str = line.partition('=')
                key = key.strip()
                val_str = val_str.strip()
                if val_str in ('True', 'true'):
                    config[key] = True
                elif val_str in ('False', 'false'):
                    config[key] = False
                else:
                    try:
                        config[key] = float(val_str)
                    except ValueError:
                        raise ValidationError(
                            f"config {filepath} line {lineno}: "
                            f"cannot parse value for '{key}': '{val_str}'"
                        )
    except FileNotFoundError:
        raise ValidationError(f"Config file not found: {filepath}")
    return config


def read_env_matrix(filepath: str) -> Tuple[pd.DataFrame, np.ndarray]:
    """
    Read environmental matrix with validation.

    Expected format (CSV):
      pixel_id, x, y, temperature_K, K_plant_0, K_plant_1, ...

    Returns:
        env_df: DataFrame with pixel_id index and columns x, y, temperature_K, K_*
        temp_array: array of temperatures (N_cells,)
    """
    # Print diagnostic header
    print("\n[ENV_MAT] Reading environmental matrix...")

    # Try to read CSV file with first column as index
    try:
        env_df = pd.read_csv(filepath, index_col=0)
    except FileNotFoundError:
        # Raise if file doesn't exist
        raise ValidationError(f"Environment file not found: {filepath}")
    except Exception as e:
        # Raise if CSV parsing fails
        raise ValidationError(f"Failed to parse {filepath}: {e}")

    # Check for required spatial coordinate columns
    for coord in ('x', 'y'):
        if coord not in env_df.columns:
            raise ValidationError(
                f"Missing required column '{coord}' in {filepath}.\n"
                f"  Found columns: {list(env_df.columns)}\n"
                f"  env_mat.txt must have columns: pixel_id, x, y, temperature_K, ..."
            )

    # Validate x/y: must be non-negative integers with no duplicate (x, y) pairs
    if (env_df['x'] < 0).any() or (env_df['y'] < 0).any():
        raise ValidationError("'x' and 'y' values must be non-negative integers.")
    dupes = env_df.duplicated(subset=['x', 'y'])
    if dupes.any():
        raise ValidationError(
            f"Duplicate (x, y) pairs found in {filepath}:\n"
            f"  {env_df[dupes][['x', 'y']].values.tolist()}"
        )

    n_x = int(env_df['x'].max()) + 1
    n_y = int(env_df['y'].max()) + 1
    # Print grid dimensions so the user can verify the spatial domain loaded as expected
    # before committing to a potentially long simulation run
    print(f"  ✓ Loaded {len(env_df)} cells on a {n_x}×{n_y} grid")

    # Check for required column: temperature (must exist)
    if 'temperature_K' not in env_df.columns:
        raise ValidationError(
            f"Missing required column 'temperature_K' in {filepath}.\n"
            f"  Found columns: {list(env_df.columns)}"
        )

    # Extract temperature range and check for realistic values
    temp_min = env_df['temperature_K'].min()
    temp_max = env_df['temperature_K'].max()
    # Reject if outside physiological range (260-330 K = -13 to +57°C)
    if temp_min < 260 or temp_max > 330:
        raise ValidationError(
            f"Temperature out of realistic range: {temp_min:.2f}–{temp_max:.2f} K.\n"
            f"  Expected ~260–330 K (−13°C to 57°C)."
        )

    # Print temperature range in both Kelvin and Celsius so the user can immediately spot
    # unit errors (e.g. Celsius supplied instead of Kelvin) or unrealistic outlier cells
    # before they silently distort all allometric and metabolic rates in the simulation
    print(f"  ✓ Temperature range: {temp_min:.2f}–{temp_max:.2f} K ({temp_min-273.15:.1f}–{temp_max-273.15:.1f}°C)")

    # Check for required NPP column (single value per cell, no monthly breakdown)
    if 'NPP' not in env_df.columns:
        raise ValidationError(
            f"Missing required column 'NPP' in {filepath}.\n"
            f"  env_mat.txt must contain an 'NPP' column with one value per cell.\n"
            f"  Found columns: {list(env_df.columns)}"
        )
    # Validate that all NPP values are strictly positive
    npp_vals = env_df['NPP'].values
    if np.any(npp_vals <= 0):
        raise ValidationError(
            f"All NPP values must be positive. "
            f"Found {np.sum(npp_vals <= 0)} non-positive values: "
            f"{npp_vals[npp_vals <= 0].tolist()}"
        )
    # Print NPP range so the user can catch order-of-magnitude unit errors
    # (e.g. kg C vs g C) before they propagate into nonsensical basal biomass dynamics
    print(f"  ✓ NPP range: {npp_vals.min():.2e}–{npp_vals.max():.2e}")

    # Check for missing values (NaN) in any column
    if env_df.isnull().any().any():
        n_nan = env_df.isnull().sum().sum()  # count total NaN values
        raise ValidationError(
            f"Found {n_nan} NaN values in environment matrix.\n"
            f"  NaN columns: {list(env_df.columns[env_df.isnull().any()])}"
        )

    # Return DataFrame and extracted temperature array
    return env_df, env_df['temperature_K'].values


def read_adjacency_matrix(filepath: str, expected_n_species: Optional[int] = None) -> np.ndarray:
    """
    Read adjacency/food-web matrix with validation.
    
    Format: space or comma-separated binary matrix (rows=resources, cols=consumers).
    
    Parameters:
        filepath: path to adjacency matrix
        expected_n_species: if provided, check that matrix size matches
    
    Returns:
        adj_mat: (N_species, N_species) binary matrix
    """
    print("\n[ADJ_MAT] Reading adjacency matrix...")
    
    try:
        adj_mat = np.loadtxt(filepath, dtype=int)
    except FileNotFoundError:
        raise ValidationError(f"Adjacency matrix file not found: {filepath}")
    except Exception as e:
        raise ValidationError(f"Failed to parse adjacency matrix: {e}")
    
    # Check shape
    if adj_mat.ndim != 2:
        raise ValidationError(
            f"Adjacency matrix must be 2D, got {adj_mat.ndim}D with shape {adj_mat.shape}"
        )
    
    if adj_mat.shape[0] != adj_mat.shape[1]:
        raise ValidationError(
            f"Adjacency matrix must be square, got {adj_mat.shape[0]}×{adj_mat.shape[1]}"
        )
    
    n_spp = adj_mat.shape[0]
    print(f"  ✓ Loaded {n_spp}×{n_spp} adjacency matrix")
    
    # Check binary
    unique_vals = np.unique(adj_mat)
    if not np.all(np.isin(unique_vals, [0, 1])):
        raise ValidationError(
            f"Adjacency matrix must be binary (0/1), found values: {unique_vals}"
        )
    
    # Check no self-loops (species cannot eat themselves)
    n_selfloops = np.sum(np.diag(adj_mat))
    if n_selfloops > 0:
        raise ValidationError(
            f"Found {n_selfloops} self-loops (diagonal entries). Species cannot eat themselves.\n"
            f"  Diagonal entries: {np.diag(adj_mat)}"
        )
    
    # Check connectivity
    n_links = np.sum(adj_mat)
    print(f"  ✓ {n_links} feeding links (density: {n_links / (n_spp**2):.3f})")
    
    # Check for disconnected consumers (carnivores with no prey)
    consumers_with_prey = np.sum(adj_mat, axis=0)  # sum over resources
    n_orphan_consumers = np.sum(consumers_with_prey == 0)
    if n_orphan_consumers > 0:
        orphans = np.where(consumers_with_prey == 0)[0]
        print(f"  ⚠ {n_orphan_consumers} species have no resources:")
        for i in orphans[:5]:
            print(f"    Species {i}")
        if len(orphans) > 5:
            print(f"    ... and {len(orphans)-5} more")
    
    # Check for species with no consumers
    resources_with_consumers = np.sum(adj_mat, axis=1)  # sum over consumers
    n_orphan_resources = np.sum(resources_with_consumers == 0)
    if n_orphan_resources > 0:
        orphans = np.where(resources_with_consumers == 0)[0]
        print(f"  ℹ {n_orphan_resources} species have no consumers (may be expected):")
        for i in orphans[:5]:
            print(f"    Species {i}")
        if len(orphans) > 5:
            print(f"    ... and {len(orphans)-5} more")
    
    if expected_n_species is not None and n_spp != expected_n_species:
        raise ValidationError(
            f"Species count mismatch: adjacency matrix is {n_spp}×{n_spp} "
            f"but expected {expected_n_species} species"
        )
    
    return adj_mat


def read_traits(filepath: str, expected_n_species: Optional[int] = None) -> pd.DataFrame:
    """
    Read species traits with comprehensive validation.
    
    Expected format (CSV):
      species_id, body_mass_g, is_basal, initial_biomass_g_per_m2, 
      metabolic_rate_base, metabolic_rate_exponent, assimilation_plant, 
      assimilation_animal, hill_exponent, ...
    
    Parameters:
        filepath: path to traits file
        expected_n_species: if provided, check that count matches
    
    Returns:
        traits_df: DataFrame indexed by species_id
    """
    print("\n[TRAITS] Reading species traits...")
    
    try:
        traits_df = pd.read_csv(filepath, index_col=0)
    except FileNotFoundError:
        raise ValidationError(f"Traits file not found: {filepath}")
    except Exception as e:
        raise ValidationError(f"Failed to parse traits file: {e}")
    
    n_spp = len(traits_df)
    print(f"  ✓ Loaded {n_spp} species")
    
    # Check for required columns
    required = ['body_mass_g', 'is_basal', 'initial_biomass_g_per_m2']
    missing = [col for col in required if col not in traits_df.columns]
    if missing:
        raise ValidationError(
            f"Missing required columns in {filepath}:\n"
            f"  Required: {required}\n"
            f"  Missing: {missing}\n"
            f"  Found: {list(traits_df.columns)}"
        )
    
    # Check for required columns: vegetation_type must be present
    if 'vegetation_type' not in traits_df.columns:
        raise ValidationError(
            f"Missing required column 'vegetation_type' in {filepath}.\n"
            f"  Each basal species must have vegetation_type set to 'herb' or 'tree'.\n"
            f"  Found columns: {list(traits_df.columns)}"
        )
    # Validate that vegetation_type is only set on basal species and has valid values
    basal_mask = traits_df['is_basal'] == 1
    veg_types = traits_df.loc[basal_mask, 'vegetation_type']
    invalid_types = veg_types[~veg_types.isin(['herb', 'tree'])]
    if len(invalid_types) > 0:
        raise ValidationError(
            f"vegetation_type must be 'herb' or 'tree' for basal species. "
            f"Found invalid values: {invalid_types.tolist()}"
        )
    # Non-basal species should not have vegetation_type set
    consumer_veg = traits_df.loc[~basal_mask, 'vegetation_type'] if 'vegetation_type' in traits_df.columns else pd.Series(dtype=str)
    non_null_consumers = consumer_veg.dropna()
    if len(non_null_consumers) > 0:
        print(f"  ⚠ vegetation_type set on {len(non_null_consumers)} consumer species (will be ignored)")
    n_herbs = (veg_types == 'herb').sum()
    n_trees = (veg_types == 'tree').sum()
    print(f"  ✓ vegetation_type: {n_herbs} herb species, {n_trees} tree species")

    # Check NaN only in required numeric columns (optional columns like f_struct may have NaN)
    required_cols_for_nan = ['body_mass_g', 'is_basal', 'initial_biomass_g_per_m2']
    for col in required_cols_for_nan:
        if col in traits_df.columns and traits_df[col].isnull().any():
            raise ValidationError(
                f"Found NaN values in required column '{col}' in {filepath}."
            )

    # Report optional vegetation trait columns if present
    if 'f_struct' in traits_df.columns:
        # f_struct may have NaN (will use config default for those species)
        n_present = traits_df['f_struct'].notna().sum()
        print(f"  ✓ f_struct: {n_present}/{n_spp} species have per-species values (rest use config default)")
    
    # Check body mass
    M = traits_df['body_mass_g'].values
    if np.any(M <= 0):
        raise ValidationError(
            f"Body mass must be positive. Found {np.sum(M <= 0)} non-positive values."
        )
    print(f"  ✓ Body mass range: {M.min():.2e}–{M.max():.2e} g (ratio: {M.max()/M.min():.2e})")
    
    # Check is_basal
    if 'is_basal' in traits_df.columns:
        n_basal = (traits_df['is_basal'] == 1).sum()
        n_consumer = (traits_df['is_basal'] == 0).sum()
        if n_basal == 0:
            raise ValidationError("No basal species (is_basal=1) found. At least 1 required.")
        print(f"  ✓ {n_basal} basal species, {n_consumer} consumer species")
        
        # Check for invalid is_basal values
        invalid_basal = np.sum(~traits_df['is_basal'].isin([0, 1]))
        if invalid_basal > 0:
            raise ValidationError(
                f"is_basal must be 0 or 1. Found {invalid_basal} other values: "
                f"{np.unique(traits_df['is_basal'])}"
            )
    
    # Check initial biomass
    B0 = traits_df['initial_biomass_g_per_m2'].values
    if np.any(B0 < 0):
        raise ValidationError(
            f"Initial biomass must be non-negative. Found {np.sum(B0 < 0)} negative values."
        )
    if np.all(B0 == 0):
        raise ValidationError("All initial biomasses are zero; simulation will be trivial.")
    print(f"  ✓ Initial biomass range: {B0[B0>0].min():.2e}–{B0.max():.2e} g/m²")
    
    # Report optional per-species metabolic columns if present
    if 'metabolic_rate_base' in traits_df.columns:
        # Count non-NaN entries so the user knows how many species use per-species values
        n_present = traits_df['metabolic_rate_base'].notna().sum()
        X0_vals = traits_df['metabolic_rate_base'].dropna().values
        if np.any(X0_vals < 0):
            raise ValidationError(
                f"metabolic_rate_base must be non-negative. Found negative values."
            )
        print(f"  ✓ metabolic_rate_base: {n_present}/{n_spp} species have per-species values "
              f"(range: {X0_vals.min():.2e}–{X0_vals.max():.2e}; rest use config X0)")

    if 'metabolic_rate_exponent' in traits_df.columns:
        # Count non-NaN entries; exponent is unconstrained but should be ~-0.25
        n_present = traits_df['metabolic_rate_exponent'].notna().sum()
        bX_vals = traits_df['metabolic_rate_exponent'].dropna().values
        print(f"  ✓ metabolic_rate_exponent: {n_present}/{n_spp} species have per-species values "
              f"(range: {bX_vals.min():.3f}–{bX_vals.max():.3f}; rest use config b_X)")

    # Check assimilation efficiencies if present; report coverage
    for col in ['assimilation_plant', 'assimilation_animal']:
        if col in traits_df.columns:
            n_present = traits_df[col].notna().sum()
            e_vals = traits_df[col].dropna().values
            if np.any((e_vals < 0) | (e_vals > 1)):
                raise ValidationError(
                    f"'{col}' must be in [0, 1]. Found range: {e_vals.min():.3f}–{e_vals.max():.3f}"
                )
            cfg_key = 'e_plant' if col == 'assimilation_plant' else 'e_animal'
            print(f"  ✓ {col}: {n_present}/{n_spp} species have per-species values "
                  f"(range: {e_vals.min():.3f}–{e_vals.max():.3f}; rest use config {cfg_key})")

    # Check Hill exponent if present
    if 'hill_exponent' in traits_df.columns:
        q = traits_df['hill_exponent'].values
        if np.any(q <= 0):
            raise ValidationError(
                f"Hill exponent must be positive. Found: {np.unique(q)}"
            )
        print(f"  ✓ Hill exponent: {q.min():.2f}–{q.max():.2f}")
    
    if expected_n_species is not None and n_spp != expected_n_species:
        raise ValidationError(
            f"Species count mismatch: traits has {n_spp} species "
            f"but expected {expected_n_species}"
        )
    
    return traits_df


def validate_inputs(env_df: pd.DataFrame, adj_mat: np.ndarray, 
                    traits_df: pd.DataFrame) -> bool:
    """
    Comprehensive cross-file validation.
    
    Checks:
      - Species count consistency
      - Basal species count
      - Reasonable parameter ranges
      - Logical food-web structure
    
    Returns:
        True if all checks pass
    """
    print("\n[VALIDATION] Cross-file consistency checks...")
    
    n_species_adj = adj_mat.shape[0]
    n_species_traits = len(traits_df)
    n_cells = len(env_df)
    
    # Check species counts match
    if n_species_adj != n_species_traits:
        raise ValidationError(
            f"Species count mismatch:\n"
            f"  Adjacency matrix: {n_species_adj}×{n_species_adj}\n"
            f"  Traits: {n_species_traits} rows\n"
            f"  These must match."
        )
    
    # Check basal species are leaf nodes in food web
    n_basal = (traits_df['is_basal'] == 1).sum()
    basal_idx = np.where(traits_df['is_basal'] == 1)[0]

    # Basal species should only be consumed, not consume.
    # Check columns (consumer axis): basal species must not appear as consumers.
    basal_as_consumers = np.sum(adj_mat[:, basal_idx])
    if basal_as_consumers > 0:
        raise ValidationError(
            f"Basal species should not consume. Found {basal_as_consumers} links from basal spp.\n"
            f"  (Basal indices: {basal_idx})"
        )
    
    print(f"  ✓ Species count consistent: {n_species_adj} species, {n_basal} basal")
    
    # Check cells and environment
    if n_cells < 1:
        raise ValidationError(f"Must have at least 1 cell, found {n_cells}")
    print(f"  ✓ Environment: {n_cells} cell(s)")
    
    # Check biomass > 0 for at least one species
    B0_total = traits_df['initial_biomass_g_per_m2'].sum()
    if B0_total <= 0:
        raise ValidationError("Total initial biomass is zero; simulation will not run.")
    print(f"  ✓ Total initial biomass: {B0_total:.2e} g/m²")
    
    # Check trophic structure
    n_links = np.sum(adj_mat)
    if n_links == 0:
        raise ValidationError("No feeding links in adjacency matrix; trivial simulation.")
    
    # Compute food-web properties
    mean_degree_in = np.mean(np.sum(adj_mat, axis=0))  # in-degree (# prey)
    mean_degree_out = np.mean(np.sum(adj_mat, axis=1))  # out-degree (# predators)
    
    print(f"  ✓ Food web: {n_links} links, mean prey/species: {mean_degree_in:.2f}, "
          f"mean predators/species: {mean_degree_out:.2f}")
    
    # Warning: very weak connectivity
    connectance = n_links / (n_species_adj ** 2)
    if connectance < 0.01:
        print(f"  ⚠ Very low connectance: {connectance:.4f} "
              f"(sparse food web; may exclude many species)")
    
    print("  ✓ All validation checks passed!")
    return True


def check_parameter_completeness(config: Dict) -> bool:
    """
    Verify that required configuration parameters are present.
    
    Parameters:
        config: configuration dictionary loaded by read_config()
    
    Returns:
        True if all required parameters present
    """
    print("\n[CONFIG] Checking parameter completeness...")
    
    required_params = {
        # Vegetation growth (NPP-driven, vegetation.md equation)
        'psi':                 'carbon to wet matter conversion factor (g wet / g C)',
        'f_struct_default':    'default fractional NPP allocation to structural tissue',
        'alpha_herbs_default': 'half-saturation constant for herb/tree competitive partition (g/m2)',
        # Allometric — consumers and metabolic loss only (no basal logistic growth)
        'X0': 'metabolic rate normalization',
        'b_X': 'metabolic exponent',
        'a0': 'attack rate normalization',
        'b_a_prey': 'attack rate prey exponent',
        'b_a_pred': 'attack rate predator exponent',
        'h0': 'handling time normalization',
        'b_h_prey': 'handling time prey exponent',
        'b_h_pred': 'handling time predator exponent',
        # Functional response
        'q_hill': 'Hill exponent',
        'interference': 'consumer interference',
        # Efficiency
        'e_plant': 'plant assimilation efficiency',
        'e_animal': 'animal assimilation efficiency',
        # Temperature
        'use_temperature': 'use temperature dependence',
        'T0_K': 'reference temperature',
        'k_B': 'Boltzmann constant',
        'E_a': 'activation energy',
        # Numerical
        'ext_threshold': 'extinction threshold',
        'extinction_timescale': 'extinction timescale',
    }
    
    missing = []
    for param, desc in required_params.items():
        if param not in config:
            missing.append(f"  {param}: {desc}")
    
    if missing:
        raise ValidationError(
            f"Missing configuration parameters:\n" + "\n".join(missing)
        )
    
    # Check parameter ranges
    warnings = []
    
    if config['X0'] <= 0:
        warnings.append(f"X0={config['X0']} should be positive (metabolic rate)")
    
    if config['a0'] <= 0:
        warnings.append(f"a0={config['a0']} should be positive (attack rate)")
    
    if not (-1 < config['b_X'] < 0):
        warnings.append(f"b_X={config['b_X']} is outside typical range [-1, 0]")
    
    if config['q_hill'] <= 0:
        warnings.append(f"q_hill={config['q_hill']} should be positive")
    
    if not (0 < config['e_plant'] < 1):
        warnings.append(f"e_plant={config['e_plant']} should be in (0, 1)")
    
    if not (0 < config['e_animal'] < 1):
        warnings.append(f"e_animal={config['e_animal']} should be in (0, 1)")
    
    if config['e_animal'] <= config['e_plant']:
        warnings.append(f"e_animal ({config['e_animal']}) should be > e_plant ({config['e_plant']})")
    
    if config['T0_K'] < 260 or config['T0_K'] > 330:
        warnings.append(f"T0_K={config['T0_K']} is outside typical range [260, 330] K")
    
    if config['k_B'] <= 0:
        warnings.append(f"k_B={config['k_B']} (Boltzmann constant) should be positive")
    
    if config['E_a'] < 0:
        warnings.append(f"E_a={config['E_a']} (activation energy) should be non-negative")
    
    if config['ext_threshold'] < 0 or config['ext_threshold'] > 1:
        warnings.append(f"ext_threshold={config['ext_threshold']} should be in [0, 1]")
    
    print(f"  ✓ All {len(required_params)} required parameters present")
    
    if warnings:
        print("\n  Parameter range warnings:")
        for w in warnings:
            print(f"    ⚠ {w}")
    
    return True


def print_summary(env_df: pd.DataFrame, adj_mat: np.ndarray, 
                  traits_df: pd.DataFrame, config: Dict):
    """Print a summary of loaded data."""
    print("\n" + "="*70)
    print("INPUT SUMMARY")
    print("="*70)
    
    n_spp = len(traits_df)
    n_cells = len(env_df)
    n_basal = (traits_df['is_basal'] == 1).sum()
    
    print(f"\nStructure:")
    print(f"  Species: {n_spp} ({n_basal} basal, {n_spp - n_basal} consumers)")
    print(f"  Cells: {n_cells}")
    print(f"  Links: {np.sum(adj_mat)} (connectance: {np.sum(adj_mat)/(n_spp**2):.3f})")
    
    print(f"\nSpecies range:")
    M = traits_df['body_mass_g'].values
    print(f"  Body mass: {M.min():.2e}–{M.max():.2e} g")
    B0 = traits_df['initial_biomass_g_per_m2'].values
    print(f"  Initial biomass: {B0[B0>0].min():.2e}–{B0.max():.2e} g/m²")
    
    print(f"\nEnvironment:")
    T = env_df['temperature_K'].values
    print(f"  Temperature: {T.min():.2f}–{T.max():.2f} K ({T.min()-273.15:.1f}–{T.max()-273.15:.1f}°C)")
    
    print(f"\nModel settings:")
    print(f"  Temp dependence: {config['use_temperature']}")
    print(f"  q_hill: {config['q_hill']}")
    print(f"  e_plant: {config['e_plant']}, e_animal: {config['e_animal']}")
    
    print("="*70 + "\n")
