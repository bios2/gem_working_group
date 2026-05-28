"""
Input loading and parameter preparation for the ATN model.

This layer does ONE job: read the 4 files (config, env, adj, traits), validate
them, and precompute every array that is constant in time (base_a, base_h,
base_X, assimilation matrix E, temperatures, basal/consumer indices).

The result is:
  - params : dict consumed directly by gem_working_group.src.gem.atn_processes.derivatives()
  - vegetation : a PlantVegetationModel instance shared across the run

No ODE integration happens here. No process computation happens here.
"""
from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray
import pandas as pd

from atn_io import (
    read_config, read_env_matrix, read_adjacency_matrix, read_traits,
    validate_inputs, check_parameter_completeness,
)
from gem_working_group.src.gem.atn import (
    attack_rate_matrix, handling_time_matrix,
    metabolic_base_rate, assimilation_matrix,
)
from vegetation_model import PlantVegetationModel


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #
def _per_species_param(traits_df: pd.DataFrame, col_name: str,
                       default_value: float, n_species: int) -> NDArray[np.float64]:
    """
    Read a per-species column from traits_df, with a global default fallback.

    If the column is absent, return an (S,) array filled with default_value.
    If it exists but has NaNs, replace NaNs with default_value.
    """
    if col_name in traits_df.columns:
        return traits_df[col_name].fillna(default_value).values.astype(float)
    return np.full(n_species, default_value, dtype=float)


# --------------------------------------------------------------------------- #
# 1. Raw file reading + validation                                             #
# --------------------------------------------------------------------------- #
def load_inputs(env_file: str, adj_file: str, traits_file: str, config_file: str
                ) -> Tuple[Dict, pd.DataFrame, NDArray[np.float64], pd.DataFrame]:
    """
    Read the 4 input files and apply the cross-file validation from atn_io.

    Returns: (config, env_df, adj_mat, traits_df)
    """
    config    = read_config(config_file)
    env_df, _ = read_env_matrix(env_file)
    adj_mat   = read_adjacency_matrix(adj_file)
    traits_df = read_traits(traits_file)

    validate_inputs(env_df, adj_mat, traits_df)
    check_parameter_completeness(config)
    return config, env_df, adj_mat, traits_df


# --------------------------------------------------------------------------- #
# 2. Build the precomputed-parameter dict                                      #
# --------------------------------------------------------------------------- #
def build_atn_params(adj_mat: NDArray[np.float64], traits_df: pd.DataFrame,
                     env_df: pd.DataFrame, config: Dict
                     ) -> Tuple[Dict, PlantVegetationModel]:
    """
    Build the `params` dict consumed by the process functions.

    Precomputes everything that is independent of time:
      - base_a, base_h, base_X (allometric matrices/vectors)
      - E (assimilation matrix)
      - T_K (per-cell temperatures)
      - basal_idx / consumer_idx indices, is_basal mask

    Returns: (params, vegetation)
    """
    n_species = len(traits_df)
    n_cells   = len(env_df)

    # --- Per-type indices/masks ---
    is_basal     = (traits_df['is_basal'].values == 1)
    basal_idx    = np.where(is_basal)[0]
    consumer_idx = np.where(~is_basal)[0]

    # --- Per-species parameters (with fallback to global defaults) ---
    X0_spp       = _per_species_param(traits_df, 'metabolic_rate_base',     config['X0'],       n_species)
    bX_spp       = _per_species_param(traits_df, 'metabolic_rate_exponent', config['b_X'],      n_species)
    e_plant_spp  = _per_species_param(traits_df, 'assimilation_plant',      config['e_plant'],  n_species)
    e_animal_spp = _per_species_param(traits_df, 'assimilation_animal',     config['e_animal'], n_species)

    # --- Allometric precomputations (constant in time) ---
    M = traits_df['body_mass_g'].values.astype(float)
    base_a = attack_rate_matrix(M, adj_mat,
                                config['a0'], config['b_a_prey'], config['b_a_pred'])
    base_h = handling_time_matrix(M, adj_mat,
                                  config['h0'], config['b_h_prey'], config['b_h_pred'])
    base_X = metabolic_base_rate(X0_spp, M, bX_spp)
    E      = assimilation_matrix(is_basal, e_plant_spp, e_animal_spp)

    # --- Per-cell temperatures ---
    T_K = env_df['temperature_K'].values.astype(float)

    # --- Vegetation sub-model (reused as-is) ---
    vegetation = PlantVegetationModel(traits_df, env_df, config)

    params: Dict = {
        # Dimensions
        'n_species'   : n_species,
        'n_cells'     : n_cells,
        # Indices / masks
        'is_basal'    : is_basal,
        'basal_idx'   : basal_idx,
        'consumer_idx': consumer_idx,
        # Per-species (kept around for diagnostics)
        'M'           : M,
        'X0_spp'      : X0_spp,
        'bX_spp'      : bX_spp,
        'e_plant_spp' : e_plant_spp,
        'e_animal_spp': e_animal_spp,
        # Precomputed arrays (consumed by derivatives)
        'base_a'      : base_a,
        'base_h'      : base_h,
        'base_X'      : base_X,
        'E'           : E,
        'T_K'         : T_K,
        # Scalar parameters (copied from config so params is self-contained)
        'use_temperature'     : bool(config['use_temperature']),
        'T0_K'                : float(config['T0_K']),
        'k_B'                 : float(config['k_B']),
        'E_a'                 : float(config['E_a']),
        'q_hill'              : float(config['q_hill']),
        'interference'        : float(config['interference']),
        'ext_threshold'       : float(config['ext_threshold']),
        'extinction_timescale': float(config['extinction_timescale']),
    }

    print(f"ATN params built: {n_species} species, {n_cells} cells.")
    return params, vegetation


# --------------------------------------------------------------------------- #
# 3. Initial conditions                                                         #
# --------------------------------------------------------------------------- #
def initial_biomass(traits_df: pd.DataFrame, n_cells: int,
                    noise_frac: float = 0.01) -> NDArray[np.float64]:
    """
    Build the initial-biomass matrix (n_cells, n_species).

    Reads the 'initial_biomass_g_per_m2' column of traits_df, replicates it
    across all cells, then adds a relative gaussian noise of `noise_frac` to
    break inter-cell symmetry (a numpy seed must be set upstream).

    noise_frac=0.01 → gaussian noise with amplitude 1 % of the initial value.
    """
    n_species = len(traits_df)
    B0 = np.zeros((n_cells, n_species))
    B0[:, :] = traits_df['initial_biomass_g_per_m2'].values.astype(float)
    B0 += np.random.normal(0.0, noise_frac * np.maximum(B0, 1e-6), B0.shape)
    return np.maximum(B0, 0.0)
