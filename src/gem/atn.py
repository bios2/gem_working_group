"""
ATN model process functions — modular version.

All functions here are PURE:
  - they depend only on numpy;
  - they read no files and write no files;
  - they hold no class state and no global state.

Each function receives all the arrays/scalars it needs as explicit arguments
and returns an array. This makes every process independently testable
(see docs/processes_implementation_specification.md §5.3).

Shape conventions:
  S       = number of species
  n_cells = number of spatial cells
  B       : (n_cells, S)            biomass per cell and species (g/m²)
  a_ij,   : (n_cells, S, S) or (S, S) — axis 1 = prey, axis 2 = predator
  h_ij    : same
  E       : (S, S) — assimilation efficiency (axis 0 = prey, axis 1 = predator)
  F       : (n_cells, S, S) — functional response
"""
from typing import Tuple

import numpy as np
from numpy.typing import NDArray


# --------------------------------------------------------------------------- #
# 1. Temperature factor (Boltzmann–Arrhenius)                                  #
# --------------------------------------------------------------------------- #
def temperature_factor(T_K: NDArray[np.float64],
                       T0_K: float, E_a: float, k_B: float) -> NDArray[np.float64]:
    """
    Boltzmann–Arrhenius multiplicative factor, one value per cell.

    f(T) = exp(-E_a * (T0 - T) / (k_B * T * T0))

    Inputs:
        T_K  : per-cell temperatures in Kelvin, shape (n_cells,)
        T0_K : reference temperature (K)
        E_a  : activation energy (eV)
        k_B  : Boltzmann constant (eV/K)

    Returns:
        tf : shape (n_cells,)
    """
    return np.exp(-E_a * (T0_K - T_K) / (k_B * T_K * T0_K))


# --------------------------------------------------------------------------- #
# 2. Allometric matrices (constant in time)                                    #
# --------------------------------------------------------------------------- #
def attack_rate_matrix(M: NDArray[np.float64], adj_mat: NDArray[np.float64],
                       a0: float, b_a_prey: float, b_a_pred: float) -> NDArray[np.float64]:
    """
    Base attack-rate matrix: a_ij = a0 * M_prey^b_prey * M_pred^b_pred * adj_mat.

    M       : (S,) — body masses
    adj_mat : (S, S) — adj_mat[prey, predator] = 1 if trophic link exists

    Returns: (S, S)
    """
    assert M.ndim == 1 and adj_mat.ndim == 2
    assert adj_mat.shape == (M.size, M.size)
    M_prey = M[:, np.newaxis]    # (S, 1)
    M_pred = M[np.newaxis, :]    # (1, S)
    return a0 * np.power(M_prey, b_a_prey) * np.power(M_pred, b_a_pred) * adj_mat


def handling_time_matrix(M: NDArray[np.float64], adj_mat: NDArray[np.float64],
                         h0: float, b_h_prey: float, b_h_pred: float) -> NDArray[np.float64]:
    """
    Base handling-time matrix h_ij. Same structure as attack_rate_matrix.
    """
    assert M.ndim == 1 and adj_mat.ndim == 2
    assert adj_mat.shape == (M.size, M.size)
    M_prey = M[:, np.newaxis]
    M_pred = M[np.newaxis, :]
    return h0 * np.power(M_prey, b_h_prey) * np.power(M_pred, b_h_pred) * adj_mat


def metabolic_base_rate(X0_spp: NDArray[np.float64],
                        M: NDArray[np.float64],
                        bX_spp: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Allometric base of the metabolic rate: X_i = X0_i * M_i^bX_i  (no temperature).

    Returns: (S,)
    """
    assert X0_spp.shape == M.shape == bX_spp.shape
    return X0_spp * np.power(M, bX_spp)


def assimilation_matrix(is_basal: NDArray[np.bool_],
                        e_plant_spp: NDArray[np.float64],
                        e_animal_spp: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Assimilation efficiency matrix E[prey, predator].

    E[i, j] = e_plant_spp[j]  if prey i is basal
            = e_animal_spp[j] otherwise
    """
    assert is_basal.shape == e_plant_spp.shape == e_animal_spp.shape
    return np.where(is_basal[:, np.newaxis],
                    e_plant_spp[np.newaxis, :],
                    e_animal_spp[np.newaxis, :]).astype(float)


# --------------------------------------------------------------------------- #
# 3. Functional response, feeding gain, predation loss                         #
# --------------------------------------------------------------------------- #
def functional_response(B: NDArray[np.float64],
                        a_ij: NDArray[np.float64],
                        h_ij: NDArray[np.float64],
                        q_hill: float,
                        interference: float,
                        ext_threshold: float
                        ) -> Tuple[NDArray[np.float64], NDArray[np.bool_]]:
    """
    Generalized Holling type II functional response, vectorised over cells.

        F[c, i, j] = a_ij * B_i^q  /  (1 + c * B_j + Σ_k h_kj * a_kj * B_k^q)

    Inputs:
        B    : (n_cells, S) — biomass (will be clamped to 0)
        a_ij : (n_cells, S, S) or (1, S, S)
        h_ij : same shape as a_ij
        q_hill, interference, ext_threshold : scalars

    Returns:
        F       : (n_cells, S, S)
        extinct : (n_cells, S) boolean — species below the extinction threshold
    """
    assert B.ndim == 2
    B_safe = np.maximum(B, 0.0)
    B_q = np.power(B_safe, q_hill)                                  # (n_cells, S)

    numerator = a_ij * B_q[:, :, np.newaxis]                        # (n_cells, S, S)
    handling_sum = (h_ij * a_ij * B_q[:, :, np.newaxis]).sum(axis=1)  # (n_cells, S)
    denominator = 1.0 + interference * B_safe + handling_sum        # (n_cells, S)

    F = numerator / denominator[:, np.newaxis, :]                   # (n_cells, S, S)

    # An extinct prey contributes to no functional response (row i of F = 0).
    extinct = B_safe < ext_threshold
    F = F * (~extinct)[:, :, np.newaxis]
    return F, extinct


def feeding_gain(B: NDArray[np.float64],
                 F: NDArray[np.float64],
                 E: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Assimilated energy gain for each consumer:
        gain[c, j] = B[c, j] * Σ_prey E[prey, j] * F[c, prey, j]

    B : (n_cells, S),  F : (n_cells, S, S),  E : (S, S)
    Returns: (n_cells, S)
    """
    assert B.ndim == 2 and F.ndim == 3 and E.ndim == 2
    return B * (E[np.newaxis, :, :] * F).sum(axis=1)


def predation_loss(B: NDArray[np.float64],
                   F: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Biomass loss of each prey due to predation:
        loss[c, i] = Σ_predator B[c, predator] * F[c, i, predator]

    Returns: (n_cells, S)
    """
    assert B.ndim == 2 and F.ndim == 3
    return (F * B[:, np.newaxis, :]).sum(axis=-1)


# --------------------------------------------------------------------------- #
# 4. Total derivative dB/dt — process orchestration                            #
# --------------------------------------------------------------------------- #
def derivatives(B: NDArray[np.float64],
                t: float,
                params: dict,
                vegetation) -> NDArray[np.float64]:
    """
    Compute dB/dt for all cells and all species simultaneously.

    Composition:
        dB/dt = feeding_gain + vegetation_growth - metabolism*B - predation_loss

    For basal species:    feeding_gain = 0 (no prey columns via adj_mat).
    For consumer species: vegetation_growth = 0 (vegetation model returns 0).

    params : dict built by atn_data.build_atn_params(), containing:
        base_a, base_h, base_X, E, T_K,
        use_temperature, T0_K, k_B, E_a,
        q_hill, interference, ext_threshold, extinction_timescale
    vegetation : PlantVegetationModel instance (interface: growth_all_cells(B_safe))
    """
    assert B.ndim == 2, "B must be (n_cells, S)"
    B_safe = np.maximum(B, 0.0)

    # --- 1. Temperature scaling (vectorised over cells) ---
    if params['use_temperature']:
        tf = temperature_factor(params['T_K'],
                                params['T0_K'], params['E_a'], params['k_B'])
        tf_3d = tf[:, np.newaxis, np.newaxis]                       # (n_cells, 1, 1)
        a_ij = params['base_a'][np.newaxis, :, :] * tf_3d           # (n_cells, S, S)
        h_ij = params['base_h'][np.newaxis, :, :] * tf_3d
        X = params['base_X'][np.newaxis, :] * tf[:, np.newaxis]     # (n_cells, S)
    else:
        a_ij = params['base_a'][np.newaxis, :, :]                   # (1, S, S)
        h_ij = params['base_h'][np.newaxis, :, :]
        X = params['base_X'][np.newaxis, :]                         # (1, S)

    # --- 2. Functional response + extinction mask ---
    F, extinct = functional_response(
        B_safe, a_ij, h_ij,
        params['q_hill'], params['interference'], params['ext_threshold']
    )

    # --- 3. Energy budget ---
    gain = feeding_gain(B_safe, F, params['E'])
    loss = predation_loss(B_safe, F)
    G    = vegetation.growth_all_cells(B_safe)                      # (n_cells, S), zero for non-basal

    dBdt = gain + G - X * B_safe - loss

    # --- 4. Smooth decay for species below the extinction threshold ---
    dBdt[extinct] = -B_safe[extinct] / params['extinction_timescale']
    return dBdt
