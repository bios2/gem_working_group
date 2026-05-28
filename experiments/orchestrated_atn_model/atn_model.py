"""
Spatially explicit unscaled ATN model (following Binzer et al. 2016).
Implementation of Section 8 from ATN_model_spatiotemporal_formulas_parameters.Rmd

Option B time integration (processes_implementation_specification.md §7):
derivatives() returns dB/dt for all cells at once; run_all_cells() advances the
state with a fixed-step RK4 loop — no scipy.odeint, no per-cell Python loop.

All per-species loops have been replaced with vectorised numpy matrix operations.
"""
from pathlib import Path

import numpy as np
from numpy.typing import NDArray
import pandas as pd
from typing import Dict
from vegetation_model import PlantVegetationModel



class ATNModel:
    """
    Unscaled ATN model with spatial grid cells and optional temperature dependence.


    All cells are integrated simultaneously. derivatives() accepts B: (n_cells, S)
    and returns dBdt: (n_cells, S) without any Python loops over cells or species.
    """

    def __init__(self, adj_mat: NDArray[np.float64], traits_df: pd.DataFrame,
                 env_df: pd.DataFrame, config: Dict):
        """
        Initialize ATN model and precompute temperature-independent allometric arrays.

        Parameters:
            adj_mat:    (S, S) binary matrix, adj_mat[prey, consumer] = 1
            traits_df:  species traits table
            env_df:     per-cell environment (temperature_K, NPP, x, y)
            config:     model parameters dict
        """
        self.adj_mat = adj_mat
        self.traits  = traits_df
        self.env     = env_df
        self.cfg     = config

        self.n_species = len(traits_df)
        self.n_cells   = len(env_df)

        self.basal_idx    = np.where(traits_df['is_basal'] == 1)[0]
        self.consumer_idx = np.where(traits_df['is_basal'] == 0)[0]


        self.vegetation = PlantVegetationModel(traits_df, env_df, config)
        self.herb_idx = self.vegetation.herb_idx
        self.tree_idx = self.vegetation.tree_idx

        # Per-species metabolic parameters (fall back to config global if absent/NaN)
        x0_default = config['X0']
        bx_default = config['b_X']
        if 'metabolic_rate_base' in traits_df.columns:
            self.X0_spp = traits_df['metabolic_rate_base'].fillna(x0_default).values.astype(float)
        else:
            self.X0_spp = np.full(self.n_species, x0_default)
        if 'metabolic_rate_exponent' in traits_df.columns:
            self.bX_spp = traits_df['metabolic_rate_exponent'].fillna(bx_default).values.astype(float)
        else:
            self.bX_spp = np.full(self.n_species, bx_default)

        # Per-species assimilation efficiencies
        e_plant_default  = config['e_plant']
        e_animal_default = config['e_animal']
        if 'assimilation_plant' in traits_df.columns:
            self.e_plant_spp = traits_df['assimilation_plant'].fillna(e_plant_default).values.astype(float)
        else:
            self.e_plant_spp = np.full(self.n_species, e_plant_default)
        if 'assimilation_animal' in traits_df.columns:
            self.e_animal_spp = traits_df['assimilation_animal'].fillna(e_animal_default).values.astype(float)
        else:
            self.e_animal_spp = np.full(self.n_species, e_animal_default)

        # Temperature parameters kept as scalars for the Boltzmann-Arrhenius formula
        self.T0_K = config['T0_K']
        self.k_B  = config['k_B']
        self.E_a  = config['E_a']
        self.q_hill = config['q_hill']
        self.ext_threshold = config['ext_threshold']

        # Pre-extract temperatures for vectorised derivatives()
        self.T_K = env_df['temperature_K'].values.astype(float)  # (n_cells,)

        # ------------------------------------------------------------------ #
        # Precompute temperature-independent allometric bases (shape S × S)  #
        # Temperature scaling is applied at runtime in derivatives()          #
        # ------------------------------------------------------------------ #
        M = traits_df['body_mass_g'].values.astype(float)   # (S,)
        self.M = M

        M_prey = M[:, np.newaxis]   # (S, 1) — prey body mass, varies by row
        M_pred = M[np.newaxis, :]   # (1, S) — predator body mass, varies by column

        # base_a[prey, consumer] = a0 * M_prey^b_prey * M_pred^b_pred * adj_mat
        self.base_a = (config['a0']
                       * np.power(M_prey, config['b_a_prey'])
                       * np.power(M_pred, config['b_a_pred'])
                       * adj_mat).astype(float)   # (S, S)

        # base_h[prey, consumer] = h0 * M_prey^b_prey * M_pred^b_pred * adj_mat
        self.base_h = (config['h0']
                       * np.power(M_prey, config['b_h_prey'])
                       * np.power(M_pred, config['b_h_pred'])
                       * adj_mat).astype(float)   # (S, S)

        # base_X[species] = X0_i * M_i^bX_i  (no temperature yet)
        self.base_X = (self.X0_spp * np.power(M, self.bX_spp)).astype(float)  # (S,)

        # ------------------------------------------------------------------ #
        # Precompute assimilation efficiency matrix E[prey, consumer]         #
        # E[i, j] = efficiency of consumer j eating prey i                    #
        # ------------------------------------------------------------------ #
        is_basal = (traits_df['is_basal'].values == 1)   # (S,) bool
        # np.where broadcasts: is_basal[:, newaxis] selects row-wise (prey type)
        self.E = np.where(
            is_basal[:, np.newaxis],               # (S, 1) — True if prey i is basal
            self.e_plant_spp[np.newaxis, :],       # (1, S) — plant assimilation for consumer j
            self.e_animal_spp[np.newaxis, :]       # (1, S) — animal assimilation for consumer j
        ).astype(float)                            # (S, S)

        print(f"ATN Model initialized: {self.n_species} species, {self.n_cells} cells")

    def derivatives(self, B: NDArray[np.float64], t: float) -> NDArray[np.float64]:
        """
        Compute dB/dt for all species in all cells simultaneously. No Python loops.

        B: (n_cells, S)   — current biomass, all cells
        t: float          — current time (reserved; vegetation growth is time-invariant)

        Returns dBdt: (n_cells, S)

        Equation structure (same formula for basal and consumers — the zeros fall
        out of the math because adj_mat and the vegetation model enforce them):
            dB/dt = feeding_gain + vegetation_growth - metabolic_loss - predation_loss
        For basal species: feeding_gain = 0 (no adj_mat columns)
        For consumers:     vegetation_growth = 0 (PlantVegetationModel returns 0)
        """
        assert B.shape == (self.n_cells, self.n_species)
        B_safe = np.maximum(B, 0.0)   # (n_cells, S)

        # ---- Temperature scaling (vectorised over cells) ----
        if self.cfg['use_temperature']:
            # Boltzmann-Arrhenius factor, one value per cell
            temp_factor = np.exp(
                -self.E_a * (self.T0_K - self.T_K)
                / (self.k_B * self.T_K * self.T0_K)
            )                                                    # (n_cells,)
            tf = temp_factor[:, np.newaxis, np.newaxis]          # (n_cells, 1, 1)
            a_ij = self.base_a[np.newaxis, :, :] * tf           # (n_cells, S, S)
            h_ij = self.base_h[np.newaxis, :, :] * tf           # (n_cells, S, S)
            X    = self.base_X[np.newaxis, :] * temp_factor[:, np.newaxis]  # (n_cells, S)
        else:
            # No temperature: broadcast (1, S, S) and (1, S) over cells at zero cost
            a_ij = self.base_a[np.newaxis, :, :]   # (1, S, S)
            h_ij = self.base_h[np.newaxis, :, :]   # (1, S, S)
            X    = self.base_X[np.newaxis, :]       # (1, S)

        # ---- Functional response F[cell, prey, consumer] ----
        # F[c, i, j] = a_ij * B_i^q / (1 + c_j*B_j + Σ_k h_kj*a_kj*B_k^q)
        B_q = np.power(B_safe, self.q_hill)                      # (n_cells, S)

        # B_q[:, :, newaxis] broadcasts prey biomass over the consumer axis
        numerator = a_ij * B_q[:, :, np.newaxis]                 # (n_cells, S, S)

        # handling_sum[cell, consumer] = Σ_prey h*a*B_prey^q
        # axis=1 sums over the prey axis of the (n_cells, S_prey, S_consumer) tensor
        handling_sum = (h_ij * a_ij * B_q[:, :, np.newaxis]).sum(axis=1)  # (n_cells, S)

        denominator = (1.0 + self.cfg['interference'] * B_safe
                       + handling_sum)                            # (n_cells, S) — per consumer

        # denominator[:, newaxis, :] broadcasts the consumer denominator over the prey axis
        F = numerator / denominator[:, np.newaxis, :]            # (n_cells, S, S)

        # Zero out prey rows where species is locally extinct
        extinct = B_safe < self.ext_threshold                    # (n_cells, S)
        F *= (~extinct)[:, :, np.newaxis]                        # broadcast over consumer axis

        # ---- Predation / herbivory loss ----
        # predation_loss[cell, prey] = Σ_consumer B[consumer] * F[prey, consumer]
        # B_safe[:, newaxis, :] → (n_cells, 1, S): consumer biomass broadcast over prey axis
        # .sum(axis=-1) sums out the consumer axis
        predation_loss = (F * B_safe[:, np.newaxis, :]).sum(axis=-1)   # (n_cells, S)

        # ---- Assimilated feeding gain ----
        # feeding_gain[cell, j] = B[j] * Σ_prey E[prey,j] * F[prey,j]
        # (E * F) shape: (n_cells, S_prey, S_consumer); .sum(axis=1) sums out prey
        feeding_gain = B_safe * (self.E[np.newaxis, :, :] * F).sum(axis=1)  # (n_cells, S)

        # ---- Vegetation growth (basal species only, zero for consumers) ----
        G = self.vegetation.growth_all_cells(B_safe)             # (n_cells, S)

        # ---- dB/dt ----
        dBdt = feeding_gain + G - X * B_safe - predation_loss    # (n_cells, S)

        # Smooth extinction: drive sub-threshold species toward zero
        dBdt[extinct] = -B_safe[extinct] / self.cfg['extinction_timescale']

        return dBdt

    def run_all_cells(self, B_initial: NDArray[np.float64],
                      t_eval: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        Option B: fixed-step RK4 over all cells simultaneously. No per-cell loop.

        B_initial: (n_cells, S)  — initial biomass
        t_eval:    (n_tp,)       — output time points (uniform spacing assumed)

        Returns B_traj: (n_tp, n_cells, S)

        RK4 takes four derivative evaluations per step, each covering all cells
        at once. Each step costs ~4 × one vectorised derivatives() call instead of
        n_cells × (odeint overhead + adaptive substeps).

        Note: accuracy depends on the step size dt = t_eval[1] - t_eval[0]. ATN
        dynamics can be moderately stiff; verify against a reference run on the
        test data before trusting production results.
        """
        n_tp  = len(t_eval)
        B_traj = np.zeros((n_tp, self.n_cells, self.n_species))
        B = np.maximum(B_initial.astype(float), 0.0)
        B_traj[0] = B

        for i in range(1, n_tp):
            dt = float(t_eval[i] - t_eval[i - 1])
            t  = float(t_eval[i - 1])

            k1 = self.derivatives(B,                                    t)
            k2 = self.derivatives(np.maximum(B + 0.5 * dt * k1, 0.0), t + 0.5 * dt)
            k3 = self.derivatives(np.maximum(B + 0.5 * dt * k2, 0.0), t + 0.5 * dt)
            k4 = self.derivatives(np.maximum(B + dt * k3,        0.0), t + dt)

            B = np.maximum(B + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4), 0.0)
            B_traj[i] = B
            print(f"  t = {t_eval[i]:.1f} / {t_eval[-1]:.1f} days", end='\r')

        print(f"\n✓ Completed {n_tp - 1} RK4 steps ({self.n_cells} cells).")
        return B_traj

    def save_consumer_output(self, B_traj: NDArray[np.float64],
                             t_eval: NDArray[np.float64], output_dir) -> None:
        """
        Save instantaneous dB/dt for all consumer species to atn_model.txt.

        Calls derivatives() once per time point (all cells in one call) — no inner
        cell loop.


        Columns: pixel_id, x, y, time, species, delta_biomass
          delta_biomass = dB_j/dt (g/m²/day) — full consumer dynamics
        """
        n_tp       = len(t_eval)

        n_consumers = len(self.consumer_idx)

        dBdt_arr = np.zeros((n_tp, self.n_cells, n_consumers))
        for t_idx in range(n_tp):

            dBdt = self.derivatives(B_traj[t_idx], t_eval[t_idx])  # (n_cells, S)
            dBdt_arr[t_idx] = dBdt[:, self.consumer_idx]


        cell_x = self.env['x'].values.astype(int)
        cell_y = self.env['y'].values.astype(int)

        t_rep    = np.repeat(t_eval, self.n_cells * n_consumers)
        cell_rep = np.tile(np.repeat(np.arange(self.n_cells), n_consumers), n_tp)
        x_rep    = np.tile(np.repeat(cell_x, n_consumers), n_tp)
        y_rep    = np.tile(np.repeat(cell_y, n_consumers), n_tp)
        sp_rep   = np.tile(self.consumer_idx, n_tp * self.n_cells)
        d_rep    = dBdt_arr.ravel()

        table = np.column_stack([cell_rep, x_rep, y_rep, t_rep, sp_rep, d_rep])
        np.savetxt(
            Path(output_dir) / 'atn_model.txt', table,
            fmt=['%d', '%d', '%d', '%.4f', '%d', '%.6e'],
            header='pixel_id x y time species delta_biomass',
            comments=''
        )
