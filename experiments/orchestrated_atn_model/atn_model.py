"""
Spatially explicit unscaled ATN model (following Binzer et al. 2016).
Implementation of Section 8 from ATN_model_spatiotemporal_formulas_parameters.Rmd

Each consumer species eats from resources following Holling Type II functional response.
Basal species (plants) grow logistically and are consumed by animals.
All rates scale allometrically with body mass and optionally with temperature.
"""
# Import required libraries
from pathlib import Path

import numpy as np  # numerical arrays and mathematics
from scipy.integrate import odeint  # ODE solver (4th-order Runge-Kutta adaptive)
import pandas as pd  # tabular data handling
from typing import Tuple, Dict  # type hints for function signatures
from vegetation_model import PlantVegetationModel  # NPP-driven plant growth model

class ATNModel:
    """
    Unscaled ATN model with spatial grid cells and optional temperature dependence.
    
    Each grid cell has independent local dynamics (no dispersal by default).
    Species biomass B_{i,g}(t) changes due to:
    - For basal: logistic growth, metabolic loss, herbivory
    - For consumers: feeding gain, metabolic loss, predation
    """
    
    def __init__(self, adj_mat: np.ndarray, traits_df: pd.DataFrame, 
                 env_df: pd.DataFrame, config: Dict):
        """
        Initialize ATN model.
        
        Parameters:
            adj_mat: (n_spp, n_spp) binary adjacency matrix (rows=resources, cols=consumers)
            traits_df: DataFrame with body_mass, is_basal, initial_biomass, metabolic params, etc.
            env_df: DataFrame with temperature and carrying capacities per cell
            config: Dict with parameters (see config.py for defaults)
        """
        # Store input data
        self.adj_mat = adj_mat  # food-web topology
        self.traits = traits_df  # species traits table
        self.env = env_df  # environmental conditions per cell
        self.cfg = config  # model parameters
        
        # Get dimensions
        self.n_species = len(traits_df)  # total number of species
        self.n_cells = len(env_df)  # total number of spatial cells
        
        # Identify basal (plant) and consumer (animal) species indices
        self.basal_idx = np.where(traits_df['is_basal'] == 1)[0]  # array of basal species indices
        self.consumer_idx = np.where(traits_df['is_basal'] == 0)[0]  # array of consumer species indices

        # Plant vegetation model owns herb/tree bookkeeping and NPP-driven growth.
        self.vegetation = PlantVegetationModel(traits_df, env_df, config)
        self.herb_idx = self.vegetation.herb_idx
        self.tree_idx = self.vegetation.tree_idx

        # Per-species metabolic normalization X0_i [day^-1] with global config fallback
        # Allows ectotherm/endotherm distinction (e.g. endotherms have higher X0)
        x0_default = config['X0']  # global fallback from config
        if 'metabolic_rate_base' in traits_df.columns:
            # Use per-species values; fill any NaN with global default
            self.X0_spp = traits_df['metabolic_rate_base'].fillna(x0_default).values.astype(float)
        else:
            # Column absent: uniform global default for all species
            self.X0_spp = np.full(self.n_species, x0_default)

        # Per-species metabolic exponent b_X_i [dimensionless] with global config fallback
        # Negative: larger species have lower mass-specific metabolic cost
        bx_default = config['b_X']  # global fallback from config
        if 'metabolic_rate_exponent' in traits_df.columns:
            # Use per-species values; fill any NaN with global default
            self.bX_spp = traits_df['metabolic_rate_exponent'].fillna(bx_default).values.astype(float)
        else:
            # Column absent: uniform global default for all species
            self.bX_spp = np.full(self.n_species, bx_default)

        # Per-species assimilation efficiency for plant prey [0-1] with global config fallback
        # Consumer j's efficiency when eating a basal (plant) resource
        e_plant_default = config['e_plant']  # global fallback from config
        if 'assimilation_plant' in traits_df.columns:
            # Use per-species values; fill any NaN with global default
            self.e_plant_spp = traits_df['assimilation_plant'].fillna(e_plant_default).values.astype(float)
        else:
            # Column absent: uniform global default for all species
            self.e_plant_spp = np.full(self.n_species, e_plant_default)

        # Per-species assimilation efficiency for animal prey [0-1] with global config fallback
        # Consumer j's efficiency when eating another consumer (animal) resource
        e_animal_default = config['e_animal']  # global fallback from config
        if 'assimilation_animal' in traits_df.columns:
            # Use per-species values; fill any NaN with global default
            self.e_animal_spp = traits_df['assimilation_animal'].fillna(e_animal_default).values.astype(float)
        else:
            # Column absent: uniform global default for all species
            self.e_animal_spp = np.full(self.n_species, e_animal_default)

        # Extract allometric rate constants from config (kept for reference; X0/b_X
        # overridden per-species by X0_spp/bX_spp above)
        self.X0 = config['X0']  # metabolic loss rate normalization constant
        self.bX = config['b_X']  # metabolic loss mass exponent (typically -0.25)
        self.a0 = config['a0']  # attack rate normalization constant
        self.ba_prey = config['b_a_prey']  # attack rate prey mass exponent
        self.ba_pred = config['b_a_pred']  # attack rate predator mass exponent
        self.h0 = config['h0']  # handling time normalization constant
        self.bh_prey = config['b_h_prey']  # handling time prey mass exponent
        self.bh_pred = config['b_h_pred']  # handling time predator mass exponent
        
        # Extract functional response parameters
        self.q_hill = config['q_hill']  # Hill exponent for functional response (2 ≈ Type II)
        # Temperature parameters (for Boltzmann-Arrhenius model)
        self.T0_K = config['T0_K']  # reference temperature in Kelvin (20°C = 293.15 K)
        self.k_B = config['k_B']  # Boltzmann constant in eV/K
        self.E_a = config['E_a']  # activation energy in eV
        
        # Extinction parameters
        self.ext_threshold = config['ext_threshold']  # biomass threshold for local extinction
        
        # Print initialization message
        print(f"ATN Model initialized: {self.n_species} species, {self.n_cells} cells")
    
    def _allometric_rate(self, rate_type: str, M_prey: np.ndarray, 
                         M_pred: np.ndarray, T_K: float) -> np.ndarray:
        """
        Compute allometric rate with optional temperature dependence.
        
        Formula: p(M_i, M_j, T) = p0 * M_prey^b_prey * M_pred^b_pred * 
                                  exp(-E(T0-T)/(k*T*T0))
        
        This follows Boltzmann-Arrhenius kinetics: rates increase exponentially with temp.
        """
        # Select parameters based on rate type requested
        if rate_type == 'attack':
            # Attack rate scales with predator mass^0.5 and prey mass^-0.5
            p0, b_prey, b_pred = self.a0, self.ba_prey, self.ba_pred
        elif rate_type == 'handling':
            # Handling time scales inversely (smaller predators slower at handling)
            p0, b_prey, b_pred = self.h0, self.bh_prey, self.bh_pred
        else:
            # Raise error if unknown rate type requested
            raise ValueError(f"Unknown rate type: {rate_type}")
        
        # Compute allometric part: p0 * M_prey^b_prey * M_pred^b_pred
        p_allom = p0 * np.power(M_prey, b_prey) * np.power(M_pred, b_pred)
        
        # Add temperature dependence if enabled in config
        if self.cfg['use_temperature']:
            # Boltzmann-Arrhenius factor: exp(-E(T0-T)/(k*T*T0))
            # Higher temp → larger exponential → faster rates
            temp_factor = np.exp(-self.E_a * (self.T0_K - T_K) / 
                                (self.k_B * T_K * self.T0_K))
            p_allom *= temp_factor
        
        return p_allom
    
    def _metabolic_rate(self, M: np.ndarray, T_K: float) -> np.ndarray:
        """Metabolic loss rate X_i = X0_i * M_i^(b_X_i) * temp_factor"""
        # Compute per-species metabolic rate using per-species X0 and b_X
        # X0_spp and bX_spp are either per-species from traits or the global config default
        X = self.X0_spp * np.power(M, self.bX_spp)
        # Apply temperature dependence if enabled
        if self.cfg['use_temperature']:
            temp_factor = np.exp(-self.E_a * (self.T0_K - T_K) /
                                (self.k_B * T_K * self.T0_K))
            X *= temp_factor
        return X
    
    def _functional_response(self, B: np.ndarray, j_idx: int) -> np.ndarray:
        """
        Unscaled Holling type II functional response for consumer j.
        
        F_ij = a_ij * B_i^q / (1 + c_j * B_j + sum_k h_kj * a_kj * B_k^q)
        
        Returns feeding rate of consumer j on each resource i.
        """
        # Get consumer j's biomass (used for interference term)
        B_j = B[j_idx]
        
        # Get attack rates a_ij for this consumer's resources
        a_ij = self.a_ij[:, j_idx]  # column j of attack rate matrix
        # Get handling times h_kj for this consumer's resources
        h_kj = self.h_ij[:, j_idx]  # column j of handling time matrix
        
        # Numerator: a_ij * B_i^q for each resource i
        numerator = a_ij * np.power(B, self.q_hill)
        
        # Denominator: 1 + interference term + handling time terms
        # Interference term: c_j * B_j (intraspecific competition among predators)
        # Handling time term: sum_k h_kj * a_kj * B_k^q (time spent handling all prey types)
        denominator = 1.0 + self.cfg['interference'] * B_j + \
                      np.sum(h_kj * a_ij * np.power(B, self.q_hill))
        
        # Compute functional response for each resource
        F_ij = numerator / denominator
        # Set feeding rate to 0 for extinct species (below extinction threshold)
        F_ij[B < self.ext_threshold] = 0
        
        return F_ij
    
    def derivatives(self, y: np.ndarray, t: float, cell_idx: int) -> np.ndarray:
        """
        Compute dB/dt for all species in a given cell.
        
        State vector: y = [B_1, B_2, ..., B_n] (biomass density of each species)
        Returns: dydt = [dB_1/dt, dB_2/dt, ..., dB_n/dt]
        """
        # Clamp biomass to non-negative (numerical solver might produce small negatives)
        B = np.maximum(y, 0)
        # Initialize derivative array
        dydt = np.zeros_like(B)
        
        # Get environmental conditions for this cell
        T_K = self.env.loc[cell_idx, 'temperature_K']  # temperature in Kelvin
        # Get body masses for all species
        M = self.traits['body_mass_g'].values
        
        # Precompute allometric rates for this cell (constant within one time step)
        # These matrices have shape (n_species, n_species)
        # Convention: rows = resource/prey (i), columns = consumer/predator (j)
        M_prey = M[:, np.newaxis]  # resource/prey mass varies by row
        M_pred = M[np.newaxis, :]  # consumer/predator mass varies by column

        # Compute attack rate matrix a_ij
        self.a_ij = self._allometric_rate('attack', M_prey, M_pred, T_K) * self.adj_mat
        # Compute handling time matrix h_ij
        self.h_ij = self._allometric_rate('handling', M_prey, M_pred, T_K) * self.adj_mat

        # Compute metabolic loss rates X_i for all species
        X = self._metabolic_rate(M, T_K)

        # Compute NPP-driven vegetation growth rates for all basal species in this cell.
        G = self.vegetation.growth(B, cell_idx, t)

        # ===== EQUATIONS FOR BASAL SPECIES (PLANTS) =====
        for i in self.basal_idx:
            # Growth from NPP-driven vegetation equation (G_i already computed above)
            growth = G[i]

            # Loss to consumers: sum_j B[j] * F[i,j]
            # F[i,j] uses B[i]^q (resource biomass) in the numerator and includes
            # the interference term; delegate to _functional_response to avoid
            # duplicating that logic incorrectly here.
            loss_to_consumers = sum(B[j] * self._functional_response(B, j)[i]
                                    for j in self.consumer_idx)

            # dB_i/dt = NPP_growth - metabolic_loss - loss_to_consumers
            dydt[i] = growth - X[i] * B[i] - loss_to_consumers
        
        # ===== EQUATIONS FOR CONSUMER SPECIES (ANIMALS) =====
        for j in self.consumer_idx:
            # Compute functional response F_ij for this consumer on all resources
            F_ij = self._functional_response(B, j)
            
            # Get assimilation efficiencies for this consumer j eating each resource i
            # e_plant_spp[j]: consumer j's efficiency for plant (basal) prey
            # e_animal_spp[j]: consumer j's efficiency for animal (consumer) prey
            e_i = np.array([self.e_plant_spp[j] if self.traits.iloc[i]['is_basal']
                           else self.e_animal_spp[j] for i in range(self.n_species)])
            
            # Consumption gain: B_j * sum_i (e_i * F_ij)
            # This is the biomass gained from eating all available resources
            gain = B[j] * np.sum(e_i * F_ij)
            
            # Loss to predators: sum_{j'} B[j'] * F[j, j']
            # F[j, j'] = feeding rate of predator j' on focal species j.
            loss_to_predators = sum(B[jp] * self._functional_response(B, jp)[j]
                                    for jp in self.consumer_idx)
            loss = X[j] * B[j] + loss_to_predators
            
            # dB_j/dt = consumption_gain - metabolic_loss - loss_to_predators
            dydt[j] = gain - loss
        
        # Drive extinct species smoothly to zero (optional extinction dynamics)
        dydt[B < self.ext_threshold] = -B[B < self.ext_threshold] / self.cfg['extinction_timescale']
        
        return dydt
    
    def run_cell(self, B_initial: np.ndarray, cell_idx: int, 
                 t_eval: np.ndarray) -> np.ndarray:
        """
        Run ODE for a single cell using scipy's odeint solver.
        
        Parameters:
            B_initial: initial biomass vector (n_species,) in g/m²
            cell_idx: which cell to simulate
            t_eval: time points to evaluate (days)
        
        Returns:
            B_traj: (len(t_eval), n_species) biomass trajectory
        """
        # Define derivative function (closure over cell_idx)
        def deriv(y, t):
            return self.derivatives(y, t, cell_idx)
        
        # Solve ODE system using 4th-order Runge-Kutta with adaptive timesteps
        # rtol: relative tolerance for error control
        # atol: absolute tolerance for error control
        B_traj = odeint(deriv, B_initial, t_eval, full_output=False, rtol=1e-6, atol=1e-8)
        return B_traj
    
    def run_all_cells(self, B_initial: np.ndarray, t_eval: np.ndarray) -> np.ndarray:
        """
        Run ODE for all cells independently (no dispersal).
        
        Parameters:
            B_initial: (n_cells, n_species) initial biomass matrix
            t_eval: time points to evaluate
        
        Returns:
            B_traj: (len(t_eval), n_cells, n_species) full trajectory
        """
        # Allocate output array for trajectory
        B_traj = np.zeros((len(t_eval), self.n_cells, self.n_species))
        
        # Loop over cells and integrate each independently
        for cell_idx in range(self.n_cells):
            # Print progress indicator
            print(f"  Running cell {cell_idx} / {self.n_cells} ...", end='\r')
            # Integrate this cell's dynamics
            B_traj[:, cell_idx, :] = self.run_cell(B_initial[cell_idx, :], cell_idx, t_eval)
        
        # Print completion message
        print(f"✓ Completed all {self.n_cells} cells.              ")
        return B_traj

    def save_consumer_output(self, B_traj: np.ndarray, t_eval: np.ndarray, output_dir) -> None:
        """
        Save instantaneous dB/dt for all consumer species to atn_model.txt.

        Re-evaluates derivatives() at every recorded time point and cell using the
        saved biomass trajectory, then writes a long-format table.

        Columns: pixel_id, x, y, time, species, delta_biomass
          delta_biomass = dB_j/dt (g/m²/day) — full consumer dynamics
        """
        n_tp = len(t_eval)
        n_consumers = len(self.consumer_idx)

        dBdt_arr = np.zeros((n_tp, self.n_cells, n_consumers))
        for t_idx in range(n_tp):
            for cell_idx in range(self.n_cells):
                dydt = self.derivatives(B_traj[t_idx, cell_idx, :], t_eval[t_idx], cell_idx)
                dBdt_arr[t_idx, cell_idx, :] = dydt[self.consumer_idx]

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
