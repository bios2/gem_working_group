"""
Spatially explicit unscaled ATN model (following Binzer et al. 2016).
Implementation of Section 8 from ATN_model_spatiotemporal_formulas_parameters.Rmd

Each consumer species eats from resources following Holling Type II functional response.
Basal species (plants) grow logistically and are consumed by animals.
All rates scale allometrically with body mass and optionally with temperature.
"""
# Import required libraries
import numpy as np  # numerical arrays and mathematics
from scipy.integrate import odeint  # ODE solver (4th-order Runge-Kutta adaptive)
import pandas as pd  # tabular data handling
from typing import Tuple, Dict  # type hints for function signatures

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
        
        # Extract allometric rate constants from config
        self.r0 = config['r0']  # basal growth rate normalization constant
        self.br = config['b_r']  # basal growth mass exponent (typically -0.25)
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
        self.Ropt = config['R_opt']  # optimal consumer/resource body mass ratio for feeding
        self.gamma = config['gamma']  # L-matrix sharpness (how peaked is the feeding kernel)
        self.link_threshold = config['link_threshold']  # minimum link strength to keep
        
        # Temperature parameters (for Boltzmann-Arrhenius model)
        self.T0_K = config['T0_K']  # reference temperature in Kelvin (20°C = 293.15 K)
        self.k_B = config['k_B']  # Boltzmann constant in eV/K
        self.E_a = config['E_a']  # activation energy in eV
        
        # Extinction parameters
        self.ext_threshold = config['ext_threshold']  # biomass threshold for local extinction
        
        # Print initialization message
        print(f"ATN Model initialized: {self.n_species} species, {self.n_cells} cells")
    
    def _L_matrix(self, cell_idx: int) -> np.ndarray:
        """
        Compute the body-size matching matrix L_{ij}.
        
        The L-matrix represents the feeding efficiency based on body size ratio.
        L_{ij} = [z * exp(1 - z)]^gamma, where z = M_j / (M_i * R_opt)
        
        This is bell-shaped: peaks at z=1 (optimal ratio), decays on both sides.
        """
        # Extract body masses for all species
        M = self.traits['body_mass_g'].values
        
        # Create predator body mass column vector (for broadcasting)
        M_pred = M[:, np.newaxis]  # shape (n_species, 1)
        # Create prey body mass row vector (for broadcasting)
        M_prey = M[np.newaxis, :]   # shape (1, n_species)
        
        # Compute body size ratio z = M_predator / (M_prey * R_opt)
        z = M_pred / (M_prey * self.Ropt)
        # Compute L-matrix: bell-shaped curve with sharpness controlled by gamma
        L = np.power(z * np.exp(1.0 - z), self.gamma)
        
        # Remove weak links (below threshold)
        L[L < self.link_threshold] = 0
        # Basal species cannot be consumers, so remove all outgoing links from basal species
        L[self.basal_idx, :] = 0
        
        # Return L-matrix masked by adjacency (only keep links where adj_mat = 1)
        return L * self.adj_mat
    
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
            raise ValueError(f\"Unknown rate type: {rate_type}\")
        
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
        """Metabolic loss rate X_i = X0 * M^(-0.25) * temp_factor"""
        # Compute metabolic rate: X_i = X0 * M^b_X
        X = self.X0 * np.power(M, self.bX)
        # Apply temperature dependence if enabled
        if self.cfg['use_temperature']:
            temp_factor = np.exp(-self.E_a * (self.T0_K - T_K) / 
                                (self.k_B * T_K * self.T0_K))
            X *= temp_factor
        return X
    
    def _basal_growth_rate(self, M: np.ndarray, T_K: float) -> np.ndarray:
        """Maximum basal (plant) growth rate r_i = r0 * M^(-0.25) * temp_factor"""
        # Compute basal growth rate: r_i = r0 * M^b_r
        r = self.r0 * np.power(M, self.br)
        # Apply temperature dependence if enabled
        if self.cfg['use_temperature']:
            temp_factor = np.exp(-self.E_a * (self.T0_K - T_K) / 
                                (self.k_B * T_K * self.T0_K))
            r *= temp_factor
        return r
    
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
        M_prey = M[np.newaxis, :]  # row vector for broadcasting
        M_pred = M[:, np.newaxis]  # column vector for broadcasting
        
        # Compute attack rate matrix a_ij
        self.a_ij = self._allometric_rate('attack', M_prey, M_pred, T_K) * self.adj_mat
        # Compute handling time matrix h_ij
        self.h_ij = self._allometric_rate('handling', M_prey, M_pred, T_K) * self.adj_mat
        
        # Compute metabolic loss rates X_i for all species
        X = self._metabolic_rate(M, T_K)
        # Compute basal growth rates r_i for all species
        r = self._basal_growth_rate(M, T_K)
        
        # ===== EQUATIONS FOR BASAL SPECIES (PLANTS) =====
        for i in self.basal_idx:
            # Get carrying capacity for this species in this cell
            K_i = self.env.loc[cell_idx, f'K_plant_{i}'] if f'K_plant_{i}' in self.env.columns else self.cfg['K_default']
            
            # Logistic growth term: r_i * B_i * (1 - B_i/K_i)
            # This reduces growth when biomass approaches carrying capacity
            growth = r[i] * B[i] * (1.0 - B[i] / K_i) if K_i > 0 else 0
            
            # Loss to consumers: sum over all consumers j eating species i
            # This requires computing F_ij (feeding rate) for this species' consumers
            loss_to_consumers = np.sum(B[self.consumer_idx] * self.a_ij[i, self.consumer_idx] * 
                                      np.power(B[self.consumer_idx], self.q_hill) / 
                                      (1.0 + np.sum(self.h_ij[:, self.consumer_idx] * 
                                                   self.a_ij[:, self.consumer_idx] * 
                                                   np.power(B[:, np.newaxis], self.q_hill), axis=0)))
            
            # dB_i/dt = growth - metabolic_loss - loss_to_consumers
            dydt[i] = growth - X[i] * B[i] - loss_to_consumers
        
        # ===== EQUATIONS FOR CONSUMER SPECIES (ANIMALS) =====
        for j in self.consumer_idx:
            # Compute functional response F_ij for this consumer on all resources
            F_ij = self._functional_response(B, j)
            
            # Get assimilation efficiencies e_i for all resources
            # Plants have lower efficiency than animal prey
            e_i = np.array([self.cfg['e_plant'] if self.traits.iloc[i]['is_basal'] 
                           else self.cfg['e_animal'] for i in range(self.n_species)])
            
            # Consumption gain: B_j * sum_i (e_i * F_ij)
            # This is the biomass gained from eating all available resources
            gain = B[j] * np.sum(e_i * F_ij)
            
            # Loss to consumers: sum over all j' that eat consumer j
            # This is more complex as it involves predation on this consumer
            loss = X[j] * B[j] + np.sum(B[self.consumer_idx] * F_ij[j] * 
                                       (1.0 - np.eye(len(self.consumer_idx))[np.where(self.consumer_idx == j)[0][0]] 
                                        if j in self.consumer_idx else 0))
            
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
            print(f\"  Running cell {cell_idx} / {self.n_cells} ...\", end='\r')
            # Integrate this cell's dynamics
            B_traj[:, cell_idx, :] = self.run_cell(B_initial[cell_idx, :], cell_idx, t_eval)
        
        # Print completion message
        print(f\"✓ Completed all {self.n_cells} cells.              \")
        return B_traj
