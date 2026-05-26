"""
Default configuration and allometric parameters for ATN model.
Based on Binzer et al. (2016) and literature averages.

All parameters are organized by function to aid in understanding the model structure.
"""

# Configuration dictionary: all model parameters in one place
CONFIG = {
    # ===== ALLOMETRIC RATE CONSTANTS =====
    # These control how biological rates scale with body mass
    
    'X0': 0.5,           # Metabolic loss rate normalization constant (per day)
                         # Actual loss: X_i = X0 * M_i^(b_X)
    'b_X': -0.25,        # Metabolic loss mass exponent (~-0.25, same as growth)
                         # Negative exponent: metabolic cost scales sublinearly with mass
    
    'a0': 0.001,         # Attack/clearance rate normalization (per day)
                         # Actual attack rate: a_ij = a0 * M_prey^(b_a_prey) * M_pred^(b_a_pred)
    'b_a_prey': -0.5,    # Attack rate prey mass exponent
                         # Negative: predators are more efficient at eating small prey
    'b_a_pred': 0.5,     # Attack rate predator mass exponent
                         # Positive: larger predators attack faster
    
    'h0': 0.01,          # Handling time normalization (days)
                         # Actual handling time: h_ij = h0 * M_prey^(b_h_prey) * M_pred^(b_h_pred)
    'b_h_prey': 0.5,     # Handling time prey mass exponent
                         # Positive: large prey take longer to handle
    'b_h_pred': -0.5,    # Handling time predator mass exponent
                         # Negative: large predators process prey faster
    
    # ===== FUNCTIONAL RESPONSE PARAMETERS =====
    # These control feeding behavior (Holling Type II / Holling Type III)
    
    'q_hill': 2.0,       # Hill exponent for functional response
                         # q=1: linear (Type I), q=2: Type II (sigmoid), q>2: super-sigmoid (Type III)
                         # q=2 is default for many consumer-resource models
    
    'interference': 0.0, # Intraspecific consumer interference coefficient
                         # 0 = no interference, >0 = cannibalism/competition among predators
    
    'R_opt': 100.0,      # Optimal consumer/resource body mass ratio
                         # Predators are most efficient at consuming prey M_prey * R_opt ≈ M_pred
                         # E.g., R_opt=100 means predator is 100x larger than optimal prey
    
    'gamma': 2.0,        # L-matrix sharpness (body-size matching curve width)
                         # Higher γ: more specific feeding (peaked bell curve)
                         # Lower γ: more generalist feeding (flatter bell curve)
    
    'link_threshold': 0.01,  # Minimum link strength to keep in L-matrix
                             # Links with L_ij < threshold are removed (set to 0)
    
    # ===== ASSIMILATION EFFICIENCIES =====
    # Fraction of consumed biomass that gets incorporated into predator biomass
    
    'e_plant': 0.45,     # Assimilation efficiency when eating plants
                         # Lower than animal prey due to cell walls, lignin, etc.
                         # 45% of plant biomass consumed converts to predator biomass
    
    'e_animal': 0.85,    # Assimilation efficiency when eating animal prey
                         # Higher than plants due to similar biochemistry
                         # 85% of animal biomass consumed converts to predator biomass
    
    # ===== VEGETATION GROWTH =====
    # Parameters for NPP-driven basal species growth (vegetation.md equation)

    'psi': 9.813,            # Carbon to wet matter conversion factor (g wet matter / g C)
                             # From Kattge et al. (2011) via Harfoot et al. (2014) Text S1

    'f_struct_default': 0.3, # Default fractional allocation of NPP to structural tissue
                             # Used when f_struct is absent from traits.txt
                             # ~30% structural allocation (De Kauwe et al. 2014)

    'alpha_herbs_default': 1.0,  # Default half-saturation constant for herb competitive partition
                                 # C_herb = alpha / (alpha + B_trees); units match B_trees (g/m²)
                                 # Used when alpha_herbs is absent from traits.txt


    # ===== TEMPERATURE DEPENDENCE =====
    # These parameters control temperature-dependent rate modifications
    # Following Boltzmann-Arrhenius kinetics: rate ∝ exp(-E/k_B/T)
    
    'use_temperature': True,  # Whether to include temperature effects on rates
                              # If False, all rates are constant (isothermal)
    
    'T0_K': 293.15,      # Reference temperature in Kelvin (20°C)
                         # Rates are standardized to this temperature
    
    'k_B': 8.617e-5,     # Boltzmann constant (eV/K)
                         # Fundamental constant in statistical mechanics
    
    'E_a': 0.65,         # Activation energy (eV)
                         # Higher E_a: rates more sensitive to temperature changes
                         # Typical range: 0.4-1.0 eV for biological processes
    
    # ===== NUMERICAL / SOLVER PARAMETERS =====
    # These control the numerical integration and stopping criteria
    
    'ext_threshold': 1e-6,      # Extinction threshold (g/m²)
                                # Species with B_i < threshold are considered locally extinct
                                # Set to zero in dynamics to prevent numerical errors
    
    'extinction_timescale': 0.1,  # Timescale for extinct species to decay (days)
                                  # If B_i < ext_threshold, dB_i/dt = -B_i / extinction_timescale
                                  # Faster decay prevents oscillations near zero
}
