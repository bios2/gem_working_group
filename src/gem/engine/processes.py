from gem.metabolism import calculate_massspecific_metabolism, calculate_biomass_metabolism

from .environment_state import EnvironmentState
import numpy as np

from .ecosystem_grid_state import EcosystemGridState
from ..vegetation import logistic_growth_delta

# ============================================================================
# PROCESS ARCHITECTURE: Pure Science + Thin Adapter
# ============================================================================
#
# Every process has two parts:
#
# 1. PURE SCIENCE FUNCTION (e.g., src/gem/vegetation.py)
#    - Input: explicit numpy arrays (biomass, parameters, environment data)
#    - Output: delta array (change to be integrated)
#    - No grid, no engine, no side effects
#    - Works on any shape: (S,), (Y, S), or (X, Y, S) through broadcasting
#    - Easy to test with hand-built data
#    - Easy to reason about for ecology contributors
#    - Example signature:
#
#      def logistic_growth_delta(biomass, growth_rate, carrying_capacity, dt):
#          """Logistic growth: dB/dt = r*B*(1 - B/K)"""
#          assert biomass.shape == growth_rate.shape == carrying_capacity.shape
#          return dt * growth_rate * biomass * (1.0 - biomass / carrying_capacity)
#
# 2. ADAPTER FUNCTION (e.g., apply_vegetation_growth below)
#    - Input: grid and env (the engine's shared state)
#    - Pulls the slices and parameters the science function needs
#    - Calls the pure science function
#    - Writes the delta back to the correct delta layer
#    - Example pattern:
#
#      def apply_vegetation_growth(grid, env):
#          r = grid.registry.get_group_parameter("plants", "base_growth_rate")
#          K = env.get_layer("carrying_capacity")[..., np.newaxis]
#          B = grid.get_layer_view("biomass", "plants")
#          
#          delta = logistic_growth_delta(B, r, K, dt=1.0)
#          
#          with grid.edit_group_data("vegetation_delta", "plants") as d:
#              d += delta
#
# WHY THIS PATTERN?
# - The science function signature documents what the process depends on
#   (no hidden dependencies buried in grid/env reads)
# - Non-programmer ecologists can read, prototype, and test the science in a notebook
#   without touching engine code
# - Adapters are boilerplate; the biology is in the science function
# - Delta layers make composition order-independent: every process sees frozen biomass_t
#
# ============================================================================


# --- Mock Processes (What your colleagues will write) ---

def apply_vegetation_growth(grid: EcosystemGridState, env: EnvironmentState):
    """Adapter: vegetation growth process.
    
    Pulls data from grid/env, calls the pure science function, writes delta back.
    The actual biology is in src.gem.vegetation.logistic_growth_delta().
    """
    # Extract parameters for the plants group
    growth_rate = grid.registry.get_group_parameter("plants", "base_growth_rate")
    carrying_capacity = env.get_layer("carrying_capacity")[..., np.newaxis]
    biomass = grid.get_layer_view("biomass", "plants")
    
    # Broadcast growth_rate (scalar) and carrying_capacity to match biomass shape
    growth_rate_broadcast = np.full_like(biomass, growth_rate)
    
    # Call pure science function
    delta = logistic_growth_delta(
        biomass=biomass,
        growth_rate=growth_rate_broadcast,
        carrying_capacity=carrying_capacity,
        dt=1.0,  # one time step
    )
    
    # Write delta back to the grid
    with grid.edit_group_data("vegetation_delta", "plants") as plants_delta:
        plants_delta += delta
     
    
def apply_metabolism_mass_specific(grid: EcosystemGridState, env: EnvironmentState):
    """
    Dependency adapter: Calculates W/g and stores it in a shared layer.
    Reused by ATN (loss) and Dispersal (trigger).
    """
    group = "vertebrates"

    # 1. Fetch per-species traits (1D arrays: (S_group,))
    # Note: mass_g is individual mass, not population biomass.
    m_indiv = grid.registry.get_group_parameter(group, "mass_g")
    c_int = grid.registry.get_group_parameter(group, "c_int")
    b = grid.registry.get_group_parameter(group, "b")
    
    # np.nan for all t_reg values
    t_reg = np.full(m_indiv.shape, np.nan)

    # 2. Fetch ambient temperature (2D array: (Y, X))
    ambient_temp = env.get_layer("temperature")

    # 3. Calculate mass-specific metabolic rate (W/g).
    # Traits are (S,) and ambient_temp is given a trailing axis to be (Y, X, 1).
    # NumPy automatically broadcasts these to (Y, X, S) inside the science function.
    ms_rate = calculate_massspecific_metabolism(
        mass_g=m_indiv,
        body_temp_C=t_reg,
        c_int=c_int,
        b=b,
        ambient_temp_C=ambient_temp[..., np.newaxis],
    )

    # We use [:] to perform an in-place update on the yielded view.
    with grid.edit_group_data("metabolism_mass_specific", group) as rate_view:
        rate_view[:] = ms_rate




def apply_metabolism_biomass(grid: EcosystemGridState, env: EnvironmentState): 
    """
    Biomass-modifying adapter: Applies the rate to population biomass over dt.
    """
    group = "vertebrates"
    biomass = grid.get_layer_view("biomass", group)
    ms_rate = grid.get_layer_view("metabolic_rate", group)

    delta = calculate_biomass_metabolism(
        initial_biomass_g=biomass,
        mass_specific_metabolic_rate=ms_rate,
        dt=1.0 # 1 day
    )

    with grid.edit_group_data("vegetation_delta", group) as d:
        d += delta
        

def apply_atn_step(grid: EcosystemGridState, env: EnvironmentState):
    """ATN handles feeding, metabolism, and sets up growth rates for dispersal.
    
    Reads biomass and feeding links, writes trophic dynamics delta to biomass_delta.
    """
    # The ATN process reads biomass to determine trophic interactions,
    # computes feeding and metabolism fluxes,
    # and accumulates the net trophic delta into biomass_delta.
    # Also populates intermediate matrices like grid.net_growth_rate for other processes.
    pass

def apply_dispersal(grid: EcosystemGridState, env: EnvironmentState):
    """Dispersal moves biomass between adjacent grid cells.
    
    Reads current biomass and net growth rates, writes dispersal delta to biomass_delta.
    """
    # Dispersal reads grid.biomass (current state at start of step)
    # and grid.net_growth_rate (computed by ATN),
    # calculates diffusion and movement fluxes,
    # and accumulates the dispersal delta into biomass_delta.
    pass


# ============================================================================
# CONTRIBUTING A NEW PROCESS
# ============================================================================
#
# Step 1: Write the pure science function in a new module (e.g., src/gem/mortality.py):
#
#   from numpy.typing import NDArray
#   import numpy as np
#
#   def mortality_delta(
#       biomass: NDArray[np.float64],
#       mortality_rate: NDArray[np.float64],
#       dt: float,
#   ) -> NDArray[np.float64]:
#       """Loss due to background mortality.
#       
#       Args:
#           biomass: Current biomass. Same shape for all inputs.
#           mortality_rate: Per-capita mortality rate. Same shape as biomass.
#           dt: Time step.
#       
#       Returns:
#           Biomass delta (negative, since mortality removes biomass).
#       """
#       assert biomass.shape == mortality_rate.shape
#       return -dt * mortality_rate * biomass
#
# Step 2: Write tests in tests/test_mortality.py:
#   - Use hand-built arrays, no engine
#   - Test edge cases: zero mortality, all dead, etc.
#   - Verify shape invariance: (S,), (Y, S), (X, Y, S)
#
# Step 3: Write the adapter here (e.g., apply_mortality below):
#   - Register your delta layer (if new)
#   - Extract parameters from grid/env
#   - Call the science function
#   - Write delta back
#
# Step 4: Register the adapter with the engine:
#   model.add_process(apply_mortality)
#
# ============================================================================

# --- Placeholder adapters for ATN and dispersal ---
# (To be implemented by the respective teams; same pattern as apply_vegetation_growth)
