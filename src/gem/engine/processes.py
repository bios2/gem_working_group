from .environment_state import EnvironmentState
import numpy as np

from .ecosystem_grid_state import EcosystemGridState
from ..vegetation import logistic_growth_delta
from ..dispersal import compute_disperse_delta

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
    """Adapter: density-dependent dispersal (Ryser et al. 2021).

    Reads biomass + per-species rate layers from grid (ATN populates the rate
    layers), reads boundary_number from env, calls compute_disperse_delta,
    writes the dispersal delta into the registered dispersal_delta layer.
    """
    biomass    = grid.get_layer_view("biomass", "all")
    net_growth = grid.get_layer_view("net_growth_rate", "all")
    metabolism = grid.get_layer_view("metabolism_rate", "all")

    # Species-level parameter (shape (S,)) -> broadcast to (X, Y, S) so the
    # science function's shape assert passes.
    max_disp_rate = grid.registry.get_group_parameter("all", "max_dispersal_rate")
    max_disp_rate = np.broadcast_to(max_disp_rate, biomass.shape)

    # boundary_number is a 2D geometric field on env; add a trailing axis so it
    # broadcasts over species inside compute_disperse_delta.
    boundary_number = env.get_layer("boundary_number")[..., np.newaxis].astype(np.int64)

    delta = compute_disperse_delta(
        biomass=biomass,
        net_growth=net_growth,
        metabolism=metabolism,
        max_disp_rate=max_disp_rate,
        boundary_number=boundary_number,
        b=10.0,
        dt=1.0,
    )

    with grid.edit_group_data("dispersal_delta", "all") as d:
        d += delta


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
