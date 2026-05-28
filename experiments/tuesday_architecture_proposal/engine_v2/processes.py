from .EnvironmentState import EnvironmentState
import numpy as np

from .EcosystemGridState import EcosystemGridState

# --- Mock Processes (What your colleagues will write) ---

def apply_vegetation_growth(grid: EcosystemGridState, env: EnvironmentState):
    """Vegetation updates the basal biomass layer."""
    
    r = grid.registry.get_group_parameter("plants", "base_growth_rate")
    #broadcasting the carrying capacity across the species dimension for plants. 
    #inital shape of cc is (X, Y) and becomes (X, Y, 1) after adding the new axis, allowing it to broadcast correctly with the plants layer which has shape (X, Y, N_plants).
    cc = env.get_layer("carrying_capacity")[..., np.newaxis]
    
    # The 'with' block handles the extraction and the write-back automatically!
    with grid.edit_group_data("biomass", "plants") as plants:
        # We do the math in-place on the yielded object
        plants += (r * plants * (1 - (plants / cc)))
        
        

def apply_atn_step(grid: EcosystemGridState, env: EnvironmentState):
    """ATN handles feeding, metabolism, and sets up growth rates for dispersal."""
    # The ATN process modifies biomass based on trophic interactions
    # and populates intermediate matrices like grid.net_growth_rate
    pass

def apply_dispersal(grid: EcosystemGridState, env: EnvironmentState):
    """Dispersal moves biomass between adjacent grid cells."""
    # Dispersal reads grid.biomass and grid.net_growth_rate, 
    # calculates diffusion, and shifts biomass left/right/up/down.
    pass