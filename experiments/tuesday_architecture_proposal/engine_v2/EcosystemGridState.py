import numpy as np
from .SpeciesRegistry import SpeciesRegistry
from contextlib import contextmanager


class EcosystemGridState:
    """
    Holds dynamic 3D biological states [X, Y, Species_ID].
    Modules can dynamically add their own tracking layers.
    """
    def __init__(self, shape: tuple, registry: SpeciesRegistry):
        self.shape = shape
        self.registry = registry
        
        # Dictionary to hold all 3D [X, Y, Species] matrices
        self.layers = {}
        
        # The core currency requested by the ATN specs is always initialized
        self.add_layer("biomass")

    def add_layer(self, layer_name: str, initial_value: float = 0.0):
        """Creates a new 3D tracking matrix for all cells and species."""
        if layer_name in self.layers:
            raise ValueError(f"Layer '{layer_name}' already exists.")
            
        grid_x, grid_y = self.shape
        num_species = self.registry.num_species
        
        self.layers[layer_name] = np.full(
            (grid_x, grid_y, num_species), 
            initial_value, 
            dtype=np.float32
        )

    def get_layer_view(self, layer_name: str, group_name: str) -> np.ndarray:
        """
        Returns a mutable view of a specific layer, restricted to a species group.
        Shape returned: (X, Y, N_species_in_group).
        """
        if layer_name not in self.layers:
            raise KeyError(f"Layer '{layer_name}' has not been registered.")
            
        indices = self.registry.get_group_indices(group_name)
        if len(indices) == 0:
            raise ValueError(f"No species found for group '{group_name}'")
            
        return self.layers[layer_name][:, :, indices]
    
    @contextmanager
    def edit_group_data(self, layer_name: str, group_name: str):
        """Context manager that auto-saves data back to the master grid."""
        idx = self.registry.get_group_indices(group_name)
        
        # 1. Extract the working copy
        working_copy = self.layers[layer_name][:, :, idx].copy()
        
        # 2. Yield it to the user's 'with' block
        yield working_copy
        
        # 3. Auto-save it back when the block finishes
        self.layers[layer_name][:, :, idx] = working_copy