import numpy as np
from .species_registry import SpeciesRegistry
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
        
        # Registry: maps source_layer_name -> list of delta_layer_names.
        # Used by the engine to integrate all deltas into their source layers at step end.
        self.delta_layers = {}
        
        # The core currency requested by the ATN specs is always initialized
        self.add_layer("biomass")
        # Register a default delta for biomass; processes can add more.
        self.add_delta_layer("vegetation_delta", source_layer="biomass")

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
    
    def add_delta_layer(self, delta_name: str, source_layer: str):
        """Register a delta layer that will update a source layer at step end.
        
        Multiple deltas can be registered for the same source layer.
        The engine integrates all registered deltas into the source layer
        and resets them to zero after each step.
        
        Args:
            delta_name: Name of the delta layer (will be created if it doesn't exist).
            source_layer: Name of the source layer this delta updates.
        """
        # Create the delta layer if needed
        if delta_name not in self.layers:
            self.add_layer(delta_name)
        
        # Register in the delta registry
        if source_layer not in self.delta_layers:
            self.delta_layers[source_layer] = []
        
        if delta_name not in self.delta_layers[source_layer]:
            self.delta_layers[source_layer].append(delta_name)

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