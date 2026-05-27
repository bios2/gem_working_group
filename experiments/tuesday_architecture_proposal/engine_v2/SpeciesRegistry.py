import numpy as np

class SpeciesRegistry:
    def __init__(self, num_species: int):
        self.num_species = num_species
        self.groups = {"all": np.arange(num_species)}
        
        # Still stores the 1D arrays under the hood for fast matrix math
        self.params = {}
        self.adj_mat = np.zeros((num_species, num_species), dtype=np.int8)
        
        
        
        
        
        
    def add_species_to_group(self, group_name: str, species_indices: list, params: dict = None):
        """Assigns indices to a group and optionally sets their parameters."""
        if group_name not in self.groups:
            self.groups[group_name] = []
        self.groups[group_name].extend(species_indices)
        
        # Automatically unpack the dictionary into our fast arrays
        if params:
            for param_name, value in params.items():
                # If we've never seen this parameter before, create a zeroed array for it
                if param_name not in self.params:
                    self.params[param_name] = np.zeros(self.num_species, dtype=np.float32)
                
                # Assign the parameter value only to the indices in this group
                self.params[param_name][species_indices] = value

    def get_group_indices(self, group_name: str) -> np.ndarray:
        return np.array(self.groups.get(group_name, []), dtype=int)

    def get_group_parameter(self, group_name: str, param_name: str) -> np.ndarray:
        """Returns the parameter array sliced specifically for the requested group."""
        if param_name not in self.params:
            raise KeyError(f"Parameter '{param_name}' does not exist.")
        indices = self.get_group_indices(group_name)
        return self.params[param_name][indices]
    
    def add_feeding_link(self, resource_group: str, consumer_group: str):
        """Sets the adjacency matrix so the consumer group eats the resource group."""
        res_idx = self.get_group_indices(resource_group)
        con_idx = self.get_group_indices(consumer_group)
        
        # NumPy meshgrid logic to set the intersections to 1
        for r in res_idx:
            for c in con_idx:
                self.adj_mat[r, c] = 1
    
    