import numpy as np

class EnvironmentState:
    """
    Represents the physical world grid.
    Passed through our modular functions to provide environmental context.
    """
    def __init__(self, shape: tuple):
        # shape could be (latitude_bins, longitude_bins) e.g., (180, 360) for 1-degree resolution
        self.shape = shape 
        
        # Dictionary to hold our vectorized spatial data (e.g., temperature, land/water mask)
        self.layers = {}

    def add_layer(self, name: str, data: np.ndarray):
        """Adds a new environmental variable layer."""
        if data.shape != self.shape:
            raise ValueError(f"Layer shape {data.shape} does not match world shape {self.shape}")
        self.layers[name] = data
        
    def get_layer(self, name: str) -> np.ndarray:
        return self.layers[name]