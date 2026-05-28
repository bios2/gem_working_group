import numpy as np

class EnvironmentState:
    """
    Represents the physical world grid and its environmental layers (e.g., Temperature, NPP).
    
    Following the project's geographic guidelines, this class supports equal-area 
    projections (e.g., Albers North America) rather than assuming lat-lon. 
    This ensures that one cell equals one comparable unit of surface area.
    """
    def __init__(
        self, 
        shape: tuple, 
        crs: str = "ESRI:102008",
        proj_string: str = "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs",
        cell_size: float = 100000.0,
        origin: tuple = (-7000000.0, -2000000.0)
    ):
        """
        Initializes the environment state with metadata for an equal-area grid.
        """
        self.shape = shape  # (grid_x, grid_y)
        self.crs = crs
        self.proj_string = proj_string
        self.cell_size = cell_size
        self.origin = origin

        # Dictionary to hold our vectorized spatial data (e.g., temperature, land/water mask)
        self.layers = {}

    def add_layer(self, name: str, data: np.ndarray):
        """Adds a new environmental variable layer."""
        if data.shape != self.shape:
            raise ValueError(f"Layer shape {data.shape} does not match world shape {self.shape}")
        self.layers[name] = data
        
    def get_layer(self, name: str) -> np.ndarray:
        return self.layers[name]

    @property
    def cell_area(self) -> float:
        """Returns the surface area of a single cell in square meters."""
        return float(self.cell_size ** 2)