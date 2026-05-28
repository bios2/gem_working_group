from .environment_state import EnvironmentState
from .ecosystem_grid_state import EcosystemGridState

class EcosystemEngine:
    def __init__(self, grid_state: EcosystemGridState, env_state: EnvironmentState):
        self.grid = grid_state
        self.env = env_state
        self.processes = []
        
    def add_process(self, process_func):
        """Registers a process to the simulation pipeline."""
        self.processes.append(process_func)
        
    def step(self):
        """
        Executes one time step:
        1. All processes accumulate their contributions into registered delta layers.
        2. Integrate all deltas: for each (source_layer, delta_layers) pair,
           source_layer += Σ(delta_layers).
        3. Zero all delta layers for the next step.
        
        This makes computation order-independent: all processes see the same state_t,
        and multiple processes can contribute to the same source layer through separate deltas.
        """
        # Run all processes; each accumulates into registered delta layers
        for process in self.processes:
            process(self.grid, self.env)
        
        # Integrate all registered deltas into their source layers
        for source_layer, delta_layer_list in self.grid.delta_layers.items():
            for delta_layer in delta_layer_list:
                self.grid.layers[source_layer] += self.grid.layers[delta_layer]
                # Reset delta for next step
                self.grid.layers[delta_layer][:] = 0.0