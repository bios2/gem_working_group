from .EnvironmentState import EnvironmentState
from .EcosystemGridState import EcosystemGridState
class EcosystemEngine:
    def __init__(self, grid_state: EcosystemGridState, env_state: EnvironmentState):
        self.grid = grid_state
        self.env = env_state
        self.processes = []
        
    def add_process(self, process_func):
        """Registers a process to the simulation pipeline."""
        self.processes.append(process_func)
        
    def step(self):
        """Executes one time step by running all processes in order."""
        for process in self.processes:
            process(self.grid, self.env)