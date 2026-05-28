"""
Test 1 - Sublinear invasion front under the "b = 0" equivalent.

"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from gem.dispersal import enforce_boundary_conditions
from gem.engine.ecosystem_engine import EcosystemEngine
from gem.engine.ecosystem_grid_state import EcosystemGridState
from gem.engine.environment_state import EnvironmentState
from gem.engine.processes import apply_dispersal
from gem.engine.species_registry import SpeciesRegistry

# --- Parameters --------------------------------------------------------------
GRID_SHAPE = (51, 51)
CENTER = (25, 25)
N_STEPS = 50
THRESHOLD = 1e-6
MAX_DISP_RATE = 0.5
METABOLISM = 0.5
NET_GROWTH = 0.5  # equal to metabolism -> sigmoid degenerates to 0.5 (b = 0 equiv)

# --- Engine setup ------------------------------------------------------------
# Initialize species list (only one species to test the species dispersal)
registry = SpeciesRegistry(num_species=1)
# Initialize maximum dispersal rate
registry.params["max_dispersal_rate"] = np.array([MAX_DISP_RATE], dtype=np.float32)

# Create a grid object. We add layers to it iteratively ; these layers will
# be able to change over time as the engine runs
grid = EcosystemGridState(shape=GRID_SHAPE, registry=registry)
grid.add_layer("net_growth_rate")
grid.add_layer("metabolism_rate")
# This is a special type of layer, that gets updated only by the dispersal process,
# at the end of each time step
grid.add_delta_layer("dispersal_delta", source_layer="biomass")

# Initialize starting values for net growth rate and metabolism
grid.layers["net_growth_rate"][:] = np.float32(NET_GROWTH)
grid.layers["metabolism_rate"][:] = np.float32(METABOLISM)

# Initialize starting biomass distribution 
grid.layers["biomass"][CENTER[0], CENTER[1], 0] = 1.0  # seed central cell

# EnvironmentState objects encode abiotic characteristics. We use the same
# data structure to specify boundary conditions 
env = EnvironmentState(shape=GRID_SHAPE)
# enforce_boundary_conditions returns (n_row, n_col, 1); env.add_layer wants 2D.
env.add_layer(
    "boundary_number",
    enforce_boundary_conditions(*GRID_SHAPE)[..., 0].astype(np.int64), #ensures the boundary conditions matrix has the same dimensions as the grid, and is of correct integer type
)

# EcosystemEngine objects run the simulation. Processes that we want to
# consider in the simulation are added one at the time (here we only care about dispersal)
engine = EcosystemEngine(grid, env)
engine.add_process(apply_dispersal)

# --- Invasion front tracking -------------------------------------------------
# Compute distances from the center for all cells. This will be used downstream
# to track invasion front position over time
rows, cols = np.indices(GRID_SHAPE)
dist_from_center = np.sqrt((rows - CENTER[0]) ** 2 + (cols - CENTER[1]) ** 2)

# Initialize empty array to store results
front = np.zeros(N_STEPS + 1)
# t = 0: biomass is at the central cell only, so front = 0.
for t in range(1, N_STEPS + 1):
    engine.step() # Core function that implements one step of all processes
    biomass_2d = grid.layers["biomass"][:, :, 0]
    above = biomass_2d > THRESHOLD
    #Compute the front as the maximum distance from the distance
    front[t] = float(dist_from_center[above].max()) if above.any() else 0.0

# --- Plot --------------------------------------------------------------------
t_axis = np.arange(N_STEPS + 1)
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_axis, front, "o-", color="steelblue", markersize=4, label="observed front")
# Reference curves rescaled to the final value, for visual comparison.
if front[-1] > 0:
    ax.plot(t_axis, front[-1] * np.sqrt(t_axis / N_STEPS), "--", color="gray",
            label=r"$\sqrt{t}$ reference")
ax.plot(t_axis, front[-1] * t_axis / N_STEPS, ":", color="black",
        label="linear reference")
ax.set_xlabel("Time step")
ax.set_ylabel("Invasion front distance (cells)")
ax.set_title("Test 1 - Sublinear front (b = 0 equivalent)\n"
             f"net_growth = metabolism = {METABOLISM}, max_disp_rate = {MAX_DISP_RATE}")
ax.grid(alpha=0.3)
ax.legend()

out_path = Path(__file__).parent / "test_1.png"
plt.tight_layout()
plt.savefig(out_path, dpi=120)
plt.close(fig)

# --- Verdict -----------------------------------------------------------------
# Compare sum of squared residuals for a linear fit vs a sqrt(t) fit. The smaller
# SSR is the better-matching shape; sublinear means sqrt(t) wins.
linear_fit = np.poly1d(np.polyfit(t_axis, front, 1))(t_axis)
sqrt_fit = np.poly1d(np.polyfit(np.sqrt(t_axis), front, 1))(np.sqrt(t_axis))
ssr_lin = float(((front - linear_fit) ** 2).sum())
ssr_sqrt = float(((front - sqrt_fit) ** 2).sum())
verdict = "PASS" if ssr_sqrt < ssr_lin else "FAIL"

print(f"Test 1 -> {out_path.name}")
print(f"  SSR linear fit:  {ssr_lin:.3f}")
print(f"  SSR sqrt(t) fit: {ssr_sqrt:.3f}")
print(f"  Hypothesis (sublinear): {verdict}")
