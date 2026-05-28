"""
Test 2 - Metabolism contrast.

Hypothesis (tests.md): same net_growth and same max_dispersal_rate, a species
with high metabolism (1.0) should disperse more than a species with low
metabolism (0.1).

Mechanism: per-capita dispersal rate is
    d = max_disp_rate / (1 + exp(-b * (metabolism - net_growth)))
With net_growth = 0 and b = 10:
    high-met (1.0): d ~= 0.500 (sigmoid saturated)
    low-met  (0.1): d ~= 0.366
So the high-metabolism species advances further per step.

Both species are run inside the same engine; their fronts are tracked side by
side.

Output: test_2.png saved next to this script.
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
METABOLISM_HIGH = 1.0
METABOLISM_LOW = 0.1
NET_GROWTH = 0.0

# --- Engine setup ------------------------------------------------------------
registry = SpeciesRegistry(num_species=2)
registry.params["max_dispersal_rate"] = np.array(
    [MAX_DISP_RATE, MAX_DISP_RATE], dtype=np.float32
)

grid = EcosystemGridState(shape=GRID_SHAPE, registry=registry)
grid.add_layer("net_growth_rate")
grid.add_layer("metabolism_rate")
grid.add_delta_layer("dispersal_delta", source_layer="biomass")

grid.layers["net_growth_rate"][:] = np.float32(NET_GROWTH)
# Per-species metabolism: species 0 high, species 1 low. (S,) broadcasts to (X, Y, S).
grid.layers["metabolism_rate"][:] = np.array(
    [METABOLISM_HIGH, METABOLISM_LOW], dtype=np.float32
)
grid.layers["biomass"][CENTER[0], CENTER[1], :] = 1.0  # seed both species at center

env = EnvironmentState(shape=GRID_SHAPE)
env.add_layer(
    "boundary_number",
    enforce_boundary_conditions(*GRID_SHAPE)[..., 0].astype(np.int64),
)

engine = EcosystemEngine(grid, env)
engine.add_process(apply_dispersal)

# --- Invasion front tracking -------------------------------------------------
rows, cols = np.indices(GRID_SHAPE)
dist_from_center = np.sqrt((rows - CENTER[0]) ** 2 + (cols - CENTER[1]) ** 2)

front_high = np.zeros(N_STEPS + 1)
front_low = np.zeros(N_STEPS + 1)
for t in range(1, N_STEPS + 1):
    engine.step()
    above_high = grid.layers["biomass"][:, :, 0] > THRESHOLD
    above_low = grid.layers["biomass"][:, :, 1] > THRESHOLD
    front_high[t] = float(dist_from_center[above_high].max()) if above_high.any() else 0.0
    front_low[t] = float(dist_from_center[above_low].max()) if above_low.any() else 0.0

# --- Plot --------------------------------------------------------------------
t_axis = np.arange(N_STEPS + 1)
fig, ax = plt.subplots(figsize=(8, 5))
ax.plot(t_axis, front_high, "o-", color="tomato", markersize=4,
        label=f"High metabolism ({METABOLISM_HIGH})")
ax.plot(t_axis, front_low, "x-", color="steelblue", markersize=4,
        label=f"Low metabolism ({METABOLISM_LOW})")
ax.set_xlabel("Time step")
ax.set_ylabel("Invasion front distance (cells)")
ax.set_title("Test 2 - Metabolism contrast\n"
             f"net_growth = {NET_GROWTH}, max_disp_rate = {MAX_DISP_RATE}, b = 10")
ax.legend()
ax.grid(alpha=0.3)

out_path = Path(__file__).parent / "test_2.png"
plt.tight_layout()
plt.savefig(out_path, dpi=120)
plt.close(fig)

# --- Verdict -----------------------------------------------------------------
passes = bool((front_high >= front_low).all()) and bool(front_high[-1] > front_low[-1])
verdict = "PASS" if passes else "FAIL"

print(f"Test 2 -> {out_path.name}")
print(f"  Final front high metabolism: {front_high[-1]:.2f} cells")
print(f"  Final front low  metabolism: {front_low[-1]:.2f} cells")
print(f"  Hypothesis (high >= low at every step, strict at end): {verdict}")
