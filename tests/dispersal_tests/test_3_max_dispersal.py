"""
Test 3 - Maximum dispersal rate contrast.

Hypothesis (tests.md): same net_growth and same metabolism, a species with high
max_dispersal_rate (a = 0.8) should reach farther than a species with low
max_dispersal_rate (a = 0.2).

Mechanism: per-capita dispersal rate is
    d = max_disp_rate / (1 + exp(-b * (metabolism - net_growth)))
With net_growth = 0, metabolism = 0.5, b = 10:
    sigmoid factor = 1 / (1 + exp(-5)) ~= 0.993
So d ~= 0.79 (high) vs 0.20 (low) - per-capita rate scales directly with
max_dispersal_rate, and the high-a species advances further per step.

Both species are run inside the same engine; their fronts are tracked side by
side.

Output: test_3.png saved next to this script.
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
MAX_DISP_RATE_HIGH = 0.8
MAX_DISP_RATE_LOW = 0.2
METABOLISM = 0.5
NET_GROWTH = 0.0

# --- Engine setup ------------------------------------------------------------
registry = SpeciesRegistry(num_species=2)
registry.params["max_dispersal_rate"] = np.array(
    [MAX_DISP_RATE_HIGH, MAX_DISP_RATE_LOW], dtype=np.float32
)

grid = EcosystemGridState(shape=GRID_SHAPE, registry=registry)
grid.add_layer("net_growth_rate")
grid.add_layer("metabolism_rate")
grid.add_delta_layer("dispersal_delta", source_layer="biomass")

grid.layers["net_growth_rate"][:] = np.float32(NET_GROWTH)
grid.layers["metabolism_rate"][:] = np.float32(METABOLISM)  # same for both species
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
ax.plot(t_axis, front_high, "o-", color="steelblue", markersize=4,
        label=f"High max dispersal (a = {MAX_DISP_RATE_HIGH})")
ax.plot(t_axis, front_low, "x-", color="tomato", markersize=4,
        label=f"Low max dispersal (a = {MAX_DISP_RATE_LOW})")
ax.set_xlabel("Time step")
ax.set_ylabel("Invasion front distance (cells)")
ax.set_title("Test 3 - Maximum dispersal rate contrast\n"
             f"net_growth = {NET_GROWTH}, metabolism = {METABOLISM}, b = 10")
ax.legend()
ax.grid(alpha=0.3)

out_path = Path(__file__).parent / "test_3.png"
plt.tight_layout()
plt.savefig(out_path, dpi=120)
plt.close(fig)

# --- Verdict -----------------------------------------------------------------
passes = bool((front_high >= front_low).all()) and bool(front_high[-1] > front_low[-1])
verdict = "PASS" if passes else "FAIL"

print(f"Test 3 -> {out_path.name}")
print(f"  Final front high max_disp_rate: {front_high[-1]:.2f} cells")
print(f"  Final front low  max_disp_rate: {front_low[-1]:.2f} cells")
print(f"  Hypothesis (high >= low at every step, strict at end): {verdict}")
