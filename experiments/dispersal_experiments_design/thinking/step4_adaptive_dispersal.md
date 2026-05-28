# Step 4 — Adaptive (density-dependent) dispersal

## Goal

Upgrade the density-independent diffusion function to include the dispersal
trigger from Ryser et al. (2021): per-capita emigration rate responds to local
trophic conditions rather than being a fixed constant.

---

## Key idea from Ryser et al. (2021)

In the density-independent version, every cell loses the same fraction of its
biomass per step regardless of how well the population is doing locally. In
Ryser, individuals are more likely to leave when local conditions are bad
(low food, high predation pressure) and less likely to leave when conditions
are good. This is implemented as a sigmoid function of the **net growth rate**
$\nu_{i,z}$:

$$d_{i,z} = \frac{a}{1 + e^{-b(x_i - \nu_{i,z})}}$$

- $a$ — maximum dispersal rate (ceiling)
- $b = 10$ — steepness; controls how sharply dispersal switches around the inflection point
- $x_i$ — metabolic rate of the species; acts as the inflection point
- $\nu_{i,z}$ — net growth rate at cell $z$: feeding gains minus predation losses minus metabolism, per unit biomass

When $\nu_{i,z} > x_i$ (the population is growing faster than it costs to
maintain itself), the sigmoid argument is negative and $d_{i,z}$ is close to
zero — individuals stay put. When $\nu_{i,z} < x_i$ (the population is energy
limited or over-exploited), the argument is positive and $d_{i,z}$ rises
toward $a$ — individuals disperse.

---

## Design decisions

### What changes vs. the density-independent function

| Component | Density-independent | Adaptive (this step) |
|---|---|---|
| Per-capita dispersal rate | fixed scalar `disp_rate` | cell-specific matrix `d_grid` from sigmoid |
| Extra inputs needed | none | `net_growth_grid`, `metabolic_rate`, `max_disp_rate`, `b` |
| Survival during dispersal | 1 (no matrix mortality) | 1 (no matrix mortality, same) |
| Routing to neighbours | equal split via `/4` | equal split via `/4` (unchanged) |
| Boundary handling | `nbr_count` matrix | `nbr_count` matrix (unchanged) |

The only structural change is replacing the scalar `disp_rate / 4` with a
matrix `d_grid / 4`, where each cell has its own rate.

### Why keep survival = 1

Setting the survival factor to 1 (no matrix mortality during dispersal)
decouples the effect of the dispersal trigger from the effect of landscape
hostility. This lets us first validate that the sigmoid trigger produces the
expected spatial patterns before adding a second source of complexity.

### Why create a new file rather than overwrite

AGENTS.md step 6 requires a direct comparison between density-independent and
adaptive dispersal. Both functions must coexist so they can be called
side-by-side in a script.

---

## Net growth rate — what the caller must supply

The function takes `net_growth_grid` as an input rather than computing it
internally. This is intentional: the net growth rate depends on the full ATN
(feeding rates, predation, metabolism), which lives outside the dispersal
function. The dispersal function should stay focused on movement only.

For a single isolated species (no predators, no prey), the net growth rate
simplifies to just the intrinsic growth minus metabolism. For a full food web
the caller computes it from the ATN state before calling this function.

---

## Expected behaviour

- Cells with low biomass or poor trophic conditions → high $d_{i,z}$ → strong
  net emigration (drainage-like effect on struggling patches).
- Cells with high biomass and good conditions → low $d_{i,z}$ → population
  stays put and accumulates.
- At equilibrium, dispersal flows should be smaller overall than in the
  density-independent case when conditions are uniformly good, but larger from
  patches under stress.
