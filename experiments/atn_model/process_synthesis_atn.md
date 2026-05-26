# Process synthesis: spatial ATN compartment

This document summarizes the information needed to treat the spatialized ATN as a modular compartment that can later be merged with vegetation, movement, and broader spatial ecosystem models. The working unit is one local ATN per spatial cell. At present, cells are independent; future versions can couple them through movement, dispersal, vegetation inputs, or other landscape processes.

The central state variable is local biomass:

```text
B[i, g, t] = biomass density of species i in cell g at time t
```

where biomass is currently expressed in `g/m2`, time in days, body mass in grams, and temperature in Kelvin.

## 1. Inputs

### State variables

These are dynamic variables that the ATN owns or updates during simulation.

| Variable | Meaning | Current source |
|---|---|---|
| `B_i,g(t)` | Biomass density of species `i` in cell `g` | Initialized from `traits.txt`, updated by ODE integration |
| `presence_i,g(t)` | Local presence/absence, currently implicit through biomass and extinction threshold | Derived from `B_i,g(t) > ext_threshold` |

### Spatial inputs

| Variable | Meaning | Current source |
|---|---|---|
| `pixel_id` | Cell identifier | `env_mat.txt` |
| `x`, `y` | Grid coordinates | `env_mat.txt` |
| `temperature_K` | Local cell temperature | `env_mat.txt` |
| `K_plant_i` | Optional local carrying capacity for basal species `i` | `env_mat.txt`; otherwise `K_default` |

### Species and food-web inputs

| Variable | Meaning | Current source |
|---|---|---|
| `species_id` | Species or trophic node identifier | `traits.txt` |
| `body_mass_g` | Species body mass | `traits.txt` |
| `is_basal` | Basal producer flag: `1` for plants/basal nodes, `0` for consumers | `traits.txt` |
| `initial_biomass_g_per_m2` | Initial biomass density | `traits.txt` |
| `adj_mat[i, j]` | Feeding link: consumer `j` eats resource `i` | `adj_mat.txt` |

### Model parameters

| Parameter | Meaning |
|---|---|
| `r0`, `b_r` | Basal growth normalization and body-mass exponent |
| `X0`, `b_X` | Metabolic loss normalization and body-mass exponent |
| `a0`, `b_a_prey`, `b_a_pred` | Attack-rate normalization and prey/predator mass exponents |
| `h0`, `b_h_prey`, `b_h_pred` | Handling-time normalization and prey/predator mass exponents |
| `q_hill` | Functional response Hill exponent |
| `interference` | Consumer interference coefficient |
| `R_opt`, `gamma`, `link_threshold` | Body-size feeding-kernel parameters |
| `e_plant`, `e_animal` | Assimilation efficiencies for plant and animal resources |
| `K_default` | Default plant carrying capacity when cell-specific values are absent |
| `use_temperature`, `T0_K`, `k_B`, `E_a` | Temperature-scaling switch and Boltzmann-Arrhenius constants |
| `ext_threshold`, `extinction_timescale` | Local extinction threshold and decay timescale |

## 2. Outputs

### Primary output

| Output | Meaning | Current file |
|---|---|---|
| `biomass.txt` | Long-format biomass trajectory by cell, time, and species | `atn_output/<timestamp>/biomass.txt` |

Current columns:

```text
pixel_id, x, y, time_step, species_id, biomass
```

### Summary output

| Output | Meaning | Current file |
|---|---|---|
| Species traits used in the run | Body mass, basal/consumer type, initial biomass | `simulation_summary.txt` |
| Model constants | Full configuration used for the run | `simulation_summary.txt` |
| Run dimensions | Number of species, cells, time steps, duration, grid dimensions | `simulation_summary.txt` |
| Persistence diagnostics | Fraction of cells where each species persists above threshold | Printed to console |

### Useful future outputs for module integration

These would make the ATN easier to merge with vegetation and movement modules:

| Output | Meaning |
|---|---|
| `consumption_flux[i, j, g, t]` | Biomass of resource `i` removed by consumer `j` in cell `g` |
| `assimilation_flux[i, j, g, t]` | Biomass gained by consumer `j` from resource `i` after assimilation efficiency |
| `metabolic_loss[i, g, t]` | Biomass lost to metabolism |
| `basal_production[i, g, t]` | Realized plant/basal production before herbivory |
| `local_extinction[i, g, t]` | Indicator for local extinction events |
| `movement_demand_or_pressure[i, g, t]` | Optional signal passed to a future movement module |

## 3. External Data

The ATN needs environmental and biological data that may eventually be provided by other modules rather than static files.

| External data | Role in ATN | Possible future provider |
|---|---|---|
| Temperature | Controls temperature-dependent growth, metabolism, and attack rates | Climate or microclimate module |
| Plant carrying capacity or resource availability | Limits basal biomass and controls bottom-up forcing | Vegetation, productivity, nutrient, or land-cover module |
| Local plant biomass | Could replace or constrain ATN basal biomass | Dynamic vegetation module |
| Species body masses and trophic guilds | Define allometric rates and food-web role | Trait database or regional species pool |
| Feeding links or diet matrix | Defines who can eat whom | Food-web generator, empirical diet database, or body-size kernel |
| Cell geometry and neighbor graph | Needed for dispersal/movement coupling | Spatial grid or landscape module |
| Habitat suitability and barriers | Could mask local presence or modify movement | Land-cover or movement module |

The cleanest future interface is for the ATN to accept a per-cell environment table, a regional species/trait table, and either a fixed food-web matrix or a dynamically generated local interaction matrix.

## 4. Dependencies

### Internal biological dependencies

The ATN couples biomass dynamics through trophic and allometric processes:

| Process | Depends on | Affects |
|---|---|---|
| Basal growth | `B_i,g`, `r_i,g`, `K_i,g`, temperature, optional competition matrix | Basal biomass |
| Metabolism | Body mass, temperature, metabolic parameters | All species biomass |
| Attack rate | Resource mass, consumer mass, temperature, food-web link | Feeding rates |
| Handling time | Resource mass, consumer mass, food-web link | Feeding-rate saturation |
| Functional response | Resource biomass, consumer biomass, attack rate, handling time, Hill exponent, interference | Consumer gains and resource losses |
| Assimilation | Resource type, assimilation efficiency | Consumer biomass gain |
| Predation/herbivory | Local food-web links and biomasses | Resource biomass loss |
| Extinction handling | Local biomass threshold | Local persistence |

### Cross-module dependencies

These are the main seams for a common architecture:

| Dependency | What ATN receives | What ATN returns |
|---|---|---|
| Vegetation/productivity | Plant carrying capacity, plant biomass, basal growth constraints, nutrient limitation | Herbivory pressure, plant consumption fluxes, basal biomass trajectories |
| Movement/dispersal | Immigration/emigration rates, neighbor graph, movement rules, habitat suitability | Local biomass or abundance gradients, extinction pressure, food availability |
| Climate/environment | Temperature and possibly moisture or seasonality | Temperature-sensitive biomass and trophic flux responses |
| Species pool/traits | Body mass, guild, diet constraints, metabolism type, reproduction type | Persistence and biomass outcomes under those traits |

The current model has no explicit movement term:

```text
D_i,g(t) = 0
```

For integration with a movement module, the biomass equation can later become:

```text
dB_i,g/dt = local ATN gains - local ATN losses + D_i,g(t)
```

where `D_i,g(t)` is net immigration minus emigration.

## 5. Expected Behaviors and Ecological Examples

These are target ecological behaviors or scenario classes that the ATN should be able to reproduce or explore.

| Behavior | Brief description |
|---|---|
| Bottom-up control | Higher plant carrying capacity or productivity supports more consumer biomass, up to stability limits. |
| Top-down control | Predators suppress herbivores, indirectly releasing basal species from herbivory. |
| Trophic cascades | Changes at high trophic levels propagate downward through herbivores to plants. |
| Paradox of enrichment | Increasing basal carrying capacity can destabilize consumer-resource dynamics, producing oscillations or extinctions. |
| Body-size filtering | Feeding links and interaction strengths depend on predator/resource body-mass ratios. |
| Temperature sensitivity | Warmer cells can change growth, metabolism, attack rates, and persistence through Boltzmann-Arrhenius scaling. |
| Local extinction and persistence | Species can persist in some cells and disappear in others depending on food availability and environmental context. |
| Spatial rescue effects, future version | Movement among cells can rescue populations that would otherwise go locally extinct. |

## 6. Complexity Roadmap

The goal is to increase complexity in controlled steps while keeping a common interface for the ATN compartment.

### Stage 1: Minimal vertical food web

Start with one simple vertical chain:

```text
plant -> herbivore -> predator
```

Purpose:

- verify biomass accounting;
- test bottom-up and top-down effects;
- reproduce basic consumer-resource oscillations and extinction behavior;
- make the ATN easy to couple with vegetation and movement prototypes.

### Stage 2: Horizontal diversity within guilds

Increase from one species per trophic role to multiple species per guild:

```text
n plants -> n herbivores -> n predators
```

Purpose:

- represent food-web richness;
- test redundancy and competition within trophic levels;
- explore how connectance and diet breadth affect stability;
- support guild-level comparison with vegetation and movement models.

### Stage 3: Metabolic diversity

Add explicit metabolic categories, especially:

```text
ectotherms vs endotherms
```

Purpose:

- allow different metabolic normalizations, exponents, and temperature sensitivities;
- represent vertebrate groups more realistically;
- test how warming affects taxa with different energetic strategies.

### Stage 4: Reproductive strategy diversity

Add life-history categories:

```text
semelparous species = reproduce once
iteroparous species = reproduce multiple times
```

Purpose:

- explore different population recovery and persistence patterns;
- connect biomass dynamics to demographic events;
- allow seasonal or age/life-stage structure if needed later.

### Stage 5: Vertebrate food webs as the core module

The first full ecological target is a vertebrate food-web module. This keeps the initial scope manageable while still allowing:

- body-size structured diets;
- herbivore-carnivore-omnivore structure;
- endotherm/ectotherm contrasts;
- explicit coupling to vegetation as basal resources;
- explicit coupling to movement as spatial redistribution.

### Stage 6: Invertebrate extension

After the vertebrate module is stable, add invertebrates as an additional complexity axis. Possible groups:

| Invertebrate group | Role |
|---|---|
| Leaf-eaters | Additional herbivory pressure on vegetation |
| Detritivores, future option | Link dead organic matter to food webs if detritus is added |
| Animal prey | Prey base for vertebrates and larger invertebrates |
| Pollinators or mutualists, future option | Possible link to plant reproduction if vegetation dynamics require it |

This extension should probably be modular rather than forced into the first vertebrate ATN, because invertebrates may need different body-size ranges, generation times, feeding rules, and reproductive assumptions.

## Interface Summary for Common Architecture

At each coupling step, the ATN compartment can be thought of as:

```text
ATN step(
  cell_environment,
  species_traits,
  local_or_regional_food_web,
  biomass_state,
  optional_movement_flux,
  optional_vegetation_constraints
) -> updated_biomass_state, trophic_fluxes, persistence_status
```

The minimum stable contract is:

- inputs: local environment, species traits, food-web topology, initial/current biomass, model parameters;
- outputs: biomass trajectories, persistence status, and eventually trophic fluxes;
- external coupling: vegetation modifies basal resources, movement modifies spatial redistribution, climate modifies biological rates.
