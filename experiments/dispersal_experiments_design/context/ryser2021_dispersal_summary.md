# Dispersal implementation in Ryser et al. 2021

**Source:** Ryser et al. (2021) *Nature Communications* — "Landscape heterogeneity buffers biodiversity of simulated meta-food-webs under global change through rescue and drainage effects." DOI: 10.1038/s41467-021-24877-0

---

## Core principle

The ATN (allometric trophic network) runs independently on each habitat patch. Spatialization is achieved by adding emigration and immigration terms to the biomass ODE of each animal species. Plants do not disperse.

---

## Biomass dynamics with dispersal

For animals on patch *z*:

$$\frac{dB_{i,z}}{dt} = \underbrace{B_{i,z}\sum_j e_j F_{ij,z} - \sum_j B_{j,z}F_{ji,z} - x_i B_{i,z}}_{\text{local ATN}} \underbrace{- E_{i,z} + I_{i,z}}_{\text{dispersal}}$$

For plants (no dispersal):

$$\frac{dB_{i,z}}{dt} = r_i G_i B_{i,z} - \sum_j B_{j,z}F_{ji,z} - x_i B_{i,z}$$

Food-web topology (who eats whom) and all trophic parameters are identical across all patches.

---

## Emigration — adaptive, density-dependent

Total biomass leaving patch *z* for species *i*:

$$E_{i,z} = d_{i,z} \cdot B_{i,z}$$

The per-capita dispersal rate is a **sigmoid function of local net growth rate**:

$$d_{i,z} = \frac{a}{1 + e^{-b(x_i - \nu_{i,z})}}$$

| Parameter | Value | Role |
|---|---|---|
| $a$ | max dispersal rate (scenario-specific) | ceiling on emigration |
| $b = 10$ | steepness constant | sharpness of the switch |
| $x_i = x_A m_i^{-0.305}$ | metabolic rate of species *i* | inflection point |
| $\nu_{i,z}$ | net growth rate on patch *z* | local condition signal |

The net growth rate is:

$$\nu_{i,z} = \frac{B_{i,z}\sum_j e_j F_{ij,z} - \sum_j B_{j,z}F_{ji,z} - x_i B_{i,z}}{B_{i,z}}$$

**Biological interpretation:** when local conditions are good ($\nu_{i,z} > x_i$), emigration is low. When conditions deteriorate (resource scarcity, heavy predation), $\nu_{i,z}$ drops below $x_i$ and emigration increases sharply. Dispersal is thus an adaptive response to local trophic context.

---

## Immigration — distance-weighted with matrix mortality

Biomass arriving at patch *z* for species *i*:

$$I_{i,z} = \sum_{n \in N_z} E_{i,n} \cdot \underbrace{\max(1-\delta_{i,nz},\, 0)}_{\text{survival}} \cdot \underbrace{\frac{\max(1-\delta_{i,nz},\, 0)}{\sum_{m \in N_n}\max(1-\delta_{i,nm},\, 0)}}_{\text{routing toward } z}$$

where $N_z$ is the set of patches within dispersal range of species *i* relative to patch *z*, and $\delta_{i,nz}$ is the distance between patches *n* and *z* normalised by species *i*'s maximum dispersal distance.

Two mechanisms operate simultaneously:

1. **Matrix mortality:** fraction surviving the journey from *n* to *z* is $(1 - \delta_{i,nz})$. Loss is proportional to distance. If $\delta_{i,nz} \geq 1$, no biomass arrives.
2. **Routing:** emigrating biomass from patch *n* is partitioned across all reachable destinations proportionally to proximity. Closer patches attract a larger share, reflecting that dispersing organisms encounter nearby patches first.

---

## Maximum dispersal distance — body mass scaling

$$\delta_i = \delta_0 \cdot m_i^\epsilon \qquad (\delta_0 = 0.1256,\ \epsilon = 0.05)$$

Larger-bodied species disperse farther. The exponent is small (0.05), so the scaling is weak but consistent — top predators always have access to a broader spatial network than small herbivores.

---

## Landscape structure

Patches are placed as a **random geometric graph**: positions are drawn from a uniform distribution on a unit square $[0,1]^2$. Each patch has its own nutrient supply concentration $S$ (the eutrophication axis). Two patches are connected for species *i* if and only if their distance is within $\delta_i$.

---

## What differs in our implementation

Per AGENTS.md, our spatial ATN differs from Ryser in two ways:

1. **Density-independent dispersal:** emigration rate $d_{i,z}$ is a fixed constant rather than a sigmoid function of net growth rate.
2. **Grid topology with von Neumann neighbourhood:** dispersal only occurs between the four cardinal neighbours of each cell, not across a continuous distance-weighted network.

This simplification removes the adaptive feedback between local trophic state and movement, and restricts spatial connectivity to a regular lattice.
