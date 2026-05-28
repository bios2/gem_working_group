---
title: "Spatially explicit Allometric Trophic Network (ATN) model: formulas, parameters, and implementation notes"
author: "Working notes"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    number_sections: true
  word_document:
    toc: true
---

# Purpose

This document rewrites the Allometric Trophic Network (ATN) model in a **time- and space-explicit notation**, because the objective is to run one local ATN in each spatial grid cell of a spatially explicit ecosystem model.

The key change relative to a single local ATN is that every local state variable and every local rate can now depend on:

- species or node identity, indexed by $i$, $j$, or $k$;
- grid cell, indexed by $g$;
- time, indexed continuously by $t$.

The central state variable is therefore not only species biomass $B_i(t)$, but local species biomass:

$$B_{i,g}(t)$$

where $B_{i,g}(t)$ is the biomass density of species or node $i$ in grid cell $g$ at time $t$.

This document covers the three ATN model families used in `ATNr`:

- **Unscaled ATN**, following Binzer et al. (2016);
- **Scaled ATN**, following Delmas et al. (2017);
- **Unscaled ATN with nutrients**, following Schneider et al. (2016).

The notation below is written for a spatial model in which each grid cell has its own local food web and local biomass dynamics. A movement/dispersal term is included as optional, because a first implementation can run independent local ATNs, while a later implementation can couple cells through dispersal.

# 1. Spatial and temporal indexing

## 1.1 Indices

| Symbol | Meaning |
|----|----|
| $i$ | focal species or resource species |
| $j$ | consumer species |
| $k$ | resource species summed over in functional responses |
| $n$ | nutrient pool |
| $g$ | spatial grid cell |
| $t$ | time |
| $S$ | total number of species or model nodes in the regional species pool |
| $S_g(t)$ | number of species present in cell $g$ at time $t$ |
| $n_b$ | number of basal species in the regional model |
| $n_{b,g}(t)$ | number of basal species present in cell $g$ at time $t$ |
| $n_n$ | number of nutrient pools |

## 1.2 Local biomass state variable

The biomass of species $i$ in cell $g$ at time $t$ is:

$$B_{i,g}(t)$$

A species can be absent from a grid cell. There are two equivalent ways to represent absence.

### Option 1: zero biomass

$$B_{i,g}(t) = 0$$

### Option 2: presence mask

Define a presence indicator:

$$I_{i,g}(t)=
\begin{cases}
1, & \text{if species } i \text{ is present in cell } g \text{ at time } t \\
0, & \text{otherwise}
\end{cases}$$

Then local interactions are masked by $I_{i,g}(t)I_{j,g}(t)$.

In code, the simplest approach is usually to keep a global species vector of length $S$ and set absent species to zero biomass in a given cell.

## 1.3 Body mass

In standard ATN models, body mass is a fixed species trait:

$$M_i$$

For a spatially explicit model, the default assumption should be:

$$M_{i,g}(t) = M_i$$

This means that body mass does not change across cells or through time. However, the more general notation is:

$$M_{i,g}(t)$$

This allows future extensions where body mass changes with local temperature, ontogeny, or local trait variation.

Important: standard ATN models do **not** include ontogeny. Therefore, unless explicitly added, $M_i$ is fixed and only biomass $B_{i,g}(t)$ changes.

# 2. Local food-web convention

ATN models use a resource-by-consumer matrix convention.

The local food-web matrix in grid cell $g$ is:

$$fw_{ij,g}(t) =
\begin{cases}
1, & \text{if consumer } j \text{ eats resource } i \text{ in cell } g \text{ at time } t \\
0, & \text{otherwise}
\end{cases}$$

Thus:

- rows are resources;
- columns are consumers;
- $F_{ij,g}(t)$ is the feeding rate of consumer $j$ on resource $i$ in cell $g$ at time $t$;
- $B_{j,g}(t)F_{ij,g}(t)$ is the biomass of resource $i$ removed by consumer $j$ in cell $g$ per unit time.

If the food-web topology is fixed within each cell, write:

$$fw_{ij,g}$$

If species ranges, local co-occurrence, or body sizes change through time, write:

$$fw_{ij,g}(t)$$

# 3. Spatially explicit L matrix

The L matrix gives the probability that consumer $j$ successfully attacks resource $i$ based on body-size matching.

For each cell $g$ and time $t$, define:

$$z_{ij,g}(t) = \frac{M_{j,g}(t)}{M_{i,g}(t)R_{opt}}$$

where:

- $M_{i,g}(t)$ is the body mass of resource $i$ in cell $g$;
- $M_{j,g}(t)$ is the body mass of consumer $j$ in cell $g$;
- $R_{opt}$ is the optimal consumer/resource body-mass ratio.

The feeding kernel is:

$$L_{ij,g}(t) = \left[z_{ij,g}(t)\exp\left(1-z_{ij,g}(t)\right)\right]^\gamma$$

Links below the threshold are removed:

$$L_{ij,g}(t) = 0 \quad \text{if} \quad L_{ij,g}(t) < th$$

Basal species are not consumers, so:

$$L_{ij,g}(t) = 0 \quad \text{for all } j \leq n_b$$

If species are absent from the cell, links are also removed:

$$L_{ij,g}(t) = L_{ij,g}(t) I_{i,g}(t)I_{j,g}(t)$$

The local binary food web is:

$$fw_{ij,g}(t) = \mathbb{1}\left(L_{ij,g}(t)>0\right)$$

## Interpretation

The L-matrix function is maximal when:

$$z_{ij,g}(t)=1$$

which implies:

$$M_{j,g}(t)=R_{opt}M_{i,g}(t)$$

Thus consumer $j$ is most efficient at consuming resources approximately $R_{opt}$ times smaller than itself.

# 4. Local allometric and temperature-dependent rates

Many ATN parameters are calculated using allometric relationships with body mass and temperature.

For a pairwise parameter $p_{ij,g}(t)$, such as attack rate, capture coefficient, or handling time:

$$p_{ij,g}(t)
=
p_0 M_{i,g}(t)^{b_{prey}}M_{j,g}(t)^{b_{pred}}
\exp\left(
\frac{-E\left[T_0-T_{K,g}(t)\right]}{kT_{K,g}(t)T_0}
\right)$$

where:

- $p_0$ is a normalization constant;
- $M_{i,g}(t)$ is resource body mass;
- $M_{j,g}(t)$ is consumer body mass;
- $b_{prey}$ is the resource-mass exponent;
- $b_{pred}$ is the consumer-mass exponent;
- $E$ is activation energy;
- $T_{K,g}(t)$ is local temperature in Kelvin;
- $T_0$ is the reference temperature in Kelvin;
- $k$ is Boltzmann's constant.

If temperature is constant over time in a given cell, write $T_{K,g}$. If temperature is the same everywhere, write $T_K$.

# 5. Local basal growth and metabolic rates

The maximum growth rate of basal species $i$ in cell $g$ can be written as:

$$r_{i,g}(t) = r_0 M_{i,g}(t)^{b_r}
\exp\left(
\frac{-E_r\left[T_0-T_{K,g}(t)\right]}{kT_{K,g}(t)T_0}
\right)$$

The metabolic loss rate of species $i$ in cell $g$ is:

$$X_{i,g}(t) = X_0 M_{i,g}(t)^{b_X}
\exp\left(
\frac{-E_X\left[T_0-T_{K,g}(t)\right]}{kT_{K,g}(t)T_0}
\right)$$

If the model does not include temperature dependence, these simplify to:

$$r_i = r_0M_i^{b_r}$$

and:

$$X_i = X_0M_i^{b_X}$$

In many ATN applications:

$$b_r \approx b_X \approx -0.25$$

# 6. General spatial ATN biomass equations

The local dynamics of each species in each cell are given by:

$$\frac{dB_{i,g}(t)}{dt}
=
\text{local gains}
-
\text{local losses}
+
\mathcal{D}_{i,g}(t)$$

where $\mathcal{D}_{i,g}(t)$ is an optional dispersal or movement term.

For a first implementation with independent local ATNs:

$$\mathcal{D}_{i,g}(t)=0$$

## 6.1 Basal species / plants

For basal species $i$ in cell $g$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)r_{i,g}(t)G_{i,g}(t)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_j B_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

where:

- $B_{i,g}(t)r_{i,g}(t)G_{i,g}(t)$ is realized local basal production;
- $X_{i,g}(t)B_{i,g}(t)$ is local metabolic loss;
- $\sum_j B_{j,g}(t)F_{ij,g}(t)$ is local biomass lost to consumers;
- $\mathcal{D}_{i,g}(t)$ is optional dispersal, redistribution, or external biomass input.

## 6.2 Non-basal species / consumers

For consumer species $i$ in cell $g$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)\sum_j e_{j,g}(t)F_{ji,g}(t)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_j B_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

where:

- $F_{ji,g}(t)$ is the feeding rate of consumer $i$ on resource $j$ in cell $g$;
- $e_{j,g}(t)$ is the assimilation efficiency associated with resource $j$;
- $B_{i,g}(t)\sum_j e_{j,g}(t)F_{ji,g}(t)$ is biomass gained by consumer $i$;
- $X_{i,g}(t)B_{i,g}(t)$ is metabolic loss;
- $\sum_j B_{j,g}(t)F_{ij,g}(t)$ is biomass lost because species $i$ is eaten by its consumers.

# 7. Optional spatial coupling by dispersal

If cells are independent, set:

$$\mathcal{D}_{i,g}(t)=0$$

If cells are coupled by movement, a simple biomass-dispersal term is:

$$\mathcal{D}_{i,g}(t)
=
\sum_{g' \neq g} m_{i,g'\rightarrow g}(t)B_{i,g'}(t)
-
\sum_{g' \neq g} m_{i,g\rightarrow g'}(t)B_{i,g}(t)$$

where:

- $m_{i,g'\rightarrow g}(t)$ is the rate at which biomass of species $i$ moves from cell $g'$ into cell $g$;
- $m_{i,g\rightarrow g'}(t)$ is the rate at which biomass of species $i$ leaves cell $g$ for cell $g'$.

This term is not part of standard local ATNr, but it is the natural place to add spatial coupling in a spatial ATN.

# 8. Unscaled ATN model in space and time (WE WILL USE THIS ONE)

The unscaled model uses biological rates in their original units.

## 8.1 Functional response

For consumer $j$ feeding on resource $i$ in cell $g$:

$$F_{ij,g}(t)
=
\frac{a_{ij,g}(t)B_{i,g}(t)^{q_j}}
{1+c_{j,g}(t)B_{j,g}(t)+\sum_k h_{kj,g}(t)a_{kj,g}(t)B_{k,g}(t)^{q_j}}$$

where:

- $a_{ij,g}(t)$ is the local clearance/attack rate of consumer $j$ on resource $i$;
- $B_{i,g}(t)$ is local resource biomass;
- $q_j$ is the Hill exponent of consumer $j$;
- $c_{j,g}(t)$ is local intraspecific consumer interference;
- $B_{j,g}(t)$ is local consumer biomass;
- $h_{kj,g}(t)$ is handling time of consumer $j$ on resource $k$;
- the sum over $k$ includes all local resources of consumer $j$.

If the Hill exponent is the same for all consumers, write $q$ instead of $q_j$.

## 8.2 Local logistic growth limitation for basal species

For basal species $i$ in cell $g$:

$$G_{i,g}(t)
=
1-
\frac{\sum_{\ell \in basal}\alpha_{i\ell,g}(t)B_{\ell,g}(t)}{K_{i,g}(t)}$$

where:

- $K_{i,g}(t)$ is the local carrying capacity of basal species $i$;
- $\alpha_{i\ell,g}(t)$ is the competitive effect of basal species $\ell$ on basal species $i$ in cell $g$.

## 8.3 Full unscaled basal equation

For basal species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)r_{i,g}(t)
\left(
1-
\frac{\sum_{\ell \in basal}\alpha_{i\ell,g}(t)B_{\ell,g}(t)}{K_{i,g}(t)}
\right)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_jB_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

## 8.4 Full unscaled consumer equation

For consumer species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)\sum_j e_{j,g}(t)F_{ji,g}(t)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_jB_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

# 9. Scaled ATN model in space and time

The scaled model rescales biological rates relative to the growth rate of the smallest basal species.

## 9.1 Local scaling

If scaling is done separately in each cell, define the reference basal species in cell $g$ as the smallest basal species present in that cell. Let its growth rate be:

$$r_{ref,g}(t)$$

Then:

$$\tilde{r}_{i,g}(t)=\frac{r_{i,g}(t)}{r_{ref,g}(t)}$$

and:

$$\tilde{x}_{i,g}(t)=\frac{X_{i,g}(t)}{r_{ref,g}(t)}$$

If scaling is done globally, the same reference species and same $r_{ref}$ are used in all cells.

## 9.2 Scaled functional response

For consumer $j$ feeding on resource $i$ in cell $g$:

$$F_{ij,g}(t)
=
\frac{w_{ij,g}(t)B_{i,g}(t)^{q_j}}
{B_{0j,g}(t)^{q_j}+c_{j,g}(t)B_{j,g}(t)+\sum_k w_{kj,g}(t)B_{k,g}(t)^{q_j}}$$

where:

- $w_{ij,g}(t)$ is the local relative consumption weight;
- $B_{0j,g}(t)$ is the local half-saturation biomass of consumer $j$;
- $c_{j,g}(t)$ is local consumer interference;
- $q_j$ is the Hill exponent.

Diet weights usually satisfy:

$$\sum_i w_{ij,g}(t)=1$$

for each consumer $j$ in cell $g$.

## 9.3 Local logistic growth limitation

For basal species:

$$G_{i,g}(t)
=
1-
\frac{\sum_{\ell \in basal}\alpha_{i\ell,g}(t)B_{\ell,g}(t)}{K_{i,g}(t)}$$

## 9.4 Full scaled basal equation

For basal species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)\tilde{r}_{i,g}(t)G_{i,g}(t)
-
\tilde{x}_{i,g}(t)B_{i,g}(t)
-
\sum_j\frac{\tilde{x}_{j,g}(t)}{y_{j,g}(t)}B_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

## 9.5 Full scaled consumer equation

For consumer species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
\frac{\tilde{x}_{i,g}(t)}{y_{i,g}(t)}B_{i,g}(t)
\sum_j e_{j,g}(t)F_{ji,g}(t)
-
\tilde{x}_{i,g}(t)B_{i,g}(t)
-
\sum_j\frac{\tilde{x}_{j,g}(t)}{y_{j,g}(t)}B_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

where $y_{i,g}(t)$ is the maximum feeding rate relative to metabolic rate.

# 10. Unscaled ATN with nutrients in space and time

The nutrient model uses explicit nutrient pools in each grid cell.

## 10.1 State vector per grid cell

For cell $g$, the local state vector is:

$$y_g(t)=\left(N_{1,g}(t),\dots,N_{n_n,g}(t),B_{1,g}(t),\dots,B_{S,g}(t)\right)$$

where:

- $N_{n,g}(t)$ is the concentration of nutrient $n$ in cell $g$;
- $B_{i,g}(t)$ is species biomass in cell $g$.

## 10.2 Nutrient-model functional response

For consumer $j$ feeding on resource $i$ in cell $g$:

$$F_{ij,g}(t)
=
\frac{w_{ij,g}(t)b_{ij,g}(t)B_{i,g}(t)^{q_j}}
{M_{j,g}(t)\left[1+c_{j,g}(t)B_{j,g}(t)+\sum_k w_{kj,g}(t)h_{kj,g}(t)b_{kj,g}(t)B_{k,g}(t)^{q_j}\right]}$$

where:

- $w_{ij,g}(t)$ is local relative consumption weight;
- $b_{ij,g}(t)$ is local capture coefficient;
- $M_{j,g}(t)$ is consumer body mass;
- $h_{kj,g}(t)$ is handling time of consumer $j$ on resource $k$;
- $c_{j,g}(t)$ is local consumer interference.

## 10.3 Nutrient-limited basal growth

For basal species $i$ in cell $g$:

$$G_{i,g}(t)
=
\min_{n\in nutrients}
\left(
\frac{N_{n,g}(t)}{K_{ni,g}(t)+N_{n,g}(t)}
\right)$$

where:

- $N_{n,g}(t)$ is the concentration of nutrient $n$ in cell $g$;
- $K_{ni,g}(t)$ is the local half-saturation or uptake-efficiency parameter for nutrient $n$ by basal species $i$.

## 10.4 Nutrient dynamics

For nutrient $n$ in cell $g$:

$$\frac{dN_{n,g}(t)}{dt}
=
D_{n,g}(t)\left[S_{n,g}(t)-N_{n,g}(t)\right]
-
\sum_{i\in basal}v_{ni,g}(t)r_{i,g}(t)G_{i,g}(t)B_{i,g}(t)
+
\mathcal{D}^{N}_{n,g}(t)$$

where:

- $D_{n,g}(t)$ is the local nutrient turnover rate;
- $S_{n,g}(t)$ is the local nutrient supply concentration;
- $v_{ni,g}(t)$ is the relative content of nutrient $n$ in basal species $i$;
- $r_{i,g}(t)G_{i,g}(t)B_{i,g}(t)$ is realised local production of basal species $i$;
- $\mathcal{D}^{N}_{n,g}(t)$ is optional nutrient movement or external nutrient input.

If nutrients do not move among cells:

$$\mathcal{D}^{N}_{n,g}(t)=0$$

## 10.5 Full nutrient-model basal equation

For basal species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)r_{i,g}(t)G_{i,g}(t)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_jB_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

with:

$$G_{i,g}(t)=\min_n\left(\frac{N_{n,g}(t)}{K_{ni,g}(t)+N_{n,g}(t)}\right)$$

## 10.6 Full nutrient-model consumer equation

For consumer species $i$:

$$\frac{dB_{i,g}(t)}{dt}
=
B_{i,g}(t)\sum_j e_{j,g}(t)F_{ji,g}(t)
-
X_{i,g}(t)B_{i,g}(t)
-
\sum_jB_{j,g}(t)F_{ij,g}(t)
+
\mathcal{D}_{i,g}(t)$$

# 11. Local diet weights

Diet weights distribute consumer feeding effort among locally available resources.

A simple local normalization is:

$$w_{ij,g}(t)=
\frac{fw_{ij,g}(t)}{\sum_k fw_{kj,g}(t)}$$

if all resources are equally preferred.

In the Schneider-style model, weights can be based on the L matrix:

$$w_{ij,g}(t)=
\frac{L_{ij,g}(t)}{\sum_k L_{kj,g}(t)}$$

for each consumer $j$ in cell $g$.

The local normalization condition is:

$$\sum_i w_{ij,g}(t)=1$$

for every consumer $j$ with at least one local resource.

If consumer $j$ has no local resources, all its weights should be set to zero:

$$w_{ij,g}(t)=0 \quad \text{for all } i$$

# 12. Local assimilation efficiencies

Assimilation efficiency can depend on the resource type and can also vary by cell and time:

$$e_{i,g}(t)=
\begin{cases}
e_{P,g}(t), & \text{if resource } i \text{ is basal} \\
e_{A,g}(t), & \text{if resource } i \text{ is animal}
\end{cases}$$

Usually:

$$e_{P,g}(t)<e_{A,g}(t)$$

because plant biomass is less efficiently converted into consumer biomass than animal prey.

# 13. Local biomass flow interpretation

The local feeding rate $F_{ij,g}(t)$ is a per-biomass feeding rate of consumer $j$ on resource $i$.

The local biomass removed from resource $i$ by consumer $j$ is:

$$\Phi_{ij,g}^{loss}(t)=B_{j,g}(t)F_{ij,g}(t)$$

The local biomass gained by consumer $j$ from resource $i$ is:

$$\Phi_{ij,g}^{gain}(t)=e_{i,g}(t)B_{j,g}(t)F_{ij,g}(t)$$

Total local loss of species $i$ to consumers is:

$$\Phi_{i,g}^{consumed}(t)=\sum_jB_{j,g}(t)F_{ij,g}(t)$$

Total local consumption gain of consumer $j$ is:

$$\Phi_{j,g}^{food}(t)=B_{j,g}(t)\sum_i e_{i,g}(t)F_{ij,g}(t)$$

# 14. Local extinction threshold

A local population can be treated as extinct in cell $g$ if:

$$B_{i,g}(t)<ext$$

Then set:

$$B_{i,g}(t)=0$$

and usually:

$$I_{i,g}(t)=0$$

This is a **local extinction**, not necessarily a global extinction. Species $i$ may still persist in other cells.

# 15. Spatial solver structure

## 15.1 Independent local ATNs

If each cell is independent, solve one ODE system per cell:

$$\frac{dy_g(t)}{dt}=f_g\left(y_g(t),t,\theta_g(t)\right)$$

where:

- $y_g(t)$ is the state vector for cell $g$;
- $\theta_g(t)$ is the local parameter set;
- $f_g$ is the ATN derivative function for cell $g$.

This is equivalent to running one ATN per spatial cell.

## 15.2 Coupled spatial ATNs

If cells are coupled by dispersal, the full state vector is:

$$y(t)=\left(y_1(t),y_2(t),\dots,y_G(t)\right)$$

and the full system is:

$$\frac{dy(t)}{dt}=f\left(y(t),t,\theta(t)\right)$$

where the derivative for each cell includes local ATN dynamics plus spatial movement terms.

# 16. Parameter tables for the spatial ATN

## 16.1 Core state variables

| Symbol | Meaning | Dimension |
|----|----|----|
| $B_{i,g}(t)$ | biomass of species $i$ in cell $g$ | species $\times$ cells $\times$ time |
| $N_{n,g}(t)$ | nutrient concentration of nutrient $n$ in cell $g$ | nutrients $\times$ cells $\times$ time |
| $I_{i,g}(t)$ | local presence indicator | species $\times$ cells $\times$ time |
| $y_g(t)$ | local state vector | variable by model |
| $\mathcal{D}_{i,g}(t)$ | movement/dispersal term for species $i$ in cell $g$ | species $\times$ cells $\times$ time |

## 16.2 Structural parameters

| Symbol | Meaning | Spatial/time version |
|----|----|----|
| $S$ | total regional species pool size | usually global scalar |
| $S_g(t)$ | local richness in cell $g$ | cell and time specific |
| $n_b$ | number of basal species in regional pool | usually global scalar |
| $n_{b,g}(t)$ | number of basal species in cell $g$ | cell and time specific |
| $n_n$ | number of nutrient pools | usually global scalar |
| $M_i$ | species body mass | fixed by default |
| $M_{i,g}(t)$ | local/time-specific body mass | optional extension |
| $fw_{ij,g}(t)$ | local food-web adjacency | resource $\times$ consumer $\times$ cell $\times$ time |
| $L_{ij,g}(t)$ | local L-matrix value | resource $\times$ consumer $\times$ cell $\times$ time |
| $ext$ | extinction threshold | scalar or cell specific |

## 16.3 Environmental parameters

| Symbol | Meaning | Spatial/time version |
|----|----|----|
| $T_{K,g}(t)$ | temperature in Kelvin | cell and time specific |
| $T_0$ | reference temperature | global scalar |
| $k$ | Boltzmann constant | global scalar |
| $S_{n,g}(t)$ | nutrient supply concentration | nutrient, cell, and time specific |
| $D_{n,g}(t)$ | nutrient turnover rate | nutrient, cell, and time specific |
| $K_{i,g}(t)$ | local carrying capacity for basal species | basal species, cell, and time specific |

## 16.4 Biological-rate parameters

| Symbol | Meaning | Spatial/time version |
|----|----|----|
| $r_{i,g}(t)$ | basal maximum growth rate | basal species, cell, and time specific |
| $X_{i,g}(t)$ | metabolic loss rate | species, cell, and time specific |
| $a_{ij,g}(t)$ | attack/clearance rate | resource, consumer, cell, and time specific |
| $b_{ij,g}(t)$ | capture coefficient | resource, consumer, cell, and time specific |
| $h_{ij,g}(t)$ | handling time | resource, consumer, cell, and time specific |
| $c_{j,g}(t)$ | consumer interference | consumer, cell, and time specific |
| $q_j$ | Hill exponent | consumer specific, often fixed |
| $e_{i,g}(t)$ | assimilation efficiency of resource $i$ | resource, cell, and time specific |
| $w_{ij,g}(t)$ | diet weight | resource, consumer, cell, and time specific |
| $B_{0j,g}(t)$ | half-saturation biomass in scaled model | consumer, cell, and time specific |
| $y_{j,g}(t)$ | maximum feeding rate relative to metabolism | consumer, cell, and time specific |

## 16.5 Nutrient-model parameters

| Symbol | Meaning | Spatial/time version |
|----|----|----|
| $K_{ni,g}(t)$ | nutrient half-saturation parameter | nutrient, basal species, cell, and time specific |
| $v_{ni,g}(t)$ | nutrient content of basal species | nutrient, basal species, cell, and time specific |
| $N_{n,g}(t)$ | nutrient concentration | nutrient, cell, and time specific |
| $S_{n,g}(t)$ | nutrient supply concentration | nutrient, cell, and time specific |
| $D_{n,g}(t)$ | nutrient turnover rate | nutrient, cell, and time specific |

# 17. What changes relative to a non-spatial ATN?

The non-spatial ATN has:

$$B_i(t), F_{ij}(t), r_i, X_i, fw_{ij}$$

The spatial ATN has:

$$B_{i,g}(t), F_{ij,g}(t), r_{i,g}(t), X_{i,g}(t), fw_{ij,g}(t)$$

The biological meaning is the same, but every biomass, rate, interaction, and environmental driver can now be local to a grid cell.

# 18. Implementation checklist for one ATN per grid cell

For each grid cell $g$:

1.  Define the local species pool using species ranges or local occurrence data.
2.  Set initial local biomasses $B_{i,g}(0)$.
3.  Set local environmental variables such as $T_{K,g}(t)$, $K_{i,g}(t)$, $S_{n,g}(t)$, or NPP-derived parameters.
4.  Build the local L matrix $L_{ij,g}(t)$ or use a fixed local food web $fw_{ij,g}$.
5.  Calculate local allometric rates: $r_{i,g}(t)$, $X_{i,g}(t)$, $a_{ij,g}(t)$, $b_{ij,g}(t)$, $h_{ij,g}(t)$.
6.  Normalize local diet weights $w_{ij,g}(t)$.
7.  Compute the local feeding matrix $F_{ij,g}(t)$.
8.  Compute local biomass derivatives $dB_{i,g}(t)/dt$.
9.  Compute local nutrient derivatives $dN_{n,g}(t)/dt$ if using the nutrient model.
10. Apply local extinction threshold.
11. If spatial coupling is included, add dispersal terms $\mathcal{D}_{i,g}(t)$.

# 19. Recommended first implementation

The recommended first spatial implementation is:

$$\mathcal{D}_{i,g}(t)=0$$

That is, run one independent ATN per grid cell. This lets you check that the local ATN equations work before adding spatial movement.

Then, after the local model is validated, add dispersal:

$$\mathcal{D}_{i,g}(t)
=
\sum_{g' \neq g} m_{i,g'\rightarrow g}(t)B_{i,g'}(t)
-
\sum_{g' \neq g} m_{i,g\rightarrow g'}(t)B_{i,g}(t)$$

This creates a true spatially coupled ATN model.

# 20. Validation strategy against ATNr

ATNr is non-spatial by default. To compare a spatial implementation against ATNr, validate one cell at a time.

For a chosen grid cell $g$:

1.  Extract the local food web $fw_{ij,g}$.
2.  Extract local body masses $M_{i,g}$ or fixed $M_i$.
3.  Extract local starting biomasses $B_{i,g}(0)$.
4.  Extract local environmental parameters such as $T_{K,g}$.
5.  Run ATNr as if this grid cell were a non-spatial food web.
6.  Run the custom spatial code for the same grid cell with $\mathcal{D}_{i,g}=0$.
7.  Compare one-step derivatives:

$$f_{custom,g}\left(y_g(0),t_0\right)
\quad \text{vs.} \quad
f_{ATNr}\left(y_g(0),t_0\right)$$

8.  If derivatives match, compare trajectories:

$$B_{i,g}^{custom}(t)
\quad \text{vs.} \quad
B_i^{ATNr}(t)$$

# 21. Common mistakes in the spatial version

1.  **Forgetting the cell index.** Every local biomass and local rate should carry $g$.
2.  **Mixing regional and local species pools.** A species absent from a cell should have $B_{i,g}=0$ and should not feed locally.
3.  **Using the wrong food-web orientation.** Rows are resources and columns are consumers.
4.  **Normalizing diet weights globally instead of locally.** $w_{ij,g}$ must be normalized across resources available to consumer $j$ in cell $g$.
5.  **Using global temperature instead of local temperature.** Temperature-dependent rates should use $T_{K,g}(t)$.
6.  **Treating local extinction as global extinction.** $B_{i,g}<ext$ means species $i$ is extinct in cell $g$, not necessarily everywhere.
7.  **Adding dispersal before validating local dynamics.** First validate local ATN equations with $\mathcal{D}_{i,g}=0$.

# 22. References

- Binzer, A., Guill, C., Rall, B. C., & Brose, U. (2016). Interactive effects of warming, eutrophication and size structure: impacts on biodiversity and food-web structure. *Global Change Biology*, 22, 220-227.
- Delmas, E., Brose, U., Gravel, D., Stouffer, D. B., & Poisot, T. (2017). Simulations of biomass dynamics in community food webs. *Methods in Ecology and Evolution*, 8, 881-886.
- Gauzens, B., Brose, U., Delmas, E., & Berti, E. (2023). ATNr: Allometric Trophic Network models in R. *Methods in Ecology and Evolution*, 14, 2766-2773.
- Schneider, F. D., Brose, U., Rall, B. C., & Guill, C. (2016). Animal diversity and ecosystem functioning in dynamic food webs. *Nature Communications*, 7, 12718.
- Williams, R. J., & Martinez, N. D. (2000). Simple rules yield complex food webs. *Nature*, 404, 180-183.
- Yodzis, P., & Innes, S. (1992). Body size and consumer-resource dynamics. *The American Naturalist*, 139, 1151-1175.
