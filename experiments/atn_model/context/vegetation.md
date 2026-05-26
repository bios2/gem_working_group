$G_i = NPP_m^T \cdot \psi \cdot \delta t_{NPP} \cdot (1 - f_{\text{struct.}}) \cdot C_i$

$A_{\text{cell}}$: dropped because we have a regional model for $A$

$C_i$: for herbs:

```
1 |\.
  |  \
  |   \___________
  +-----------------> B_trees
```

$$C_{\text{herb}} = \frac{\alpha_{\text{herbs}}}{\alpha_{\text{herbs}} + B_{\text{trees}}}$$

![herb curve](plot_herb.png)

and for trees &nbsp;&nbsp;&nbsp;&nbsp; $(1 - C_{\text{herb}})$

```
  |        ________
  |    ___/
  | __/
  |/
  +-----------------> B_trees
```

$$C_{\text{trees}} = 1 - C_{\text{herb}} = \frac{B_{\text{trees}}}{\alpha_{\text{herbs}} + B_{\text{trees}}}$$

![trees curve](plot_trees.png)

parameter &nbsp;&nbsp;&nbsp;&nbsp; $\alpha_{\text{herbs}}$ shapes the curve.

---

## Equation explanation

*All definitions below are drawn exclusively from Harfoot et al. (2014) Text S1 (pbio.1001841.s018.docx), Terrestrial plant model section.*

$G_i$ is the growth rate of leaf biomass for vegetation type $i$ (evergreen or deciduous) over one model timestep. All biomass units are grams of wet biomass (Text S1, opening note). It appears in the terrestrial leaf biomass difference equation:

$$B_{l,(t+\Delta t)} = B_{l,t} + \delta B_{l,\text{Growth}} - \delta B_{l,\text{Mort}}$$

where $\delta B_{l,\text{Growth}} = G_{\text{ever}} + G_{\text{decid}}$ is the sum of the growth contributions from evergreen and deciduous stocks.

### Terms

**$NPP_m^T$** — monthly terrestrial NPP in kg C m⁻² month⁻¹ (Text S1, Terrestrial plant model). It is derived in two steps. First, yearly NPP is modelled following the Miami model (Lieth 1975), where $NPP_{\text{max}}$ is the maximum possible net primary production, $c_p$ and $m_p$ are coefficients relating NPP to temperature, and $\rho$ relates NPP to total annual precipitation $P$. Monthly NPP is then estimated by multiplying the annual value by $\omega$, "an estimate of the contribution that the NPP in month $m$ makes to yearly total NPP", a term introduced "to capture seasonal patterns of productivity" and estimated from Terra/MODIS net primary productivity data (Text S1, Terrestrial plant model).

**$\psi$** — "a conversion factor from carbon to wet matter" (Text S1, Terrestrial plant model). Converts $NPP_m^T$ from carbon units into grams of wet biomass, which are the units used throughout the model.

**$\delta t_{NPP}$** — "a scalar to convert [NPP] from its monthly value to the model time step" (Text S1, Terrestrial plant model). Scales the monthly NPP value to match the duration of one model timestep.

**$(1 - f_{\text{struct.}})$** — leaf allocation fraction. $f_{\text{struct.}}$ is "the fractional allocation of primary production to structural tissue" (Text S1, Terrestrial plant model). The complement $(1 - f_{\text{struct.}})$ is therefore the fraction of NPP directed to leaves rather than to woody or structural biomass. Text S1 defines $f_{\text{struct.}}$ as a derived quantity based on the fine root mortality rate $\mu_{\text{FineRoot}}$ and monthly average temperature.

**$C_i$** — competitive partitioning coefficient introduced in these notes to replace the original Madingley terms. In the original Text S1 equations, the partition between evergreen and deciduous growth uses two terms that $C_i$ absorbs:

- $f_{\text{ever}}$: "the proportion of NPP produced by evergreen leaves at a particular location" (Text S1). In the original Madingley formulation this is derived from frost frequency and controls how NPP is split between evergreen and deciduous stocks.
- $f_{\text{LeafMort}}$: "the proportion of total mortality that is leaf mortality" (Text S1). A scalar that scales the growth term so that the leaf allocation matches the leaf component of total plant mortality.

The $C_i$ term here replaces both $f_{\text{ever}} \cdot f_{\text{LeafMort}}$ (for evergreen) and $(1-f_{\text{ever}}) \cdot f_{\text{LeafMort}}$ (for deciduous) with a biomass-based competitive partition between herbs and trees.

**$A_{\text{cell}}$** — grid-cell area in m². Present in the original Text S1 growth equations, where it scales the per-unit-area NPP to the total biomass increment for the whole grid cell. Dropped here because the regional model supplies $NPP_m^T$ already as a grid-cell total rather than an areal density (noted in the handwritten original).