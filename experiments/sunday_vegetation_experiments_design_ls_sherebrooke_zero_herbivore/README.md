# Vegetation model design notes

## Problems to address

The current vegetation representation relies on broad functional groups.

Vegetation growth appears to be independent of existing biomass, or at least assumes a full vegetation stock in a way that resembles a simple MTE-style formulation.

There is no explicit competition between functional groups for space, light, water, nutrients, or biomass allocation.

The model should represent sensitivity to disturbances.

Tree sensitivity to droughts should be represented, possibly through a distribution of sensitivities rather than a single fixed response.

Seed consumption may need to be represented because it can affect vegetation recruitment and regeneration.

## Expected model behaviour

Vegetation growth should follow a logistic or biomass-limited growth formulation rather than unlimited growth.

Herbivory should reduce vegetation biomass and growth.

The model should be able to reproduce a Whittaker-style biome or vegetation distribution plot.

Vegetation should be able to escape herbivory through height or structural development.

Seasonality should affect vegetation biomass and growth.

The temperature-growth relationship should stop growth at or below 0 °C.

## First exercise

Extract the existing vegetation model and rewrite it as a function with temperature and precipitation as arguments.

Run time series simulations for different vegetation functional groups.

Reinterpret the current functional groups as herbs, coniferous trees, and deciduous trees.

Add a non-structured deciduous functional group if needed.

## Possible model structure

A possible solution is to use a hierarchical vegetation model.

At the regional scale, the model could represent disturbance dynamics, especially fire disturbance, which would drive stand age and disturbance-frequency distributions.

At the local scale, the model could represent biomass dynamics following a Harfoot/Madingley-like vegetation formulation, but with biomass as a time-dependent state variable.

A limitation of using only a prescribed biomass trajectory, `B(t)`, is that biomass becomes difficult to compute consistently because parameters and environmental drivers are time-dependent.

An alternative is to represent vegetation using cohorts or structured state variables.

For example, vegetation biomass could be represented as `B(a, t)`, where `a` is cohort age or structural stage and `t` is time.

The model could track separate biomass dynamics for herbs, coniferous trees, and deciduous trees.

Possible state variables include `dH/dt`, `dC/dt`, and `dT/dt`, where `H` represents herbs, `C` represents coniferous trees, and `T` represents deciduous trees.

This structure creates two main challenges: coexistence among vegetation groups and transitions between vegetation states or cohorts.

The spatial unit could be approximately 5 m, which may allow the model to neglect mixed cases at first.