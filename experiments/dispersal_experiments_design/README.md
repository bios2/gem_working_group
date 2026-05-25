# Dispersal experiments

Monday : The objective for today is to build a dispersal model on a grid that mimics a diffusion process. To this end, we take inspiration from the Ryzer paper (2021) - Landscape heterogeneity buffers biodiversity of
simulated meta-food-webs under global change through rescue and drainage effects.

The Ryzer paper is slightly more complicated than a standard dispersal on a grid, as it accounts for emigration triggers,
such as resource availability, predation pressure, and inter- and intra-specific
competition.

Another difference is that the Ryzer model works with patch coordinates (dispersal can reach any patch, but it decreases with distance), while dispersal on a grid will always occur between adjacent cells.

For now, the dispersal function will take as an input a matrix of biomass, and will output an updated matrix of biomass after dispersal has occured. Once the functions for the dispersal on a grid are done, we will need to validate the model. To this end, we will initiate an homogenous grid with a peak of biomass in the center and we will verify that biomass spreads outwards over time.