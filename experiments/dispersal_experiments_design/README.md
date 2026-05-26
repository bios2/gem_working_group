# Dispersal experiments

## Monday
The objective for today is to build a dispersal model on a grid that mimics a diffusion process. To this end, we take inspiration from the Ryzer paper (2021) - Landscape heterogeneity buffers biodiversity of
simulated meta-food-webs under global change through rescue and drainage effects.

The Ryzer paper is slightly more complicated than a standard dispersal on a grid, as it accounts for emigration triggers,
such as resource availability, predation pressure, and inter- and intra-specific
competition.

Another difference is that the Ryzer model works with patch coordinates (dispersal can reach any patch, but it decreases with distance), while dispersal on a grid will always occur between adjacent cells.

For now, the dispersal function will take as an input a matrix of biomass, and will output an updated matrix of biomass after dispersal has occured. Once the functions for the dispersal on a grid are done, we will need to validate the model. To this end, we will initiate an homogenous grid with a peak of biomass in the center and we will verify that biomass spreads outwards over time.

## Tuesday 
The objective for today is to implement dispersal using the density-dependant formula from the Ryser model (2021), and to scale the model to take into account multiple species that interact.

To facilitate integration with other teams, here is our code architecture

### Function 1 : Emigration rate
#### Inputs 
$a$ : Maximal dispersal rate (parameter, depends on external trait data)
$b$ : Dispersal sensitivity to environmental conditions, can also be seen as dispersal plasticity (paramater, arbitrary value)
$x_i$ : Metabolism of species i (output of the ATN)
$v_{i,z}$ : Net population growth rate (output of the ATN)
species : Species identity
biomass_matrix

#### Outputs 
$E_{i,z}$ : matrix of emigration rates for species i

#### External data needs
Species list
Trait data (body mass, other dispersal related traits)
Allometric coefficients

#### Dependencies
Function to compute metabolism ($x_i$)
Function to compute net population growth rate ($v_{i,z}$)

#### Test cases
Compare dispersal distance over time of a bird (high $a$) vs a lizard (low $a$)
Compare dispersal distane over time high metabolism species, small species (e.g. hummingbird) vs a low metabolism, large species (e.g. black bear)

#### Integration tests
Compare static resources (baseline scenario) with resources shifting northwards (climate change scenario)

Test model behaviour (stability, chaotic dynamic) when an herbivore and its consumer interact

#### Minimal implementation
Density-independant dispersal on 1-Dimensional homogenous grid

### Function 2 : Immigration rate
#### Inputs 
biomass_matrix
$E_{i,z} : emigration_rates_matrix (output of emigration rate function)

#### Outputs 
$I_{i,z}$ : matrix of immigration rates for species i

#### External data needs
None

#### Dependencies
Function to compute emigration rates

#### Test cases
Same as emigration ?

#### Integration tests
Same as emigration ?

#### Minimal implementation
Same as emigration ?