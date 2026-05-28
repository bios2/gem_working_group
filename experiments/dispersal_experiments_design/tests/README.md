# This folder contains expected behaviour of the dispersal functions

All tests must be represented as a plot with time in x axis and distance  of the invasion front (from the origin) in the y axis. Biomass is initiated within only one cell (the central cell). The invasion front can be calculated as the farthest cell above a threshold 

The outputs of tests must be written in the figures folder as .png figure

Here are listed the different tests and the expected behaviours (these hypotheses may be proven wrong).

Tests are done in a 51x51 grid and 50 time steps.

## Test 1 
H : The function diffuse (density-independent dispersal) should output a sublinear relationship between dispersal distance and time

O : Expected behaviour confirmed


## Test 2 
H : The functions diffuse and diffuse_density_dependent should have the same output if b = 0 and max_disp_rate/2 = disp_rate

O : Expected behaviour confirmed

## Test 3 
H : Density-dependent dispersal should lead to a faster initial dispersal than density-independent dispersal under the same conditions 

O : Contrary to expected behaviour, density-dependent dispersal and independent dispersal have the same pace initially, but density dependent dispersal is better at having sustained dispersal. The fact that they are both identical is actually an artifact of the way we compute invasion front distances (since some biomass always move at each time step, we have to set a threshold). The density-dependent formulation in fact does have a faster initial dispersal than density-independent dispersal. 

## Test 4 
H : With the same resource grid (growth_rate grid), density-dependent dispersal should be higher for species with high metabolism (= 1.0) than low metabolism  (=0.1) (hummingbird vs black bear scenario), all else being equal

O : Expected behaviour confirmed

## Test 5 : 
H : All else being equal, species with high maximum dispersal (a = 0.8) should disperse further than species with low dispersal (a = 0.2)

O : Expected behaviour confirmed

