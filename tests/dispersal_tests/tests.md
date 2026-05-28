# This file contains expected behaviour of the dispersal functions

All tests must be represented as a plot with time in x axis and distance  of the invasion front (from the origin) in the y axis. Biomass is initiated within only one cell (the central cell). The invasion front can be calculated as the farthest cell above a threshold 

Here are listed the different tests and the expected behaviours (these hypotheses may be proven wrong).

Tests are done in a 51x51 grid and 50 time steps.

Each test must be written in a different file and output a .png figure called test_<number>.png

## Test 1 
H : When b = 0, the relationship between invasion front (1E-06 threshold) and time should be sublinear

## Test 2
H : With the same resource grid (growth_rate grid), density-dependent dispersal should be higher for species with high metabolism (= 1.0) than low metabolism  (=0.1) (hummingbird vs black bear scenario), all else being equal

## Test 3
H : All else being equal, species with high maximum dispersal (a = 0.8) should disperse further than species with low dispersal (a = 0.2)


