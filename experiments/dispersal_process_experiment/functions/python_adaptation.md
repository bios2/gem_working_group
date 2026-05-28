# Python adaptation of R functions

This file specifies what has changed for the dispersal functions from the R prototypes

## 1 - Function names
Function names within the dispersal.py that describe change in biomasses must have the following naming convention : "<process>_delta".

diffuse_density_dependent --> disperse_delta

We removed the diffuse (density-independent) function, because it could be retrieved from the more general diffuse_density_dependent function

## 2 - Time scale resolution

We added the time resolution (dt) as an input in the disperse_delta function, and ensured that the biomass flux would be scaled to match the time scale resolution

## 3 - Numpy array manipulation
The R functions were using a biomass matrix, and species identity was not defined. However, within the ATN, many species can interact together. We changed the data format from [X, Y] to [X, Y, S] to account for species identity. This also necessitated changes to parts of the code that were determining neighbours for each cell, in order to compute diffusion. Three changes were made : 

*enforce_boundary_conditions* now has a dummy trailing edge to be able to match other functions format ([X, Y, 1] format)

*disperse_delta* now has two assert statement, to ensure that it has the same spatial dimensions as the boundary_conditions matrix, and that it is in the same format as the net_growth matrix

*run_diffusion_simulation* now calls the updated *enforce_boundary_conditions* function

## 4 - Parameter format input (maximum dispersal distance and metabolism)

The R functions were considering that both metabolism and maximum dispersal distance were homogenous across space and species. However, metabolism will vary across species (different body masses) and space (different temperatures), and maximum dispersal distance will vary across species.

This will not pose any problem, since the adapter script will ensure that all of them are in the right format. The signatures and docstring of dispersal function were updated to take into account the new data format.


## 5 - Integrate *enforce_boundary_conditions* within the *disperse_delta* function ?

Right now, the boundary conditions are enforced within the function that implements a test simulation run, but this simulation function will not be called in the engine. 

After discussing with Claude, I realized this will not be a necessary change. Computing *enforce_boundary_conditions* within *disperse_delta* would necessitate computing it at every time step. This would lead to unecessary computations, since it does not change from one time step to the next.

*enforce_boundary_conditions* instead work as an helper that will be called once in the processes.py file (for initiation before starting the simulation)




