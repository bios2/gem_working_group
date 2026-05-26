# Function descriptions

This file contains function descriptions for the implementation of diffusion in the Allometric Trophic Network model.

## enforce_boundary_conditions.R
This function ensures that boundary are reflective (individuals on the borders cannot flow out of the system). It takes into input the grid matrix and outputs the number of boundaries for eah cell.

## diffuse.R
This function implements one step of density-independent dispersal. It takes into input the current state of the system (biomass matrix), the dispersal coefficient (species-specific) and the output of enforce_boundary_conditions.R. It outputs the net migration for each cell.

## diffuse_density_dependent.R
This function implements one step of density-dependent dispersal. It takes into input the current state of the system (biomass matrix), the maximum dispersal rate (a), the environmental sensitivity of dispersal (b), the metabolism (species-specific, output of the ATN) and the net growth rate (species and cell-specific, output of the ATN). It outputs the net migration for each cell. 

Note : this gives the same result as the diffuse.R function, if b=0. (disp_rate = a/2)

## run_diffusion_simulation.R
This function runs a simulation for one species for density dependent dispersal. It takes as an input an initial distribution for a species (biomass matrix), the same set of parameters as diffuse_density_dependent.R and number of time steps required for the simulation. It outputs a list of matrices, where each list index corresponds to the system state for one species at a specific time step.