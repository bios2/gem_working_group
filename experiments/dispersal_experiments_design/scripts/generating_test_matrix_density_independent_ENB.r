####################### 
#Generate a very basic diffusion matrix, 
#where species travel like ink in water through the matrix 

#Created: 25 May 2026 
#Created by: ENB
#Last updated: 


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

#Part 0. Script setup
#clear R's brain
rm(list=ls())

#load relevant libraries for script
pkgs <- c("tidyverse", "sp", "sf", "terra", "raster",  "parallel", 
          "lubridate", "elevatr", 
          'vcfR', 'robust', 'vegan', 'qvalue', 'data.table', 'maps',
          'doParallel')
#install.packages(pkgs)
lapply(pkgs, library, character.only = TRUE)
rm(pkgs)

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
 #Create a diffusion matrix of ~100 empty cells 
 test_rast <- rast(nrow = 10, ncol = 10, vals = NA)

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 1. Density-independent diffusion simulation

source("experiments/dispersal_experiments_design/functions/build_rook_nbr_count_matrix.r")
source("experiments/dispersal_experiments_design/functions/diffuse_step_vectorised.r")
source("experiments/dispersal_experiments_design/functions/run_diffusion_simulation.r")

# Simulation parameters
n_time    <- 100
disp_rate <- 0.05
n_row     <- nrow(test_rast)
n_col     <- ncol(test_rast)
n_cells   <- n_row * n_col

# Initial population: two origin cells at full density
init_pop                                        <- rep(0, n_cells)
origin_cells                                    <- c(cellFromRowCol(test_rast, 3, 3),
                                                     cellFromRowCol(test_rast, 8, 8))
init_pop[origin_cells]                          <- 100

# Run simulation
sim_results <- run_diffusion_simulation(
  init_pop  = init_pop,
  n_row     = n_row,
  n_col     = n_col,
  disp_rate = disp_rate,
  n_time    = n_time
)

# Convert each time-step matrix to a terra SpatRaster
# c(t(m)) converts row-major matrix back to terra cell-order vector
pop_rasters <- lapply(sim_results, function(grid) {
  r <- test_rast
  values(r) <- c(t(grid))
  r
})

# Sanity check: total population should be conserved across all time steps
total_pop <- sapply(pop_rasters, function(r) sum(values(r), na.rm = TRUE))
cat("Population range across time steps (should be constant):",
    round(range(total_pop), 4), "\n")

# Plot a selection of time steps with a shared colour scale
plot_steps   <- c(1, 2, 10, 25, 50, 100)
all_vals     <- unlist(lapply(plot_steps, function(t) values(pop_rasters[[t]])))
global_range <- range(all_vals, na.rm = TRUE)

par(mfrow = c(1, length(plot_steps)), mar = c(2, 2, 2, 1))
for (t in plot_steps) {
  plot(pop_rasters[[t]], main = paste("t =", t - 1),
       col   = hcl.colors(50, "YlOrRd", rev = TRUE),
       range = global_range)
}

