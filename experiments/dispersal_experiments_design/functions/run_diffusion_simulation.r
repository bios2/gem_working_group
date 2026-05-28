# run_diffusion_simulation.r
# Runs a density-dependent diffusion simulation over n_time steps.
#
# Arguments:
#   init_pop          - biomass_matrix
#   n_row             - number of raster rows
#   n_col             - number of raster columns
#   max_disp_rate     - maximum per-capita dispersal rate (parameter a in Ryser et al. 2021)
#   n_time            - number of time steps to simulate
#   net_growth_matrix - n_row x n_col matrix of per-capita net growth rates (nu_{i,z});
#                       output of the ATN; held constant across steps when ATN is not coupled
#   metabolism        - scalar; species metabolic rate x_i; sigmoid inflection point
#   b                 - sigmoid steepness (default 10, as in Ryser et al. 2021)
#
# Returns a list of (n_time + 1) matrices (n_row x n_col), where index 1 is
# the initial state. Convert a matrix m back to terra cell-order values with
# c(t(m)).

run_diffusion_simulation <- function(biomass_matrix,
                                     max_disp_rate, n_time,
                                     net_growth_matrix, metabolism, b = 10) {
  # pre-compute neighbour counts once; reused at every time step
  boundary_number <- enforce_boundary_conditions(nrow(biomass_matrix), ncol(biomass_matrix))

  # pre-allocate list to store all snapshots (t=0 through t=n_time)
  results      <- vector("list", n_time + 1)
  results[[1]] <- biomass_matrix

  for (t in seq_len(n_time)) {
    # advance grid by one dispersal step and add net migration to current biomass
    biomass_matrix         <- biomass_matrix + diffuse_density_dependent(biomass_matrix, net_growth_matrix,
                                                             metabolism, max_disp_rate,
                                                             boundary_number, b)
    results[[t + 1]] <- biomass_matrix
  }

  return(results)
}
