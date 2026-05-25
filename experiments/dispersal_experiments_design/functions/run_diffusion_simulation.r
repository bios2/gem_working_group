# run_diffusion_simulation.r
# Runs a density-independent diffusion simulation over n_time steps.
#
# Arguments:
#   init_pop  - numeric vector of initial populations in terra cell order
#               (row-major: left to right, top to bottom); length = n_row * n_col
#   n_row     - number of raster rows
#   n_col     - number of raster columns
#   disp_rate - fraction of individuals dispersing per step
#   n_time    - number of time steps to simulate
#
# Returns a list of (n_time + 1) matrices (n_row x n_col), where index 1 is
# the initial state. Convert a matrix m back to terra cell-order values with
# c(t(m)).

run_diffusion_simulation <- function(init_pop, n_row, n_col,
                                     disp_rate, n_time) {
  # pre-compute neighbour counts once; reused at every time step
  nbr_count <- build_rook_nbr_count_matrix(n_row, n_col)

  # terra cell order is row-major; byrow = TRUE preserves that in the matrix
  pop_grid <- matrix(init_pop, nrow = n_row, ncol = n_col, byrow = TRUE)

  # pre-allocate list to store all snapshots (t=0 through t=n_time)
  results      <- vector("list", n_time + 1)
  results[[1]] <- pop_grid

  for (t in seq_len(n_time)) {
    # advance grid by one dispersal step and store the new state
    pop_grid         <- diffuse_step_vectorised(pop_grid, disp_rate, nbr_count)
    results[[t + 1]] <- pop_grid
  }

  return(results)
}
