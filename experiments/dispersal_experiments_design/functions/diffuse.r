# diffuse.r
# Advances the population grid by one density-independent diffusion step.
#
# Each cell loses (disp_rate / 4) * n_nbrs of its population, where n_nbrs is
# the number of rook neighbours (hard-wall boundary: no flux through edges
# without an adjacent cell). Individuals disperse equally to all neighbours.
#
# Arguments:
#   biomass_matrix  - n_row x n_col numeric matrix of current populations
#   disp_rate - Dispersal coefficient
#   boundary_number - n_row x n_col integer matrix from build_rook_nbr_count_matrix()
#
# Returns an n_row x n_col numeric matrix of updated biomass immigration - emigration.
#
# Implementation uses 2D matrix shifts instead of per-cell loops, making it
# efficient for large rasters.

diffuse <- function(biomass_matrix, disp_rate, boundary_number) {
  n_row <- nrow(biomass_matrix)
  n_col <- ncol(biomass_matrix)

  # biomass each cell sends toward each one of its four potential neighbours
  flux_per_edge <- biomass_matrix * (disp_rate / 4)

  # total biomass leaving a cell equals flux_per_edge times the number of
  # actual neighbours (edge/corner cells have fewer, so they lose less)
  flux_out <- flux_per_edge * boundary_number

  # zero-padding strips used to enforce reflective boundary conditions during shifts
  zero_row <- matrix(0, nrow = 1,     ncol = n_col)
  zero_col <- matrix(0, nrow = n_row, ncol = 1)

  # shift flux_per_edge in each cardinal direction so that cell (r, c) receives
  # the flux emitted toward it by its neighbour in that direction
  from_north <- rbind(zero_row,               flux_per_edge[-n_row, ])  # row above sends down
  from_south <- rbind(flux_per_edge[-1, ],    zero_row)                 # row below sends up
  from_west  <- cbind(zero_col,               flux_per_edge[, -n_col])  # col to the left sends right
  from_east  <- cbind(flux_per_edge[, -1],    zero_col)                 # col to the right sends left

  # net biomass: subtract total outflow, add inflow from all four directions
  return(- flux_out + from_north + from_south + from_west + from_east)
}
