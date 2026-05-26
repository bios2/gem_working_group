# diffuse_density_dependent.r
# Advances the population grid by one density-dependent diffusion step.
#
# Per-capita emigration rate is a sigmoid function of the local net growth rate
# (Ryser et al. 2021, Eq. 10): high when conditions are bad (net growth below
# metabolic demands), low when conditions are good. Survival during dispersal
# is set to 1 (no matrix mortality). Dispersal occurs only between the four
# adjacent neighbours (reflective boundary conditions).
#
# Arguments:
#   biomass_matrix     - n_row x n_col numeric matrix of current biomass densities
#   net_growth_matrix  - n_row x n_col numeric matrix of per-capita net growth
#                      rates (nu_{i,z}): feeding gains minus predation losses
#                      minus metabolism, divided by biomass; computed from ATN
#                      state before calling this function
#   metabolism         - scalar; species metabolic rate x_i; acts as the sigmoid
#                      inflection point (dispersal switches around this value)
#   max_disp_rate    - scalar; maximum per-capita dispersal rate (parameter a)
#   boundary_number  - n_row x n_col integer matrix from build_rook_nbr_count_matrix()
#   b                - sigmoid steepness (default 10, as in Ryser et al. 2021)
#
# Returns an n_row x n_col numeric matrix of updated biomass densities.

diffuse_density_dependent <- function(biomass_matrix, net_growth_matrix, metabolism,
                                   max_disp_rate, boundary_number, b = 10) {
  n_row <- nrow(biomass_matrix)
  n_col <- ncol(biomass_matrix)

  # per-capita dispersal rate: sigmoid of (metabolism - net_growth_rate)
  # rises toward max_disp_rate when net growth falls below metabolic demands,
  # falls toward 0 when the population is growing well
  d_grid <- max_disp_rate / (1 + exp(-b * (metabolism - net_growth_matrix)))

  # biomass each cell sends toward each one of its four potential neighbours
  # (survival = 1, so no mortality term; equal split across directions)
  flux_per_edge <- biomass_matrix * (d_grid / 4)

  # total biomass leaving a cell; edge/corner cells have fewer neighbours so
  # they lose less (boundary_number is 2, 3, or 4 depending on position)
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
