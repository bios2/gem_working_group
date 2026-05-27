# dispersal_metabolism.r
# Advances a species' biomass grid by one dispersal step, with density-dependent
# emigration driven by metabolic demand (Ryser et al. 2021, Eq. 10).
#
# Computes the per-capita metabolic rate from species traits (one scalar per
# species), then derives local net growth at each grid cell, and delegates to
# diffuse_density_dependent(). All matrix operations are vectorised; no loops
# are needed, so the function scales to arbitrarily large grids.
#
# Arguments:
#   biomass_matrix         - n_row x n_col numeric matrix of biomass at time t (g)
#   delta_biomass_external - n_row x n_col numeric matrix of biomass change from
#                            feeding and predation (ATN dynamics), excluding
#                            metabolism, at time t+1 (g/day)
#   alpha                  - scalar; maximum per-capita dispersal rate
#   mass_g                 - scalar; individual body mass (g)
#   temp_C                 - scalar; ambient temperature (°C); NULL or NA for endotherms
#   group                  - "ectotherm" or "endotherm"
#   endotherm_group        - "bird" or "mammal"; NULL or NA for ectotherms
#   b                      - sigmoid steepness (default 10, as in Ryser et al. 2021)
#
# Returns an n_row x n_col numeric matrix of net dispersal flux (change in biomass
# due to dispersal only); add to biomass_matrix to advance one time step.

source("experiments/experiment_diversity_metabolism/calculate_metabolism.r")
source(
  "experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r"
)
source(
  "experiments/dispersal_experiments_design/functions/diffuse_density_dependent.r"
)

dispersal_from_growth_metabolism <- function(biomass_matrix, delta_biomass_external,
                                 alpha, mass_g,
                                 temp_C = NULL,
                                 group = c("ectotherm", "endotherm"),
                                 endotherm_group = NULL,
                                 b = 10) {
  # Metabolic rate in W/g (scalar, species-level); converted to fractional
  # biomass consumed per day: W/g * 86400 s/day / 7000 J/g = g/g/day
  fmr_w_per_g     <- calculate_metabolism(mass_g, temp_C, group, endotherm_group)
  metabolic_rate <- fmr_w_per_g * 86400 / 7000

  boundary_number   <- enforce_boundary_conditions(
    nrow(biomass_matrix), ncol(biomass_matrix)
  )

  # Per-capita net growth at each cell (1/day):
  #   (external flux - metabolic loss) / local biomass
  # Vectorised over the matrix; metabolism_rate is a scalar so no loop needed.
  net_growth_matrix <- delta_biomass_external / biomass_matrix

  diffuse_density_dependent(
    biomass_matrix    = biomass_matrix,
    net_growth_matrix = net_growth_matrix,
    metabolism        = metabolic_rate,
    max_disp_rate     = alpha,
    boundary_number   = boundary_number,
    b                 = b
  )
}
