# Calculates daily wet biomass loss (g/day) for a population due to
# field metabolic rate, using the Blyth et al. (2026) allometric equation.
#
# Energy conversion: 1 kg wet mass = 7,000,000 J  →  7,000 J/g
# FMR (W/g) * 86,400 s/day / 7,000 J/g = g wet mass consumed per g per day
#
# To apply across a data frame of mixed endo/ectotherm species, use mapply.
# Set temp_C = NA for endotherms and endotherm_group = NA for ectotherms:
#
#   species$delta_biomass_g <- with(species, mapply(
#     calculate_biomass_metabolism,
#     initial_biomass_g = initial_biomass_g,
#     mass_g               = mass_g,
#     temp_C               = temp_C,
#     group                = group,
#     endotherm_group      = endotherm_group
#   ))

source("experiments/experiment_diversity_metabolism/calculate_metabolism.r")

######################################
# function parameters:
# initial_biomass_g = initial wet biomass of the population (g)
# mass_g, temp_C, group, endotherm_group = see calculate_metabolism()
#   use temp_C = NA for endotherms; endotherm_group = NA for ectotherms
#output: change in biomass of tthe population due to metabolism (deltaB)

calculate_biomass_metabolism <- function(initial_biomass_g, mass_g,
                                         temp_C = NULL,
                                         group = c("ectotherm", "endotherm"),
                                         endotherm_group = NULL) {
  # mapply passes NA from data frame columns rather than NULL; convert here
  # so calculate_metabolism() is.null() checks work correctly
  if (!is.null(temp_C) && length(temp_C) == 1 && is.na(temp_C)) {
    temp_C <- NULL
  }
  if (!is.null(endotherm_group) &&
      length(endotherm_group) == 1 && is.na(endotherm_group)) {
    endotherm_group <- NULL
  }

  fmr_w_per_g <- calculate_metabolism(mass_g, temp_C, group, endotherm_group)

  # FMR (W/g) * 86400 s/day / 7000 J/g = g wet biomass consumed per day
  (initial_biomass_g * fmr_w_per_g * 86400) / 7000
}
