# Test script for calculate_biomass_metabolism()
# Run from project root (gem_working_group/)

source("experiments/experiment_diversity_metabolism/calculate_biomass_metabolism.r")

# Example species data frame
# mass_g          = individual body mass (assumed; adjust as needed)
# temp_C          = ambient temperature for ectotherms; NA for endotherms
# initial_biomass_g = initial wet biomass of the population (g)
species <- data.frame(
  species_name         = c("example_bird", "example_mammal", "example_lizard"),
  group                = c("endotherm", "endotherm", "ectotherm"),
  endotherm_group      = c("bird", "mammal", NA),
  mass_g               = c(500, 5000, 200),  # g: ~pigeon, ~cat, ~lizard
  temp_C               = c(NA, NA, 20),       # 20 C ambient for lizard
  initial_biomass_g = c(1000, 1000, 1000),   # 1000 kg, 10000 kg, 1000 kg
  stringsAsFactors     = FALSE
)

species$delta_biomass_g <- with(species, mapply(
  calculate_biomass_metabolism,
  initial_biomass_g = initial_biomass_g,
  mass_g               = mass_g,
  temp_C               = temp_C,
  group                = group,
  endotherm_group      = endotherm_group
))

species$next_biomass_g <- species$initial_biomass_g - species$delta_biomass_g

print(species[, c("species_name", "initial_biomass_g",
                  "delta_biomass_g", "next_biomass_g")])
