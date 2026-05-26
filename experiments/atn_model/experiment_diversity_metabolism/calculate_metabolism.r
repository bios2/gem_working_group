# Functions for calculating mass-specific metabolism of endotherms and ectotherms
#
# Source: Blyth et al. (2026). The Critical Role of Coefficients: Updating Allometric
#   Normalisation Constants for Modern Ecology and Modelling. Ecology Letters, 29, e70330.
#   https://doi.org/10.1111/ele.70330
#
# Metabolic rate scales with body mass and temperature via the allometric equation:
#   X [W]      = exp(C) * M^b * exp(-Ea / kT)          (whole-body)
#   x [W/g]    = exp(C) * M^(b-1) * exp(-Ea / kT)      (mass-specific)
#
# C is the PGLS intercept (ln W) from Blyth et al. Table 1:
#   ectothermic vertebrates: C = 17.4  [95% CRI: 16.46–18.28]
#   endothermic vertebrates: C = 19.53 [95% CRI: 18.93–20.13]
# M is body mass in grams, b is the allometric exponent (taken from Table S3 in Blyth et al.):
# Ea = 0.63 eV (activation energy), k = 8.617e-5 eV/K (Boltzmann constant),
# T is body temperature in Kelvin for endotherms or environmental temperature for ectotherms.
# Multiply resting rates by 3 to obtain field metabolic rate (FMR; after Nagy 1987 and Brose et al. 2008)


calculate_metabolism <- function(mass_g, temp_C, group = c("ectotherm", "endotherm"), fmr = FALSE ) {
  # Convert temperature to Kelvin
  temp_K <- temp_C + 273.15
  
  # Set coefficients based on group
  if (group == "ectotherm") {
    C <- 17.4
    b <- 0.84
  } else if (group == "endotherm") {
    C <- 19.53
    b <- 0.73
  } else {
    stop("Group must be either 'ectotherm' or 'endotherm'")
  }
  
  # Calculate mass-specific metabolic rate (W/g)
  x <- exp(C) * mass_g^(b - 1) * exp(-0.63 / (8.617e-5 * temp_K))
  
  # If FMR is requested, multiply by 3
  if (fmr) {
    x <- x * 3
  }
  
  return(x)
}

