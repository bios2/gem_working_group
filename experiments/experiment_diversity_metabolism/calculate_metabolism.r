# Functions for calculating mass-specific metabolism of endotherms and ectotherms
#
# Source: Blyth et al. (2026). The Critical Role of Coefficients: Updating Allometric
#   Normalisation Constants for Modern Ecology and Modelling. Ecology Letters, 29, e70330.
#   https://doi.org/10.1111/ele.70330
#Gillooly, J. F., Gomez, J. P., & Mavrodiev, E. V. (2017). 
#A broad-scale comparison of aerobic activity levels in vertebrates: Endotherms versus ectotherms.
# Proceedings of the Royal Society B: Biological Sciences, 284(1849), 20162328.
# https://doi.org/10.1098/rspb.2016.2328

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
# Multiply resting rates by 3 to obtain field metabolic rate (FMR; after Brose et al. 2008)

######################################
#function specific parameters: 
#mass_g = average body mass of individual of species (g)
#temp_C = body temperature in Celsius/ambient temperature (accepted for ectotherms only)
#group = "ectotherm" or "endotherm"
#endotherm_group = "bird" or "mammal", sets body temperature (only used if group = "endotherm")
  #temperatures sourced from Gillooly et al. (2017; 10.1098/rspb.2016.2328)


calculate_metabolism <- function(mass_g, temp_C = NULL,
                                 group = c("ectotherm", "endotherm"),
                                 endotherm_group = NULL) {
  # Set coefficients and temperature based on group
  if (group == "ectotherm") {
    if (is.null(temp_C)) stop("temp_C must be provided for ectotherms")
    c_int  <- 17.4
    b      <- 0.84
    temp_k <- temp_C + 273.15
  } else if (group == "endotherm") {
    if (is.null(endotherm_group)) {
      stop("endotherm_group must be 'bird' or 'mammal'")
    }
    c_int  <- 19.53
    b      <- 0.73
    temp_k <- switch(endotherm_group,
      bird   = 41.5 + 273.15,
      mammal = 36.5 + 273.15,
      stop("endotherm_group must be 'bird' or 'mammal'")
    )
  } else {
    stop("group must be 'ectotherm' or 'endotherm'")
  }

  # Calculate mass-specific metabolic rate (W/g)
  x <- exp(c_int) * mass_g^(b - 1) * exp(-0.63 / (8.617e-5 * temp_k))
  #convert to field metabolic rate
  x <- x * 3


  return(x)
}
