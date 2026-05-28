# Extract and compute Madingley vegetation parameters for Sherbrooke, Canada
# Location: lon = -71.9, lat = 45.4
# Climate data extracted from MadingleyR spatial rasters using GDAL.
# Vegetation model parameters follow Harfoot et al. 2014 (TerrestrialCarbon.cpp).

# ---- Madingley parameter constants (from TerrestrialCarbon.cpp) ----
mad <- list(
  MaxNPP          = 0.961644704,
  T1NPP           = 0.237468183,
  T2NPP           = 0.100597089,
  PNPP            = 0.001184101,
  FracStructScalar = 7.154615419,
  MaxFracStruct    = 0.362742634,
  A_FracEver       = 1.270782192,
  B_FracEver       = -1.828591558,
  C_FracEver       = 0.844864063,
  MEG_LeafMort     = 0.040273936,
  CEG_LeafMort     = 1.013070062,
  MD_LeafMort      = 0.020575964,
  CD_LeafMort      = -1.195235464,
  MFRoot_Mort      = 0.04309283,
  CFRoot_Mort      = -1.478393163,
  ErMin = 0.01,  ErMax = 24.0,
  DrMin = 0.01,  DrMax = 24.0,
  FrmMin = 0.01, FrmMax = 12.0,
  P2StMort = 0.139462774,
  P1StMort = -4.395910091,
  StmMin = 0.001, StmMax = 1.0,
  MassCperDry   = 0.476,
  MassDryperWet = 0.213,
  # Fire disturbance (TerrestrialCarbon.cpp)
  NPPScalarFire   = 8.419,
  NPPHalfSatFire  = 1.149,
  LFSScalarFire   = 19.984,
  LFSHalfSatFire  = 0.388,
  BaseScalarFire  = 2.0,
  MinFireRate     = 2.26e-6   # minimum annual disturbance rate (yr^-1)
)

# ---- Helper functions mirroring TerrestrialCarbon.cpp ----

calc_miami_npp <- function(temp, precip) {
  npp_t <- mad$MaxNPP / (1 + exp(mad$T1NPP - mad$T2NPP * temp))
  npp_p <- mad$MaxNPP * (1 - exp(-mad$PNPP * precip))
  min(npp_t, npp_p)
}

calc_frac_struct <- function(npp) {
  min_fs <- 0.01
  fs <- min_fs * (exp(mad$FracStructScalar * npp) /
        (1 + min_fs * (exp(mad$FracStructScalar * npp) - 1)))
  fs <- min(fs, 1 - min_fs)
  fs * mad$MaxFracStruct
}

calc_frac_struct_herb <- function(npp) {
  # Herbs have negligible structural tissue (no woody stems)
  min_fs <- 0.01
  max_fs_herb <- 0.05  # much lower structural allocation than trees
  fs <- min_fs * (exp(mad$FracStructScalar * npp) /
        (1 + min_fs * (exp(mad$FracStructScalar * npp) - 1)))
  fs <- min(fs, 1 - min_fs)
  fs * max_fs_herb
}

calc_frac_evergreen <- function(ndf) {
  frac <- mad$A_FracEver * ndf^2 + mad$B_FracEver * ndf + mad$C_FracEver
  max(0, min(1, frac))
}

calc_eg_leaf_mort <- function(temp) {
  r <- exp(mad$MEG_LeafMort * temp - mad$CEG_LeafMort)
  max(mad$ErMin, min(mad$ErMax, r))
}

calc_decid_leaf_mort <- function(temp) {
  r <- exp(-(mad$MD_LeafMort * temp + mad$CD_LeafMort))
  max(mad$DrMin, min(mad$DrMax, r))
}

calc_froot_mort <- function(temp) {
  r <- exp(mad$MFRoot_Mort * temp + mad$CFRoot_Mort)
  max(mad$FrmMin, min(mad$FrmMax, r))
}

calc_leaf_frac_alloc <- function(leaf_mort, froot_mort) {
  leaf_mort / (leaf_mort + froot_mort)
}

calc_struct_mort <- function(aet) {
  r <- exp(mad$P2StMort * aet / 1000 + mad$P1StMort)
  max(mad$StmMin, min(mad$StmMax, r))
}

# Annual fire-driven disturbance rate (yr^-1), from Madingley TerrestrialCarbon.cpp.
# npp : annual NPP (kg C/m2/yr); lfs : fraction of year with dry soil (0-1).
calc_fire_rate <- function(npp, lfs) {
  f_npp <- 1 / (1 + exp(-mad$NPPScalarFire  * (npp - mad$NPPHalfSatFire)))
  f_lfs <- 1 / (1 + exp(-mad$LFSScalarFire  * (lfs - mad$LFSHalfSatFire)))
  rate  <- mad$BaseScalarFire * f_npp * f_lfs
  max(mad$MinFireRate, min(1.0, rate))
}

# ---- Sherbrooke climate (extracted from MadingleyR rasters) ----

clim <- list(
  temps   = c(-11.73, -10.35, -4.62, 3.13, 10.46, 15.50, 17.99,
               16.88,  12.31,  6.10, -0.35, -8.02),  # monthly °C
  precips = c(92.36, 77.94, 88.51, 86.67, 97.06, 123.32, 124.97,
              139.18, 110.56, 109.04, 105.17, 101.65),  # monthly mm
  # Raw monthly NPP from raster (kg C/m2/month, used for seasonality)
  npp_raw = c(0, 0.003937, 0.122047, 1.244094, 2.070866, 4.196850,
              3.488189, 2.395669, 2.779528, 0.505906, 0.240158, 0),
  # Monthly AET (mm) from Madingley soil-water balance simulation
  aet_monthly = c(10.43, 11.65, 18.63, 32.26, 50.93, 69.27, 79.43,
                  75.21, 59.07, 40.52, 25.49, 15.06),
  ndf = 0.5168,  # fraction year frost (NDF)
  lfs = 0.0000   # fraction year fire (LFS)
)

# ---- Derive parameters ----

mean_temp   <- mean(clim$temps)
total_precip <- sum(clim$precips)
total_aet    <- sum(clim$aet_monthly)

# Annual NPP via Miami model (kg C/m2/year)
annual_npp <- calc_miami_npp(mean_temp, total_precip)

# Monthly NPP: annual NPP distributed by seasonality from raster
npp_pos      <- pmax(0, clim$npp_raw)
seasonality  <- npp_pos / sum(npp_pos)
npp_monthly  <- annual_npp * seasonality  # kg C/m2/month (sums to annual_npp)

# Fractional allocation to structural tissue (uses annual NPP)
f_struct_ever  <- calc_frac_struct(annual_npp)
f_struct_decid <- calc_frac_struct(annual_npp)  # same formula, same annual NPP
f_struct_herb  <- calc_frac_struct_herb(annual_npp)

# Annual leaf mortality rates (yr^-1)
mu_ever  <- calc_eg_leaf_mort(mean_temp)
mu_decid <- calc_decid_leaf_mort(mean_temp)
mu_herb  <- calc_decid_leaf_mort(mean_temp)  # herbs: same seasonal turnover as deciduous

# Fine root mortality (yr^-1)
froot_mort <- calc_froot_mort(mean_temp)

# Fractional leaf allocation (fraction of non-structural NPP going to leaves)
f_leafmort_ever  <- calc_leaf_frac_alloc(mu_ever,  froot_mort)
f_leafmort_decid <- calc_leaf_frac_alloc(mu_decid, froot_mort)
f_leafmort_herb  <- calc_leaf_frac_alloc(mu_herb,  froot_mort)

# Structural mortality (yr^-1)
struct_mort <- calc_struct_mort(total_aet)

# Fraction year evergreen (used for reference, not directly in our competition formulation)
frac_ever <- calc_frac_evergreen(clim$ndf)

# Cell area for 1-degree cell at lat 45.4° (m^2)
A_cell <- cos(45.4 * pi / 180) * (111000)^2  # ~8.65e9 m^2

# ---- Build params list ----

get_sherbrooke_params <- function(
  alpha        = 5.0,   # Beer-Lambert shading coefficient (m2 / kg C)
  beta_ever    = 0.01,  # annual herbivory rate on evergreen (yr^-1)
  beta_decid   = 0.01,  # annual herbivory rate on deciduous (yr^-1)
  beta_herb    = 0.05,  # annual herbivory rate on herbs (yr^-1)
  k_fire_ever  = 1.5,   # fire susceptibility multiplier for evergreen (> 1 = more fire-prone)
  k_fire_decid = 0.67,  # fire susceptibility multiplier for deciduous (< 1 = less fire-prone)
  a_max        = 600L   # maximum patch age (months, ~50 years)
) {
  # Climate-derived base fire rate, then scaled by type
  e_base    <- calc_fire_rate(annual_npp, clim$lfs)
  e_ever_yr  <- min(1.0, e_base * k_fire_ever)
  e_decid_yr <- min(1.0, e_base * k_fire_decid)

  list(
    # Monthly NPP (kg C/m2/month)
    npp_monthly     = npp_monthly,

    # Structural allocation (dimensionless fraction)
    f_struct_ever   = f_struct_ever,
    f_struct_decid  = f_struct_decid,
    f_struct_herb   = f_struct_herb,

    # Fractional leaf allocation (dimensionless)
    f_leafmort_ever  = f_leafmort_ever,
    f_leafmort_decid = f_leafmort_decid,
    f_leafmort_herb  = f_leafmort_herb,

    # Annual leaf mortality rates (yr^-1), divided by 12 in model for monthly step
    mu_ever  = mu_ever,
    mu_decid = mu_decid,
    mu_herb  = mu_herb,

    # Competition coefficient
    alpha = alpha,

    # Herbivory rates (yr^-1), converted to monthly in model
    beta_ever  = beta_ever,
    beta_decid = beta_decid,
    beta_herb  = beta_herb,

    # Monthly disturbance rates (convert from annual)
    e_monthly = c(1 - (1 - e_ever_yr)^(1/12),
                  1 - (1 - e_decid_yr)^(1/12)),

    # Cell area
    A_cell = A_cell,

    # Maximum patch age (months)
    a_max = a_max,

    # Diagnostics
    annual_npp   = annual_npp,
    frac_ever    = frac_ever,
    struct_mort  = struct_mort,
    mean_temp    = mean_temp,
    total_precip = total_precip,
    total_aet    = total_aet,
    e_ever_yr    = e_ever_yr,
    e_decid_yr   = e_decid_yr
  )
}

# ---- Default initial conditions ----

make_initial_conditions <- function(params) {
  a_max <- params$a_max

  # Occupancy: distribute uniformly across ages and tree types, proportional to frac_ever
  p_init <- matrix(0, 2, a_max)
  # Steady-state approximation: exponential age distribution
  ages   <- 1:a_max
  e      <- params$e_monthly
  p_init[1, ] <- params$frac_ever     * exp(-e[1] * ages) * e[1] / (1 - exp(-e[1]))
  p_init[2, ] <- (1 - params$frac_ever) * exp(-e[2] * ages) * e[2] / (1 - exp(-e[2]))

  # Normalize to sum to 1
  p_init <- p_init / sum(p_init)

  # Biomass: small positive values proportional to age (patch age as proxy for succession)
  b_tree_init <- matrix(0, 2, a_max)
  b_herb_init <- matrix(0, 2, a_max)

  # Equilibrium tree biomass approximation (ages linearly toward a carrying capacity)
  b_max_tree <- 0.3  # kg C/m2 (rough equilibrium)
  b_max_herb <- 0.05 # kg C/m2
  age_half   <- 120  # months to half-saturation (~10 years)

  for (a in 1:a_max) {
    frac_mature <- a / (a + age_half)
    b_tree_init[1, a] <- b_max_tree * frac_mature
    b_tree_init[2, a] <- b_max_tree * frac_mature
    b_herb_init[1, a] <- b_max_herb * (1 - frac_mature)
    b_herb_init[2, a] <- b_max_herb * (1 - frac_mature)
  }

  list(p = p_init, b_tree = b_tree_init, b_herb = b_herb_init)
}
