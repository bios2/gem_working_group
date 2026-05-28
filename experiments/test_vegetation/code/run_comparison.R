# Comparative simulation for Phoenix, Iqaluit, Kelowna, and Sherbrooke
source("code/extract_parameters.R")
source("code/vegetation_model.R")

# ---- Generic parameter function ----
# Takes pre-computed climate scalars and returns a params list
# (same structure as get_sherbrooke_params, but location-agnostic)

get_params <- function(temps, precips, npp_raw, aet_monthly, frac_frost, lfs,
                       lat,
                       alpha        = 5.0,
                       beta_ever    = 0.01,
                       beta_decid   = 0.01,
                       beta_herb    = 0.05,
                       k_fire_ever  = 1.5,   # evergreen more fire-prone
                       k_fire_decid = 0.67,  # deciduous less fire-prone
                       a_max        = 600L) {

  mean_temp    <- mean(temps)
  total_precip <- sum(precips)
  total_aet    <- sum(aet_monthly)

  annual_npp   <- calc_miami_npp(mean_temp, total_precip)

  npp_pos     <- pmax(0, npp_raw)
  seasonality <- if (sum(npp_pos) > 0) npp_pos / sum(npp_pos) else rep(1/12, 12)
  npp_monthly <- annual_npp * seasonality

  f_struct_ev  <- calc_frac_struct(annual_npp)
  f_struct_dc  <- calc_frac_struct(annual_npp)
  f_struct_hb  <- calc_frac_struct_herb(annual_npp)

  mu_ev  <- calc_eg_leaf_mort(mean_temp)
  mu_dc  <- calc_decid_leaf_mort(mean_temp)
  mu_hb  <- calc_decid_leaf_mort(mean_temp)
  froot  <- calc_froot_mort(mean_temp)

  f_lm_ev <- calc_leaf_frac_alloc(mu_ev, froot)
  f_lm_dc <- calc_leaf_frac_alloc(mu_dc, froot)
  f_lm_hb <- calc_leaf_frac_alloc(mu_hb, froot)

  frac_ev <- calc_frac_evergreen(frac_frost)

  A_cell <- cos(lat * pi / 180) * (111000)^2

  # Climate-derived fire rates; evergreen patches burn more readily than deciduous
  e_base    <- calc_fire_rate(annual_npp, lfs)
  e_ever_yr  <- min(1.0, e_base * k_fire_ever)
  e_decid_yr <- min(1.0, e_base * k_fire_decid)

  list(
    npp_monthly     = npp_monthly,
    f_struct_ever   = f_struct_ev,
    f_struct_decid  = f_struct_dc,
    f_struct_herb   = f_struct_hb,
    f_leafmort_ever  = f_lm_ev,
    f_leafmort_decid = f_lm_dc,
    f_leafmort_herb  = f_lm_hb,
    mu_ever  = mu_ev,
    mu_decid = mu_dc,
    mu_herb  = mu_hb,
    alpha    = alpha,
    beta_ever = beta_ever, beta_decid = beta_decid, beta_herb = beta_herb,
    e_monthly = c(1 - (1 - e_ever_yr)^(1/12),
                  1 - (1 - e_decid_yr)^(1/12)),
    A_cell    = A_cell,
    a_max     = a_max,
    annual_npp   = annual_npp,
    frac_ever    = frac_ev,
    mean_temp    = mean_temp,
    total_precip = total_precip,
    total_aet    = total_aet,
    e_ever_yr    = e_ever_yr,
    e_decid_yr   = e_decid_yr
  )
}

# ---- Climate data (extracted from MadingleyR rasters via Python/GDAL) ----

clim_data <- list(

  Sherbrooke = list(
    lat   = 45.4,
    temps   = c(-11.73,-10.35,-4.62,3.13,10.46,15.50,17.99,16.88,12.31,6.10,-0.35,-8.02),
    precips = c(92.36,77.94,88.51,86.67,97.06,123.32,124.97,139.18,110.56,109.04,105.17,101.65),
    npp_raw = c(0,0.003937,0.122047,1.244094,2.070866,4.196850,3.488189,2.395669,2.779528,0.505906,0.240158,0),
    aet_monthly = c(10.43,11.65,18.63,32.26,50.93,69.27,79.43,75.21,59.07,40.52,25.49,15.06),
    frac_frost = 0.5168, lfs = 0.0
  ),

  Phoenix = list(
    lat   = 33.45,
    temps   = c(10.13,12.28,14.94,18.87,23.86,28.90,32.40,31.14,27.86,21.46,14.65,10.25),
    precips = c(30.16,29.50,31.77,10.54,3.84,3.56,31.29,41.22,23.82,19.49,22.08,29.81),
    npp_raw = c(0.683071,1.125984,1.421260,0.860236,0.240157,0,0,0.092520,0,0.387795,0.358268,0.240157),
    aet_monthly = c(18.05,19.63,22.81,13.18,5.51,3.95,25.10,30.44,20.04,16.88,16.99,18.71),
    frac_frost = 0.090, lfs = 1.0
  ),

  Iqaluit = list(
    lat   = 63.75,
    temps   = c(-25.70,-26.62,-23.21,-14.91,-4.67,2.29,6.73,5.80,1.31,-5.30,-12.85,-20.84),
    precips = c(20.79,20.55,20.76,26.81,26.08,39.11,55.40,63.88,49.51,41.29,32.95,22.32),
    npp_raw = c(0,0,0,0,0,0.102053,0.280356,0.193539,0.037849,0,0,0),
    aet_monthly = c(0.00,0.00,0.00,0.00,2.15,19.88,43.02,39.64,17.34,5.21,0.35,0.00),
    frac_frost = 0.817, lfs = 0.0
  ),

  Kelowna = list(
    lat   = 49.9,
    temps   = c(-7.38,-5.09,-2.68,1.86,6.71,10.21,13.45,12.92,9.04,3.53,-3.20,-6.88),
    precips = c(77.79,58.92,58.74,47.42,61.27,74.17,45.54,45.10,42.96,43.24,63.77,77.66),
    npp_raw = c(0.003937,0.092520,0.476378,1.864173,2.897638,3.842520,3.990157,3.576772,1.805118,0.712598,0.033465,0),
    aet_monthly = c(4.14,5.79,11.35,20.18,32.98,43.36,40.85,39.57,29.62,17.47,7.06,4.28),
    frac_frost = 0.651, lfs = 0.0
  )
)

# ---- Run simulations ----

n_years <- 200
n_steps <- n_years * 12

results_all <- list()

for (loc in names(clim_data)) {
  cd <- clim_data[[loc]]
  params <- get_params(
    temps       = cd$temps,
    precips     = cd$precips,
    npp_raw     = cd$npp_raw,
    aet_monthly = cd$aet_monthly,
    frac_frost  = cd$frac_frost,
    lfs         = cd$lfs,
    lat         = cd$lat
  )

  cat(sprintf("%-12s  NPP=%.3f  f_ever=%.2f  T=%.1fC  P=%.0fmm  e_ever=%.2e  e_decid=%.2e\n",
              loc, params$annual_npp, params$frac_ever,
              params$mean_temp, params$total_precip,
              params$e_ever_yr, params$e_decid_yr))

  init <- make_initial_conditions(params)

  res <- run_vegetation(init$p, init$b_tree, init$b_herb, params, n_steps)
  res$location <- loc
  res$year     <- res$step / 12

  # Annual means (kg C → Tg C)
  res_yr <- aggregate(cbind(B_ever, B_decid, B_herb) ~ floor(year),
                      data = res, FUN = mean)
  names(res_yr)[1] <- "year"
  res_yr$B_total   <- rowSums(res_yr[, c("B_ever","B_decid","B_herb")]) / 1e9
  res_yr$B_ever    <- res_yr$B_ever  / 1e9
  res_yr$B_decid   <- res_yr$B_decid / 1e9
  res_yr$B_herb    <- res_yr$B_herb  / 1e9
  res_yr$location  <- loc

  results_all[[loc]] <- res_yr
}

cat("\n")

# ---- Figure ----

all_data <- do.call(rbind, results_all)
write.csv(all_data, "results/comparison_simulation.csv", row.names = FALSE)

locations <- names(clim_data)
col_ever  <- "#2E8B57"
col_decid <- "#DAA520"
col_herb  <- "#8B4513"
col_total <- "#1C1C1C"

png("figures/comparison_biomass.png", width = 2400, height = 2000, res = 150)
par(mfrow = c(2, 2), mar = c(4.5, 5, 3.5, 1.5), oma = c(0, 0, 2, 0))

for (loc in locations) {
  d    <- results_all[[loc]]
  ymax <- max(d$B_total) * 1.15
  if (ymax == 0) ymax <- 0.01

  plot(d$year, d$B_total,
       type = "l", lwd = 2.5, col = col_total,
       ylim = c(0, ymax),
       xlab = "Année", ylab = "Biomasse régionale (Tg C)",
       main = loc, cex.main = 1.4, cex.lab = 1.2, cex.axis = 1.1)

  lines(d$year, d$B_ever,  col = col_ever,  lwd = 2)
  lines(d$year, d$B_decid, col = col_decid, lwd = 2)
  lines(d$year, d$B_herb,  col = col_herb,  lwd = 2)

  # Final-year values in legend
  last <- tail(d, 1)
  leg_txt <- c(
    sprintf("Total   %.3f", last$B_total),
    sprintf("Résineux %.3f", last$B_ever),
    sprintf("Décidus  %.3f", last$B_decid),
    sprintf("Herbes   %.3f", last$B_herb)
  )
  legend("topright", legend = leg_txt,
         col = c(col_total, col_ever, col_decid, col_herb),
         lwd = c(2.5, 2, 2, 2), bty = "n", cex = 0.95)
}

mtext("Biomasse végétale au fil du temps — 4 localités canadiennes et américaines",
      outer = TRUE, cex = 1.3, font = 2)

dev.off()
cat("Figure sauvegardée : figures/comparison_biomass.png\n")
