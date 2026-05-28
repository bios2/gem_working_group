# Example simulation for Sherbrooke, Canada
# Runs the vegetation model and saves total regional biomass figure

source("code/sherbrooke_parameters.R")
source("code/vegetation_model.R")

# ---- Parameters and initial conditions ----

params <- get_sherbrooke_params()

cat("=== Sherbrooke parameters ===\n")
cat(sprintf("  Annual NPP (Miami):    %.4f kg C/m2/yr\n", params$annual_npp))
cat(sprintf("  Fraction evergreen:    %.3f\n",             params$frac_ever))
cat(sprintf("  Mean temperature:      %.2f C\n",           params$mean_temp))
cat(sprintf("  Total precipitation:   %.1f mm\n",          params$total_precip))
cat(sprintf("  Total AET:             %.1f mm\n",          params$total_aet))
cat(sprintf("  Mu evergreen (yr-1):   %.4f\n",             params$mu_ever))
cat(sprintf("  Mu deciduous (yr-1):   %.4f\n",             params$mu_decid))
cat(sprintf("  f_struct evergreen:    %.4f\n",             params$f_struct_ever))
cat(sprintf("  f_leafmort evergreen:  %.4f\n",             params$f_leafmort_ever))
cat(sprintf("  f_leafmort deciduous:  %.4f\n",             params$f_leafmort_decid))
cat(sprintf("  e_ever (monthly):      %.6f\n",             params$e_monthly[1]))
cat(sprintf("  e_decid (monthly):     %.6f\n",             params$e_monthly[2]))
cat("\n")

init <- make_initial_conditions(params)

# ---- Run simulation ----

n_years <- 200
n_steps <- n_years * 12

cat(sprintf("Running %d-year simulation (%d monthly steps)...\n", n_years, n_steps))
t0 <- proc.time()
result <- run_vegetation(init$p, init$b_tree, init$b_herb, params, n_steps)
cat(sprintf("Done in %.1f seconds.\n\n", (proc.time() - t0)[3]))

# Convert biomass from kg C to Tg C for the cell (1 Tg = 10^12 g = 10^9 kg)
result$B_ever_Tg  <- result$B_ever  / 1e9
result$B_decid_Tg <- result$B_decid / 1e9
result$B_herb_Tg  <- result$B_herb  / 1e9
result$year       <- result$step / 12

# ---- Save results ----

write.csv(result, "results/sherbrooke_simulation.csv", row.names = FALSE)
cat("Results saved to results/sherbrooke_simulation.csv\n")

# ---- Summary statistics (last year) ----

last_yr <- result[result$step > (n_steps - 12), ]
cat("=== Final year means ===\n")
cat(sprintf("  Evergreen biomass:  %.4f Tg C\n", mean(last_yr$B_ever_Tg)))
cat(sprintf("  Deciduous biomass:  %.4f Tg C\n", mean(last_yr$B_decid_Tg)))
cat(sprintf("  Herb biomass:       %.4f Tg C\n", mean(last_yr$B_herb_Tg)))
total <- mean(last_yr$B_ever_Tg + last_yr$B_decid_Tg + last_yr$B_herb_Tg)
cat(sprintf("  Total biomass:      %.4f Tg C\n", total))
cat("\n")

# ---- Figure ----

png("figures/sherbrooke_biomass.png", width = 1800, height = 1100, res = 150)

# Colour palette
col_ever  <- "#2E8B57"  # sea green
col_decid <- "#DAA520"  # goldenrod
col_herb  <- "#8B4513"  # saddlebrown
col_total <- "#1C1C1C"  # near-black

# Annual means for smoother line
result_yr <- aggregate(
  cbind(B_ever_Tg, B_decid_Tg, B_herb_Tg) ~ floor(year),
  data = result,
  FUN  = mean
)
names(result_yr)[1] <- "year"
result_yr$B_total_Tg <- rowSums(result_yr[, c("B_ever_Tg", "B_decid_Tg", "B_herb_Tg")])

ymax <- max(result_yr$B_total_Tg) * 1.1

par(mar = c(5, 5, 4, 2))
plot(result_yr$year, result_yr$B_total_Tg,
     type = "l", lwd = 2.5, col = col_total,
     ylim = c(0, ymax),
     xlab = "Year",
     ylab = "Total regional biomass (Tg C)",
     main = "Vegetation model — Sherbrooke, Canada (45.4°N, 71.9°W)",
     cex.axis = 1.1, cex.lab = 1.3, cex.main = 1.3)

lines(result_yr$year, result_yr$B_ever_Tg,  col = col_ever,  lwd = 2)
lines(result_yr$year, result_yr$B_decid_Tg, col = col_decid, lwd = 2)
lines(result_yr$year, result_yr$B_herb_Tg,  col = col_herb,  lwd = 2)

legend("topright",
       legend = c("Total", "Evergreen", "Deciduous", "Herbs"),
       col    = c(col_total, col_ever, col_decid, col_herb),
       lwd    = c(2.5, 2, 2, 2),
       bty    = "n", cex = 1.2)

dev.off()
cat("Figure saved to figures/sherbrooke_biomass.png\n")
