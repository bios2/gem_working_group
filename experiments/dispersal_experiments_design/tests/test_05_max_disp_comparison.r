# test_05_max_disp_comparison.r
# Test 5: all else being equal, species with a higher maximum dispersal rate
# should disperse further than species with a lower maximum dispersal rate.
#
# Scenario: net_growth = 0, metabolism = 0.5, b = 10. Two species differ only
# in max_disp_rate (0.8 vs 0.2). The sigmoid output is the same fraction for
# both, so the absolute per-capita dispersal rate scales directly with max_disp_rate.
# Expected: high max_disp_rate front lies above low max_disp_rate front.

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse_density_dependent.r")

# --- Parameters ---
n_row             <- 51
n_col             <- 51
n_time            <- 50
max_disp_rate_high <- 0.8
max_disp_rate_low  <- 0.2
metabolism         <- 0.5
b                  <- 10
threshold          <- 1e-6
center_r           <- 26
center_c           <- 26

# --- Initialization ---
init <- matrix(0, nrow = n_row, ncol = n_col)
init[center_r, center_c] <- 1

biomass_high <- init
biomass_low  <- init

boundary_number   <- enforce_boundary_conditions(n_row, n_col)
net_growth_matrix <- matrix(0, nrow = n_row, ncol = n_col)

dist_from_center <- sqrt((row(init) - center_r)^2 + (col(init) - center_c)^2)

# Expected per-capita dispersal rates (for reference)
sigmoid_val <- 1 / (1 + exp(-b * (metabolism - 0)))
message(sprintf("Sigmoid value (metabolism = %.1f, net_growth = 0): %.4f", metabolism, sigmoid_val))
message(sprintf("Expected d (high max_disp_rate = %.1f): %.4f", max_disp_rate_high, max_disp_rate_high * sigmoid_val))
message(sprintf("Expected d (low  max_disp_rate = %.1f): %.4f", max_disp_rate_low,  max_disp_rate_low  * sigmoid_val))

# --- Simulation loops ---
front_high <- numeric(n_time + 1)
front_low  <- numeric(n_time + 1)
front_high[1] <- 0
front_low[1]  <- 0

for (t in seq_len(n_time)) {
  biomass_high <- biomass_high + diffuse_density_dependent(
    biomass_high, net_growth_matrix, metabolism, max_disp_rate_high, boundary_number, b
  )
  biomass_low <- biomass_low + diffuse_density_dependent(
    biomass_low, net_growth_matrix, metabolism, max_disp_rate_low, boundary_number, b
  )

  above_high <- biomass_high > threshold
  above_low  <- biomass_low  > threshold
  front_high[t + 1] <- if (any(above_high)) max(dist_from_center[above_high]) else 0
  front_low[t + 1]  <- if (any(above_low))  max(dist_from_center[above_low])  else 0
}

# --- Plot ---
png("experiments/dispersal_experiments_design/figures/test_05_max_disp_comparison.png",
    width = 800, height = 600)

ylim_range <- range(c(front_high, front_low))

plot(
  0:n_time, front_high,
  type = "b", pch = 19, cex = 0.6, col = "steelblue",
  xlab = "Time step",
  ylab = "Invasion front distance (cells)",
  main = "Test 5 — Maximum dispersal rate effect\n(net growth = 0, metabolism = 0.5, b = 10)",
  ylim = ylim_range
)

lines(0:n_time, front_low, type = "b", pch = 4, cex = 0.6, col = "tomato", lty = 2)

legend(
  "topleft",
  legend = c(
    sprintf("High max dispersal  [max_disp_rate = %.1f]", max_disp_rate_high),
    sprintf("Low  max dispersal  [max_disp_rate = %.1f]", max_disp_rate_low)
  ),
  col    = c("steelblue", "tomato"),
  pch    = c(19, 4),
  lty    = c(1, 2),
  pt.cex = 0.8
)

dev.off()

message("Test 5 complete — figure saved to experiments/dispersal_experiments_design/figures/test_05_max_disp_comparison.png")
