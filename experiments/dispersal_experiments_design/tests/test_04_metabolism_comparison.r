# test_04_metabolism_comparison.r
# Test 4: with the same resource grid, density-dependent dispersal should be
# higher for species with high metabolism than low metabolism.
#
# Scenario: net_growth = 0 everywhere, two species differ only in metabolism.
# High metabolism (hummingbird, x = 1.0): sigmoid is close to max_disp_rate
#   because (metabolism - net_growth) = 1.0 is large and positive.
# Low metabolism (black bear, x = 0.1): sigmoid gives a lower dispersal rate
#   because (metabolism - net_growth) = 0.1 is only slightly positive.
# Expected: hummingbird front lies above bear front.

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse_density_dependent.r")

# --- Parameters ---
n_row            <- 51
n_col            <- 51
n_time           <- 50
max_disp_rate    <- 0.4
metabolism_high  <- 1.0   # hummingbird
metabolism_low   <- 0.1   # black bear
b                <- 10
threshold        <- 1e-6
center_r         <- 26
center_c         <- 26

# --- Initialization ---
init <- matrix(0, nrow = n_row, ncol = n_col)
init[center_r, center_c] <- 1

biomass_high <- init
biomass_low  <- init

boundary_number   <- enforce_boundary_conditions(n_row, n_col)
net_growth_matrix <- matrix(0, nrow = n_row, ncol = n_col)

dist_from_center <- sqrt((row(init) - center_r)^2 + (col(init) - center_c)^2)

# Expected per-capita dispersal rates (for reference)
d_high <- max_disp_rate / (1 + exp(-b * (metabolism_high - 0)))
d_low  <- max_disp_rate / (1 + exp(-b * (metabolism_low  - 0)))
message(sprintf("Expected d (hummingbird, metabolism = %.1f): %.4f", metabolism_high, d_high))
message(sprintf("Expected d (black bear,  metabolism = %.1f): %.4f", metabolism_low,  d_low))

# --- Simulation loops ---
front_high <- numeric(n_time + 1)
front_low  <- numeric(n_time + 1)
front_high[1] <- 0
front_low[1]  <- 0

for (t in seq_len(n_time)) {
  biomass_high <- biomass_high + diffuse_density_dependent(
    biomass_high, net_growth_matrix, metabolism_high, max_disp_rate, boundary_number, b
  )
  biomass_low <- biomass_low + diffuse_density_dependent(
    biomass_low, net_growth_matrix, metabolism_low, max_disp_rate, boundary_number, b
  )

  above_high <- biomass_high > threshold
  above_low  <- biomass_low  > threshold
  front_high[t + 1] <- if (any(above_high)) max(dist_from_center[above_high]) else 0
  front_low[t + 1]  <- if (any(above_low))  max(dist_from_center[above_low])  else 0
}

# --- Plot ---
png("experiments/dispersal_experiments_design/figures/test_04_metabolism_comparison.png",
    width = 800, height = 600)

ylim_range <- range(c(front_high, front_low))

plot(
  0:n_time, front_high,
  type = "b", pch = 19, cex = 0.6, col = "steelblue",
  xlab = "Time step",
  ylab = "Invasion front distance (cells)",
  main = "Test 4 — Metabolism effect on density-dependent dispersal\n(net growth = 0, max_disp_rate = 0.4, b = 10)",
  ylim = ylim_range
)

lines(0:n_time, front_low, type = "b", pch = 4, cex = 0.6, col = "tomato", lty = 2)

legend(
  "topleft",
  legend = c(
    sprintf("High metabolism — hummingbird  [x = %.1f]", metabolism_high),
    sprintf("Low metabolism  — black bear   [x = %.1f]", metabolism_low)
  ),
  col    = c("steelblue", "tomato"),
  pch    = c(19, 4),
  lty    = c(1, 2),
  pt.cex = 0.8
)

dev.off()

message("Test 4 complete — figure saved to experiments/dispersal_experiments_design/figures/test_04_metabolism_comparison.png")
