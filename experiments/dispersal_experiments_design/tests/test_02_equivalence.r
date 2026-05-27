# test_02_equivalence.r
# Test 2: diffuse() and diffuse_density_dependent() should produce identical
# invasion front trajectories when b = 0 and max_disp_rate = 2 * disp_rate.
#
# When b = 0, the sigmoid collapses to 0.5 regardless of metabolism or net
# growth, so d_grid = max_disp_rate / 2 everywhere — matching the constant
# disp_rate used by diffuse(). The net_growth_matrix and metabolism values
# are therefore irrelevant and set to 0.
#
# Expected: both curves overlap exactly on the plot.

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse.r")
source("experiments/dispersal_experiments_design/functions/diffuse_density_dependent.r")

# --- Parameters ---
n_row         <- 51
n_col         <- 51
n_time        <- 50
disp_rate     <- 0.4
max_disp_rate <- 2 * disp_rate  # equivalence condition: max_disp_rate / 2 = disp_rate
b             <- 0
threshold     <- 1e-6
center_r      <- 26
center_c      <- 26

# --- Initialization (identical starting state for both runs) ---
init <- matrix(0, nrow = n_row, ncol = n_col)
init[center_r, center_c] <- 1

biomass_ind <- init
biomass_dep <- init

boundary_number   <- enforce_boundary_conditions(n_row, n_col)
net_growth_matrix <- matrix(0, nrow = n_row, ncol = n_col)  # irrelevant when b = 0
metabolism        <- 0                                        # irrelevant when b = 0

dist_from_center <- sqrt((row(init) - center_r)^2 + (col(init) - center_c)^2)

# --- Simulation loops ---
front_ind <- numeric(n_time + 1)
front_dep <- numeric(n_time + 1)
front_ind[1] <- 0
front_dep[1] <- 0

for (t in seq_len(n_time)) {
  biomass_ind <- biomass_ind + diffuse(biomass_ind, disp_rate, boundary_number)
  biomass_dep <- biomass_dep + diffuse_density_dependent(
    biomass_dep, net_growth_matrix, metabolism, max_disp_rate, boundary_number, b
  )

  above_ind <- biomass_ind > threshold
  above_dep <- biomass_dep > threshold
  front_ind[t + 1] <- if (any(above_ind)) max(dist_from_center[above_ind]) else 0
  front_dep[t + 1] <- if (any(above_dep)) max(dist_from_center[above_dep]) else 0
}

# Numerical check: report max absolute difference across all time steps
max_front_diff <- max(abs(front_ind - front_dep))
message(sprintf("Max absolute difference in front distance: %.2e", max_front_diff))

# --- Plot ---
png("experiments/dispersal_experiments_design/figures/test_02_equivalence.png", width = 800, height = 600)

plot(
  0:n_time, front_ind,
  type = "b", pch = 19, cex = 0.6, col = "steelblue",
  xlab = "Time step",
  ylab = "Invasion front distance (cells)",
  main = "Test 2 — Equivalence: diffuse vs. diffuse_density_dependent (b = 0)",
  ylim = range(c(front_ind, front_dep))
)

lines(0:n_time, front_dep, type = "b", pch = 4, cex = 0.6, col = "tomato", lty = 2)

legend(
  "topleft",
  legend = c(
    "diffuse()  [disp_rate = 0.4]",
    "diffuse_density_dependent()  [b = 0, max_disp_rate = 0.8]"
  ),
  col    = c("steelblue", "tomato"),
  pch    = c(19, 4),
  lty    = c(1, 2),
  pt.cex = 0.8
)

dev.off()

message("Test 2 complete — figure saved to figures/test_02_equivalence.png")
