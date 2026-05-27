# test_03_density_dependent_faster.r
# Test 3: density-dependent dispersal should produce faster initial spread than
# density-independent dispersal under the same conditions.
#
# Scenario: net_growth = 0 everywhere (world at carrying capacity), metabolism > 0.
# This places every cell in "bad" conditions (net growth < metabolic demand), driving
# the sigmoid close to max_disp_rate. The density-independent run uses disp_rate =
# max_disp_rate / 2 as the neutral baseline (the rate both functions share at b = 0).
# Expected: density-dependent front lies above density-independent front.

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse.r")
source("experiments/dispersal_experiments_design/functions/diffuse_density_dependent.r")

# --- Parameters ---
n_row         <- 51
n_col         <- 51
n_time        <- 50
max_disp_rate <- 0.4
disp_rate     <- max_disp_rate / 2  # neutral baseline for density-independent
metabolism    <- 0.5                 # inflection point; net_growth = 0 < metabolism → bad conditions
b             <- 10
threshold     <- 1e-6
center_r      <- 26
center_c      <- 26

# --- Initialization ---
init <- matrix(0, nrow = n_row, ncol = n_col)
init[center_r, center_c] <- 1

biomass_ind <- init
biomass_dep <- init

boundary_number   <- enforce_boundary_conditions(n_row, n_col)
net_growth_matrix <- matrix(0, nrow = n_row, ncol = n_col)  # net growth = 0 everywhere

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

# --- Plot ---
png("experiments/dispersal_experiments_design/figures/test_03_density_dependent_faster.png",
    width = 800, height = 600)

ylim_range <- range(c(front_ind, front_dep))

plot(
  0:n_time, front_dep,
  type = "b", pch = 19, cex = 0.6, col = "steelblue",
  xlab = "Time step",
  ylab = "Invasion front distance (cells)",
  main = "Test 3 — Density-dependent vs. density-independent dispersal\n(net growth = 0, metabolism = 0.5)",
  ylim = ylim_range
)

lines(0:n_time, front_ind, type = "b", pch = 4, cex = 0.6, col = "tomato", lty = 2)

legend(
  "topleft",
  legend = c(
    sprintf("diffuse_density_dependent()  [max_disp_rate = %.1f, b = %g, metabolism = %.1f]",
            max_disp_rate, b, metabolism),
    sprintf("diffuse()  [disp_rate = %.1f = max_disp_rate / 2]", disp_rate)
  ),
  col    = c("steelblue", "tomato"),
  pch    = c(19, 4),
  lty    = c(1, 2),
  pt.cex = 0.8
)

dev.off()

message("Test 3 complete — figure saved to experiments/dispersal_experiments_design/figures/test_03_density_dependent_faster.png")
