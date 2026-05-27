# test_01_sublinear_dispersal.r
# Test 1: density-independent dispersal should produce a sublinear relationship
# between invasion front distance and time.
#
# Setup: biomass = 1 at the central cell of a 51x51 grid, 0 elsewhere.
# At each time step, the invasion front is defined as the farthest cell from
# the center with biomass above the threshold (1e-6).
# Expected: distance grows as sqrt(t), i.e., concave-down curve.

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse.r")

# --- Parameters ---
n_row      <- 51
n_col      <- 51
n_time     <- 50
disp_rate  <- 0.4
threshold  <- 1e-6
center_r   <- 26
center_c   <- 26

# --- Initialization ---
biomass <- matrix(0, nrow = n_row, ncol = n_col)
biomass[center_r, center_c] <- 1

boundary_number <- enforce_boundary_conditions(n_row, n_col)

# Euclidean distance from each cell to the center cell
dist_from_center <- sqrt((row(biomass) - center_r)^2 + (col(biomass) - center_c)^2)

# --- Simulation loop ---
front_distance <- numeric(n_time + 1)
front_distance[1] <- 0  # only the center cell is occupied at t = 0

for (t in seq_len(n_time)) {
  biomass <- biomass + diffuse(biomass, disp_rate, boundary_number)
  above   <- biomass > threshold
  front_distance[t + 1] <- if (any(above)) max(dist_from_center[above]) else 0
}

# --- Plot ---
png("experiments/dispersal_experiments_design/figures/test_01_sublinear_dispersal.png", width = 800, height = 600)

plot(
  0:n_time, front_distance,
  type = "b", pch = 19, cex = 0.5,
  xlab = "Time step",
  ylab = "Invasion front distance (cells)",
  main = "Test 1 — Density-independent dispersal: front distance vs. time"
)

# Reference: sqrt(t) curve scaled to the data to guide the eye
t_seq   <- seq(0, n_time, length.out = 200)
scale   <- front_distance[n_time + 1] / sqrt(n_time)
lines(t_seq, scale * sqrt(t_seq), col = "tomato", lty = 2, lwd = 1.5)

legend(
  "topleft",
  legend = c("observed front", expression(paste("reference: ", sqrt(t)))),
  col    = c("black", "tomato"),
  pch    = c(19, NA),
  lty    = c(1, 2),
  pt.cex = 0.8
)

dev.off()

message("Test 1 complete — figure saved to figures/test_01_sublinear_dispersal.png")
