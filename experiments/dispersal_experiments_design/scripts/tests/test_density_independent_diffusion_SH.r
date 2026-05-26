getwd()
# test_density_independent_diffusion_SH.r
#
# Tests diffuse_step_vectorised on a 51 x 51 grid with all biomass placed at
# the center cell. Outputs a time series of the maximum distance from the
# center reached by any non-zero biomass -- tracking the leading dispersal
# front over time.
#
# Expected behaviour: maximum distance starts at 0 and increases over time
# as the dispersal front propagates outward from the center.
#
# Run this script with the working directory set to:
#   experiments/dispersal_experiments_design/

# --- Source functions ---------------------------------------------------------

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse_step_vectorised.r")
source("experiments/dispersal_experiments_design/functions/run_diffusion_simulation.r")

# --- Parameters ---------------------------------------------------------------

n_row     <- 51      # odd so the center cell is exact
n_col     <- 51
disp_rate <- 0.1     # fraction of biomass dispersing per step
n_time    <- 300     # number of time steps

total_biomass <- 1000   # arbitrary biomass units placed at center

# --- Initial conditions -------------------------------------------------------

# place all biomass at the single center cell; all other cells start at 0
center_row <- ceiling(n_row / 2)   # row 26
center_col <- ceiling(n_col / 2)   # col 26

init_pop             <- numeric(n_row * n_col)
center_idx           <- (center_row - 1) * n_col + center_col  # row-major index
init_pop[center_idx] <- total_biomass

# --- Run simulation -----------------------------------------------------------

results <- run_diffusion_simulation(init_pop  = init_pop,
                                    n_row     = n_row,
                                    n_col     = n_col,
                                    disp_rate = disp_rate,
                                    n_time    = n_time)

# --- Distance matrix (precomputed, reused at each time step) ------------------

# Euclidean distance of each cell from the center cell (in cell units)
row_idx          <- matrix(rep(seq_len(n_row), times = n_col),
                           nrow = n_row, ncol = n_col)
col_idx          <- matrix(rep(seq_len(n_col), each  = n_row),
                           nrow = n_row, ncol = n_col)
dist_from_center <- sqrt((row_idx - center_row)^2 + (col_idx - center_col)^2)

# --- Compute maximum dispersal distance at each time step ---------------------

# furthest cell from center that contains any non-negligible biomass
biomass_threshold <- 1e-10   # below this value a cell is considered empty

max_dist <- vapply(results, function(grid) {
  occupied <- dist_from_center[grid > biomass_threshold]
  if (length(occupied) == 0) return(0)
  max(occupied)
}, numeric(1))

time_steps <- seq(0, n_time)

# --- Sanity checks ------------------------------------------------------------

# total biomass must be conserved at every step (hard-wall boundary: no flux out)
total_per_step    <- vapply(results, sum, numeric(1))
max_biomass_error <- max(abs(total_per_step - total_biomass))

cat(sprintf("Max biomass conservation error : %.2e\n", max_biomass_error))
cat(sprintf("Max distance at t=0            : %.4f (should be 0)\n", max_dist[1]))
cat(sprintf("Max distance at t=%d           : %.4f cells\n", n_time, max_dist[n_time + 1]))

# --- Plot time series ---------------------------------------------------------

figures_dir <- "figures"
if (!dir.exists(figures_dir)) dir.create(figures_dir)

png(file.path(figures_dir, "test_density_independent_max_distance.png"),
    width = 800, height = 500, res = 120)

plot(time_steps, max_dist,
     type = "l",
     lwd  = 2,
     col  = "steelblue",
     xlab = "Time step",
     ylab = "Maximum dispersal distance from center (cells)",
     main = "Density-independent diffusion — leading dispersal front over time",
     las  = 1)

grid(col = "grey85", lty = 1)

dev.off()

cat("Figure saved to figures/test_density_independent_max_distance.png\n")
