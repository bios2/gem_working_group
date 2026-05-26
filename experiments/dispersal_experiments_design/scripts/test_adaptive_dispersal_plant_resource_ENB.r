#######################
# Test diffuse_step_adaptive() with plant-resource driven population dynamics.
# Species growth per individual is proportional to plant material per individual
# (density-dependent via shared resource). High max_disp_rate drives rapid
# spread from an initial single-cell cluster.
#
# Consumer-resource model:
#   per_cap_growth  = alpha * (plant / pop)     [decreases as pop grows]
#   net_growth      = per_cap_growth - metabolic_rate
#   d(pop)/dt       = net_growth * pop
#   d(plant)/dt     = r_plant * (K_plant - plant) - feed_rate * pop
#
# Created: 26 May 2026
# Created by: ENB
# Last updated:


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

# Part 0. Script setup
rm(list = ls())

source("experiments/dispersal_experiments_design/functions/enforce_boundary_conditions.r")
source("experiments/dispersal_experiments_design/functions/diffuse_step_adaptive.r")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 1. Grid and plant resource setup

n_row <- 50
n_col <- 50

# Plant carrying capacity: heterogeneous grid with a resource-rich patch
# in the bottom-right quadrant to create an attractive sink for dispersers
K_plant     <- matrix(5,  nrow = n_row, ncol = n_col)
K_plant[6:10, 6:10] <- 15   # resource-rich 5x5 patch

# Initialise plant biomass at carrying capacity (before any grazing)
plant_grid <- K_plant

# Fixed neighbour-count matrix: 2 (corner), 3 (edge), or 4 (interior)
nbr_count <- enforce_boundary_conditions(n_row, n_col)


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 2. Species and simulation parameters

metabolic_rate <- 0.1    # constant per-capita metabolic cost (sigmoid inflection)
max_disp_rate  <- 0.8    # high: species disperses aggressively when growth is poor
alpha          <- 0.5    # per_cap_growth = alpha * (plant / pop); scales plant-to-growth conversion
r_plant        <- 0.3    # logistic plant regeneration rate
feed_rate      <- 0.05   # per-capita plant consumption rate

dt     <- 0.1    # Euler integration time step
n_time <- 300    # number of steps (= 30 time units)


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 3. Initial conditions: 100 individuals concentrated in the top-left cell

pop_grid       <- matrix(0, nrow = n_row, ncol = n_col)
pop_grid[1, 1] <- 100


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 4. Storage

pop_history   <- vector("list", n_time + 1)
plant_history <- vector("list", n_time + 1)
pop_history[[1]]   <- pop_grid
plant_history[[1]] <- plant_grid


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 5. Simulation loop

for (t in seq_len(n_time)) {

  # Per-capita growth proportional to plant material per individual.
  # pmax() prevents division-by-zero in unpopulated cells.
  safe_pop       <- pmax(pop_grid, 1e-6)
  per_cap_growth <- alpha * (plant_grid / safe_pop)

  # Net per-capita rate: growth gain minus constant metabolic cost.
  # Negative when the cell is overcrowded relative to available plant.
  # Passed to diffuse_step_adaptive as the signal driving dispersal intensity.
  net_growth_grid <- per_cap_growth - metabolic_rate

  # Population update (Euler): grows when net_growth > 0, declines otherwise
  pop_grid <- pop_grid + pop_grid * net_growth_grid * dt
  pop_grid <- pmax(pop_grid, 0)

  # Plant update: logistic regeneration offset by grazing from the population
  plant_grid <- plant_grid +
    r_plant * (K_plant - plant_grid) * dt -
    feed_rate * pop_grid * dt
  plant_grid <- pmax(plant_grid, 0)

  # Density-dependent dispersal: high emigration from cells where
  # net_growth < metabolic_rate (plant per individual is too low)
  pop_grid <- diffuse_step_adaptive(
    pop_grid        = pop_grid,
    net_growth_grid = net_growth_grid,
    metabolic_rate  = metabolic_rate,
    max_disp_rate   = max_disp_rate,
    nbr_count       = nbr_count
  )

  pop_history[[t + 1]]   <- pop_grid
  plant_history[[t + 1]] <- plant_grid
}


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 6. Diagnostics

total_pop <- sapply(pop_history, sum)

cat("--- Population summary ---\n")
cat("Initial total:  ", round(total_pop[1], 2), "\n")
cat("Final total:    ", round(total_pop[n_time + 1], 2), "\n")
cat("Min total:      ", round(min(total_pop), 2), "\n")
cat("Max total:      ", round(max(total_pop), 2), "\n")

# Analytical equilibrium: alpha * plant_eq / pop_eq = metabolic_rate
#   and  plant_eq = K - feed_rate * pop_eq / r_plant
# => pop_eq = K / (metabolic_rate / alpha + feed_rate / r_plant)
pop_eq_poor <- K_plant[1, 1] / (metabolic_rate / alpha + feed_rate / r_plant)
pop_eq_rich <- K_plant[10, 10] / (metabolic_rate / alpha + feed_rate / r_plant)
expected_eq  <- 75 * pop_eq_poor + 25 * pop_eq_rich

cat("\n--- Expected equilibrium (analytical) ---\n")
cat("Per-cell pop (K=5 patch): ", round(pop_eq_poor, 2), "\n")
cat("Per-cell pop (K=15 patch):", round(pop_eq_rich, 2), "\n")
cat("Grid total:               ", round(expected_eq, 2), "\n")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 7. Plots

plot_steps <- c(1, 5, 20, 75, 150, 300)

# Helper: display a matrix with row 1 at the top (north) using image()
plot_grid <- function(mat, main, col_palette, zlim = NULL) {
  z <- t(mat)[, nrow(mat):1]
  if (is.null(zlim)) zlim <- range(z, na.rm = TRUE)
  image(z, main = main, col = col_palette, zlim = zlim,
        xaxt = "n", yaxt = "n", asp = 1)
}

# --- 7a. Total population trajectory ---
par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))
plot(seq(0, n_time) * dt, total_pop, type = "l", lwd = 2,
     xlab = "Time", ylab = "Total population",
     main = "Total population over time")
abline(h = expected_eq, lty = 2, col = "grey50")
legend("bottomright", legend = "Analytical equilibrium",
       lty = 2, col = "grey50", bty = "n")

# --- 7b. Population snapshots ---
all_pop_vals <- unlist(lapply(pop_history[plot_steps + 1], as.vector))
pop_zlim     <- range(all_pop_vals, na.rm = TRUE)
pop_pal      <- hcl.colors(50, "YlOrRd", rev = TRUE)

par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (t in plot_steps) {
  plot_grid(pop_history[[t + 1]],
            main      = paste("Pop  t =", t * dt),
            col_palette = pop_pal,
            zlim      = pop_zlim)
}

# --- 7c. Plant resource snapshots ---
all_plant_vals <- unlist(lapply(plant_history[plot_steps + 1], as.vector))
plant_zlim     <- range(all_plant_vals, na.rm = TRUE)
plant_pal      <- hcl.colors(50, "Greens", rev = TRUE)

par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (t in plot_steps) {
  plot_grid(plant_history[[t + 1]],
            main      = paste("Plant  t =", t * dt),
            col_palette = plant_pal,
            zlim      = plant_zlim)
}
 






#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 2. Shared biological parameters (identical for both species)

metabolic_rate     <- 0.1    # constant per-capita metabolic cost
alpha              <- 0.5    # per_cap_growth = alpha * (plant / pop)
r_plant            <- 0.3    # plant logistic regeneration rate
feed_rate          <- 0.05   # per-capita plant consumption rate

# Species-specific dispersal rates
max_disp_fast   <- 0.80   # fast disperser (original species)
max_disp_lizard <- 0.05   # lizard: same biology, far lower mobility

dt     <- 0.1
n_time <- 300


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 3. Helper: run one simulation

run_sim <- function(max_disp_rate) {
  pop_grid   <- matrix(0, nrow = n_row, ncol = n_col)
  pop_grid[1, 1] <- 100
  plant_grid <- K_plant

  pop_hist   <- vector("list", n_time + 1)
  plant_hist <- vector("list", n_time + 1)
  pop_hist[[1]]   <- pop_grid
  plant_hist[[1]] <- plant_grid

  for (t in seq_len(n_time)) {
    safe_pop        <- pmax(pop_grid, 1e-6)
    net_growth_grid <- alpha * (plant_grid / safe_pop) - metabolic_rate

    pop_grid   <- pmax(pop_grid + pop_grid * net_growth_grid * dt, 0)
    plant_grid <- pmax(
      plant_grid + r_plant * (K_plant - plant_grid) * dt - feed_rate * pop_grid * dt,
      0
    )

    pop_grid <- diffuse_step_adaptive(
      pop_grid        = pop_grid,
      net_growth_grid = net_growth_grid,
      metabolic_rate  = metabolic_rate,
      max_disp_rate   = max_disp_rate,
      nbr_count       = nbr_count
    )

    pop_hist[[t + 1]]   <- pop_grid
    plant_hist[[t + 1]] <- plant_grid
  }

  list(pop = pop_hist, plant = plant_hist)
}


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 4. Run both simulations

cat("Running fast-disperser simulation...\n")
sim_fast   <- run_sim(max_disp_fast)

cat("Running lizard simulation...\n")
sim_lizard <- run_sim(max_disp_lizard)


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 5. Diagnostics

total_fast   <- sapply(sim_fast$pop,   sum)
total_lizard <- sapply(sim_lizard$pop, sum)

# Analytical equilibrium (same for both species since biology is identical)
pop_eq_poor <- K_plant[1, 1]    / (metabolic_rate / alpha + feed_rate / r_plant)
pop_eq_rich <- K_plant[10, 10]  / (metabolic_rate / alpha + feed_rate / r_plant)
expected_eq <- 75 * pop_eq_poor + 25 * pop_eq_rich

cat("\n--- Fast disperser ---\n")
cat("Final total:", round(total_fast[n_time + 1], 2),
    "  (expected:", round(expected_eq, 2), ")\n")

cat("\n--- Lizard ---\n")
cat("Final total:", round(total_lizard[n_time + 1], 2),
    "  (expected:", round(expected_eq, 2), ")\n")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 6. Plots

plot_steps <- c(1, 5, 20, 75, 150, 300)
time_axis  <- seq(0, n_time) * dt

# Shared colour palettes
pop_pal  <- hcl.colors(50, "YlOrRd",  rev = TRUE)
liz_pal  <- hcl.colors(50, "BluYl",   rev = TRUE)
diff_pal <- hcl.colors(50, "Blue-Red 3")   # diverging: blue = lizard ahead, red = fast ahead

# Shared population zlim across both species for fair spatial comparison
all_pop_vals <- c(
  unlist(lapply(sim_fast$pop[plot_steps + 1],   as.vector)),
  unlist(lapply(sim_lizard$pop[plot_steps + 1], as.vector))
)
pop_zlim <- range(all_pop_vals, na.rm = TRUE)

# Helper: display a matrix with row 1 at the top (north)
plot_grid <- function(mat, main, col_palette, zlim = NULL) {
  z <- t(mat)[, nrow(mat):1]
  if (is.null(zlim)) zlim <- range(z, na.rm = TRUE)
  image(z, main = main, col = col_palette, zlim = zlim,
        xaxt = "n", yaxt = "n", asp = 1)
}

# --- 6a. Total population trajectories ---
par(mfrow = c(1, 1), mar = c(4, 4, 2, 4))
plot(time_axis, total_fast, type = "l", lwd = 2, col = "#D62728",
     xlab = "Time", ylab = "Total population",
     main = "Total population: fast disperser vs. lizard",
     ylim = range(c(total_fast, total_lizard)))
lines(time_axis, total_lizard, lwd = 2, col = "#1F77B4")
abline(h = expected_eq, lty = 2, col = "grey50")
legend("bottomright",
       legend = c(paste0("Fast disperser  (d_max = ", max_disp_fast, ")"),
                  paste0("Lizard          (d_max = ", max_disp_lizard, ")"),
                  "Analytical equilibrium"),
       col  = c("#D62728", "#1F77B4", "grey50"),
       lty  = c(1, 1, 2), lwd = c(2, 2, 1), bty = "n")

# --- 6b. Fast-disperser population snapshots ---
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (t in plot_steps) {
  plot_grid(sim_fast$pop[[t + 1]],
            main        = paste("Fast  t =", t * dt),
            col_palette = pop_pal,
            zlim        = pop_zlim)
}

# --- 6c. Lizard population snapshots (same colour scale) ---
par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (t in plot_steps) {
  plot_grid(sim_lizard$pop[[t + 1]],
            main        = paste("Lizard  t =", t * dt),
            col_palette = pop_pal,
            zlim        = pop_zlim)
}

# --- 6d. Difference maps: fast disperser minus lizard ---
# Red: fast disperser has more individuals; blue: lizard has more
diff_vals <- lapply(plot_steps, function(t)
  sim_fast$pop[[t + 1]] - sim_lizard$pop[[t + 1]]
)
max_abs_diff <- max(abs(unlist(diff_vals)), na.rm = TRUE)
diff_zlim    <- c(-max_abs_diff, max_abs_diff)

par(mfrow = c(2, 3), mar = c(1, 1, 2, 1))
for (i in seq_along(plot_steps)) {
  t <- plot_steps[i]
  plot_grid(diff_vals[[i]],
            main        = paste("Diff (fast - lizard)  t =", t * dt),
            col_palette = diff_pal,
            zlim        = diff_zlim)
}


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_
# Part 7. Empirical speed of spread: distance vs. time

# Euclidean distance of each cell from the origin cell [1, 1]
row_idx   <- matrix(rep(1:n_row, n_col),       nrow = n_row, ncol = n_col)
col_idx   <- matrix(rep(1:n_col, each = n_row), nrow = n_row, ncol = n_col)
cell_dist <- sqrt((row_idx - 1)^2 + (col_idx - 1)^2)  # in grid cells

# Centroid distance: population-weighted mean distance from [1, 1].
# Tracks where the average individual is over time.
centroid_dist_ts <- function(pop_hist) {
  sapply(pop_hist, function(pop) {
    total <- sum(pop)
    if (total < 1e-10) return(0)
    sum(cell_dist * pop) / total
  })
}

# Leading-edge distance: furthest cell from [1, 1] with any detectable population.
# Tracks the colonisation frontier.
leading_edge_ts <- function(pop_hist, threshold = 0.5) {
  sapply(pop_hist, function(pop) {
    if (!any(pop > threshold)) return(0)
    max(cell_dist[pop > threshold])
  })
}

cd_fast   <- centroid_dist_ts(sim_fast$pop)
cd_lizard <- centroid_dist_ts(sim_lizard$pop)
le_fast   <- leading_edge_ts(sim_fast$pop)
le_lizard <- leading_edge_ts(sim_lizard$pop)

# Maximum possible distance on the grid (origin to opposite corner)
max_dist <- sqrt((n_row - 1)^2 + (n_col - 1)^2)

# --- 7a. Centroid displacement vs. time ---
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

plot(time_axis, cd_fast, type = "l", lwd = 2, col = "#D62728",
     ylim = c(0, max_dist),
     xlab = "Time", ylab = "Distance from origin (grid cells)",
     main = "Centroid displacement")
lines(time_axis, cd_lizard, lwd = 2, col = "#1F77B4")
abline(h = max_dist, lty = 3, col = "grey70")
legend("topleft",
       legend = c(paste0("Fast  (d_max = ", max_disp_fast, ")"),
                  paste0("Lizard (d_max = ", max_disp_lizard, ")"),
                  "Grid diagonal"),
       col = c("#D62728", "#1F77B4", "grey70"),
       lty = c(1, 1, 3), lwd = c(2, 2, 1), bty = "n")

# --- 7b. Leading-edge distance vs. time ---
plot(time_axis, le_fast, type = "l", lwd = 2, col = "#D62728",
     ylim = c(0, max_dist),
     xlab = "Time", ylab = "Distance from origin (grid cells)",
     main = "Leading-edge distance")
lines(time_axis, le_lizard, lwd = 2, col = "#1F77B4")
abline(h = max_dist, lty = 3, col = "grey70")
legend("topleft",
       legend = c(paste0("Fast  (d_max = ", max_disp_fast, ")"),
                  paste0("Lizard (d_max = ", max_disp_lizard, ")"),
                  "Grid diagonal"),
       col = c("#D62728", "#1F77B4", "grey70"),
       lty = c(1, 1, 3), lwd = c(2, 2, 1), bty = "n")
