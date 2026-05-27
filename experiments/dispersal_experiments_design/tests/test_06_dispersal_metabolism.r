# test_06_dispersal_metabolism.r
# Tests dispersal_metabolism() over three time steps.
#
# Scenario: Canada lynx (~10 kg mammal) starts in the left half (cols 1-5) of
# a 10x10 grid. At t=1, growth is just below metabolic demands, triggering
# high dispersal into the right half. At t=2, conditions improve to slightly
# above metabolic demands, reducing dispersal. t=3 shows the resulting spread.

source("experiments/dispersal_experiments_design/functions/dispersal_metabolism.r")

# --- Species parameters: Canada lynx ---
mass_g_lynx <- 10000  # ~10 kg
alpha_lynx  <- 0.3    # moderate max per-capita dispersal rate

# Metabolic rate scalar (g biomass / g biomass / day) — same unit conversion
# used inside dispersal_metabolism(), reproduced here for delta_ext setup
fmr_w_per_g     <- calculate_metabolism(mass_g_lynx, group = "endotherm",
                                        endotherm_group = "mammal")
metabolism_rate <- fmr_w_per_g * 86400 / 7000

message(sprintf("Lynx metabolism_rate: %.5f g/g/day", metabolism_rate))

# --- t=1: initial state — lynx occupy left half (columns 1-5) uniformly ---
# Unoccupied cells hold a near-zero floor to avoid division by zero in
# net_growth = delta_ext / biomass - metabolism_rate inside the function
biomass_t1          <- matrix(1e-10, nrow = 10, ncol = 10)
biomass_t1[, 1:5]   <- 100  # g per cell

# --- t=1 -> t=2: net growth just below metabolic demands (high dispersal) ---
# Target: net_growth = 0.99 * metabolism_rate
# => delta_ext / biomass = metabolism_rate + 0.99 * metabolism_rate
# => delta_ext = biomass * 1.99 * metabolism_rate
delta_ext_t1 <- biomass_t1 * 1.99 * metabolism_rate

dispersal_flux_t2 <- dispersal_metabolism(
  biomass_matrix         = biomass_t1,
  delta_biomass_external = delta_ext_t1,
  alpha                  = alpha_lynx,
  mass_g                 = mass_g_lynx,
  group                  = "endotherm",
  endotherm_group        = "mammal"
)

biomass_t2 <- biomass_t1 + dispersal_flux_t2
biomass_t2[biomass_t2 < 1e-10] <- 1e-10

# --- t=2 -> t=3: net growth slightly above metabolic demands (low dispersal) ---
# Target: net_growth = 1.1 * metabolism_rate
# => delta_ext = biomass * 2.1 * metabolism_rate
delta_ext_t2 <- biomass_t2 * 2.1 * metabolism_rate

dispersal_flux_t3 <- dispersal_metabolism(
  biomass_matrix         = biomass_t2,
  delta_biomass_external = delta_ext_t2,
  alpha                  = alpha_lynx,
  mass_g                 = mass_g_lynx,
  group                  = "endotherm",
  endotherm_group        = "mammal"
)

biomass_t3 <- biomass_t2 + dispersal_flux_t3
biomass_t3[biomass_t3 < 1e-10] <- 1e-10

# --- Plot ---
png(
  "experiments/dispersal_experiments_design/figures/test_06_dispersal_metabolism.png",
  width = 1200, height = 480
)

par(mfrow = c(1, 3), mar = c(4, 4, 4, 1), oma = c(0, 0, 3, 0))

# Shared color scale; treat near-zero cells as NA (plotted as white background)
zlim_max <- max(biomass_t1[biomass_t1 >= 1], biomass_t2[biomass_t2 >= 1],
                biomass_t3[biomass_t3 >= 1])
zlim     <- c(0, zlim_max)
cols     <- hcl.colors(64, "YlOrRd", rev = TRUE)

plot_grid <- function(mat, title) {
  display              <- mat
  display[display < 1] <- NA  # mask near-zero as absent
  # image() expects z with nrow = length(x), ncol = length(y);
  # t(mat) maps [col, row] -> (x=col, y=row), row 1 at bottom
  image(
    x    = 1:ncol(mat),
    y    = 1:nrow(mat),
    z    = t(display),
    zlim = zlim,
    col  = cols,
    xlab = "Column",
    ylab = "Row",
    main = title,
    axes = FALSE
  )
  axis(1, at = 1:10)
  axis(2, at = 1:10, las = 1)
  abline(h = seq(0.5, 10.5, 1), col = "grey70", lwd = 0.5)
  abline(v = seq(0.5, 10.5, 1), col = "grey70", lwd = 0.5)
  box()
}

plot_grid(biomass_t1, "t = 1   (initial)")
plot_grid(biomass_t2, "t = 2   (net growth just below metabolism)")
plot_grid(biomass_t3, "t = 3   (net growth slightly above metabolism)")

mtext(
  sprintf(
    "Lynx biomass dispersal  |  mass = %g g,  alpha = %.1f,  metabolism_rate = %.4f /day",
    mass_g_lynx, alpha_lynx, metabolism_rate
  ),
  outer = TRUE, cex = 1.0, font = 2
)

dev.off()

message("Test 6 complete — figure saved to experiments/dispersal_experiments_design/figures/test_06_dispersal_metabolism.png")
