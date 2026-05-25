
suppressPackageStartupMessages(library(ATNr))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_csv <- function(x, name, row.names = FALSE) {
  utils::write.csv(x, file.path(out_dir, name), row.names = row.names)
}

# ---------------- Figure 3: temperature versus extinctions ----------------
set.seed(12)
n_species <- 50
n_basal <- 20
n_nut <- 2
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
L <- create_Lmatrix(masses, n_basal, Ropt = 50, gamma = 2, th = 0.01)
fw <- L
fw[fw > 0] <- 1
model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
temperatures <- seq(4, 22, by = 2)
biomasses <- runif(n_species + n_nut, 2, 3)
times <- seq(0, 100000, 100)
extinctions <- rep(NA, length(temperatures))

for (i in seq_along(temperatures)) {
  temp <- temperatures[[i]]
  model <- initialise_default_Unscaled_nuts(model, L, temperature = temp)
  model$q <- rep(1.4, n_species - n_basal)
  model$S <- rep(10, n_nut)
  sol <- lsoda_wrapper(times, biomasses, model, verbose = FALSE)
  extinctions[[i]] <- sum(sol[nrow(sol), 4:ncol(sol)] < 1e-6)
}

write_csv(data.frame(species = seq_len(n_species) - 1,
                     body_mass_g = masses,
                     is_basal = as.integer(seq_len(n_species) <= n_basal)),
          "figure3_species.csv")
write_csv(fw, "figure3_food_web.csv", row.names = FALSE)
write_csv(data.frame(state = c(paste0("nutrient_", seq_len(n_nut) - 1),
                               paste0("species_", seq_len(n_species) - 1)),
                     initial = biomasses),
          "figure3_initial_state.csv")
write_csv(data.frame(temperature_C = temperatures,
                     extinctions_atnr = extinctions),
          "figure3_atnr_extinctions.csv")

# ---------------- Figure 4: K=1 versus K=10 enrichment scenario ------------
set.seed(1234)
S <- 10
fw4 <- create_niche_model(S, C = .15)
TL <- TroLev(fw4)
masses4 <- as.numeric(0.01 * 100 ^ (TL - 1))
n_basal4 <- sum(colSums(fw4) == 0)
mod <- create_model_Scaled(nb_s = S, nb_b = n_basal4, BM = masses4, fw = fw4)
mod <- initialise_default_Scaled(mod)
times4 <- seq(0, 300, by = 2)
biomasses4 <- runif(S, 2, 3)

mod$K <- 1
sol1 <- lsoda_wrapper(times4, biomasses4, mod, verbose = FALSE)
mod$K <- 10
sol10 <- lsoda_wrapper(times4, biomasses4, mod, verbose = FALSE)

write_csv(data.frame(species = seq_len(S) - 1,
                     body_mass_g = masses4,
                     trophic_level = as.numeric(TL),
                     is_basal = as.integer(seq_len(S) <= n_basal4)),
          "figure4_species.csv")
write_csv(fw4, "figure4_food_web.csv", row.names = FALSE)
write_csv(data.frame(species = seq_len(S) - 1,
                     initial_biomass = biomasses4),
          "figure4_initial_biomass.csv")
write_csv(as.data.frame(sol1), "figure4_atnr_K1_timeseries.csv")
write_csv(as.data.frame(sol10), "figure4_atnr_K10_timeseries.csv")
