source("code/sherbrooke_parameters.R")  # pour calc_miami_npp

# ---- Grille T x P ----
temps  <- seq(-15, 35, length.out = 300)
precips <- seq(0, 3000, length.out = 300)
grid   <- expand.grid(T = temps, P = precips)
grid$NPP <- mapply(calc_miami_npp, grid$T, grid$P)
npp_mat  <- matrix(grid$NPP, nrow = length(temps), ncol = length(precips))

# ---- Villes ----
cities <- data.frame(
  name   = c("Sherbrooke", "Phoenix", "Iqaluit", "Kelowna", "Vancouver", "Panama"),
  T      = c(3.94, 20.56, -9.83, 2.71, 5.28, 25.91),
  P      = c(1256, 277, 419, 697, 2293, 2498),
  NPP    = c(0.519, 0.269, 0.218, 0.489, 0.551, 0.879),
  stringsAsFactors = FALSE
)

# ---- Palette ----
npp_breaks <- seq(0, 0.962, length.out = 101)
pal <- colorRampPalette(c("#2c1654", "#3b4f9e", "#2e8b8b",
                           "#5ab552", "#d4e44a", "#f5c542", "#e8512a"))(100)

png("figures/npp_climate_space.png", width = 1800, height = 1400, res = 150)
par(mar = c(5, 5.5, 3, 6))

image(temps, precips, npp_mat,
      col   = pal,
      breaks = npp_breaks,
      xlab  = "Température annuelle moyenne (°C)",
      ylab  = "Précipitations totales annuelles (mm)",
      main  = "NPP annuelle selon le modèle Miami (Harfoot et al. 2014)",
      cex.lab = 1.3, cex.axis = 1.1, cex.main = 1.3)

# Contours
contour(temps, precips, npp_mat,
        levels = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
        labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9"),
        col = "white", lwd = 0.8, labcex = 0.85, add = TRUE)

# Points villes
points(cities$T, cities$P, pch = 21, bg = "white", col = "black",
       cex = 2.2, lwd = 2)

# Étiquettes (ajustées pour éviter les chevauchements)
offsets <- data.frame(
  dx = c( 0,   0,   0,   1.2,  0,    0),
  dy = c(120, 120, 120, 120, -180, 120),
  adj_x = c(0.5, 0.5, 0.5, 0, 0.5, 0.5)
)
for (i in seq_len(nrow(cities))) {
  text(cities$T[i] + offsets$dx[i],
       cities$P[i] + offsets$dy[i],
       labels = sprintf("%s\n%.3f", cities$name[i], cities$NPP[i]),
       cex = 0.95, font = 2, adj = c(offsets$adj_x[i], 0))
}

# Barre de couleur manuelle
usr <- par("usr")
x1 <- usr[2] + (usr[2] - usr[1]) * 0.02
x2 <- x1     + (usr[2] - usr[1]) * 0.04
y_seq <- seq(usr[3], usr[4], length.out = 101)
for (k in 1:100) {
  rect(x1, y_seq[k], x2, y_seq[k+1],
       col = pal[k], border = NA, xpd = TRUE)
}
rect(x1, usr[3], x2, usr[4], border = "black", lwd = 0.8, xpd = TRUE)
npp_ticks <- c(0, 0.2, 0.4, 0.6, 0.8, 0.962)
for (v in npp_ticks) {
  y_pos <- usr[3] + (v / 0.962) * (usr[4] - usr[3])
  text(x2 + (usr[2]-usr[1])*0.01, y_pos,
       labels = sprintf("%.2f", v), cex = 0.9, adj = 0, xpd = TRUE)
}
text(x1 + (x2-x1)/2, usr[4] + (usr[4]-usr[3])*0.03,
     "NPP\n(kg C m⁻² an⁻¹)", cex = 0.95, adj = 0.5, xpd = TRUE)

dev.off()
cat("Figure sauvegardée : figures/npp_climate_space.png\n")
