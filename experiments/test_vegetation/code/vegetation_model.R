# Vegetation simulation model
# Hierarchical patch-occupancy model with local biomass dynamics
# Based on Harfoot et al. 2014 (Madingley), modified with functional groups
# and local/regional dynamics distinction.
#
# State variables (all matrices: 2 rows = [evergreen, deciduous], a_max cols = age classes):
#   p       : fraction of cell occupied by patches of type j at age a
#   b_tree  : tree biomass (kg C/m2) in patches of type j at age a
#   b_herb  : herb biomass (kg C/m2) in patches of type j at age a
#
# Functional groups: 1=evergreen, 2=deciduous (trees); herbs present in all patches.
# Time step: monthly.

run_vegetation <- function(p_init, b_tree_init, b_herb_init, params, n_steps) {
  # Arguments:
  #   p_init      : 2 x a_max matrix, initial patch occupancy (must sum to 1)
  #   b_tree_init : 2 x a_max matrix, initial tree biomass (kg C/m2 per patch)
  #   b_herb_init : 2 x a_max matrix, initial herb biomass (kg C/m2 per patch)
  #   params      : named list (see extract_parameters.R for structure)
  #   n_steps     : number of monthly time steps

  a_max  <- params$a_max
  A_cell <- params$A_cell

  p      <- p_init
  b_tree <- b_tree_init
  b_herb <- b_herb_init

  # Output storage: total biomass (kg C) in the cell per group per time step
  out <- data.frame(
    step   = 1:n_steps,
    month  = ((1:n_steps - 1) %% 12) + 1,
    B_ever  = NA_real_,
    B_decid = NA_real_,
    B_herb  = NA_real_
  )

  for (t in 1:n_steps) {
    m <- ((t - 1) %% 12) + 1  # calendar month (1-12)

    npp_m <- params$npp_monthly[m]  # kg C/m2 for this month

    # ---- Local biomass dynamics ----
    # Competition factors (Beer-Lambert law)
    # C_herb[j,a] = exp(-alpha * b_tree[j,a])  (shading of herbs by trees)
    # C_tree[j,a] = 1 - C_herb[j,a]            (complement, asymmetric)
    C_herb <- exp(-params$alpha * b_tree)
    C_tree <- 1 - C_herb

    # Growth: delta_G = NPP_m * (1 - f_struct) * f_leafmort * C  [kg C/m2/month]
    # Rows: 1=ever, 2=decid
    dG_tree <- matrix(0, 2, a_max)
    dG_tree[1, ] <- npp_m * (1 - params$f_struct_ever)  * params$f_leafmort_ever  * C_tree[1, ]
    dG_tree[2, ] <- npp_m * (1 - params$f_struct_decid) * params$f_leafmort_decid * C_tree[2, ]

    dG_herb <- matrix(0, 2, a_max)
    dG_herb[1, ] <- npp_m * (1 - params$f_struct_herb) * params$f_leafmort_herb * C_herb[1, ]
    dG_herb[2, ] <- npp_m * (1 - params$f_struct_herb) * params$f_leafmort_herb * C_herb[2, ]

    # Mortality: leaf senescence (annual rate / 12) + herbivory
    mu_tree <- c(params$mu_ever, params$mu_decid) / 12  # monthly leaf mort rate
    mu_herb_m <- params$mu_herb / 12
    beta_tree <- c(params$beta_ever, params$beta_decid) / 12
    beta_herb_m <- params$beta_herb / 12

    dM_tree <- sweep(b_tree, 1, mu_tree + beta_tree, `*`)
    dM_herb <- (mu_herb_m + beta_herb_m) * b_herb

    # Update biomass, clipped to >= 0 (preserve matrix dimensions)
    b_tree_upd <- b_tree + dG_tree - dM_tree; b_tree_upd[b_tree_upd < 0] <- 0
    b_herb_upd <- b_herb + dG_herb - dM_herb; b_herb_upd[b_herb_upd < 0] <- 0

    # ---- Total biomass and relative biomass (trees only for colonization) ----
    B_ever  <- A_cell * sum(p[1, ] * b_tree_upd[1, ])
    B_decid <- A_cell * sum(p[2, ] * b_tree_upd[2, ])
    B_tree_total <- B_ever + B_decid

    r <- if (B_tree_total > 0) c(B_ever, B_decid) / B_tree_total else c(0.5, 0.5)

    # ---- Regional occupancy dynamics ----
    # Monthly disturbance rates
    e <- params$e_monthly  # length-2 vector: [e_ever, e_decid]

    # Total disturbance (patches returning to age 0)
    D <- sum(e[1] * p[1, ]) + sum(e[2] * p[2, ])

    # New occupancy and biomass arrays
    p_new      <- matrix(0, 2, a_max)
    b_tree_new <- matrix(0, 2, a_max)
    b_herb_new <- matrix(0, 2, a_max)

    for (j in 1:2) {
      ej <- e[j]

      # Lumped oldest class: combines surviving a_max patches + newly aged (a_max-1) patches
      p_new[j, a_max] <- (p[j, a_max] + p[j, a_max - 1]) * (1 - ej)

      wt_old  <- p[j, a_max]     * (1 - ej)
      wt_new  <- p[j, a_max - 1] * (1 - ej)
      p_tot   <- p_new[j, a_max]
      if (p_tot > 0) {
        b_tree_new[j, a_max] <- (wt_old * b_tree_upd[j, a_max] + wt_new * b_tree_upd[j, a_max - 1]) / p_tot
        b_herb_new[j, a_max] <- (wt_old * b_herb_upd[j, a_max] + wt_new * b_herb_upd[j, a_max - 1]) / p_tot
      }

      # Ages 2 to (a_max-1): shift from previous age class
      if (a_max > 2) {
        idx_new <- 2:(a_max - 1)
        idx_old <- 1:(a_max - 2)
        p_new[j, idx_new]      <- p[j, idx_old] * (1 - ej)
        b_tree_new[j, idx_new] <- b_tree_upd[j, idx_old]
        b_herb_new[j, idx_new] <- b_herb_upd[j, idx_old]
      }

      # Age 1: new patches from disturbance (fire resets biomass to 0)
      p_new[j, 1]      <- D * r[j]
      b_tree_new[j, 1] <- 0
      b_herb_new[j, 1] <- 0
    }

    # Update state
    p      <- p_new
    b_tree <- b_tree_new
    b_herb <- b_herb_new

    # Record total biomass (kg C in the cell)
    B_herb_total <- A_cell * (sum(p[1, ] * b_herb[1, ]) + sum(p[2, ] * b_herb[2, ]))
    out$B_ever[t]  <- A_cell * sum(p[1, ] * b_tree[1, ])
    out$B_decid[t] <- A_cell * sum(p[2, ] * b_tree[2, ])
    out$B_herb[t]  <- B_herb_total
  }

  return(out)
}
