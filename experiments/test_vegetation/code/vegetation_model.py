"""Vegetation simulation model (Python port of vegetation_model.R).

Hierarchical patch-occupancy model with local biomass dynamics.
Based on Harfoot et al. 2014 (Madingley), modified with functional groups
and a local/regional dynamics distinction.

State variables (all (2, a_max) arrays; columns = age classes):
    p       : fraction of cell occupied by patches of type j at age a
    b_tree  : tree biomass (kg C/m2) in patches of type j at age a
    b_herb  : herb biomass (kg C/m2) in patches of type j at age a

Functional groups (rows): row 0 = evergreen, row 1 = deciduous (trees);
herbs present in all patches. Time step: monthly.

Note on indexing vs. the R version: R is 1-based, Python is 0-based. The R
rows 1/2 (evergreen/deciduous) are rows 0/1 here, the oldest age class
`a_max` is index `a_max - 1`, and age class 1 is index 0.
"""

import numpy as np
import pandas as pd


def run_vegetation(p_init, b_tree_init, b_herb_init, params, n_steps):
    """Run the monthly vegetation simulation.

    Arguments:
        p_init      : (2, a_max) array, initial patch occupancy (must sum to 1)
        b_tree_init : (2, a_max) array, initial tree biomass (kg C/m2 per patch)
        b_herb_init : (2, a_max) array, initial herb biomass (kg C/m2 per patch)
        params      : dict (see extract_parameters for structure)
        n_steps     : number of monthly time steps

    Returns a pandas DataFrame with columns step, month, B_ever, B_decid, B_herb
    (total biomass in kg C in the cell per group per time step).
    """
    a_max = params["a_max"]
    A_cell = params["A_cell"]

    # Copy so the caller's initial-condition arrays are not mutated.
    p = np.array(p_init, dtype=float)
    b_tree = np.array(b_tree_init, dtype=float)
    b_herb = np.array(b_herb_init, dtype=float)

    npp_monthly = np.asarray(params["npp_monthly"], dtype=float)
    e = np.asarray(params["e_monthly"], dtype=float)  # [e_ever, e_decid]

    # Per-group annual -> monthly rates (row 0 = evergreen, row 1 = deciduous).
    mu_tree = np.array([params["mu_ever"], params["mu_decid"]]) / 12.0
    beta_tree = np.array([params["beta_ever"], params["beta_decid"]]) / 12.0
    mu_herb_m = params["mu_herb"] / 12.0
    beta_herb_m = params["beta_herb"] / 12.0

    B_ever_out = np.empty(n_steps)
    B_decid_out = np.empty(n_steps)
    B_herb_out = np.empty(n_steps)

    for t in range(n_steps):
        m = t % 12  # 0-based calendar month index (R used 1-based m)
        npp_m = npp_monthly[m]

        # ---- Local biomass dynamics ----
        # Competition factors (Beer-Lambert law): trees shade herbs.
        # C_herb = exp(-alpha * b_tree); C_tree = 1 - C_herb (complement, asymmetric).
        C_herb = np.exp(-params["alpha"] * b_tree)
        C_tree = 1.0 - C_herb

        # Growth: NPP_m * (1 - f_struct) * f_leafmort * C  [kg C/m2/month]
        dG_tree = np.zeros((2, a_max))
        dG_tree[0, :] = npp_m * (1 - params["f_struct_ever"]) * params["f_leafmort_ever"] * C_tree[0, :]
        dG_tree[1, :] = npp_m * (1 - params["f_struct_decid"]) * params["f_leafmort_decid"] * C_tree[1, :]

        dG_herb = np.zeros((2, a_max))
        dG_herb[0, :] = npp_m * (1 - params["f_struct_herb"]) * params["f_leafmort_herb"] * C_herb[0, :]
        dG_herb[1, :] = npp_m * (1 - params["f_struct_herb"]) * params["f_leafmort_herb"] * C_herb[1, :]

        # Mortality: leaf senescence (annual rate / 12) + herbivory.
        # The R sweep() over rows is a row-wise broadcast here: (..)[:, None].
        dM_tree = b_tree * (mu_tree + beta_tree)[:, None]
        dM_herb = (mu_herb_m + beta_herb_m) * b_herb

        # Update biomass, clipped to >= 0.
        b_tree_upd = np.maximum(b_tree + dG_tree - dM_tree, 0.0)
        b_herb_upd = np.maximum(b_herb + dG_herb - dM_herb, 0.0)

        # ---- Total biomass and relative biomass (trees only, for colonization) ----
        B_ever = A_cell * np.sum(p[0, :] * b_tree_upd[0, :])
        B_decid = A_cell * np.sum(p[1, :] * b_tree_upd[1, :])
        B_tree_total = B_ever + B_decid

        if B_tree_total > 0:
            r = np.array([B_ever, B_decid]) / B_tree_total
        else:
            r = np.array([0.5, 0.5])

        # ---- Regional occupancy dynamics ----
        # Total disturbance (patches returning to age 0).
        D = np.sum(e[0] * p[0, :]) + np.sum(e[1] * p[1, :])

        p_new = np.zeros((2, a_max))
        b_tree_new = np.zeros((2, a_max))
        b_herb_new = np.zeros((2, a_max))

        for j in range(2):
            ej = e[j]

            # Lumped oldest class: surviving oldest + newly aged second-oldest.
            p_new[j, a_max - 1] = (p[j, a_max - 1] + p[j, a_max - 2]) * (1 - ej)

            wt_old = p[j, a_max - 1] * (1 - ej)
            wt_new = p[j, a_max - 2] * (1 - ej)
            p_tot = p_new[j, a_max - 1]
            if p_tot > 0:
                b_tree_new[j, a_max - 1] = (wt_old * b_tree_upd[j, a_max - 1] + wt_new * b_tree_upd[j, a_max - 2]) / p_tot
                b_herb_new[j, a_max - 1] = (wt_old * b_herb_upd[j, a_max - 1] + wt_new * b_herb_upd[j, a_max - 2]) / p_tot

            # Middle ages (R 2..a_max-1): shift each patch up one age class.
            if a_max > 2:
                p_new[j, 1:a_max - 1] = p[j, 0:a_max - 2] * (1 - ej)
                b_tree_new[j, 1:a_max - 1] = b_tree_upd[j, 0:a_max - 2]
                b_herb_new[j, 1:a_max - 1] = b_herb_upd[j, 0:a_max - 2]

            # Youngest class: new patches from disturbance (fire resets biomass to 0).
            p_new[j, 0] = D * r[j]
            b_tree_new[j, 0] = 0.0
            b_herb_new[j, 0] = 0.0

        # Update state.
        p = p_new
        b_tree = b_tree_new
        b_herb = b_herb_new

        # Record total biomass (kg C in the cell).
        B_ever_out[t] = A_cell * np.sum(p[0, :] * b_tree[0, :])
        B_decid_out[t] = A_cell * np.sum(p[1, :] * b_tree[1, :])
        B_herb_out[t] = A_cell * (np.sum(p[0, :] * b_herb[0, :]) + np.sum(p[1, :] * b_herb[1, :]))

    steps = np.arange(1, n_steps + 1)
    return pd.DataFrame({
        "step": steps,
        "month": ((steps - 1) % 12) + 1,
        "B_ever": B_ever_out,
        "B_decid": B_decid_out,
        "B_herb": B_herb_out,
    })
