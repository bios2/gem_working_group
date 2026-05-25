"""
dispersal_3species.py
=====================
Spatial bioenergetic food web — 3-species chain (plant → herbivore → predator)
across 5 patches with Ryser-style dispersal.

Produces a 4-panel figure:
  1. Biomass time series per species, one line per patch
  2. Emigration flux over time per species (one line per patch)
  3. Immigration flux over time per species (one line per patch)
  4. Net dispersal flow between patches at the final time step (heatmap)

Requirements: numpy  scipy  matplotlib
  pip install numpy scipy matplotlib
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm

# ─────────────────────────────────────────────────────────────────────────────
# 1.  FOOD-WEB TOPOLOGY  (fixed 3-species chain)
# ─────────────────────────────────────────────────────────────────────────────
# Species index:  0 = plant (producer), 1 = herbivore, 2 = predator
# A[i, j] = 1  means species i eats species j

S = 3
species_names = ["Plant", "Herbivore", "Predator"]
species_colors = ["#2ca02c", "#ff7f0e", "#d62728"]   # green / orange / red

A = np.array([
    [0, 0, 0],   # plant eats nobody
    [1, 0, 0],   # herbivore eats plant
    [0, 1, 0],   # predator eats herbivore
], dtype=float)

is_prod = np.array([True, False, False])

# ─────────────────────────────────────────────────────────────────────────────
# 2.  ALLOMETRIC PARAMETERS  (Brose 2006)
# ─────────────────────────────────────────────────────────────────────────────
Z_ratio   = 10.0          # predator/prey body-mass ratio
M_basal   = 1.0           # plant body mass (arbitrary units)

# Body masses from trophic level (TL):  M_i = M_basal * Z^(TL_i - 1)
tl = np.array([1.0, 2.0, 3.0])
M  = M_basal * Z_ratio ** (tl - 1.0)   # [1, 10, 100]

# Metabolic rate scaling: x_i = (a_x / a_r) * (M_i / M_basal)^(-0.25)
# Producer growth rate r_i = 1  (normalised)
ax, ar = 0.88, 1.0          # ectotherm vertebrate
xi = np.where(is_prod, 0.0, (ax / ar) * (M / M_basal) ** (-0.25))
ri = np.where(is_prod, 1.0, 0.0)

# Maximum ingestion rate: y_i = a_y / a_x  (constant for ectotherm vertebrates)
yi_val = 4.0
yi = np.where(is_prod, 0.0, yi_val)

# Assimilation efficiencies  e_ij  (plant flesh easier to assimilate)
eij = np.where(A > 0, 0.45, 0.0)   # same for both trophic links here

# Predator interference (Beddington-DeAngelis c term)
c_int = 0.0

# Functional-response half-saturation & Hill exponent
B0 = 0.5
h  = 2.0   # type-III (h=2)

# Carrying capacity (producer only)
K_base = 0.5

# ─────────────────────────────────────────────────────────────────────────────
# 3.  LANDSCAPE
# ─────────────────────────────────────────────────────────────────────────────
n_patches = 5
spacing   = 1.0

# Patch positions (1-D lattice)
pos = np.arange(n_patches, dtype=float) * spacing   # [0, 1, 2, 3, 4]

# Distance matrix
d_mat = np.abs(pos[:, None] - pos[None, :])          # (Z, Z)

# K gradient: left patches oligotrophic, right patches eutrophic
K_gradient = 0.3
K_vec = K_base + K_gradient * (pos / pos.max() - 0.5)   # (Z,)
# Result: K ranges from K_base-0.15 to K_base+0.15

# ─────────────────────────────────────────────────────────────────────────────
# 4.  DISPERSAL PARAMETERS  (Ryser 2021)
# ─────────────────────────────────────────────────────────────────────────────
delta0   = 1.5    # baseline dispersal range
epsilon  = 0.05   # body-mass scaling exponent (small → weak scaling)

# Dispersal range per species: δ_i = δ_0 * M_i^ε
delta = delta0 * M ** epsilon   # (S,)

# Build dispersal-fraction tensor frac[i, n, z]
#   = fraction of species i's emigrants from patch n that arrive at patch z
def build_dispersal_fractions(delta, d_mat):
    Z = d_mat.shape[0]
    S = delta.shape[0]
    frac = np.zeros((S, Z, Z))
    for i in range(S):
        for n in range(Z):
            w = np.maximum(1.0 - d_mat[n] / delta[i], 0.0)
            w[n] = 0.0                    # no self-dispersal
            total = w.sum()
            if total > 0:
                frac[i, n] = w / total
    return frac                           # (S, Z, Z)

frac = build_dispersal_fractions(delta, d_mat)

# Drain/starvation coefficients
alpha_drain = 0.05   # emigration when r_net > 0 (productive → drain)
beta_starv  = 0.02   # emigration when r_net < 0 (starving  → rescue)

# ─────────────────────────────────────────────────────────────────────────────
# 5.  ODE RIGHT-HAND SIDE
# ─────────────────────────────────────────────────────────────────────────────
def functional_response(B_patch, h, B0, c):
    """
    Holling type-III / Beddington-DeAngelis functional response.
    B_patch : (S,) biomass in one patch.
    Returns F (S, S)  where F[i, j] = fraction of time i spends eating j.
    """
    # Prey availability weighted by adjacency
    omega = A / (A.sum(axis=1, keepdims=True) + 1e-30)  # (S, S) equal weights
    Bh   = B_patch ** h                                  # (S,)
    num  = omega * Bh[np.newaxis, :]                     # (S, S)
    denom = (B0 ** h
             + c * B_patch[:, np.newaxis] * B0 ** h
             + (omega * Bh[np.newaxis, :]).sum(axis=1, keepdims=True))
    return num / (denom + 1e-30)                         # (S, S)


EXTINCT = 1e-20   # biomass floor

def spatial_rhs(t, B_flat):
    """Full ODE for the spatial system."""
    B_mat = B_flat.reshape(n_patches, S)   # (Z, S)
    dB    = np.zeros_like(B_mat)

    # ── local bioenergetics ──────────────────────────────────────────────────
    emig = np.zeros_like(B_mat)            # emigration (Z, S)

    for z in range(n_patches):
        B = np.maximum(B_mat[z], 0.0)
        F = functional_response(B, h, B0, c_int)   # (S, S)

        # Producers
        prod_growth = ri * B * (1.0 - B / K_vec[z])   # logistic growth
        prod_loss   = np.zeros(S)
        for i in range(S):
            if not is_prod[i]:
                continue
            # Loss of producer i to consumers
            for j in range(S):
                if A[j, i] > 0:   # j eats i
                    safe_e = eij[j, i] if eij[j, i] > 0 else 1.0
                    prod_loss[i] += xi[j] * yi[j] * B[j] * F[j, i] / safe_e

        # Consumers
        cons_gain = np.zeros(S)
        cons_resp = np.zeros(S)
        for i in range(S):
            if is_prod[i]:
                continue
            cons_gain[i] = xi[i] * yi[i] * B[i] * F[i].sum()
            cons_resp[i] = xi[i] * B[i]

        dB[z] = prod_growth - prod_loss + cons_gain - cons_resp

        # ── dispersal ───────────────────────────────────────────────────────
        # Net growth rate for consumers (Ryser r_net proxy)
        r_net = np.where(~is_prod,
                         xi * (yi * F.sum(axis=1) - 1.0),
                         0.0)   # producers don't disperse

        emig[z] = (alpha_drain * np.maximum(r_net, 0.0) * B
                   + beta_starv * np.maximum(-r_net, 0.0) * B)

    # Immigration: immig[z, i] = Σ_n  frac[i, n, z] * emig[n, i]
    # frac shape: (S, Z, Z)  →  frac[i, n, z]
    immig = np.einsum("inz,ni->zi", frac, emig)   # (Z, S)

    dB += immig - emig

    # Extinction floor: don't drive extinct species further negative
    mask = B_mat < EXTINCT
    dB[mask] = np.maximum(dB[mask], 0.0)

    return dB.ravel()


# ─────────────────────────────────────────────────────────────────────────────
# 6.  INITIAL CONDITIONS  &  INTEGRATION
# ─────────────────────────────────────────────────────────────────────────────
rng = np.random.default_rng(42)

# Stagger initial biomass: only the rightmost (eutrophic) patch starts with
# all three species; left patches start with only plants.
B0_mat = np.zeros((n_patches, S))
B0_mat[:, 0] = K_vec * 0.8            # plants everywhere
B0_mat[-1, 1] = 0.1                   # herbivore only in patch 4
B0_mat[-1, 2] = 0.02                  # predator only in patch 4

t_end = 2000.0
n_t   = 500
t_eval = np.linspace(0, t_end, n_t)

print("Running simulation …")
sol = solve_ivp(
    spatial_rhs,
    [0, t_end],
    B0_mat.ravel(),
    method="RK45",
    t_eval=t_eval,
    rtol=1e-6,
    atol=1e-10,
    max_step=t_end / 200,
)
print(f"  Done. Status: {sol.message}")

# Reshape solution: (n_patches, S, n_t)
B_sol = sol.y.reshape(n_patches, S, n_t)
t_sol = sol.t

# ─────────────────────────────────────────────────────────────────────────────
# 7.  COMPUTE DISPERSAL FLUXES  (post-hoc, for plotting)
# ─────────────────────────────────────────────────────────────────────────────
emig_ts = np.zeros((n_patches, S, n_t))   # emigration time series
immig_ts = np.zeros((n_patches, S, n_t))  # immigration time series

for ti in range(n_t):
    B_mat = np.maximum(B_sol[:, :, ti], 0.0)   # (Z, S)
    emig_t = np.zeros((n_patches, S))

    for z in range(n_patches):
        B = B_mat[z]
        F = functional_response(B, h, B0, c_int)
        r_net = np.where(~is_prod,
                         xi * (yi * F.sum(axis=1) - 1.0),
                         0.0)
        emig_t[z] = (alpha_drain * np.maximum(r_net, 0.0) * B
                     + beta_starv * np.maximum(-r_net, 0.0) * B)

    immig_t = np.einsum("inz,ni->zi", frac, emig_t)
    emig_ts[:, :, ti]  = emig_t
    immig_ts[:, :, ti] = immig_t

# Net flow between every pair of patches for each species at final time step
# net_flow[i, n, z] = emig from n to z  -  emig from z to n
# Positive → net flow from n toward z
ti_final = -1
B_final = np.maximum(B_sol[:, :, ti_final], 0.0)
emig_final = emig_ts[:, :, ti_final]   # (Z, S)

# Directed flow matrix for each species: flow[n, z] = frac[i,n,z]*emig[n,i]
net_flow = {}   # species → (Z, Z) net flow matrix
for i in range(S):
    raw = frac[i] * emig_final[:, i:i+1]   # (Z, Z): row n = outflow from n
    net_flow[i] = raw - raw.T               # positive = net flow n→z

# ─────────────────────────────────────────────────────────────────────────────
# 8.  FIGURE
# ─────────────────────────────────────────────────────────────────────────────
patch_styles = ["-", "--", "-.", ":", (0, (3, 1, 1, 1))]
patch_labels = [f"Patch {z+1}" for z in range(n_patches)]

fig = plt.figure(figsize=(18, 14))
fig.suptitle("Spatial dispersal of 3-species food chain across 5 patches",
             fontsize=15, fontweight="bold", y=0.98)

gs_top = gridspec.GridSpec(1, 3, figure=fig,
                            top=0.91, bottom=0.55,
                            left=0.06, right=0.97, wspace=0.35)
gs_bot = gridspec.GridSpec(1, 3, figure=fig,
                            top=0.48, bottom=0.07,
                            left=0.06, right=0.97, wspace=0.35)

# ── Panel row 1: biomass time series ─────────────────────────────────────────
for sp in range(S):
    ax = fig.add_subplot(gs_top[0, sp])
    for z in range(n_patches):
        ax.plot(t_sol, B_sol[z, sp, :],
                color=species_colors[sp],
                linestyle=patch_styles[z],
                linewidth=1.6,
                label=patch_labels[z])
    ax.set_title(f"{species_names[sp]} — biomass per patch",
                 fontsize=11, color=species_colors[sp])
    ax.set_xlabel("Time", fontsize=9)
    ax.set_ylabel("Biomass (a.u.)", fontsize=9)
    ax.set_yscale("symlog", linthresh=1e-4)
    ax.legend(fontsize=7.5, loc="upper left")
    ax.grid(True, alpha=0.3)

# ── Panel row 2 left: emigration fluxes ──────────────────────────────────────
for sp in range(S):
    ax = fig.add_subplot(gs_bot[0, sp])

    for z in range(n_patches):
        ax.fill_between(t_sol, emig_ts[z, sp, :],
                        alpha=0.35,
                        color=species_colors[sp],
                        label=f"Emig {patch_labels[z]}" if sp == 0 else None)
        ax.plot(t_sol, immig_ts[z, sp, :],
                color=species_colors[sp],
                linestyle=patch_styles[z],
                linewidth=1.4,
                label=f"Immig {patch_labels[z]}" if sp == 0 else None)

    ax.set_title(f"{species_names[sp]} — emigration (fill) / immigration (line)",
                 fontsize=10, color=species_colors[sp])
    ax.set_xlabel("Time", fontsize=9)
    ax.set_ylabel("Flux (a.u. / time)", fontsize=9)
    ax.set_yscale("symlog", linthresh=1e-8)
    ax.grid(True, alpha=0.3)

    if sp == 0:
        ax.legend(fontsize=6.5, loc="upper left", ncol=2)

# ── Inset: net-flow heatmaps at final time (one per species) ─────────────────
# Place three small heatmaps along the right side of the bottom row
# We repurpose gs_bot columns for the net-flow heatmaps
for sp in range(S):
    # Clear the axes we just drew on and build an inset instead
    inset_ax = fig.add_axes(
        [0.06 + sp * 0.305, 0.09, 0.12, 0.12]
    )
    mat = net_flow[sp]   # (Z, Z)
    vmax = np.abs(mat).max() + 1e-30
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im = inset_ax.imshow(mat, cmap="RdBu_r", norm=norm, aspect="auto")
    inset_ax.set_title(f"Net flow\n({species_names[sp]})", fontsize=8,
                       color=species_colors[sp])
    inset_ax.set_xlabel("To patch", fontsize=7)
    inset_ax.set_ylabel("From patch", fontsize=7)
    inset_ax.set_xticks(range(n_patches))
    inset_ax.set_yticks(range(n_patches))
    inset_ax.set_xticklabels([str(z+1) for z in range(n_patches)], fontsize=7)
    inset_ax.set_yticklabels([str(z+1) for z in range(n_patches)], fontsize=7)
    plt.colorbar(im, ax=inset_ax, shrink=0.8, pad=0.05, label="net flux")

# ── Legend annotation ─────────────────────────────────────────────────────────
fig.text(0.5, 0.527,
         "Bottom panels: filled area = emigration flux, lines = immigration flux "
         "| Insets: net dispersal direction between patches at t=2000",
         ha="center", va="center", fontsize=9, color="#444444",
         style="italic")

# ─────────────────────────────────────────────────────────────────────────────
# 9.  SAVE
# ─────────────────────────────────────────────────────────────────────────────
out_path = "dispersal_3species.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Figure saved → {out_path}")
plt.show()