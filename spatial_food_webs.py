"""
Spatial Bioenergetic Food Web Model
=====================================
Brose et al. 2006 (single-patch bioenergetics + allometric scaling)
extended with spatial dispersal from Ryser et al. 2021 (meta-food-webs).

Key equations
-------------
Brose Eq. 1a  Producer:   dBi/dt = ri*(1-Bi/K)*Bi - sum_j[ xj*yj*Bj*F[j,i]/eji ]
Brose Eq. 1b  Consumer:   dBi/dt = -xi*Bi + xi*yi*Bi*sum_j[F[i,j]] - sum_j[ xj*yj*Bj*F[j,i]/eji ]
Brose Eq. 2   Func. resp: F[i,j] = wij*Bj^h / (B0^h + c*Bi*B0^h + sum_k wik*Bk^h)
Brose Eq. 4b  Allometry:  xi = (ax/ar)*(Mi/Mbasal)^(-0.25)
Brose Eq. 5   Body mass:  Mi = Mbasal * Z^(TLi - 1)
Ryser Eq. 9   Emigration: Ei,z = alpha*max(0, r_net)*Bi,z + beta*max(0,-r_net)*Bi,z
Ryser Eq. 11  Immigr.:    Ii,z = sum_n[ frac[i,n,z] * En,i ]
Ryser Eq. 12  Range:      delta_i = delta0 * Mi^epsilon

References
----------
Brose U., Williams R.J., Martinez N.D. (2006)
    Allometric scaling enhances stability in complex food webs.
    Ecology Letters 9:1228-1236.

Ryser R. et al. (2021)
    Landscape heterogeneity buffers biodiversity of simulated meta-food-webs
    under global change through rescue and drainage effects.
    Nature Communications 12:4716.
"""

from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)


# ══════════════════════════════════════════════════════════════════════════════
# §1  FOOD WEB STRUCTURE — Niche model (Williams & Martinez 2000)
# ══════════════════════════════════════════════════════════════════════════════

def niche_model(S: int, C: float, seed: int = 42):
    """
    Generate a food web using the niche model.

    Each species i has a niche value n_i and feeds on species j
    whose niche value falls in [c_i - r_i/2,  c_i + r_i/2].

    Parameters
    ----------
    S    : number of species
    C    : target connectance (fraction of realised trophic links)
    seed : RNG seed for reproducibility

    Returns
    -------
    A    : (S, S) adjacency matrix  —  A[i,j]=1 means i eats j
    nval : (S,) sorted niche values
    """
    rng = np.random.default_rng(seed)

    for _ in range(100):                        # retry until we get a valid web
        nval = np.sort(rng.uniform(0, 1, S))
        beta = max(1.0 / (2 * C) - 1.0, 0.05)  # Beta shape parameter

        # Range sizes drawn from Beta(1, beta), scaled by niche value
        r = nval * rng.beta(1.0, beta, S)

        # Feeding centre uniform in [r/2, n_i]
        c = np.array([
            rng.uniform(r[i] / 2, nval[i]) if nval[i] > r[i] / 2
            else r[i] / 2
            for i in range(S)
        ])

        # Force species 0 to be strictly basal (no prey)
        r[0] = 0.0
        c[0] = nval[0]

        # Build adjacency matrix
        A = np.zeros((S, S))
        for i in range(S):
            lo, hi = c[i] - r[i] / 2, c[i] + r[i] / 2
            for j in range(S):
                if i != j and lo <= nval[j] <= hi:
                    A[i, j] = 1.0

        # Accept only if there is at least one basal species
        if A.sum(axis=1).min() == 0:
            return A, nval

    # Fallback: simple linear chain
    A = np.zeros((S, S))
    for i in range(1, S):
        A[i, i - 1] = 1.0
    return A, np.linspace(0, 1, S)


def trophic_levels(A: np.ndarray) -> np.ndarray:
    """
    Prey-averaged trophic level (TL = 1 for basal producers).
    TL_i = 1 + mean(TL_j  for j in prey of i).
    Converged iteratively.
    """
    S  = A.shape[0]
    tl = np.ones(S)
    for _ in range(300):
        tl_new = np.ones(S)
        for i in range(S):
            prey = np.where(A[i] > 0)[0]
            if len(prey):
                tl_new[i] = 1.0 + np.mean(tl[prey])
        if np.allclose(tl, tl_new, atol=1e-9):
            break
        tl = tl_new
    return tl


def body_masses(tl: np.ndarray, Z: float, M_basal: float = 1.0) -> np.ndarray:
    """
    Body mass from trophic level and constant predator-prey ratio Z.
    Mi = M_basal * Z^(TLi - 1)   [Brose Eq. 5]
    """
    return M_basal * Z ** (tl - 1.0)


def assimilation_matrix(A: np.ndarray, is_prod: np.ndarray,
                         e_herb: float = 0.45,
                         e_carn: float = 0.85) -> np.ndarray:
    """
    eij[i,j] = assimilation efficiency of consumer i digesting prey j.
    e_herb if j is a producer (herbivory), e_carn otherwise (carnivory).
    """
    S   = A.shape[0]
    eij = np.zeros((S, S))
    for i in range(S):
        for j in range(S):
            if A[i, j] > 0:
                eij[i, j] = e_herb if is_prod[j] else e_carn
    return eij


def prey_preferences(A: np.ndarray) -> np.ndarray:
    """
    Uniform prey preference: wij = 1 / (number of prey of i).
    Zero for non-prey species.
    """
    n_prey = A.sum(axis=1, keepdims=True)           # (S, 1)
    return np.where(n_prey > 0, A / n_prey, 0.0)   # (S, S)


# ══════════════════════════════════════════════════════════════════════════════
# §2  ALLOMETRIC PARAMETERS — Brose 2006 Eqs. 3–4
# ══════════════════════════════════════════════════════════════════════════════

def allometric_params(M: np.ndarray, is_prod: np.ndarray,
                       M_basal: float = 1.0,
                       metabolic_type: str = 'ectotherm_vert'):
    """
    Compute allometric rates normalised by the basal producer growth rate.

        xi = (ax/ar) * (M/Mbasal)^(-0.25)    [Brose Eq. 4b]
        yi = ay/ax  (constant)                [Brose Eq. 4c]
        ri = 1 (producers), 0 (consumers)     [Brose Eq. 4a]

    Parameters
    ----------
    metabolic_type : 'invertebrate'  (ax=0.314, yi=8)
                   | 'ectotherm_vert' (ax=0.88,  yi=4)  ← default

    Returns
    -------
    xi, yi, ri : each (S,)
    """
    if metabolic_type == 'invertebrate':
        ax, yi_val = 0.314, 8.0
    else:                              # ectotherm vertebrate
        ax, yi_val = 0.88,  4.0

    ar = 1.0                           # normalised to unity (Brose Eq. 4a)
    xi = (ax / ar) * (M / M_basal) ** (-0.25)
    yi = np.full(len(M), yi_val)
    ri = np.where(is_prod, 1.0, 0.0)
    return xi, yi, ri


# ══════════════════════════════════════════════════════════════════════════════
# §3  FUNCTIONAL RESPONSE — Brose 2006 Eq. 2
# ══════════════════════════════════════════════════════════════════════════════

def functional_response(B: np.ndarray, omega: np.ndarray,
                         B0: float = 0.5,
                         h:  float = 2.0,
                         c:  float = 0.0) -> np.ndarray:
    """
    Generalised functional response (Brose 2006 Eq. 2):

        F[i,j] = wij * Bj^h / (B0^h + c*Bi*B0^h + sum_k wik*Bk^h)

    h=1 → Holling Type II
    h=2 → Holling Type III  (used by Brose for stability)
    c>0 → Beddington-DeAngelis predator interference

    Parameters
    ----------
    B     : (S,) biomass vector
    omega : (S,S) prey preference matrix
    B0    : half-saturation density
    h     : Hill coefficient
    c     : predator interference

    Returns
    -------
    F : (S,S)  F[i,j] = realised fraction of max consumption of i eating j
    """
    Bs   = np.maximum(B, 0.0)
    Bh   = Bs ** h                                  # (S,)
    B0h  = B0 ** h

    # numerator[i,j] = wij * Bj^h
    num   = omega * Bh[np.newaxis, :]               # (S, S)

    # denominator[i] = B0^h + c*Bi*B0^h + sum_j wij*Bj^h
    denom = B0h + c * Bs * B0h + (omega * Bh).sum(axis=1)   # (S,)

    F = np.where(denom[:, np.newaxis] > 0,
                 num / denom[:, np.newaxis], 0.0)
    return F                                        # (S, S)


# ══════════════════════════════════════════════════════════════════════════════
# §4  SINGLE-PATCH TROPHIC DYNAMICS — Brose 2006 Eqs. 1a, 1b
# ══════════════════════════════════════════════════════════════════════════════

def brose_rhs(B: np.ndarray, *, omega, xi, yi, ri, eij,
              K: float, B0: float, h: float, c: float,
              is_prod: np.ndarray) -> np.ndarray:
    """
    Rate of change dB/dt for one well-mixed patch.

    Producers  [Eq. 1a]:
        dBi/dt = ri*(1-Bi/K)*Bi  -  sum_j[ xj*yj*Bj*F[j,i]/eji ]

    Consumers  [Eq. 1b]:
        dBi/dt = -xi*Bi  +  xi*yi*Bi*sum_j[F[i,j]]  -  sum_j[ xj*yj*Bj*F[j,i]/eji ]

    Parameters
    ----------
    B      : (S,) biomass
    omega  : (S,S) prey preferences
    xi, yi, ri : (S,) allometric rates
    eij    : (S,S) assimilation efficiencies
    K      : producer carrying capacity
    B0, h, c : functional response parameters
    is_prod: (S,) boolean, True = producer

    Returns
    -------
    dBdt : (S,)
    """
    B  = np.maximum(B, 0.0)
    F  = functional_response(B, omega, B0, h, c)         # (S, S)

    # ── Predation loss for all species i (sum over predators j) ──────────────
    # flux[j, i] = xj*yj*Bj*F[j,i] / eij[j,i]
    safe_e = np.where(eij > 0, eij, 1.0)
    flux   = xi[:, None] * yi[:, None] * B[:, None] * F / safe_e   # (S, S)
    loss   = flux.sum(axis=0)                             # (S,) summed over predators

    # ── Producer logistic growth ──────────────────────────────────────────────
    growth = ri * (1.0 - B / K) * B                      # zero for consumers (ri=0)

    # ── Consumer metabolic loss and assimilation gain ─────────────────────────
    metab  = np.where(~is_prod, xi * B,              0.0)
    gain   = np.where(~is_prod, xi * yi * B * F.sum(axis=1), 0.0)

    return growth + gain - metab - loss                   # (S,)


# ══════════════════════════════════════════════════════════════════════════════
# §5  SPATIAL LANDSCAPE & DISPERSAL — Ryser 2021 Eqs. 9, 11, 12
# ══════════════════════════════════════════════════════════════════════════════

def dispersal_range(M: np.ndarray,
                    delta0: float = 0.5,
                    epsilon: float = 0.05) -> np.ndarray:
    """
    Maximum dispersal distance for each species  [Ryser Eq. 12]:
        delta_i = delta0 * Mi^epsilon

    Larger-bodied species disperse farther.
    """
    return delta0 * M ** epsilon


def patch_distances(n_patches: int, spacing: float = 1.0) -> np.ndarray:
    """
    Distance matrix for a 1-D patch chain.
    d[n, z] = |n - z| * spacing
    """
    idx = np.arange(n_patches)
    return np.abs(idx[:, None] - idx[None, :]).astype(float) * spacing


def build_dispersal_fractions(delta: np.ndarray,
                               d_mat: np.ndarray) -> np.ndarray:
    """
    Fraction of emigrants of species i leaving patch n that arrive at patch z.

    Weight = max(1 - d[n,z]/delta_i, 0)   — linear decay, zero beyond range.
    Self-dispersal (n==z) is excluded.
    Fractions are normalised so they sum to 1 across target patches.

    Returns
    -------
    frac : (S, n_patches, n_patches)
        frac[i, n, z] = fraction going from patch n to patch z for species i
    """
    S, Z = len(delta), d_mat.shape[0]

    # w[i, n, z] = max(1 - d[n,z]/delta_i, 0)
    w = np.maximum(
        1.0 - d_mat[None, :, :] / delta[:, None, None], 0.0
    )                                                       # (S, Z, Z)

    # Remove self-dispersal: w[i, n, n] = 0
    for n in range(Z):
        w[:, n, n] = 0.0

    # Normalise: each row (source patch n, species i) sums to 1
    total = w.sum(axis=2, keepdims=True)                    # (S, Z, 1)
    frac  = np.where(total > 0, w / total, 0.0)            # (S, Z, Z)
    return frac


# ══════════════════════════════════════════════════════════════════════════════
# §6  FULL SPATIAL ODE SYSTEM
# ══════════════════════════════════════════════════════════════════════════════

def spatial_rhs(t, B_flat, *, S, n_patches, omega, xi, yi, ri, eij,
                K_vec, is_prod, B0, h, c, frac,
                alpha_drain, beta_starv) -> np.ndarray:
    """
    Full dB/dt for all species × all patches.

    State vector layout:
        B_flat[ z*S : (z+1)*S ]  =  biomass of all species on patch z

    Trophic dynamics (per patch):
        Same as single-patch Brose Eqs. 1a, 1b

    Dispersal (Ryser 2021):
        r_net_i,z = xi * (yi * sum_j[F[i,j]] - 1)     [net per-capita growth]

        Emigration [Ryser Eq. 9]:
            Ei,z = alpha * max(0,  r_net) * Bi,z    [drainage: leave when well-fed]
                 + beta  * max(0, -r_net) * Bi,z    [starvation: leave when starving]

        Immigration [Ryser Eq. 11]:
            Ii,z = sum_n[ frac[i,n,z] * En,i ]

    Producers do not disperse (consistent with Ryser 2021).
    """
    B_mat = np.maximum(B_flat.reshape(n_patches, S), 0.0)  # (Z, S)
    dB    = np.zeros_like(B_mat)
    emig  = np.zeros_like(B_mat)                           # (Z, S)

    for z in range(n_patches):
        Bz = B_mat[z]

        # ── Trophic dynamics ────────────────────────────────────────────────
        dB[z] = brose_rhs(Bz, omega=omega, xi=xi, yi=yi, ri=ri, eij=eij,
                           K=K_vec[z], B0=B0, h=h, c=c, is_prod=is_prod)

        # ── Net per-capita growth rate (consumers only) ─────────────────────
        Fz    = functional_response(Bz, omega, B0, h, c)
        r_net = np.where(
            ~is_prod,
            xi * (yi * Fz.sum(axis=1) - 1.0),
            0.0
        )

        # ── Emigration: drainage + starvation ───────────────────────────────
        emig[z] = (alpha_drain * np.maximum( r_net, 0.0) * Bz
                 + beta_starv  * np.maximum(-r_net, 0.0) * Bz)

    # ── Immigration: immig[z, i] = sum_n frac[i,n,z] * emig[n,i] ────────────
    # frac shape: (S, Z_source, Z_target)   emig shape: (Z_source, S)
    immig = np.einsum('inz,ni->zi', frac, emig)            # (Z, S)

    dB += immig - emig

    # ── Extinction floor: prevent numerical negative biomass ─────────────────
    extinct = B_mat < 1e-20
    dB[extinct] = np.maximum(dB[extinct], 0.0)

    return dB.ravel()


# ══════════════════════════════════════════════════════════════════════════════
# §7  SIMULATION PARAMETERS & RUNNER
# ══════════════════════════════════════════════════════════════════════════════

@dataclass
class SimParams:
    """
    All simulation parameters in one place.
    Edit any field before calling run() to customise the simulation.
    """
    # ── Food web ─────────────────────────────────────────────────────────────
    S           : int   = 12      # species richness per patch
    C           : float = 0.15   # food web connectance
    Z_ratio     : float = 10.0   # predator-prey body mass ratio  [Brose Eq. 5]
    metabolic   : str   = 'ectotherm_vert'  # 'invertebrate' | 'ectotherm_vert'

    # ── Functional response ──────────────────────────────────────────────────
    B0          : float = 0.5    # half-saturation density
    h           : float = 2.0    # Hill coefficient (2 = Type III, most stable)
    c_int       : float = 0.0    # predator interference (0 = none)

    # ── Landscape ────────────────────────────────────────────────────────────
    n_patches   : int   = 5      # number of spatial patches
    spacing     : float = 1.0    # distance between adjacent patches
    K_base      : float = 0.5    # carrying capacity of patch 0 (oligotrophic)
    K_gradient  : float = 0.3    # K increases this amount per patch

    # ── Dispersal [Ryser 2021] ───────────────────────────────────────────────
    alpha_drain : float = 0.05   # drainage coefficient (emigrate when well-fed)
    beta_starv  : float = 0.02   # starvation coefficient (emigrate when starving)
    delta0      : float = 1.5    # dispersal range intercept  [Ryser Eq. 12]
    epsilon     : float = 0.05   # body-mass exponent for range [Ryser Eq. 12]

    # ── Integration ──────────────────────────────────────────────────────────
    t_end       : float = 3000.0 # simulation duration
    n_t         : int   = 600    # number of saved time points
    seed        : int   = 42     # RNG seed


def run(p: SimParams | None = None) -> dict:
    """
    Build the food web, assemble the spatial ODE, and integrate it.

    Parameters
    ----------
    p : SimParams  (uses defaults if None)

    Returns
    -------
    dict with keys:
        t         : (n_t,) time vector
        B         : (n_patches, S, n_t) biomass array
        A         : (S,S) adjacency matrix
        M, tl     : (S,) body masses, trophic levels
        is_prod   : (S,) boolean producer mask
        K_vec     : (n_patches,) carrying capacities
        delta     : (S,) dispersal ranges
        d_mat     : (n_patches, n_patches) distance matrix
        p         : SimParams used
        success   : bool
    """
    if p is None:
        p = SimParams()

    rng = np.random.default_rng(p.seed)

    # ── Food web ─────────────────────────────────────────────────────────────
    A, nval   = niche_model(p.S, p.C, seed=p.seed)
    is_prod   = A.sum(axis=1) == 0
    tl        = trophic_levels(A)
    M         = body_masses(tl, p.Z_ratio)
    eij       = assimilation_matrix(A, is_prod)
    omega     = prey_preferences(A)
    xi, yi, ri = allometric_params(M, is_prod, metabolic_type=p.metabolic)

    # ── Landscape ────────────────────────────────────────────────────────────
    K_vec = p.K_base + p.K_gradient * np.arange(p.n_patches)
    d_mat = patch_distances(p.n_patches, p.spacing)
    delta = dispersal_range(M, p.delta0, p.epsilon)
    frac  = build_dispersal_fractions(delta, d_mat)        # (S, Z, Z)

    # ── Initial conditions: random biomass in (0.05, 1), same across patches ─
    B0_patch = rng.uniform(0.05, 1.0, p.S)
    B_init   = np.tile(B0_patch, p.n_patches)

    # ── ODE integration ──────────────────────────────────────────────────────
    t_eval = np.linspace(0, p.t_end, p.n_t)

    sol = solve_ivp(
        lambda t, B: spatial_rhs(
            t, B,
            S=p.S, n_patches=p.n_patches,
            omega=omega, xi=xi, yi=yi, ri=ri, eij=eij,
            K_vec=K_vec, is_prod=is_prod,
            B0=p.B0, h=p.h, c=p.c_int,
            frac=frac,
            alpha_drain=p.alpha_drain,
            beta_starv=p.beta_starv,
        ),
        t_span=[0, p.t_end],
        y0=B_init,
        t_eval=t_eval,
        method='RK45',
        rtol=1e-6,
        atol=1e-10,
        max_step=p.t_end / 200,
    )

    B_sol = np.maximum(sol.y.reshape(p.n_patches, p.S, -1), 0.0)

    return dict(
        t=sol.t, B=B_sol,
        A=A, M=M, tl=tl, nval=nval, is_prod=is_prod,
        K_vec=K_vec, delta=delta, d_mat=d_mat,
        S=p.S, n_patches=p.n_patches,
        p=p, success=sol.success, message=sol.message,
    )


# ══════════════════════════════════════════════════════════════════════════════
# §8  VISUALISATION
# ══════════════════════════════════════════════════════════════════════════════

def plot_results(res: dict, save_path: str | None = None):
    """
    Five-panel figure:
      [0,0] Food web topology (nodes = species, y = trophic level)
      [0,1] Biomass time series — patch 1 (most oligotrophic)
      [0,2] Biomass time series — last patch (most eutrophic)
      [1,0:2] Heatmap: log10 final biomass across patches × species
      [1,2] Total biomass per patch over time
    """
    t, B  = res['t'], res['B']    # B: (Z, S, T)
    A, M  = res['A'], res['M']
    tl    = res['tl']
    ip    = res['is_prod']
    Z     = res['n_patches']
    S     = res['S']
    K_vec = res['K_vec']
    p     = res['p']

    # ── Colour scheme: trophic level → YlOrRd ────────────────────────────────
    cmap_tl = plt.cm.YlOrRd
    norm_tl = mcolors.Normalize(vmin=1.0, vmax=tl.max())
    sp_col  = [cmap_tl(norm_tl(tl[i])) for i in range(S)]

    fig = plt.figure(figsize=(16, 10))
    gs  = fig.add_gridspec(2, 3, hspace=0.42, wspace=0.36)
    ax_web  = fig.add_subplot(gs[0, 0])
    ax_ts0  = fig.add_subplot(gs[0, 1])
    ax_tsZ  = fig.add_subplot(gs[0, 2])
    ax_heat = fig.add_subplot(gs[1, 0:2])
    ax_tot  = fig.add_subplot(gs[1, 2])

    # ── Panel 1: Food web ────────────────────────────────────────────────────
    ax_web.set_title('Food web topology', fontsize=11, fontweight='bold')

    # x-position: spread species within each trophic level
    pos = {}
    unique_tl = sorted(set(np.round(tl, 1)))
    for i in range(S):
        same = [j for j in range(S) if abs(tl[j] - tl[i]) < 0.15]
        pos[i] = (same.index(i) / max(len(same) - 1, 1), tl[i])

    # Draw edges (prey → predator)
    for i in range(S):
        for j in range(S):
            if A[i, j] > 0:
                ax_web.annotate(
                    '', xy=pos[i], xytext=pos[j],
                    arrowprops=dict(arrowstyle='->', color='#bbb',
                                    lw=0.8, alpha=0.7,
                                    connectionstyle='arc3,rad=0.05')
                )

    # Draw nodes
    for i in range(S):
        x, y = pos[i]
        sz = 80 + 60 * np.log10(M[i] + 1)
        ax_web.scatter(x, y, s=sz, c=[sp_col[i]], zorder=5,
                       edgecolors='#333', linewidths=0.6)
        ax_web.text(x, y + 0.06, str(i), ha='center', va='bottom',
                    fontsize=7, color='#222')

    ax_web.set_xlabel('Position within trophic level', fontsize=9)
    ax_web.set_ylabel('Trophic level', fontsize=9)
    sm = plt.cm.ScalarMappable(cmap=cmap_tl, norm=norm_tl)
    plt.colorbar(sm, ax=ax_web, shrink=0.75, label='Trophic level')

    legend_el = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#888',
               markersize=5, label='Producer', markeredgecolor='k'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#c0392b',
               markersize=8, label='Consumer (large M)', markeredgecolor='k'),
    ]
    ax_web.legend(handles=legend_el, fontsize=7, loc='upper left')

    # ── Panels 2 & 3: Time series ────────────────────────────────────────────
    for ax_ts, z_idx, label in [
        (ax_ts0, 0,     f'Patch 1 — oligotrophic  K={K_vec[0]:.2f}'),
        (ax_tsZ, Z - 1, f'Patch {Z} — eutrophic  K={K_vec[-1]:.2f}'),
    ]:
        for i in range(S):
            ls = '-' if ip[i] else '--'
            ax_ts.semilogy(t, B[z_idx, i, :] + 1e-20,
                           color=sp_col[i], lw=1.3, ls=ls, alpha=0.85)
        ax_ts.set_title(label, fontsize=9.5, fontweight='bold')
        ax_ts.set_xlabel('Time', fontsize=9)
        ax_ts.set_ylabel('Biomass (log scale)', fontsize=9)
        ax_ts.set_xlim([0, t[-1]])
        ax_ts.set_ylim([1e-20, None])
        ax_ts.legend(
            handles=[
                Line2D([0], [0], color='k', ls='-',  lw=1.5, label='Producer'),
                Line2D([0], [0], color='k', ls='--', lw=1.5, label='Consumer'),
            ],
            fontsize=7, loc='best'
        )

    # ── Panel 4: Final biomass heatmap ────────────────────────────────────────
    sort_idx = np.argsort(tl)
    B_final  = np.log10(B[:, sort_idx, -1] + 1e-20)     # (Z, S) sorted by TL

    im = ax_heat.imshow(B_final.T, aspect='auto', origin='lower',
                         cmap='viridis', interpolation='nearest',
                         vmin=-15, vmax=0)
    ax_heat.set_xticks(range(Z))
    ax_heat.set_xticklabels(
        [f'P{z+1}\nK={K_vec[z]:.2f}' for z in range(Z)], fontsize=8
    )
    ax_heat.set_yticks(range(S))
    ax_heat.set_yticklabels(
        [f'sp{sort_idx[i]}  TL={tl[sort_idx[i]]:.1f}' for i in range(S)],
        fontsize=7
    )
    ax_heat.set_xlabel('Patch  (oligotrophic → eutrophic)', fontsize=9)
    ax_heat.set_ylabel('Species (sorted by trophic level)', fontsize=9)
    ax_heat.set_title('Final biomass  log₁₀(B)  across patches × species',
                       fontsize=10, fontweight='bold')
    plt.colorbar(im, ax=ax_heat, label='log₁₀ biomass', shrink=0.85)

    # Mark producers with a horizontal line
    prod_sorted = ip[sort_idx]
    for i in range(S):
        if prod_sorted[i]:
            ax_heat.axhline(i, color='lime', lw=0.6, alpha=0.5, ls=':')

    # ── Panel 5: Total biomass per patch ─────────────────────────────────────
    patch_pal = plt.cm.plasma(np.linspace(0.1, 0.9, Z))
    for z in range(Z):
        total = B[z].sum(axis=0)
        ax_tot.semilogy(t, total + 1e-20,
                         color=patch_pal[z], lw=2.0,
                         label=f'P{z+1}  K={K_vec[z]:.2f}')
    ax_tot.set_title('Total biomass per patch', fontsize=10, fontweight='bold')
    ax_tot.set_xlabel('Time', fontsize=9)
    ax_tot.set_ylabel('Σ Biomass (log scale)', fontsize=9)
    ax_tot.legend(fontsize=8, loc='best')
    ax_tot.set_xlim([0, t[-1]])

    # ── Super-title ───────────────────────────────────────────────────────────
    fig.suptitle(
        f'Spatial Brose 2006 + Ryser 2021  |  '
        f'S={p.S}  C={p.C}  Z={p.Z_ratio}  {Z} patches  '
        f'h={p.h} ({"Type III" if p.h==2 else "Type II"})\n'
        f'α_drain={p.alpha_drain}  β_starv={p.beta_starv}  '
        f'δ₀={p.delta0}  ε={p.epsilon}  '
        f'K: {K_vec[0]:.2f} → {K_vec[-1]:.2f}',
        fontsize=10, y=1.01
    )

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved → {save_path}")
    plt.show()


# ══════════════════════════════════════════════════════════════════════════════
# §9  ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    # ── Configure simulation ─────────────────────────────────────────────────
    p = SimParams(
        S           = 12,
        C           = 0.15,
        Z_ratio     = 10.0,      # Brose: Z≥10 sufficient for invertebrate stability
        metabolic   = 'ectotherm_vert',
        h           = 2.0,       # Type III — most stable (Brose 2006 Fig. 2a)
        n_patches   = 5,
        K_base      = 0.5,       # oligotrophic end
        K_gradient  = 0.3,       # each patch richer by 0.3
        alpha_drain = 0.05,      # drainage: emigrate when local growth > 0
        beta_starv  = 0.02,      # starvation: emigrate when starving
        delta0      = 1.5,       # baseline dispersal range
        epsilon     = 0.05,      # body-mass scaling of dispersal
        t_end       = 3000.0,
        seed        = 42,
    )

    print("=" * 60)
    print(" Spatial Bioenergetic Food Web Model")
    print(" Brose 2006 + Ryser 2021 dispersal")
    print("=" * 60)
    print(f" Species:      {p.S}")
    print(f" Connectance:  {p.C}")
    print(f" Z ratio:      {p.Z_ratio}")
    print(f" Patches:      {p.n_patches}")
    print(f" K gradient:   {p.K_base:.2f} → {p.K_base + p.K_gradient*(p.n_patches-1):.2f}")
    print(f" t_end:        {p.t_end:.0f}")
    print()

    # ── Run ──────────────────────────────────────────────────────────────────
    res = run(p)

    status = "✓ success" if res['success'] else f"⚠  {res['message']}"
    print(f"Integration: {status}")
    print(f"Time steps saved: {len(res['t'])}")
    print()

    # ── Summary statistics ────────────────────────────────────────────────────
    B_final = res['B'][:, :, -1]                # (Z, S)
    persist = (B_final > 1e-10).sum(axis=1)
    print("Final species persistence per patch (threshold 1e-10):")
    for z in range(p.n_patches):
        prod_count = (B_final[z] > 1e-10)[res['is_prod']].sum()
        cons_count = (B_final[z] > 1e-10)[~res['is_prod']].sum()
        print(f"  Patch {z+1}  K={res['K_vec'][z]:.2f}  "
              f"{persist[z]}/{p.S} species  "
              f"({prod_count} prod, {cons_count} cons)")

    print()
    print("Dispersal ranges (delta_i = delta0 * Mi^epsilon):")
    for i in range(p.S):
        flag = "prod" if res['is_prod'][i] else f"TL={res['tl'][i]:.1f}"
        print(f"  sp{i:2d}  M={res['M'][i]:8.3f}  delta={res['delta'][i]:.4f}  [{flag}]")

    # ── Plot ──────────────────────────────────────────────────────────────────
    out_path = "spatial_foodweb_result.png"
    plot_results(res, save_path=out_path)