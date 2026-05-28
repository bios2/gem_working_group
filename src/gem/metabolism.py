"""Mass-specific metabolic rate and the resulting biomass loss.

Four pieces, used together:

  metabolism_parameters(species, taxa, ...) -> DataFrame[species, taxa, c_int, b]
  body_temp(taxa)                            -> body temperature in degrees C
  calculate_massspecific_metabolism(mass_g, body_temp_C, c_int, b) -> W/g
  calculate_biomass_metabolism(initial_biomass_g, mass_specific_metabolic_rate, dt)
                                             -> biomass delta (g, signed)

Only ``calculate_massspecific_metabolism`` is a numpy science function in the
sense of ``docs/processes_implementation_specification.md`` (typed arrays, equal
shapes, runtime assert). ``metabolism_parameters`` and ``body_temp`` are
**parameter helpers** that build inputs for the science function from a species
list — they are intended to be called once at experiment setup, not inside the
engine's per-step loop, and ``metabolism_parameters`` may prompt the user
interactively when ``interactive=True``.

Sources:
  - Blyth et al. (2026) Ecology Letters 29: e70330 — PGLS intercepts (C) for
    endothermic vs ectothermic vertebrates and the allometric exponents (b)
    used as defaults below.
  - Gillooly, Gomez & Mavrodiev (2017) Proc. R. Soc. B 284: 20162328 — bird
    and mammal body temperatures.
  - Brose et al. (2008) — factor of 3 to convert resting rate to field
    metabolic rate (FMR).

Whole-body and mass-specific allometric form (from Blyth et al.):

    X [W]   = exp(C) * M^b     * exp(-E_a / (k_B * T))     (whole body)
    x [W/g] = exp(C) * M^(b-1) * exp(-E_a / (k_B * T))     (mass-specific)

where M is body mass in grams and T is body temperature in Kelvin.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
import pandas as pd
from numpy.typing import NDArray

# --- defaults from the literature ---------------------------------------------
# Keys here are the taxa labels the rest of the module accepts. Endotherms are
# split bird/mammal because their body temperatures differ (Gillooly 2017),
# even though they share C and b from Blyth Table 1.
DEFAULT_TAXA_PARAMS: dict[str, dict[str, float]] = {
    "ectotherm":        {"c_int": 17.4,  "b": 0.84},
    "endotherm-bird":   {"c_int": 19.53, "b": 0.73},
    "endotherm-mammal": {"c_int": 19.53, "b": 0.73},
}

# Body temperatures in degrees Celsius. Ectotherms intentionally absent — see
# body_temp() below.
DEFAULT_BODY_TEMP_C: dict[str, float] = {
    "endotherm-bird":   41.5,
    "endotherm-mammal": 36.5,
}

# Physical constants for the Boltzmann-Arrhenius term.
E_A_DEFAULT = 0.63       # activation energy, eV
K_B_DEFAULT = 8.617e-5   # Boltzmann constant, eV / K
FMR_MULTIPLIER_DEFAULT = 3.0  # Brose et al. 2008: resting -> field metabolic rate

KELVIN_OFFSET = 273.15

# Energy conversion: 1 kg wet biomass = 7,000,000 J  ->  7,000 J/g.
# Used to translate a metabolic rate in W/g (= J/s/g) into biomass loss in g/s.
ENERGY_DENSITY_J_PER_G_DEFAULT = 7000.0
SECONDS_PER_DAY = 86400.0

VALID_TAXA = tuple(DEFAULT_TAXA_PARAMS.keys())


# --- body temperature ---------------------------------------------------------
def body_temp(taxa: str | Iterable[str]) -> float | NDArray[np.float64]:
    """Body temperature in degrees Celsius for endothermic taxa.

    Accepts either a single taxa label (str) or an iterable of labels. Returns
    a single float in the first case, a 1-D ``NDArray[np.float64]`` in the
    second.

    For endotherms (``"endotherm-bird"``, ``"endotherm-mammal"``) this returns
    the species' regulated body temperature (Gillooly 2017). For ectotherms
    body temperature equals ambient temperature, which this function does not
    know about; ``"ectotherm"`` therefore returns ``np.nan`` as a sentinel,
    and ``calculate_massspecific_metabolism`` fills NaN entries in with its
    ``ambient_temp_C`` argument.

    Examples
    --------
    >>> body_temp("endotherm-bird")
    41.5
    >>> body_temp(["endotherm-bird", "ectotherm"]).tolist()
    [41.5, nan]
    """
    # ``isinstance(taxa, str)`` first because strings are also iterable in
    # Python — without this branch a single label like "endotherm-bird" would
    # be treated as a sequence of characters.
    if isinstance(taxa, str):
        return _lookup_body_temp(taxa)

    return np.array([_lookup_body_temp(t) for t in taxa], dtype=np.float64)


def _lookup_body_temp(taxa: str) -> float:
    # Ectotherms are signalled with NaN so that calculate_massspecific_metabolism
    # can fill them in from its ambient_temp_C argument. Keeping the slot in the
    # array (instead of erroring or removing the row) means downstream arrays
    # stay aligned per-species.
    if taxa == "ectotherm":
        return float("nan")
    if taxa not in DEFAULT_BODY_TEMP_C:
        raise ValueError(
            f"Unknown taxa {taxa!r}. Expected 'ectotherm' or one of "
            f"{tuple(DEFAULT_BODY_TEMP_C.keys())}."
        )
    return DEFAULT_BODY_TEMP_C[taxa]


# --- metabolism_parameters ----------------------------------------------------
def metabolism_parameters(
    species: str | Iterable[str],
    taxa: str | Iterable[str | None] | None = None,
    custom_params: dict[str, dict[str, float]] | None = None,
    interactive: bool = True,
) -> pd.DataFrame:
    """Build a per-species table of metabolism coefficients (``c_int``, ``b``).

    The table this returns is the canonical input to
    ``calculate_massspecific_metabolism``: one row per species, columns
    ``species``, ``taxa``, ``c_int``, ``b``. Pull the columns you need with
    ``.values`` and feed them in as numpy arrays.

    Parameters
    ----------
    species
        A species name or an iterable of species names. Names are free-form
        strings; they only need to be unique.
    taxa
        Either ``None`` (every species' classification is unknown), a single
        label applied to all species, or an iterable of labels — same length as
        ``species`` — where each element is one of ``VALID_TAXA`` (``"ectotherm"``,
        ``"endotherm-bird"``, ``"endotherm-mammal"``) or ``None`` to mark that
        species as unclassified.
    custom_params
        Optional override mapping ``species_name -> {"c_int": float, "b": float,
        "taxa": str}``. Use this for non-vertebrate species, or any species with
        published coefficients of its own. Custom entries skip the interactive
        prompt entirely. ``"taxa"`` may be any descriptive label — it is stored
        but never validated against ``VALID_TAXA``.
    interactive
        When ``True`` (the default), the function asks the user via ``input()``
        what to do for any species that is still unclassified after applying
        ``taxa`` and ``custom_params``. The user can either pick one of the
        default taxa or type in their own ``c_int`` and ``b``.

        When ``False``, an unclassified species raises ``ValueError`` instead.
        Use ``interactive=False`` in notebooks, tests, and batch scripts where
        an ``input()`` call would hang or be unreproducible.

    Returns
    -------
    pandas.DataFrame
        Columns: ``species`` (str), ``taxa`` (str), ``c_int`` (float),
        ``b`` (float). One row per species, in the same order as the input.

    Examples
    --------
    >>> df = metabolism_parameters(
    ...     species=["robin", "vole"],
    ...     taxa=["endotherm-bird", "endotherm-mammal"],
    ... )
    >>> df[["species", "c_int", "b"]].to_dict("records")
    [{'species': 'robin', 'c_int': 19.53, 'b': 0.73}, {'species': 'vole', 'c_int': 19.53, 'b': 0.73}]
    """
    # Normalise species to a list. ``[species]`` (rather than ``list(species)``)
    # keeps a single string from being split into characters.
    species_list = [species] if isinstance(species, str) else list(species)

    # Normalise taxa to a list of the same length, with None for unknowns.
    if taxa is None:
        taxa_list: list[str | None] = [None] * len(species_list)
    elif isinstance(taxa, str):
        taxa_list = [taxa] * len(species_list)
    else:
        taxa_list = list(taxa)
        if len(taxa_list) != len(species_list):
            raise ValueError(
                f"taxa has length {len(taxa_list)} but species has length "
                f"{len(species_list)}; they must match."
            )

    custom_params = custom_params or {}

    rows: list[dict[str, object]] = []
    for sp, tx in zip(species_list, taxa_list):
        if sp in custom_params:
            entry = custom_params[sp]
            rows.append({
                "species": sp,
                "taxa":    entry.get("taxa", "custom"),
                "c_int":   float(entry["c_int"]),
                "b":       float(entry["b"]),
            })
            continue

        if tx is None:
            if interactive:
                rows.append(_prompt_for_unknown(sp))
            else:
                raise ValueError(
                    f"Species {sp!r} has no taxa classification, no entry in "
                    f"custom_params, and interactive=False. Provide one of "
                    f"those three before calling the function."
                )
            continue

        if tx not in DEFAULT_TAXA_PARAMS:
            raise ValueError(
                f"Unknown taxa label {tx!r} for species {sp!r}. Expected one "
                f"of {VALID_TAXA}, or pass a custom_params entry."
            )
        params = DEFAULT_TAXA_PARAMS[tx]
        rows.append({
            "species": sp,
            "taxa":    tx,
            "c_int":   float(params["c_int"]),
            "b":       float(params["b"]),
        })

    return pd.DataFrame(rows, columns=["species", "taxa", "c_int", "b"])


def _prompt_for_unknown(species_name: str) -> dict[str, object]:
    """Ask the user how to classify a single species. Returns one row."""
    print(
        f"\nSpecies {species_name!r} has no taxa classification.\n"
        f"  1) ectotherm        (c_int=17.4,  b=0.84)\n"
        f"  2) endotherm-bird   (c_int=19.53, b=0.73, body_temp=41.5 C)\n"
        f"  3) endotherm-mammal (c_int=19.53, b=0.73, body_temp=36.5 C)\n"
        f"  4) custom — enter your own c_int and b\n"
    )
    choice = input(f"  Choice for {species_name} [1-4]: ").strip()

    if choice in {"1", "2", "3"}:
        tx = VALID_TAXA[int(choice) - 1]
        params = DEFAULT_TAXA_PARAMS[tx]
        return {
            "species": species_name,
            "taxa":    tx,
            "c_int":   float(params["c_int"]),
            "b":       float(params["b"]),
        }
    if choice == "4":
        c_int = float(input(f"  Enter c_int for {species_name}: ").strip())
        b = float(input(f"  Enter b for {species_name}: ").strip())
        return {
            "species": species_name,
            "taxa":    "custom",
            "c_int":   c_int,
            "b":       b,
        }
    raise ValueError(f"Invalid choice {choice!r}; expected one of 1-4.")


# --- the science function -----------------------------------------------------
def calculate_massspecific_metabolism(
    mass_g: NDArray[np.float64],
    body_temp_C: NDArray[np.float64],
    c_int: NDArray[np.float64],
    b: NDArray[np.float64],
    ambient_temp_C: float | NDArray[np.float64] | None = None,
    E_a: float = E_A_DEFAULT,
    k_B: float = K_B_DEFAULT,
    fmr_multiplier: float = FMR_MULTIPLIER_DEFAULT,
) -> NDArray[np.float64]:
    """Field-level mass-specific metabolic rate, in W/g.

    All four array inputs must have the **same shape**. The trailing axis is
    species by convention (see processes spec §2). Per-species parameters
    (``mass_g``, ``c_int``, ``b``) of shape ``(S,)`` and a per-cell ambient
    temperature of shape ``(X, Y)`` should be broadcast up to the common shape
    by the caller (e.g. via ``np.broadcast_to``) before being passed in.

    Parameters
    ----------
    mass_g
        Individual body mass in grams.
    body_temp_C
        Body temperature in degrees Celsius. Entries that are ``NaN`` are
        treated as "use the ambient temperature": this is how ectotherms come
        out of ``body_temp``, and the function fills those slots from
        ``ambient_temp_C`` before evaluating the Arrhenius term.
    c_int
        PGLS intercept (``ln W``) from ``metabolism_parameters``.
    b
        Allometric exponent from ``metabolism_parameters``.
    ambient_temp_C
        Ambient temperature in degrees Celsius. Required only if
        ``body_temp_C`` contains any ``NaN`` (i.e. any ectotherms). May be a
        single scalar (one ambient temperature for all entries) or an array
        that broadcasts against ``body_temp_C`` (e.g. a per-cell map). If
        ``body_temp_C`` has no ``NaN`` entries, this argument is ignored.
    E_a
        Activation energy in eV (default 0.63).
    k_B
        Boltzmann constant in eV/K (default 8.617e-5).
    fmr_multiplier
        Factor applied to the resting rate to obtain field metabolic rate
        (default 3.0, after Brose et al. 2008). Pass ``1.0`` to get the
        resting rate.

    Returns
    -------
    NDArray[np.float64]
        Mass-specific metabolic rate in W/g, same shape as the inputs.
    """
    # Equal-shape assert (processes spec §3, rule 5). The science function
    # consumes already-broadcast inputs; the adapter is responsible for any
    # reshaping.
    assert mass_g.shape == body_temp_C.shape == c_int.shape == b.shape, (
        f"shape mismatch: mass_g={mass_g.shape}, body_temp_C={body_temp_C.shape}, "
        f"c_int={c_int.shape}, b={b.shape}"
    )

    # Fill NaN slots in body_temp_C from ambient. NaN is the sentinel that
    # body_temp() uses for ectotherms; we resolve it here so the user sees a
    # clear error if they forgot to pass ambient_temp_C.
    needs_ambient = np.any(np.isnan(body_temp_C))
    if needs_ambient:
        if ambient_temp_C is None:
            raise ValueError(
                "body_temp_C contains NaN entries (ectotherms) but "
                "ambient_temp_C was not provided. Pass ambient_temp_C as "
                "either a scalar or an array broadcasting against body_temp_C."
            )
        effective_temp_C = np.where(
            np.isnan(body_temp_C), ambient_temp_C, body_temp_C,
        )
    else:
        effective_temp_C = body_temp_C

    temp_K = effective_temp_C + KELVIN_OFFSET
    resting = np.exp(c_int) * mass_g ** (b - 1.0) * np.exp(-E_a / (k_B * temp_K))
    return fmr_multiplier * resting


# --- biomass loss from metabolism --------------------------------------------
def calculate_biomass_metabolism(
    initial_biomass_g: NDArray[np.float64],
    mass_specific_metabolic_rate: NDArray[np.float64],
    dt: float = 1.0,
    *,
    energy_density_J_per_g: float = ENERGY_DENSITY_J_PER_G_DEFAULT,
) -> NDArray[np.float64]:
    """Biomass delta (g, **signed**) from field metabolic loss over ``dt`` days.

    This is the biomass-modifying counterpart to
    ``calculate_massspecific_metabolism``. The split is deliberate: the
    rate is a dependency quantity (no ``dt``, reusable by ATN and dispersal),
    while this function applies the rate to a population's standing biomass
    over a time step.

    The math collapses to one line:

        delta_g = - initial_biomass_g
                  * mass_specific_metabolic_rate [W/g]
                  * (dt [days] * 86400 [s/day])
                  / energy_density_J_per_g

    Sign convention — note this differs from the R version. The function
    returns a **negative** number for the loss, matching the process spec
    where biomass_delta is the signed change the engine adds to biomass
    (``B[t+1] = B[t] + sum(deltas)``). The R version returned the positive
    magnitude; multiply by -1 to compare.

    Parameters
    ----------
    initial_biomass_g
        Population wet biomass at time ``t``, in grams. From ATN's biomass
        state in the engine.
    mass_specific_metabolic_rate
        Output of ``calculate_massspecific_metabolism``, in W/g. Same shape
        as ``initial_biomass_g``.
    dt
        Time step in days (engine convention). Last positional argument per
        the process spec.
    energy_density_J_per_g
        Wet-biomass energy density, J/g. Default 7000 J/g (= 7 MJ/kg, the
        figure used in the R prototype). Keyword-only because it is a
        physical constant most callers will leave alone.

    Returns
    -------
    NDArray[np.float64]
        Biomass delta in grams, same shape as the inputs. Always <= 0 for
        non-negative inputs.
    """
    assert initial_biomass_g.shape == mass_specific_metabolic_rate.shape, (
        f"shape mismatch: initial_biomass_g={initial_biomass_g.shape}, "
        f"mass_specific_metabolic_rate={mass_specific_metabolic_rate.shape}"
    )

    seconds = dt * SECONDS_PER_DAY
    loss_g = (
        initial_biomass_g * mass_specific_metabolic_rate * seconds
        / energy_density_J_per_g
    )
    return -loss_g
