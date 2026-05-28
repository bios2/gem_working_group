"""
Configuration for the ATN model with metabolic diversity.

Most parameters are inherited conceptually from the behavior-model ATN. The
new options under "metabolic diversity" control whether consumer metabolic
loss is calculated from thermal group, body mass, and temperature.
"""

CONFIG = {
    # Basal growth
    "r0": 0.5,
    "b_r": -0.25,

    # Legacy ATN metabolic loss. Used for basal species by default, and for
    # all species when metabolism_model == "atn".
    "X0": 0.5,
    "b_X": -0.25,

    # Feeding
    "a0": 0.001,
    "b_a_prey": -0.5,
    "b_a_pred": 0.5,
    "h0": 0.01,
    "b_h_prey": 0.5,
    "b_h_pred": -0.5,
    "q_hill": 2.0,
    "interference": 0.0,
    "R_opt": 100.0,
    "gamma": 2.0,
    "link_threshold": 0.01,

    # Assimilation
    "e_plant": 0.45,
    "e_animal": 0.85,

    # Carrying capacity
    "K_default": 100.0,

    # Generic temperature scaling for basal growth and feeding rates.
    # Metabolic-diversity metabolism has its own temperature handling.
    "use_temperature": True,
    "T0_K": 293.15,
    "k_B": 8.617e-5,
    "E_a": 0.65,

    # Metabolic diversity
    # "thermal_group": consumers use Blyth/Gillooly-style endotherm and
    # ectotherm metabolism; basal species use legacy ATN metabolism.
    # "atn": all species use X0 * M^b_X with optional generic temperature.
    "metabolism_model": "thermal_group",
    "apply_thermal_metabolism_to_basal": False,
    "default_consumer_thermal_group": "ectotherm",
    "default_endotherm_group": "mammal",
    "field_metabolic_multiplier": 3.0,
    "joules_per_g_wet_biomass": 7000.0,
    "seconds_per_day": 86400.0,
    "metabolic_rate_multiplier": 1.0,

    # Numerical
    "ext_threshold": 1e-6,
    "extinction_timescale": 0.1,
}
