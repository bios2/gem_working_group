"""
Metabolic-rate helpers for endotherms and ectotherms.

This is a Python implementation of the R functions in
experiments/experiment_diversity_metabolism/calculate_metabolism.r and
calculate_biomass_metabolism.r.
"""

from __future__ import annotations

import numpy as np


BOLTZMANN_EV_PER_K = 8.617e-5
ACTIVATION_ENERGY_EV = 0.63

THERMAL_GROUP_PARAMETERS = {
    "ectotherm": {"intercept": 17.4, "mass_exponent": 0.84},
    "endotherm": {"intercept": 19.53, "mass_exponent": 0.73},
}

ENDOTHERM_BODY_TEMPERATURE_C = {
    "bird": 41.5,
    "mammal": 36.5,
}


def body_temperature_k(
    thermal_group: str,
    ambient_temperature_k: float,
    endotherm_group: str | None = None,
) -> float:
    """Return the temperature used by the metabolism equation."""
    thermal_group = str(thermal_group).lower()
    if thermal_group == "ectotherm":
        return float(ambient_temperature_k)
    if thermal_group == "endotherm":
        endotherm_group = "mammal" if endotherm_group is None else str(endotherm_group).lower()
        if endotherm_group not in ENDOTHERM_BODY_TEMPERATURE_C:
            raise ValueError("endotherm_group must be 'bird' or 'mammal'")
        return ENDOTHERM_BODY_TEMPERATURE_C[endotherm_group] + 273.15
    raise ValueError("thermal_group must be 'ectotherm' or 'endotherm'")


def calculate_metabolism_w_per_g(
    mass_g,
    ambient_temperature_k,
    thermal_group,
    endotherm_group=None,
    field_metabolic_multiplier: float = 3.0,
):
    """
    Calculate mass-specific field metabolic rate in W/g.

    Formula:
      x = exp(C) * M^(b - 1) * exp(-Ea / (kT)) * field_multiplier
    """
    mass = np.asarray(mass_g, dtype=float)
    if np.any(mass <= 0):
        raise ValueError("mass_g must be positive")

    params = THERMAL_GROUP_PARAMETERS[str(thermal_group).lower()]
    temp_k = body_temperature_k(thermal_group, ambient_temperature_k, endotherm_group)
    resting = (
        np.exp(params["intercept"])
        * np.power(mass, params["mass_exponent"] - 1.0)
        * np.exp(-ACTIVATION_ENERGY_EV / (BOLTZMANN_EV_PER_K * temp_k))
    )
    return resting * field_metabolic_multiplier


def calculate_metabolic_loss_per_day(
    mass_g,
    ambient_temperature_k,
    thermal_group,
    endotherm_group=None,
    field_metabolic_multiplier: float = 3.0,
    seconds_per_day: float = 86400.0,
    joules_per_g_wet_biomass: float = 7000.0,
):
    """
    Convert W/g metabolism to the ATN loss coefficient X_i in day^-1.

    W/g is J/s/g. Multiplying by seconds/day and dividing by J/g wet biomass
    gives g biomass lost per g biomass per day, matching the ATN term X_i B_i.
    """
    w_per_g = calculate_metabolism_w_per_g(
        mass_g,
        ambient_temperature_k,
        thermal_group,
        endotherm_group=endotherm_group,
        field_metabolic_multiplier=field_metabolic_multiplier,
    )
    return w_per_g * seconds_per_day / joules_per_g_wet_biomass


def calculate_population_biomass_loss_per_day(
    biomass_g,
    mass_g,
    ambient_temperature_k,
    thermal_group,
    endotherm_group=None,
    field_metabolic_multiplier: float = 3.0,
    seconds_per_day: float = 86400.0,
    joules_per_g_wet_biomass: float = 7000.0,
):
    """Return g biomass/day lost by a population with total biomass_g."""
    return biomass_g * calculate_metabolic_loss_per_day(
        mass_g,
        ambient_temperature_k,
        thermal_group,
        endotherm_group=endotherm_group,
        field_metabolic_multiplier=field_metabolic_multiplier,
        seconds_per_day=seconds_per_day,
        joules_per_g_wet_biomass=joules_per_g_wet_biomass,
    )
