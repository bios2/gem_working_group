# Model changes: adding metabolic diversity

## Summary

The original behavior model used one generic allometric metabolic-loss rate for all species:

$$
X_i = X_0 M_i^{b_X} f_T
$$

The updated model keeps the same ATN biomass equations and feeding response, but lets consumer metabolic loss depend on whether a species is an ectotherm or endotherm. The new loss coefficient is still used in the same ATN term:

$$
\text{metabolic loss}_i = X_i B_i
$$

So the structure of the model is unchanged; only how `X_i` is calculated changes for species with `thermal_group = ectotherm` or `thermal_group = endotherm`.

## Before

The old model calculated metabolic loss as:

$$
X_i(T_g)
=
X_0 M_i^{b_X}
\exp\left[
\frac{-E_a(T_0 - T_g)}{kT_gT_0}
\right]
$$

where:

- `X0` is a normalization constant in day^-1.
- `M_i` is individual body mass in g.
- `bX` is the mass exponent, usually near `-0.25`.
- `T_g` is local environmental temperature in K.
- `T0` is the reference temperature in K.
- `E_a` is activation energy in eV.
- `k` is Boltzmann's constant.

This means all consumers responded to environmental temperature in the same way.

## New metabolism

The new metabolism follows the R function in `experiment_diversity_metabolism/calculate_metabolism.r`.

Whole-body metabolic rate:

$$
X_{\text{whole},i}\;[\mathrm{W}]
=
\exp(C)M_i^b
\exp\left(
\frac{-E_a}{kT_i}
\right)
$$

Mass-specific metabolic rate:

$$
x_i\;[\mathrm{W\,g^{-1}}]
=
\exp(C)M_i^{b-1}
\exp\left(
\frac{-E_a}{kT_i}
\right)
$$

Field metabolic rate:

$$
\mathrm{FMR}_i\;[\mathrm{W\,g^{-1}}] = 3x_i
$$

ATN metabolic-loss coefficient:

$$
X_i\;[\mathrm{day^{-1}}]
=
\mathrm{FMR}_i
\frac{86400\;\mathrm{s\,day^{-1}}}
{7000\;\mathrm{J\,g^{-1}\;wet\;biomass}}
$$

This gives `g biomass lost per g biomass per day`, which is exactly the coefficient needed by the ATN loss term `X_i B_i`.

## Thermal groups

Ectotherms use ambient temperature:

$$
T_i = T_g
$$

$$
C_{\mathrm{ecto}} = 17.4,
\qquad
b_{\mathrm{ecto}} = 0.84
$$

Endotherms use fixed body temperature:

$$
T_i =
\begin{cases}
36.5 + 273.15, & \text{mammals}\\
41.5 + 273.15, & \text{birds}
\end{cases}
$$

$$
C_{\mathrm{endo}} = 19.53,
\qquad
b_{\mathrm{endo}} = 0.73
$$

The shared constants are:

$$
E_a = 0.63\;\mathrm{eV},
\qquad
k = 8.617\times10^{-5}\;\mathrm{eV\,K^{-1}}
$$

## Biomass equations

Basal species still follow:

$$
\frac{dB_i}{dt}
=
r_iB_i\left(1-\frac{B_i}{K_i}\right)
- X_iB_i
- \sum_j B_jF_{ij}
$$

Consumers still follow:

$$
\frac{dB_j}{dt}
=
B_j\sum_i e_iF_{ij}
- X_jB_j
- \sum_p B_pF_{jp}
$$

The functional response is unchanged:

$$
F_{ij}
=
\frac{a_{ij}B_i^q}
{1 + c_jB_j + \sum_k h_{kj}a_{kj}B_k^q}
$$

The updated part is `X_j` for consumers:

$$
X_j
=
X_{\mathrm{thermal}}\left(
M_j,\;T_g,\;\mathrm{thermal\_group}_j,\;\mathrm{endotherm\_group}_j
\right)
$$

Basal species use the original ATN metabolism by default. This is controlled by:

```text
apply_thermal_metabolism_to_basal = False
```

## Where the code changed

The metabolism equations are implemented in `metabolism.py`:

- `body_temperature_k()` chooses the temperature used by metabolism: ambient temperature for ectotherms, fixed mammal or bird body temperature for endotherms.
- `calculate_metabolism_w_per_g()` implements:

$$
x_i = \exp(C)M_i^{b-1}\exp\left(\frac{-E_a}{kT_i}\right)\times 3
$$

- `calculate_metabolic_loss_per_day()` converts `W/g` into the ATN coefficient:

$$
X_i = x_i\frac{86400}{7000}
$$

The ATN model uses those functions in `atn_model.py`, inside `_metabolic_rate()`. That function is the key integration point:

- If `metabolism_model != "thermal_group"`, it returns the old ATN metabolism:

$$
X_i = X_0M_i^{b_X}f_T
$$

- If `metabolism_model == "thermal_group"`, it first calculates the old ATN metabolism for all species, then replaces consumer `X_i` values with the endotherm/ectotherm values from `metabolism.py`.
- Basal species keep the old ATN metabolism unless `apply_thermal_metabolism_to_basal = True`.

The new metabolism enters the biomass dynamics in `derivatives()`, where the model computes:

$$
\text{loss}_j = X_jB_j + \text{loss to predators}
$$

for consumers, and:

$$
\text{loss}_i = X_iB_i + \text{loss to consumers}
$$

for basal species.

## What the change affects

The new endotherm/ectotherm implementation affects only the metabolic-loss coefficient `X_i`.

It does not directly change:

- attack rates `a_ij`;
- handling times `h_ij`;
- the functional response `F_ij`;
- assimilation efficiencies `e_i`;
- basal carrying capacity `K_i`;
- the food-web adjacency matrix.

Important nuance: the model still has the older config switch `use_temperature`. When `use_temperature = True`, the existing ATN temperature factor still affects basal growth, attack rates, handling times, and legacy ATN metabolism. The new thermal-group metabolism has its own temperature handling, so ectotherm/endotherm metabolic loss does not use the old reference-temperature correction.

In the diagnostic `temperature_response` simulation, `use_temperature = False`. That isolates the new metabolism: plant growth and feeding are fixed, and only ectotherm metabolic loss changes with ambient temperature.

## Trait-table additions

The updated model accepts these optional columns in the species traits table:

```text
thermal_group
endotherm_group
```

Examples:

```text
thermal_group,endotherm_group
ectotherm,
endotherm,mammal
endotherm,bird
atn,
```

If a consumer has no `thermal_group`, the model uses:

```text
default_consumer_thermal_group = "ectotherm"
```

## Config additions

New config fields:

```text
metabolism_model = "thermal_group"
apply_thermal_metabolism_to_basal = False
default_consumer_thermal_group = "ectotherm"
default_endotherm_group = "mammal"
field_metabolic_multiplier = 3.0
joules_per_g_wet_biomass = 7000.0
seconds_per_day = 86400.0
metabolic_rate_multiplier = 1.0
```

To recover the old model behavior, set:

```text
metabolism_model = "atn"
```

## Diagnostic simulations

The folder `simulations/` contains `run_model_checks.py`, which creates three checks:

1. `paradox_enrichment`: repeats the original enrichment experiment with the updated model. The herbivore is a mammal endotherm, so its metabolic loss comes from the new thermal-group equation.
2. `temperature_response`: turns off generic ATN temperature scaling and shows that ectotherm metabolic loss and dynamics change with ambient temperature.
3. `metabolic_scaling`: directly checks that mass-specific metabolic loss decreases with body mass, while whole-body metabolism increases with body mass.

## Why the 8 C ectotherm can rise and then crash

In the temperature-response diagnostic, the cold ectotherm has low metabolic loss:

$$
X_{\mathrm{ecto}}(8^\circ\mathrm{C})
<
X_{\mathrm{ecto}}(20^\circ\mathrm{C})
<
X_{\mathrm{ecto}}(32^\circ\mathrm{C})
$$

Because `use_temperature = False` in that diagnostic, low temperature does not slow attack rate, handling time, or plant growth. The 8 C ectotherm therefore has a special combination:

- low metabolic cost;
- normal feeding rate;
- normal plant-growth rate.

That can let the consumer increase from the initial plant resource, overshoot, and overexploit the plant. Once plant biomass is driven very low, consumer gains collapse:

$$
B_j\sum_i e_iF_{ij} \approx 0
$$

but consumer losses continue:

$$
\frac{dB_j}{dt} \approx -X_jB_j
$$

So the consumer crashes after the resource is depleted. Since `X_j` is small at 8 C, the decline can be slow rather than immediate.

This result should be interpreted as a focused metabolism test, not as a full cold-temperature ecology prediction. If we want cold ectotherms to also feed, move, digest, or attack more slowly, then attack and handling rates should also be made thermal-group-specific, or the old `use_temperature` scaling should be used carefully alongside the new metabolism.
