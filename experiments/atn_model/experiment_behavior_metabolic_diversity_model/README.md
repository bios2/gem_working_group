# ATN model with metabolic diversity

This folder is a self-contained updated version of the behavior-model ATN. It keeps the same consumer-resource equations used for the paradox-of-enrichment experiment, but replaces the single generic consumer metabolic-loss allometry with explicit thermal groups:

- `ectotherm`: body temperature equals local ambient temperature.
- `endotherm`: body temperature is fixed by `endotherm_group` (`mammal` or `bird`).
- `atn`, `legacy`, `basal`, `plant`, or missing basal values: use the original ATN metabolic loss.

Main files:

- `atn_model.py`: ODE model.
- `metabolism.py`: Python version of the R metabolism functions from `experiment_diversity_metabolism`.
- `config.py`: model parameters and metabolic-diversity switches.
- `run_atn.py`: command-line runner for custom inputs.
- `MODEL_CHANGES.md`: formulas and explanation of what changed.
- `simulations/run_model_checks.py`: diagnostic simulations and figures.

Run all diagnostic checks from the repository root:

```bash
python3 experiments/atn_model/experiment_behavior_metabolic_diversity_model/simulations/run_model_checks.py
```

Run the example model:

```bash
python3 experiments/atn_model/experiment_behavior_metabolic_diversity_model/run_atn.py \
  --env experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/env_mat.csv \
  --adj experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/adj_mat.txt \
  --traits experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/traits.csv \
  --t-max 365
```

The traits table must include:

- `body_mass_g`
- `is_basal`
- `initial_biomass_g_per_m2`

For metabolic diversity, add:

- `thermal_group`: `ectotherm`, `endotherm`, or `atn`
- `endotherm_group`: `mammal` or `bird` for endotherms
