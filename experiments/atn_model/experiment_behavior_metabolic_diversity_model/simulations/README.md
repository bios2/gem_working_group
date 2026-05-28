# Simulation diagnostics

Run all checks from the repository root:

```bash
python3 experiments/atn_model/experiment_behavior_metabolic_diversity_model/simulations/run_model_checks.py
```

The script creates a timestamped folder named `results_YYYYMMDDHHMMSS`.

## Checks

`paradox_enrichment/` repeats the original behavior-model test with the updated metabolism. The herbivore is a 100 g mammal endotherm, so consumer metabolic loss comes from the new thermal-group equation. The expected result is the familiar enrichment pattern: low enrichment cannot maintain the consumer, intermediate enrichment supports stable coexistence, and high enrichment creates large oscillations and deep troughs.

`temperature_response/` isolates ectotherm temperature dependence by turning off generic ATN temperature scaling. Feeding and basal growth are held fixed, so changes across cells come from ectotherm metabolic loss. The rate comparison figure also shows that mammal and bird endotherm losses are constant across ambient temperature in this implementation.

`metabolic_scaling/` checks the metabolism equations directly across body masses. Mass-specific biomass loss decreases with body mass, while whole-body field metabolic rate increases with body mass.

## Generated figures in the current run

The current generated run is:

```text
results_20260527110652/
```

Key figures:

```text
results_20260527110652/paradox_enrichment/paradox_enrichment_time_series.png
results_20260527110652/paradox_enrichment/paradox_enrichment_metrics.png
results_20260527110652/temperature_response/ectotherm_temperature_response.png
results_20260527110652/temperature_response/thermal_group_rate_comparison.png
results_20260527110652/metabolic_scaling/metabolic_scaling_by_body_mass.png
```
