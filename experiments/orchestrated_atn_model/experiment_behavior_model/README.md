# Behavior experiment: paradox of enrichment

This folder contains a small controlled experiment for testing whether the current spatial ATN implementation can express the paradox of enrichment.

The experiment uses a minimal food web:

```text
Plant guild -> Herbivore
```

It varies the plant guild carrying capacity `K_plant_0` across scenarios. Low or moderate enrichment should produce extinction or stable coexistence, while stronger enrichment should increase consumer-resource oscillation amplitude and create deeper biomass troughs.

Run from the repository root:

```bash
python3 experiments/atn_model/experiment_behavior_model/run_paradox_enrichment.py
```

Each run creates a timestamped folder:

```text
experiment_behavior_model/
└── paradox_enrichment_<timestamp>/
    ├── scenario_K_0010/
    │   ├── biomass.txt
    │   └── simulation_summary.txt
    ├── scenario_K_0025/
    │   ├── biomass.txt
    │   └── simulation_summary.txt
    ├── scenario_metrics.csv
    ├── paradox_enrichment_time_series.png
    ├── paradox_enrichment_phase_planes.png
    └── paradox_enrichment_metrics.png
```

The scenario-level files keep the same core output structure as `run_atn.py`: `biomass.txt` and `simulation_summary.txt`.
