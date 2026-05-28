"""
Run diagnostic simulations for the metabolic-diversity ATN model.

Outputs are written to a timestamped results folder under this directory:
  simulations/results_<timestamp>/
"""

from __future__ import annotations

import os

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), ".mplconfig"))

from datetime import datetime
from pathlib import Path
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


MODEL_DIR = Path(__file__).resolve().parents[1]
if str(MODEL_DIR) not in sys.path:
    sys.path.insert(0, str(MODEL_DIR))

from atn_model import ATNModel  # noqa: E402
from config import CONFIG  # noqa: E402
from metabolism import calculate_metabolic_loss_per_day, calculate_metabolism_w_per_g  # noqa: E402


def controlled_config(**overrides) -> dict:
    """Return a simple, deterministic parameter set for diagnostic runs."""
    cfg = CONFIG.copy()
    cfg.update(
        {
            "metabolism_model": "thermal_group",
            "use_temperature": False,
            "r0": 1.0,
            "b_r": 0.0,
            "X0": 0.2,
            "b_X": 0.0,
            "a0": 0.02,
            "b_a_prey": 0.0,
            "b_a_pred": 0.0,
            "h0": 0.5,
            "b_h_prey": 0.0,
            "b_h_pred": 0.0,
            "q_hill": 1.0,
            "interference": 0.0,
            "e_plant": 0.45,
            "e_animal": 0.85,
            "ext_threshold": 1e-9,
            "extinction_timescale": 0.1,
        }
    )
    cfg.update(overrides)
    return cfg


def build_two_species_traits(consumer_group: str, consumer_mass_g: float, endotherm_group: str = "mammal") -> pd.DataFrame:
    return pd.DataFrame(
        {
            "body_mass_g": [1.0, consumer_mass_g],
            "is_basal": [1, 0],
            "initial_biomass_g_per_m2": [8.0, 0.5],
            "thermal_group": ["atn", consumer_group],
            "endotherm_group": ["", endotherm_group],
        },
        index=pd.Index([0, 1], name="species_id"),
    )


def build_adjacency() -> np.ndarray:
    return np.array([[0, 1], [0, 0]], dtype=int)


def initial_biomass(traits: pd.DataFrame, n_cells: int) -> np.ndarray:
    b0 = np.zeros((n_cells, len(traits)))
    for i in range(len(traits)):
        b0[:, i] = traits.iloc[i]["initial_biomass_g_per_m2"]
    return b0


def run_model(adj_mat: np.ndarray, traits_df: pd.DataFrame, env_df: pd.DataFrame, cfg: dict, t_eval: np.ndarray) -> np.ndarray:
    model = ATNModel(adj_mat, traits_df, env_df, cfg)
    return np.maximum(model.run_all_cells(initial_biomass(traits_df, len(env_df)), t_eval), 0.0)


def tail_metrics(k_plant: float, t_eval: np.ndarray, b_traj: np.ndarray, cfg: dict) -> dict:
    tail = b_traj[int(0.5 * len(t_eval)) :, 0, :]
    mean = tail.mean(axis=0)
    min_ = tail.min(axis=0)
    max_ = tail.max(axis=0)
    amp = max_ - min_
    return {
        "K_plant": k_plant,
        "plant_mean_tail": mean[0],
        "consumer_mean_tail": mean[1],
        "plant_min_tail": min_[0],
        "consumer_min_tail": min_[1],
        "plant_max_tail": max_[0],
        "consumer_max_tail": max_[1],
        "plant_amplitude_tail": amp[0],
        "consumer_amplitude_tail": amp[1],
        "plant_cv_tail": tail[:, 0].std() / max(mean[0], cfg["ext_threshold"]),
        "consumer_cv_tail": tail[:, 1].std() / max(mean[1], cfg["ext_threshold"]),
        "plant_final": b_traj[-1, 0, 0],
        "consumer_final": b_traj[-1, 0, 1],
    }


def run_paradox_enrichment(run_dir: Path) -> None:
    out_dir = run_dir / "paradox_enrichment"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = controlled_config()
    traits = build_two_species_traits("endotherm", consumer_mass_g=100.0, endotherm_group="mammal")
    adj = build_adjacency()
    t_eval = np.linspace(0, 1000, 1001)
    k_scenarios = [10, 25, 50, 100, 200, 400, 800]

    results = {}
    metrics = []
    for k in k_scenarios:
        env = pd.DataFrame(
            {"x": [0], "y": [0], "temperature_K": [293.15], "K_plant_0": [float(k)]},
            index=pd.Index([0], name="pixel_id"),
        )
        b_traj = run_model(adj, traits, env, cfg, t_eval)
        results[float(k)] = b_traj
        metrics.append(tail_metrics(k, t_eval, b_traj, cfg))

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(out_dir / "scenario_metrics.csv", index=False)

    fig, axes = plt.subplots(len(results), 1, figsize=(10, 1.85 * len(results)), sharex=True)
    for ax, (k, b_traj) in zip(axes, results.items()):
        ax.plot(t_eval, b_traj[:, 0, 0], color="#2f6f4e", lw=1.2, label="Plant guild")
        ax.plot(t_eval, b_traj[:, 0, 1], color="#7048a8", lw=1.2, label="Endotherm herbivore")
        ax.set_ylabel(f"K={k:g}\nbiomass")
        ax.grid(alpha=0.25)
    axes[0].legend(loc="upper right", frameon=False, ncol=2)
    axes[-1].set_xlabel("Time (days)")
    fig.suptitle("Paradox of enrichment retained after metabolic-diversity update")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(out_dir / "paradox_enrichment_time_series.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    axes[0].plot(metrics_df["K_plant"], metrics_df["plant_amplitude_tail"], marker="o", color="#2f6f4e", label="Plant")
    axes[0].plot(metrics_df["K_plant"], metrics_df["consumer_amplitude_tail"], marker="o", color="#7048a8", label="Consumer")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Plant carrying capacity K")
    axes[0].set_ylabel("Post-burn-in amplitude")
    axes[0].set_title("Oscillation amplitude")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)
    axes[1].plot(metrics_df["K_plant"], metrics_df["consumer_min_tail"], marker="o", color="#7048a8")
    axes[1].set_xscale("log")
    axes[1].set_xlabel("Plant carrying capacity K")
    axes[1].set_ylabel("Consumer minimum biomass")
    axes[1].set_title("Deep troughs with enrichment")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "paradox_enrichment_metrics.png", dpi=180)
    plt.close(fig)

    with open(out_dir / "README.md", "w", encoding="utf-8") as f:
        f.write("# Paradox of enrichment check\n\n")
        f.write("This diagnostic repeats the original plant-guild -> herbivore enrichment test using the updated model.\n")
        f.write("The herbivore is a 100 g mammal endotherm, so its metabolic loss is calculated by the new thermal-group metabolism.\n")


def run_temperature_response(run_dir: Path) -> None:
    out_dir = run_dir / "temperature_response"
    out_dir.mkdir(parents=True, exist_ok=True)
    cfg = controlled_config(a0=0.018, h0=0.5)
    adj = build_adjacency()
    t_eval = np.linspace(0, 650, 651)
    temperatures_c = np.array([8, 14, 20, 26, 32], dtype=float)
    env = pd.DataFrame(
        {
            "x": np.arange(len(temperatures_c)),
            "y": np.zeros(len(temperatures_c), dtype=int),
            "temperature_K": temperatures_c + 273.15,
            "K_plant_0": np.repeat(100.0, len(temperatures_c)),
        },
        index=pd.Index(np.arange(len(temperatures_c)), name="pixel_id"),
    )
    traits = build_two_species_traits("ectotherm", consumer_mass_g=0.1)
    b_traj = run_model(adj, traits, env, cfg, t_eval)

    records = []
    for cell_idx, temp_c in enumerate(temperatures_c):
        tail = b_traj[int(0.5 * len(t_eval)) :, cell_idx, :]
        records.append(
            {
                "temperature_C": temp_c,
                "ectotherm_X_per_day": calculate_metabolic_loss_per_day(0.1, temp_c + 273.15, "ectotherm"),
                "plant_mean_tail": tail[:, 0].mean(),
                "ectotherm_mean_tail": tail[:, 1].mean(),
                "ectotherm_final": b_traj[-1, cell_idx, 1],
            }
        )
    metrics = pd.DataFrame(records)
    metrics.to_csv(out_dir / "temperature_response_metrics.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for cell_idx, temp_c in enumerate(temperatures_c):
        axes[0].plot(t_eval, b_traj[:, cell_idx, 1], lw=1.2, label=f"{temp_c:g} C")
    axes[0].set_xlabel("Time (days)")
    axes[0].set_ylabel("Ectotherm herbivore biomass")
    axes[0].set_title("Ambient temperature changes ectotherm dynamics")
    axes[0].grid(alpha=0.25)
    axes[0].legend(title="Ambient", frameon=False, ncol=2)
    axes[1].plot(metrics["temperature_C"], metrics["ectotherm_X_per_day"], marker="o", color="#b35c2e")
    axes[1].set_xlabel("Ambient temperature (C)")
    axes[1].set_ylabel("Metabolic loss X (day^-1)")
    axes[1].set_title("Ectotherm metabolism increases with temperature")
    axes[1].grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_dir / "ectotherm_temperature_response.png", dpi=180)
    plt.close(fig)

    endo_mammal = [
        calculate_metabolic_loss_per_day(0.1, temp_c + 273.15, "endotherm", endotherm_group="mammal")
        for temp_c in temperatures_c
    ]
    endo_bird = [
        calculate_metabolic_loss_per_day(0.1, temp_c + 273.15, "endotherm", endotherm_group="bird")
        for temp_c in temperatures_c
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    ax.plot(temperatures_c, metrics["ectotherm_X_per_day"], marker="o", color="#b35c2e", label="Ectotherm")
    ax.plot(temperatures_c, endo_mammal, marker="s", color="#3f6fb5", label="Mammal endotherm")
    ax.plot(temperatures_c, endo_bird, marker="^", color="#6f4c9b", label="Bird endotherm")
    ax.set_xlabel("Ambient temperature (C)")
    ax.set_ylabel("Metabolic loss X (day^-1)")
    ax.set_title("Endotherm loss is independent of ambient temperature in this model")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "thermal_group_rate_comparison.png", dpi=180)
    plt.close(fig)

    with open(out_dir / "README.md", "w", encoding="utf-8") as f:
        f.write("# Temperature-response check\n\n")
        f.write("This run turns off generic ATN temperature scaling, so feeding and plant growth are fixed.\n")
        f.write("Only ectotherm metabolic loss changes with ambient temperature.\n")


def run_metabolic_scaling_check(run_dir: Path) -> None:
    out_dir = run_dir / "metabolic_scaling"
    out_dir.mkdir(parents=True, exist_ok=True)
    masses = np.logspace(-2, 6, 220)
    temp_k = 293.15
    rows = []
    for mass in masses:
        for group, endo_group in [("ectotherm", None), ("endotherm", "mammal"), ("endotherm", "bird")]:
            label = group if group == "ectotherm" else f"{endo_group} endotherm"
            x_day = calculate_metabolic_loss_per_day(mass, temp_k, group, endotherm_group=endo_group)
            w_per_g = calculate_metabolism_w_per_g(mass, temp_k, group, endotherm_group=endo_group)
            rows.append(
                {
                    "mass_g": mass,
                    "group": label,
                    "metabolic_loss_per_day": x_day,
                    "fmr_w_per_g": w_per_g,
                    "whole_body_fmr_w": w_per_g * mass,
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "metabolic_scaling_rates.csv", index=False)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3))
    colors = {"ectotherm": "#b35c2e", "mammal endotherm": "#3f6fb5", "bird endotherm": "#6f4c9b"}
    for group, data in df.groupby("group"):
        axes[0].plot(data["mass_g"], data["metabolic_loss_per_day"], lw=1.8, color=colors[group], label=group)
        axes[1].plot(data["mass_g"], data["whole_body_fmr_w"], lw=1.8, color=colors[group], label=group)
    axes[0].set_xscale("log")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Individual body mass (g)")
    axes[0].set_ylabel("Mass-specific biomass loss X (day^-1)")
    axes[0].set_title("Small individuals lose more per gram")
    axes[0].grid(alpha=0.25, which="both")
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Individual body mass (g)")
    axes[1].set_ylabel("Whole-body FMR (W)")
    axes[1].set_title("Whole-body metabolism increases with size")
    axes[1].grid(alpha=0.25, which="both")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_dir / "metabolic_scaling_by_body_mass.png", dpi=180)
    plt.close(fig)

    with open(out_dir / "README.md", "w", encoding="utf-8") as f:
        f.write("# Metabolic-scaling check\n\n")
        f.write("Rates are calculated directly from the new metabolism module at 20 C ambient temperature.\n")
        f.write("Mass-specific loss decreases with body mass, while whole-body FMR increases with body mass.\n")


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    run_dir = Path(__file__).resolve().parent / f"results_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)

    run_paradox_enrichment(run_dir)
    run_temperature_response(run_dir)
    run_metabolic_scaling_check(run_dir)

    with open(run_dir / "README.md", "w", encoding="utf-8") as f:
        f.write("# Metabolic-diversity model checks\n\n")
        f.write("Generated by `simulations/run_model_checks.py`.\n\n")
        f.write("- `paradox_enrichment/`: repeats the behavior-model enrichment test with thermal-group metabolism.\n")
        f.write("- `temperature_response/`: shows ectotherm metabolism and dynamics changing with ambient temperature.\n")
        f.write("- `metabolic_scaling/`: checks expected mass-scaling patterns for endotherms and ectotherms.\n")

    print(f"Simulation checks complete: {run_dir}")


if __name__ == "__main__":
    main()
