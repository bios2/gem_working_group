"""
Paradox of enrichment experiment for the spatial ATN model.

This experiment uses the existing ATNModel implementation with a minimal
plant-guild -> herbivore food web. It varies plant carrying capacity K and
saves each scenario with the same core output files used by run_atn.py:

  biomass.txt
  simulation_summary.txt

It also writes aggregate plots and metrics in the experiment run folder.
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


ATN_DIR = Path(__file__).resolve().parents[1]
if str(ATN_DIR) not in sys.path:
    sys.path.insert(0, str(ATN_DIR))

from atn_model import ATNModel  # noqa: E402
from config import CONFIG  # noqa: E402


EXPERIMENT_NAME = "paradox_enrichment"
SPECIES_LABELS = {
    0: "Plant guild",
    1: "Herbivore",
}


def build_config() -> dict:
    """Return a controlled parameter set for a clear enrichment experiment."""
    cfg = CONFIG.copy()

    # These choices keep the experiment intentionally simple. We turn off
    # allometric and temperature variation so the only scenario driver is K.
    cfg.update(
        {
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
    return cfg


def build_traits() -> pd.DataFrame:
    """Build a minimal two-node food web: plant guild -> herbivore."""
    return pd.DataFrame(
        {
            "body_mass_g": [1.0, 100.0],
            "is_basal": [1, 0],
            "initial_biomass_g_per_m2": [8.0, 0.5],
        },
        index=pd.Index([0, 1], name="species_id"),
    )


def build_adjacency() -> np.ndarray:
    """Rows are resources and columns are consumers: herbivore eats plant."""
    return np.array(
        [
            [0, 1],
            [0, 0],
        ],
        dtype=int,
    )


def build_environment(k_plant: float) -> pd.DataFrame:
    """Create a one-cell environment with a scenario-specific plant K."""
    return pd.DataFrame(
        {
            "x": [0],
            "y": [0],
            "temperature_K": [293.15],
            "K_plant_0": [float(k_plant)],
        },
        index=pd.Index([0], name="pixel_id"),
    )


def initial_biomass(traits: pd.DataFrame, n_cells: int) -> np.ndarray:
    """Return deterministic initial biomass with shape cells x species."""
    b0 = np.zeros((n_cells, len(traits)))
    for i in range(len(traits)):
        b0[:, i] = traits.iloc[i]["initial_biomass_g_per_m2"]
    return b0


def scenario_name(k_plant: float) -> str:
    if float(k_plant).is_integer():
        return f"scenario_K_{int(k_plant):04d}"
    return f"scenario_K_{k_plant:g}".replace(".", "p")


def write_biomass(output_dir: Path, env_df: pd.DataFrame, t_eval: np.ndarray, b_traj: np.ndarray) -> None:
    """Write long-format biomass output matching run_atn.py."""
    n_tp, n_cells, n_species = b_traj.shape
    cell_x = env_df["x"].values.astype(int)
    cell_y = env_df["y"].values.astype(int)

    t_rep = np.repeat(t_eval, n_cells * n_species)
    cell_rep = np.tile(np.repeat(np.arange(n_cells), n_species), n_tp)
    x_rep = np.tile(np.repeat(cell_x, n_species), n_tp)
    y_rep = np.tile(np.repeat(cell_y, n_species), n_tp)
    sp_rep = np.tile(np.arange(n_species), n_tp * n_cells)
    bio_rep = b_traj.ravel()

    table = np.column_stack([cell_rep, x_rep, y_rep, t_rep, sp_rep, bio_rep])
    np.savetxt(
        output_dir / "biomass.txt",
        table,
        fmt=["%d", "%d", "%d", "%.4f", "%d", "%.6e"],
        header="pixel_id x y time_step species_id biomass",
        comments="",
    )


def scenario_metrics(k_plant: float, t_eval: np.ndarray, b_traj: np.ndarray, cfg: dict) -> dict:
    """Calculate simple enrichment diagnostics from the post-burn-in window."""
    burn_idx = int(0.5 * len(t_eval))
    tail = b_traj[burn_idx:, 0, :]
    tail_mean = tail.mean(axis=0)
    tail_min = tail.min(axis=0)
    tail_max = tail.max(axis=0)
    amplitude = tail_max - tail_min
    cv = tail.std(axis=0) / np.maximum(tail_mean, cfg["ext_threshold"])
    final = b_traj[-1, 0, :]

    return {
        "K_plant": k_plant,
        "plant_mean_tail": tail_mean[0],
        "herbivore_mean_tail": tail_mean[1],
        "plant_min_tail": tail_min[0],
        "herbivore_min_tail": tail_min[1],
        "plant_max_tail": tail_max[0],
        "herbivore_max_tail": tail_max[1],
        "plant_amplitude_tail": amplitude[0],
        "herbivore_amplitude_tail": amplitude[1],
        "plant_cv_tail": cv[0],
        "herbivore_cv_tail": cv[1],
        "plant_final": final[0],
        "herbivore_final": final[1],
        "plant_persisted": int(final[0] > cfg["ext_threshold"]),
        "herbivore_persisted": int(final[1] > cfg["ext_threshold"]),
    }


def write_summary(
    output_dir: Path,
    timestamp: str,
    k_plant: float,
    traits_df: pd.DataFrame,
    env_df: pd.DataFrame,
    adj_mat: np.ndarray,
    cfg: dict,
    t_eval: np.ndarray,
    metrics: dict,
) -> None:
    """Write a scenario summary with the same spirit as run_atn.py."""
    with open(output_dir / "simulation_summary.txt", "w", encoding="utf-8") as fsum:
        fsum.write("=" * 60 + "\n")
        fsum.write("SIMULATION SUMMARY\n")
        fsum.write(f"Experiment    : {EXPERIMENT_NAME}\n")
        fsum.write(f"Run timestamp : {timestamp}\n")
        fsum.write(f"Scenario      : plant carrying capacity K = {k_plant:g}\n")
        fsum.write("=" * 60 + "\n\n")

        fsum.write("SIMULATION DIMENSIONS\n")
        fsum.write("-" * 40 + "\n")
        fsum.write(f"Number of species   : {len(traits_df)}\n")
        fsum.write(f"Number of time steps: {len(t_eval)}\n")
        fsum.write(f"Simulation duration : {float(t_eval[-1]):.1f} days\n")
        fsum.write(f"Number of pixels    : {len(env_df)}\n")
        fsum.write("Grid dimensions     : 1 x value x 1 y value\n\n")

        fsum.write("EXPERIMENT DESIGN\n")
        fsum.write("-" * 40 + "\n")
        fsum.write("Food web: Plant guild -> Herbivore\n")
        fsum.write("Manipulated driver: K_plant_0, the carrying capacity of the plant guild.\n")
        fsum.write("Expected pattern: low or moderate K gives extinction or stable coexistence; high K increases oscillation amplitude and extinction risk.\n")
        fsum.write("Note: this script calls ATNModel directly because the current validator rejects this minimal resource-by-consumer matrix.\n\n")

        fsum.write("ENVIRONMENT\n")
        fsum.write("-" * 40 + "\n")
        fsum.write(env_df.to_csv())
        fsum.write("\n")

        fsum.write("ADJACENCY MATRIX (rows=resources, columns=consumers)\n")
        fsum.write("-" * 40 + "\n")
        for row in adj_mat:
            fsum.write(" ".join(str(int(x)) for x in row) + "\n")
        fsum.write("\n")

        fsum.write("SPECIES TRAITS\n")
        fsum.write("-" * 40 + "\n")
        fsum.write(traits_df.to_csv())
        fsum.write("\n")

        fsum.write("MODEL CONSTANTS (CONFIG)\n")
        fsum.write("-" * 40 + "\n")
        for key, value in cfg.items():
            fsum.write(f"  {key:<24} = {value}\n")
        fsum.write("\n")

        fsum.write("POST-BURN-IN METRICS\n")
        fsum.write("-" * 40 + "\n")
        for key, value in metrics.items():
            fsum.write(f"  {key:<24} = {value}\n")


def run_scenario(k_plant: float, run_dir: Path, timestamp: str, t_eval: np.ndarray) -> tuple[np.ndarray, dict]:
    """Run one K scenario and write scenario-level outputs."""
    cfg = build_config()
    traits_df = build_traits()
    adj_mat = build_adjacency()
    env_df = build_environment(k_plant)

    model = ATNModel(adj_mat, traits_df, env_df, cfg)
    b0 = initial_biomass(traits_df, len(env_df))
    b_traj = model.run_all_cells(b0, t_eval)
    b_traj = np.maximum(b_traj, 0.0)

    metrics = scenario_metrics(k_plant, t_eval, b_traj, cfg)
    output_dir = run_dir / scenario_name(k_plant)
    output_dir.mkdir(parents=True, exist_ok=True)
    write_biomass(output_dir, env_df, t_eval, b_traj)
    write_summary(output_dir, timestamp, k_plant, traits_df, env_df, adj_mat, cfg, t_eval, metrics)

    return b_traj, metrics


def plot_time_series(run_dir: Path, t_eval: np.ndarray, scenario_results: dict[float, np.ndarray]) -> None:
    """Plot biomass time series for all enrichment scenarios."""
    fig, axes = plt.subplots(len(scenario_results), 1, figsize=(10, 1.9 * len(scenario_results)), sharex=True)
    if len(scenario_results) == 1:
        axes = [axes]

    for ax, (k_plant, b_traj) in zip(axes, scenario_results.items()):
        ax.plot(t_eval, b_traj[:, 0, 0], color="#26734d", lw=1.3, label=SPECIES_LABELS[0])
        ax.plot(t_eval, b_traj[:, 0, 1], color="#7a3e9d", lw=1.3, label=SPECIES_LABELS[1])
        ax.set_ylabel(f"K={k_plant:g}\nbiomass")
        ax.grid(alpha=0.25)

    axes[0].legend(loc="upper right", frameon=False, ncol=2)
    axes[-1].set_xlabel("Time (days)")
    fig.suptitle("Paradox of enrichment: biomass time series across plant carrying capacity")
    fig.tight_layout(rect=[0, 0, 1, 0.98])
    fig.savefig(run_dir / "paradox_enrichment_time_series.png", dpi=180)
    plt.close(fig)


def plot_phase_planes(run_dir: Path, scenario_results: dict[float, np.ndarray]) -> None:
    """Plot plant-herbivore phase planes for the post-burn-in window."""
    n = len(scenario_results)
    n_cols = 3
    n_rows = int(np.ceil(n / n_cols))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(11, 3.2 * n_rows), squeeze=False)

    for ax in axes.ravel():
        ax.axis("off")

    for ax, (k_plant, b_traj) in zip(axes.ravel(), scenario_results.items()):
        tail = b_traj[int(0.5 * len(b_traj)) :, 0, :]
        ax.axis("on")
        ax.plot(tail[:, 0], tail[:, 1], color="#333333", lw=1.0)
        ax.scatter(tail[0, 0], tail[0, 1], s=16, color="#26734d", label="tail start")
        ax.scatter(tail[-1, 0], tail[-1, 1], s=16, color="#7a3e9d", label="final")
        ax.set_title(f"K={k_plant:g}")
        ax.set_xlabel("Plant biomass")
        ax.set_ylabel("Herbivore biomass")
        ax.grid(alpha=0.25)

    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", frameon=False)
    fig.suptitle("Post-burn-in phase planes")
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(run_dir / "paradox_enrichment_phase_planes.png", dpi=180)
    plt.close(fig)


def plot_metrics(run_dir: Path, metrics_df: pd.DataFrame) -> None:
    """Plot oscillation amplitude and persistence against enrichment."""
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))

    axes[0].plot(metrics_df["K_plant"], metrics_df["plant_amplitude_tail"], marker="o", color="#26734d", label=SPECIES_LABELS[0])
    axes[0].plot(metrics_df["K_plant"], metrics_df["herbivore_amplitude_tail"], marker="o", color="#7a3e9d", label=SPECIES_LABELS[1])
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Plant carrying capacity K")
    axes[0].set_ylabel("Post-burn-in amplitude")
    axes[0].set_title("Oscillation amplitude increases with enrichment")
    axes[0].grid(alpha=0.25)
    axes[0].legend(frameon=False)

    width = 0.35
    x = np.arange(len(metrics_df))
    axes[1].bar(x - width / 2, metrics_df["plant_min_tail"], width, color="#26734d", label=SPECIES_LABELS[0])
    axes[1].bar(x + width / 2, metrics_df["herbivore_min_tail"], width, color="#7a3e9d", label=SPECIES_LABELS[1])
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([f"{k:g}" for k in metrics_df["K_plant"]])
    axes[1].set_xlabel("Plant carrying capacity K")
    axes[1].set_ylabel("Post-burn-in minimum biomass")
    axes[1].set_title("Deep troughs indicate extinction risk")
    axes[1].grid(axis="y", alpha=0.25)
    axes[1].legend(frameon=False)

    fig.tight_layout()
    fig.savefig(run_dir / "paradox_enrichment_metrics.png", dpi=180)
    plt.close(fig)


def main() -> None:
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    run_dir = Path(__file__).resolve().parent / f"{EXPERIMENT_NAME}_{timestamp}"
    run_dir.mkdir(parents=True, exist_ok=True)

    k_scenarios = [10, 25, 50, 100, 200, 400, 800]
    t_eval = np.linspace(0, 1000, 1001)

    scenario_results = {}
    metrics = []
    for k_plant in k_scenarios:
        print(f"\nRunning {EXPERIMENT_NAME} scenario: K={k_plant:g}")
        b_traj, scenario_metric = run_scenario(k_plant, run_dir, timestamp, t_eval)
        scenario_results[float(k_plant)] = b_traj
        metrics.append(scenario_metric)

    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(run_dir / "scenario_metrics.csv", index=False)

    plot_time_series(run_dir, t_eval, scenario_results)
    plot_phase_planes(run_dir, scenario_results)
    plot_metrics(run_dir, metrics_df)

    print("\nExperiment complete.")
    print(f"Output folder: {run_dir}")
    print("Plots:")
    print(f"  {run_dir / 'paradox_enrichment_time_series.png'}")
    print(f"  {run_dir / 'paradox_enrichment_phase_planes.png'}")
    print(f"  {run_dir / 'paradox_enrichment_metrics.png'}")


if __name__ == "__main__":
    main()
