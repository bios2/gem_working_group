"""
Command-line runner for the metabolic-diversity ATN model.

Example from the repository root:
  python3 experiments/atn_model/experiment_behavior_metabolic_diversity_model/run_atn.py \
    --env experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/env_mat.csv \
    --adj experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/adj_mat.txt \
    --traits experiments/atn_model/experiment_behavior_metabolic_diversity_model/example_inputs/traits.csv \
    --t-max 365
"""

from __future__ import annotations

import argparse
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

from atn_model import ATNModel
from config import CONFIG


def read_inputs(env_file: Path, adj_file: Path, traits_file: Path) -> tuple[pd.DataFrame, np.ndarray, pd.DataFrame]:
    env_df = pd.read_csv(env_file, index_col=0)
    adj_mat = np.loadtxt(adj_file, dtype=int)
    traits_df = pd.read_csv(traits_file, index_col=0)

    required_env = {"x", "y", "temperature_K"}
    required_traits = {"body_mass_g", "is_basal", "initial_biomass_g_per_m2"}
    missing_env = required_env - set(env_df.columns)
    missing_traits = required_traits - set(traits_df.columns)
    if missing_env:
        raise ValueError(f"Environment file missing columns: {sorted(missing_env)}")
    if missing_traits:
        raise ValueError(f"Traits file missing columns: {sorted(missing_traits)}")
    if adj_mat.shape != (len(traits_df), len(traits_df)):
        raise ValueError("Adjacency matrix shape must match the number of trait rows")
    return env_df, adj_mat, traits_df


def initial_biomass(traits_df: pd.DataFrame, n_cells: int) -> np.ndarray:
    b0 = np.zeros((n_cells, len(traits_df)))
    for i in range(len(traits_df)):
        b0[:, i] = traits_df.iloc[i]["initial_biomass_g_per_m2"]
    return b0


def write_biomass(output_dir: Path, env_df: pd.DataFrame, t_eval: np.ndarray, b_traj: np.ndarray) -> None:
    n_tp, n_cells, n_species = b_traj.shape
    cell_x = env_df["x"].values.astype(int)
    cell_y = env_df["y"].values.astype(int)

    table = np.column_stack(
        [
            np.tile(np.repeat(np.arange(n_cells), n_species), n_tp),
            np.tile(np.repeat(cell_x, n_species), n_tp),
            np.tile(np.repeat(cell_y, n_species), n_tp),
            np.repeat(t_eval, n_cells * n_species),
            np.tile(np.arange(n_species), n_tp * n_cells),
            b_traj.ravel(),
        ]
    )
    np.savetxt(
        output_dir / "biomass.txt",
        table,
        fmt=["%d", "%d", "%d", "%.4f", "%d", "%.6e"],
        header="pixel_id x y time_step species_id biomass",
        comments="",
    )


def write_summary(output_dir: Path, env_df: pd.DataFrame, traits_df: pd.DataFrame, adj_mat: np.ndarray, cfg: dict, t_eval: np.ndarray) -> None:
    with open(output_dir / "simulation_summary.txt", "w", encoding="utf-8") as f:
        f.write("METABOLIC-DIVERSITY ATN SIMULATION\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Number of species: {len(traits_df)}\n")
        f.write(f"Number of cells: {len(env_df)}\n")
        f.write(f"Duration: {float(t_eval[-1]):.1f} days\n")
        f.write(f"Metabolism model: {cfg.get('metabolism_model')}\n\n")
        f.write("Traits\n")
        f.write("-" * 60 + "\n")
        f.write(traits_df.to_csv())
        f.write("\nEnvironment\n")
        f.write("-" * 60 + "\n")
        f.write(env_df.to_csv())
        f.write("\nAdjacency matrix, rows=resources, columns=consumers\n")
        f.write("-" * 60 + "\n")
        for row in adj_mat:
            f.write(" ".join(str(int(v)) for v in row) + "\n")
        f.write("\nConfig\n")
        f.write("-" * 60 + "\n")
        for key, value in cfg.items():
            f.write(f"{key}: {value}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--env", required=True, type=Path)
    parser.add_argument("--adj", required=True, type=Path)
    parser.add_argument("--traits", required=True, type=Path)
    parser.add_argument("--t-max", default=365.0, type=float)
    parser.add_argument("--output-dir", default=None, type=Path)
    args = parser.parse_args()

    env_df, adj_mat, traits_df = read_inputs(args.env, args.adj, args.traits)
    t_eval = np.linspace(0, args.t_max, int(args.t_max) + 1)
    model = ATNModel(adj_mat, traits_df, env_df, CONFIG.copy())
    b_traj = np.maximum(model.run_all_cells(initial_biomass(traits_df, len(env_df)), t_eval), 0.0)

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = Path(__file__).resolve().parent / "atn_output" / datetime.now().strftime("%Y%m%d%H%M%S")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_biomass(output_dir, env_df, t_eval, b_traj)
    write_summary(output_dir, env_df, traits_df, adj_mat, CONFIG, t_eval)
    print(f"Saved results to {output_dir}")


if __name__ == "__main__":
    main()
