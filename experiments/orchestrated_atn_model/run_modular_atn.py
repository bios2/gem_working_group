"""
ATN model runner — modular version.

This file contains NO science.
It simply chains:
  1. load inputs       (via atn_data.load_inputs)
  2. precompute params (via atn_data.build_atn_params)
  3. build B_initial   (via atn_data.initial_biomass)
  4. integrate dB/dt with RK4, calling atn_processes.derivatives
  5. save outputs

Architecture:
    atn_processes.py    — pure functions (numpy only)
    atn_data.py         — input reading + parameter precomputation
    run_modular_atn.py  — orchestration (this file)
"""
from datetime import datetime
from pathlib import Path
import subprocess
import sys

import numpy as np
from numpy.typing import NDArray

from atn_data import load_inputs, build_atn_params, initial_biomass
from atn_processes import derivatives
from atn_io import ValidationError


# --------------------------------------------------------------------------- #
# Helpers                                                                      #
# --------------------------------------------------------------------------- #
def _git_commit_hash() -> str:
    """Current HEAD commit hash (+ '(dirty)' if the repo has changes), or 'unavailable'."""
    try:
        h = subprocess.run(['git', 'rev-parse', 'HEAD'],
                           capture_output=True, text=True, check=True).stdout.strip()
        dirty = subprocess.run(['git', 'diff', '--quiet'],
                               capture_output=True).returncode != 0
        return h + (' (dirty)' if dirty else '')
    except Exception:
        return 'unavailable'


# --------------------------------------------------------------------------- #
# Vectorised RK4 integration (all cells in parallel)                           #
# --------------------------------------------------------------------------- #
def rk4_step(B: NDArray[np.float64], t: float, dt: float,
             params: dict, vegetation) -> NDArray[np.float64]:
    """One 4th-order RK step, vectorised over (n_cells, n_species)."""
    k1 = derivatives(B,                                    t,              params, vegetation)
    k2 = derivatives(np.maximum(B + 0.5 * dt * k1, 0.0),   t + 0.5 * dt,   params, vegetation)
    k3 = derivatives(np.maximum(B + 0.5 * dt * k2, 0.0),   t + 0.5 * dt,   params, vegetation)
    k4 = derivatives(np.maximum(B + dt * k3,       0.0),   t + dt,         params, vegetation)
    return np.maximum(B + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4), 0.0)


def run_simulation(B_initial: NDArray[np.float64], t_eval: NDArray[np.float64],
                   params: dict, vegetation) -> NDArray[np.float64]:
    """
    RK4 loop over the time points in t_eval.

    Returns: B_traj of shape (len(t_eval), n_cells, n_species).
    """
    n_tp = len(t_eval)
    B_traj = np.zeros((n_tp, params['n_cells'], params['n_species']))
    B = B_initial.astype(float)
    B_traj[0] = B

    for i in range(1, n_tp):
        dt = float(t_eval[i] - t_eval[i - 1])
        t  = float(t_eval[i - 1])
        B = rk4_step(B, t, dt, params, vegetation)
        B_traj[i] = B
        print(f"  t = {t_eval[i]:.1f} / {t_eval[-1]:.1f} days", end='\r')

    print(f"\n✓ {n_tp - 1} RK4 steps across {params['n_cells']} cells.")
    return B_traj


# --------------------------------------------------------------------------- #
# Outputs                                                                       #
# --------------------------------------------------------------------------- #
def save_biomass_long(B_traj: NDArray[np.float64], t_eval: NDArray[np.float64],
                      env_df, output_dir: Path) -> None:
    """Long-format table: pixel_id, x, y, time_step, species_id, biomass."""
    n_tp, n_cells, n_species = B_traj.shape
    cell_x = env_df['x'].values.astype(int)
    cell_y = env_df['y'].values.astype(int)

    t_rep    = np.repeat(t_eval, n_cells * n_species)
    cell_rep = np.tile(np.repeat(np.arange(n_cells), n_species), n_tp)
    x_rep    = np.tile(np.repeat(cell_x, n_species), n_tp)
    y_rep    = np.tile(np.repeat(cell_y, n_species), n_tp)
    sp_rep   = np.tile(np.arange(n_species), n_tp * n_cells)
    bio_rep  = B_traj.ravel()

    table = np.column_stack([cell_rep, x_rep, y_rep, t_rep, sp_rep, bio_rep])
    np.savetxt(output_dir / 'biomass.txt', table,
               fmt=['%d', '%d', '%d', '%.4f', '%d', '%.6e'],
               header='pixel_id x y time_step species_id biomass',
               comments='')


def save_consumer_dbdt(B_traj: NDArray[np.float64], t_eval: NDArray[np.float64],
                       env_df, params: dict, vegetation,
                       output_dir: Path) -> None:
    """
    Recompute dB/dt at the saved time points, for consumer species only.
    Output columns: pixel_id, x, y, time, species, delta_biomass.
    """
    n_tp    = len(t_eval)
    cons    = params['consumer_idx']
    n_cells = params['n_cells']
    n_cons  = len(cons)

    dBdt_arr = np.zeros((n_tp, n_cells, n_cons))
    for t_idx in range(n_tp):
        dBdt = derivatives(B_traj[t_idx], float(t_eval[t_idx]), params, vegetation)
        dBdt_arr[t_idx] = dBdt[:, cons]

    cell_x = env_df['x'].values.astype(int)
    cell_y = env_df['y'].values.astype(int)

    t_rep    = np.repeat(t_eval, n_cells * n_cons)
    cell_rep = np.tile(np.repeat(np.arange(n_cells), n_cons), n_tp)
    x_rep    = np.tile(np.repeat(cell_x, n_cons), n_tp)
    y_rep    = np.tile(np.repeat(cell_y, n_cons), n_tp)
    sp_rep   = np.tile(cons, n_tp * n_cells)
    d_rep    = dBdt_arr.ravel()

    table = np.column_stack([cell_rep, x_rep, y_rep, t_rep, sp_rep, d_rep])
    np.savetxt(output_dir / 'atn_model.txt', table,
               fmt=['%d', '%d', '%d', '%.4f', '%d', '%.6e'],
               header='pixel_id x y time species delta_biomass',
               comments='')


def write_summary(output_dir: Path, timestamp: str, seed: int, commit: str,
                  t_max: float, t_eval, env_df, traits_df, config) -> None:
    """Minimal text summary (seed, commit hash, dimensions, key parameters)."""
    n_cells   = len(env_df)
    n_species = len(traits_df)
    n_x = int(env_df['x'].max()) + 1
    n_y = int(env_df['y'].max()) + 1

    with open(output_dir / 'simulation_summary.txt', 'w') as f:
        f.write("=" * 60 + "\n")
        f.write("ATN MODULAR — SIMULATION SUMMARY\n")
        f.write(f"Run timestamp : {timestamp}\n")
        f.write(f"Random seed   : {seed}\n")
        f.write(f"Git commit    : {commit}\n")
        f.write("=" * 60 + "\n\n")

        f.write("DIMENSIONS\n" + "-" * 40 + "\n")
        f.write(f"Species         : {n_species}\n")
        f.write(f"Time steps      : {len(t_eval)}\n")
        f.write(f"Duration (days) : {t_max:.1f}\n")
        f.write(f"Cells           : {n_cells}\n")
        f.write(f"Grid (x × y)    : {n_x} × {n_y}\n\n")

        f.write("CONFIG\n" + "-" * 40 + "\n")
        for k, v in config.items():
            f.write(f"  {k:<22} = {v}\n")


# --------------------------------------------------------------------------- #
# Orchestration                                                                 #
# --------------------------------------------------------------------------- #
def main(env_file: str, adj_file: str, traits_file: str,
         config_file: str = 'config.txt', t_max: float = 100.0, seed: int = 42):
    """
    Full pipeline: read → precompute → integrate → save.

    Returns: (B_traj, t_eval, params, vegetation)
    """
    np.random.seed(seed)
    commit = _git_commit_hash()

    print("=" * 70)
    print("ATN MODULAR")
    print("=" * 70)

    try:
        # ------------------------------------------------------------------ #
        # 1) INPUTS                                                            #
        # ------------------------------------------------------------------ #
        print("\n[1] Reading and validating inputs...")
        config, env_df, adj_mat, traits_df = load_inputs(
            env_file, adj_file, traits_file, config_file
        )

        # ------------------------------------------------------------------ #
        # 2) PRECOMPUTED PARAMETERS                                            #
        # ------------------------------------------------------------------ #
        print("[2] Building ATN parameters...")
        params, vegetation = build_atn_params(adj_mat, traits_df, env_df, config)

        # ------------------------------------------------------------------ #
        # 3) INITIAL CONDITIONS                                                #
        # ------------------------------------------------------------------ #
        print("[3] Initial conditions...")
        B_initial = initial_biomass(traits_df, params['n_cells'])
        positive = B_initial[B_initial > 0]
        print(f"    ✓ Initial biomass: {positive.min():.2e} – {B_initial.max():.2e} g/m²")
        print(f"    ✓ Total biomass  : {B_initial.sum():.2e} g/m²")

        # ------------------------------------------------------------------ #
        # 4) INTEGRATION                                                       #
        # ------------------------------------------------------------------ #
        n_tp = int(t_max) + 1
        t_eval = np.linspace(0.0, t_max, n_tp)
        print(f"\n[4] RK4 integration over {t_max:.0f} days ({n_tp} points)...")
        B_traj = run_simulation(B_initial, t_eval, params, vegetation)

        # ------------------------------------------------------------------ #
        # 5) SAVE                                                              #
        # ------------------------------------------------------------------ #
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        out_dir   = Path('atn_output') / timestamp
        out_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n[5] Saving to {out_dir}/")

        write_summary(out_dir, timestamp, seed, commit, t_max,
                      t_eval, env_df, traits_df, config)
        print("    ✓ simulation_summary.txt")

        save_biomass_long(B_traj, t_eval, env_df, out_dir)
        print("    ✓ biomass.txt")

        vegetation.save_output(B_traj, t_eval, out_dir)
        print("    ✓ vegetation.txt")

        save_consumer_dbdt(B_traj, t_eval, env_df, params, vegetation, out_dir)
        print("    ✓ atn_model.txt")

        # ------------------------------------------------------------------ #
        # 6) QUICK STATISTICS                                                  #
        # ------------------------------------------------------------------ #
        print("\n[6] Final statistics:")
        final = B_traj[-1, :, :].mean(axis=0)
        persist = np.sum(B_traj[-1] > config['ext_threshold'], axis=0) / params['n_cells'] * 100
        for i in range(min(params['n_species'], 10)):
            label = f"Species {i} (basal)" if params['is_basal'][i] else f"Species {i}"
            print(f"    {label:<22}  B_final={final[i]:.3e}  persistence={persist[i]:5.1f}%")
        if params['n_species'] > 10:
            print(f"    ... {params['n_species'] - 10} more species")

        print("\n" + "=" * 70)
        print("✓ SIMULATION COMPLETE")
        print("=" * 70)
        return B_traj, t_eval, params, vegetation

    except ValidationError as e:
        print(f"\n✗ VALIDATION ERROR: {e}")
        sys.exit(1)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        description="Modular ATN model (processes + data + run separated)."
    )
    parser.add_argument('env_file',    nargs='?', default='env_mat.txt')
    parser.add_argument('adj_file',    nargs='?', default='adj_mat.txt')
    parser.add_argument('traits_file', nargs='?', default='traits.txt')
    parser.add_argument('config_file', nargs='?', default='config.txt')
    parser.add_argument('--seed',  type=int,   default=42,
                        help="numpy seed for initial-biomass noise (default: 42)")
    parser.add_argument('--t_max', type=float, default=100.0,
                        help="Simulation duration in days (default: 100)")
    args = parser.parse_args()

    main(args.env_file, args.adj_file, args.traits_file,
         config_file=args.config_file, t_max=args.t_max, seed=args.seed)
