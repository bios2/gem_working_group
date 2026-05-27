"""
Main runner for the spatially explicit ATN model with comprehensive validation.

This script orchestrates the entire workflow:
1. Read input files (environment, adjacency, traits)
2. Validate all inputs (20+ checks)
3. Initialize the ATN model
4. Integrate ODEs forward in time
5. Save results and print summary
"""
import numpy as np
import pandas as pd
import subprocess
import sys
from pathlib import Path
from datetime import datetime
# Import custom modules from this project
from atn_io import (
    read_config,                                           # config reader
    read_env_matrix, read_adjacency_matrix, read_traits,  # I/O functions
    validate_inputs, check_parameter_completeness, print_summary,  # validation
    ValidationError  # custom exception
)
from atn_model import ATNModel  # the ODE model

def _git_commit_hash() -> str:
    """Return the current HEAD commit hash, or 'uncommitted' if the repo is dirty."""
    try:
        hash_ = subprocess.run(
            ['git', 'rev-parse', 'HEAD'],
            capture_output=True, text=True, check=True
        ).stdout.strip()
        dirty = subprocess.run(
            ['git', 'diff', '--quiet'],
            capture_output=True
        ).returncode != 0
        return hash_ + (' (dirty)' if dirty else '')
    except Exception:
        return 'unavailable'


def main(env_file: str, adj_file: str, traits_file: str,
         config_file: str = 'config.txt', t_max: float = 100.0,
         seed: int = 42):
    """
    Run the ATN model with full validation.

    Parameters:
        env_file: path to env_mat.txt (environment: row, col, temperature, carrying capacity)
        adj_file: path to adj_mat.txt (adjacency matrix: food web links)
        traits_file: path to traits.txt (species traits: body mass, initial biomass)
        config_file: path to config.txt (model parameters)
        t_max: simulation end time (days); one output point is saved per day
        seed: random seed for numpy (controls initial-biomass noise)
    """
    np.random.seed(seed)
    commit_hash = _git_commit_hash()

    # Print header banner
    print("=" * 70)
    print("SPATIALLY EXPLICIT ATN MODEL")
    print("=" * 70)

    try:
        # ===== STEP 1: READ AND VALIDATE INPUTS =====
        print("\n[STEP 1] Reading and validating input files...")
        CONFIG = read_config(config_file)
        # Read environment matrix (returns DataFrame and temperature array)
        env_df, temps = read_env_matrix(env_file)
        # Read adjacency matrix (binary food-web: rows=prey, cols=predators)
        adj_mat = read_adjacency_matrix(adj_file)
        # Read species traits (body mass, metabolism, initial biomass)
        traits_df = read_traits(traits_file)
        
        # Cross-file validation: check consistency across all inputs
        validate_inputs(env_df, adj_mat, traits_df)
        
        # ===== STEP 2: VALIDATE CONFIGURATION =====
        print("\n[STEP 2] Validating configuration...")
        # Check that all required parameters are present in config
        check_parameter_completeness(CONFIG)
        
        # ===== STEP 3: PRINT SUMMARY =====
        # Display input summary table
        print_summary(env_df, adj_mat, traits_df, CONFIG)
        
        # ===== STEP 4: INITIALIZE MODEL =====
        print("\n[STEP 3] Initializing model...")
        # Create ATN model instance (stores all data and parameters)
        model = ATNModel(adj_mat, traits_df, env_df, CONFIG)
        
        # ===== STEP 5: SET UP INITIAL CONDITIONS =====
        print("\n[STEP 4] Setting up initial conditions...")
        # Get dimensions
        n_cells = len(env_df)  # number of spatial cells
        n_species = len(traits_df)  # number of species
        
        # Initialize biomass matrix: (n_cells, n_species)
        B_initial = np.zeros((n_cells, n_species))
        # Fill in initial biomass from traits table for each species across all cells
        for i in range(n_species):
            B_initial[:, i] = traits_df.iloc[i]['initial_biomass_g_per_m2']
        
        # Add small random noise to break symmetry and prevent exact balances
        # Noise amplitude: 1% of biomass value (or 1% of 1 if biomass is 0)
        B_initial += np.random.normal(0, 0.01 * np.maximum(B_initial, 1e-6), B_initial.shape)
        # Clamp to non-negative (remove noise artifacts below zero)
        B_initial = np.maximum(B_initial, 0)
        
        # Print initial condition statistics
        print(f"  ✓ Initial biomass range: {B_initial[B_initial > 0].min():.2e} to "
              f"{B_initial.max():.2e} g/m²")
        print(f"  ✓ Total initial biomass: {B_initial.sum():.2e} g/m²")
        
        # ===== STEP 6: TIME INTEGRATION =====
        n_timepoints = int(t_max) + 1  # one output point per day, day 0 through day t_max
        print(f"\n[STEP 5] Running simulation for {t_max} days ({n_timepoints} output points)...")
        t_eval = np.linspace(0, t_max, n_timepoints)
        
        # Run ODE integration for all cells (returns full trajectory)
        B_traj = model.run_all_cells(B_initial, t_eval)
        
        # ===== STEP 7: SAVE RESULTS =====
        # Build output folder name from current timestamp
        timestamp = datetime.now().strftime('%Y%m%d%H%M%S')
        output_dir = Path('atn_output') / timestamp
        print(f"\n[STEP 6] Saving results to {output_dir}/...")

        # Extract spatial grid dimensions from env_df
        cell_x = env_df['x'].values.astype(int)
        cell_y = env_df['y'].values.astype(int)
        n_x = cell_x.max() + 1
        n_y = cell_y.max() + 1

        # ===== WRITE SIMULATION SUMMARY =====
        output_dir.mkdir(parents=True, exist_ok=True)
        with open(output_dir / 'simulation_summary.txt', 'w') as fsum:
            fsum.write("=" * 60 + "\n")
            fsum.write("SIMULATION SUMMARY\n")
            fsum.write(f"Run timestamp : {timestamp}\n")
            fsum.write(f"Random seed   : {seed}\n")
            fsum.write(f"Git commit    : {commit_hash}\n")
            fsum.write("=" * 60 + "\n\n")

            # --- Simulation dimensions ---
            fsum.write("SIMULATION DIMENSIONS\n")
            fsum.write("-" * 40 + "\n")
            fsum.write(f"Number of species  : {n_species}\n")
            fsum.write(f"Number of time steps: {len(t_eval)}\n")
            fsum.write(f"Simulation duration : {t_max:.1f} days\n")
            fsum.write(f"Number of pixels   : {n_cells}\n")
            fsum.write(f"Grid dimensions    : {n_x} x values x {n_y} y values\n\n")

            # --- Species traits ---
            fsum.write("SPECIES TRAITS\n")
            fsum.write("-" * 40 + "\n")
            col_w = [10, 15, 8, 25]
            header = (f"{'species_id':<{col_w[0]}}"
                      f"{'body_mass_g':>{col_w[1]}}"
                      f"{'is_basal':>{col_w[2]}}"
                      f"{'initial_biomass_g_per_m2':>{col_w[3]}}\n")
            fsum.write(header)
            fsum.write("-" * sum(col_w) + "\n")
            for sp_id, row in traits_df.iterrows():
                fsum.write(
                    f"{int(sp_id):<{col_w[0]}}"
                    f"{row['body_mass_g']:>{col_w[1]}.4f}"
                    f"{int(row['is_basal']):>{col_w[2]}}"
                    f"{row['initial_biomass_g_per_m2']:>{col_w[3]}.4f}\n"
                )
            fsum.write("\n")

            # --- Model constants ---
            fsum.write("MODEL CONSTANTS (CONFIG)\n")
            fsum.write("-" * 40 + "\n")
            descriptions = {
                # Vegetation growth (NPP-driven, vegetation.md equation)
                'psi':                 'Carbon to wet matter conversion factor (g wet / g C)',
                'f_struct_default':    'Default fractional NPP allocation to structural tissue',
                'alpha_herbs_default': 'Half-saturation constant for herb/tree partition (g/m²)',
                # Allometric rates
                'X0':                  'Metabolic loss rate normalization (day^-1)',
                'b_X':                 'Metabolic loss mass exponent',
                'a0':                  'Attack rate normalization (day^-1)',
                'b_a_prey':            'Attack rate prey mass exponent',
                'b_a_pred':            'Attack rate predator mass exponent',
                'h0':                  'Handling time normalization (days)',
                'b_h_prey':            'Handling time prey mass exponent',
                'b_h_pred':            'Handling time predator mass exponent',
                # Functional response
                'q_hill':              'Hill exponent for functional response',
                'interference':        'Consumer interference coefficient',
                # Efficiency
                'e_plant':             'Assimilation efficiency — plant prey',
                'e_animal':            'Assimilation efficiency — animal prey',
                # Temperature
                'use_temperature':     'Apply temperature scaling to rates',
                'T0_K':                'Reference temperature (K)',
                'k_B':                 'Boltzmann constant (eV/K)',
                'E_a':                 'Activation energy (eV)',
                # Numerical
                'ext_threshold':       'Extinction threshold (g/m²)',
                'extinction_timescale':'Decay timescale for extinct species (days)',
            }
            for key, val in CONFIG.items():
                desc = descriptions.get(key, '')
                fsum.write(f"  {key:<22} = {str(val):<12}  {desc}\n")
            fsum.write("\n")

        print(f"  ✓ Saved simulation_summary.txt")

        model.vegetation.save_output(B_traj, t_eval, output_dir)
        print(f"  ✓ Saved vegetation.txt")
        model.save_consumer_output(B_traj, t_eval, output_dir)
        print(f"  ✓ Saved atn_model.txt")

        # ===== STEP 8: PRINT SUMMARY STATISTICS =====
        print("\n[STEP 7] Summary statistics:")
        
        # Compute final average biomass per species (averaged across all cells)
        final_biomass = B_traj[-1, :, :].mean(axis=0)
        print(f"\n  Final mean biomass per species (g/m²):")
        # Print first 10 species
        for i in range(min(n_species, 10)):
            # Label as basal or consumer
            spp_name = f"Spp {i} (basal)" if traits_df.iloc[i]['is_basal'] else f"Spp {i}"
            print(f"    {spp_name:20s}: {final_biomass[i]:.4e}")
        # Indicate if there are more species
        if n_species > 10:
            print(f"    ... and {n_species - 10} more")
        
        # Compute persistence: fraction of cells where each species survives
        # Extinction threshold from config determines "survives"
        persistence = np.sum(B_traj[-1, :, :] > CONFIG['ext_threshold'], axis=0) / n_cells * 100
        print(f"\n  Persistence (% cells where species survives > {CONFIG['ext_threshold']:.0e}):")
        for i in range(min(n_species, 10)):
            spp_name = f"Spp {i} (basal)" if traits_df.iloc[i]['is_basal'] else f"Spp {i}"
            print(f"    {spp_name:20s}: {persistence[i]:6.1f}%")
        if n_species > 10:
            print(f"    ... and {n_species - 10} more")
        
        # Count species that went locally extinct in at least one cell
        n_extinct = np.sum(persistence < 1.0)
        print(f"\n  Locally extinct in some cells: {n_extinct} species")
        
        # ===== COMPLETION =====
        print("\n" + "=" * 70)
        print("✓ SIMULATION COMPLETE")
        print("=" * 70 + "\n")
        
        # Return results for further analysis
        return B_traj, t_eval, model
    
    except ValidationError as e:
        # Handle validation errors with informative output
        print(f"\n{'='*70}")
        print("✗ VALIDATION ERROR")
        print("="*70)
        print(f"\n{e}\n")
        sys.exit(1)  # exit with error code
    
    except Exception as e:
        # Handle unexpected errors with traceback
        print(f"\n{'='*70}")
        print("✗ UNEXPECTED ERROR")
        print("="*70)
        print(f"\n{type(e).__name__}: {e}\n")
        # Print full traceback for debugging
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Run the spatially explicit ATN model.')
    parser.add_argument('env_file',    nargs='?', default='env_mat.txt')
    parser.add_argument('adj_file',    nargs='?', default='adj_mat.txt')
    parser.add_argument('traits_file', nargs='?', default='traits.txt')
    parser.add_argument('config_file', nargs='?', default='config.txt')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for initial-biomass noise (default: 42)')
    parser.add_argument('--t_max', type=float, default=100.0)
    args = parser.parse_args()

    B_traj, t_eval, model = main(
        args.env_file, args.adj_file, args.traits_file,
        config_file=args.config_file, t_max=args.t_max, seed=args.seed
    )
