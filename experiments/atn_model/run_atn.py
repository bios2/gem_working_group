"""
Main runner for the spatially explicit ATN model with comprehensive validation.

This script orchestrates the entire workflow:
1. Read input files (environment, adjacency, traits)
2. Validate all inputs (20+ checks)
3. Initialize the ATN model
4. Integrate ODEs forward in time
5. Save results and print summary
"""
# Import required libraries
import numpy as np  # numerical arrays and math
import pandas as pd  # data frame handling
import sys  # system utilities (exit codes)
from pathlib import Path  # cross-platform path handling
# Import custom modules from this project
from atn_io import (
    read_env_matrix, read_adjacency_matrix, read_traits,  # I/O functions
    validate_inputs, check_parameter_completeness, print_summary,  # validation
    ValidationError  # custom exception
)
from atn_model import ATNModel  # the ODE model
from config import CONFIG  # default parameters

def main(env_file: str, adj_file: str, traits_file: str, 
         t_max: float = 100.0, n_timesteps: int = 1000, 
         output_dir: str = './atn_output'):
    """
    Run the ATN model with full validation.
    
    Parameters:
        env_file: path to env_mat.txt (environment: temperature, carrying capacity)
        adj_file: path to adj_mat.txt (adjacency matrix: food web links)
        traits_file: path to traits.txt (species traits: body mass, initial biomass)
        t_max: simulation end time (days)
        n_timesteps: number of output timesteps (determines output resolution)
        output_dir: directory to save results
    """
    
    # Print header banner
    print("=" * 70)
    print("SPATIALLY EXPLICIT ATN MODEL")
    print("=" * 70)
    
    try:
        # ===== STEP 1: READ AND VALIDATE INPUTS =====
        print("\n[STEP 1] Reading and validating input files...")
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
        print(f\"  ✓ Initial biomass range: {B_initial[B_initial > 0].min():.2e} to \"
              f\"{B_initial.max():.2e} g/m²\")
        print(f\"  ✓ Total initial biomass: {B_initial.sum():.2e} g/m²\")
        
        # ===== STEP 6: TIME INTEGRATION =====
        print(f\"\\n[STEP 5] Running simulation for {t_max} days ({n_timesteps} timesteps)...\")
        # Create time vector: linearly spaced points from 0 to t_max
        t_eval = np.linspace(0, t_max, n_timesteps)
        
        # Run ODE integration for all cells (returns full trajectory)
        B_traj = model.run_all_cells(B_initial, t_eval)
        
        # ===== STEP 7: SAVE RESULTS =====
        print(f\"\\n[STEP 6] Saving results to {output_dir}...\")
        # Create output directory if it doesn't exist
        Path(output_dir).mkdir(exist_ok=True)
        
        # Save biomass trajectory as NumPy binary file (efficient storage)
        np.save(f\"{output_dir}/biomass_trajectory.npy\", B_traj)
        # Save time points for reference
        np.save(f\"{output_dir}/time_points.npy\", t_eval)
        
        # ===== STEP 8: PRINT SUMMARY STATISTICS =====
        print(\"\\n[STEP 7] Summary statistics:\")
        
        # Compute final average biomass per species (averaged across all cells)
        final_biomass = B_traj[-1, :, :].mean(axis=0)
        print(f\"\\n  Final mean biomass per species (g/m²):\")
        # Print first 10 species
        for i in range(min(n_species, 10)):
            # Label as basal or consumer
            spp_name = f\"Spp {i} (basal)\" if traits_df.iloc[i]['is_basal'] else f\"Spp {i}\"
            print(f\"    {spp_name:20s}: {final_biomass[i]:.4e}\")
        # Indicate if there are more species
        if n_species > 10:
            print(f\"    ... and {n_species - 10} more\")
        
        # Compute persistence: fraction of cells where each species survives
        # Extinction threshold from config determines \"survives\"
        persistence = np.sum(B_traj[-1, :, :] > CONFIG['ext_threshold'], axis=0) / n_cells * 100
        print(f\"\\n  Persistence (% cells where species survives > {CONFIG['ext_threshold']:.0e}):\")
        for i in range(min(n_species, 10)):
            spp_name = f\"Spp {i} (basal)\" if traits_df.iloc[i]['is_basal'] else f\"Spp {i}\"
            print(f\"    {spp_name:20s}: {persistence[i]:6.1f}%\")
        if n_species > 10:
            print(f\"    ... and {n_species - 10} more\")
        
        # Count species that went locally extinct in at least one cell
        n_extinct = np.sum(persistence < 1.0)
        print(f\"\\n  Locally extinct in some cells: {n_extinct} species\")
        
        # ===== COMPLETION =====
        print(\"\\n\" + \"=\" * 70)
        print(\"✓ SIMULATION COMPLETE\")
        print(\"=\" * 70 + \"\\n\")
        
        # Return results for further analysis
        return B_traj, t_eval, model
    
    except ValidationError as e:
        # Handle validation errors with informative output
        print(f\"\\n{'='*70}\")
        print(\"✗ VALIDATION ERROR\")
        print(\"=\"*70)
        print(f\"\\n{e}\\n\")
        sys.exit(1)  # exit with error code
    
    except Exception as e:
        # Handle unexpected errors with traceback
        print(f\"\\n{'='*70}\")
        print(\"✗ UNEXPECTED ERROR\")
        print(\"=\"*70)
        print(f\"\\n{type(e).__name__}: {e}\\n\")
        # Print full traceback for debugging
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    # Allow command-line arguments for filenames
    if len(sys.argv) > 1:
        # Use provided filenames
        env_f = sys.argv[1]
        adj_f = sys.argv[2]
        traits_f = sys.argv[3]
    else:
        # Use default filenames
        env_f = 'env_mat.txt'
        adj_f = 'adj_mat.txt'
        traits_f = 'traits.txt'
    
    # Run the model
    B_traj, t_eval, model = main(env_f, adj_f, traits_f, t_max=100, n_timesteps=100)
