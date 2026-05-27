"""
Quick validation test of the ATN model.

This script runs a short simulation with example input files to verify that:
1. All input files are found and readable
2. Model initializes without errors
3. ODE integration runs and completes
4. Output is generated correctly

Use this to validate your setup before running long simulations.
"""
# Import custom modules
from atn_io import (
    read_env_matrix, read_adjacency_matrix, read_traits,  # I/O functions
    validate_inputs, check_parameter_completeness  # validation
)
from atn_model import ATNModel  # the model
from config import CONFIG  # parameters
import numpy as np  # arrays
import sys  # system utilities

def test_atn():
    """Run a quick test simulation with example files."""
    
    # Print header
    print("\n" + "="*60)
    print("ATN MODEL QUICK VALIDATION TEST")
    print("="*60)
    
    try:
        # ===== TEST 1: READ INPUT FILES =====
        print("\n[TEST 1] Reading example input files...")
        # Read three example files that must exist in the current directory
        env_df, temps = read_env_matrix('example_env_mat.txt')
        print(f"  ✓ Environment matrix: {len(env_df)} cells")
        
        adj_mat = read_adjacency_matrix('example_adj_mat.txt')
        print(f"  ✓ Adjacency matrix: {adj_mat.shape} ({adj_mat.sum()} links)")
        
        traits_df = read_traits('example_traits.txt')
        print(f"  ✓ Traits table: {len(traits_df)} species")
        
        # ===== TEST 2: VALIDATE INPUTS =====
        print("\n[TEST 2] Validating cross-file consistency...")
        validate_inputs(env_df, adj_mat, traits_df)
        print("  ✓ All validation checks passed")
        
        # ===== TEST 3: VALIDATE CONFIG =====
        print("\n[TEST 3] Checking configuration...")
        check_parameter_completeness(CONFIG)
        print("  ✓ All {0} configuration parameters present".format(len(CONFIG)))
        
        # ===== TEST 4: INITIALIZE MODEL =====
        print("\n[TEST 4] Initializing ATN model...")
        model = ATNModel(adj_mat, traits_df, env_df, CONFIG)
        print(f"  ✓ Model ready: {model.n_species} species, {model.n_cells} cells")
        
        # ===== TEST 5: SET UP INITIAL CONDITIONS =====
        print("\n[TEST 5] Setting up initial biomass...")
        # Create initial biomass matrix
        B_initial = np.zeros((model.n_cells, model.n_species))
        # Fill from traits
        for i in range(model.n_species):
            B_initial[:, i] = traits_df.iloc[i]['initial_biomass_g_per_m2']
        # Add small noise
        B_initial += np.random.normal(0, 0.01 * np.maximum(B_initial, 1e-6), B_initial.shape)
        B_initial = np.maximum(B_initial, 0)
        print(f"  ✓ Initial biomass: min={B_initial[B_initial>0].min():.2e}, ", end="")
        print(f"max={B_initial.max():.2e} g/m²")
        
        # ===== TEST 6: RUN SHORT SIMULATION =====
        print("\n[TEST 6] Running 10-day simulation (light test)...")
        # Short time span for quick testing
        t_eval = np.linspace(0, 10, 11)  # 11 points from 0 to 10 days
        # Run model
        B_traj = model.run_all_cells(B_initial, t_eval)
        print(f"  ✓ Simulation complete: output shape {B_traj.shape}")
        
        # ===== TEST 7: CHECK OUTPUT =====
        print("\n[TEST 7] Verifying output...")
        # Check output shape: should be (n_timepoints, n_cells, n_species)
        assert B_traj.shape == (len(t_eval), model.n_cells, model.n_species), \
            f"Output shape mismatch: {B_traj.shape}"
        print(f"  ✓ Output shape: {B_traj.shape}")
        
        # Check for NaN values (numerical instability)
        n_nan = np.sum(np.isnan(B_traj))
        assert n_nan == 0, f"Found {n_nan} NaN values in output"
        print(f"  ✓ No NaN values detected")
        
        # Check for negative biomass (should not happen with integrator settings)
        n_neg = np.sum(B_traj < 0)
        if n_neg > 0:
            print(f"  ⚠ Warning: {n_neg} negative biomass values (expected ~0)")
        else:
            print(f"  ✓ All biomass values non-negative")
        
        # Print final biomass statistics
        final_B = B_traj[-1, :, :].mean(axis=0)
        print(f"\n  Final biomass per species:")
        print(f"    Min: {final_B[final_B > 0].min():.2e} g/m²")
        print(f"    Max: {final_B.max():.2e} g/m²")
        print(f"    Median: {np.median(final_B[final_B > 0]):.2e} g/m²")
        
        # ===== ALL TESTS PASSED =====
        print("\n" + "="*60)
        print("✓ ALL TESTS PASSED - Model is working correctly!")
        print("="*60 + "\n")
        
        return True
    
    except AssertionError as e:
        # Output validation failed
        print(f"\n✗ TEST FAILED: {e}\n")
        return False
    
    except Exception as e:
        # Other errors (file not found, etc.)
        print(f"\n✗ ERROR: {type(e).__name__}: {e}\n")
        import traceback
        traceback.print_exc()
        return False

if __name__ == '__main__':
    # Run the test and exit with appropriate code
    success = test_atn()
    sys.exit(0 if success else 1)
