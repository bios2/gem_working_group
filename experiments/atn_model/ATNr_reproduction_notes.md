# ATNr Figure 3 and Figure 4 Reproduction Notes

This workflow compares the paper examples against the current local Python
implementation without replacing the local model equations.

## What Was Reproduced

- ATNr Figure 3 recipe: temperature gradient from 4 to 22 deg C, 50 species,
  20 basal species, 2 nutrients, body-mass-based `L` matrix, and ATNr
  `create_model_Unscaled_nuts()`.
- ATNr Figure 4 recipe: 10 species, niche-model food web, body masses from
  trophic levels, and the ATNr scaled model with `K = 1` and `K = 10`.
- Python analogue: the same food webs, body masses, initial species biomasses,
  temperatures, and carrying-capacity values were passed to `ATNModel`.

## Important Scope Differences

- The current Python model has no explicit nutrient state variables, so Figure 3
  cannot be exactly reproduced yet.
- The current Python model is unscaled, while the Figure 4 example code uses
  `create_model_Scaled()`.
- The Python Figure 3 run saved here uses a shorter horizon
  (`t = 1000`) and explicit Euler integration for tractability. The ATNr
  reference values use the paper's `t = 100000`.

## Current Outputs

- `atnr_reproduction_output/figure3_temperature_extinctions_python_vs_atnr.png`
- `atnr_reproduction_output/figure4_enrichment_python_vs_atnr.png`
- `atnr_reproduction_output/figure3_python_vs_atnr_extinctions.csv`
- `atnr_reproduction_output/reproduction_notes.json`

## Formula-Alignment Issues To Check Before Claiming Reproduction

These are not fixed by the reproduction script, because the goal was to test the
current implementation rather than overwrite it with ATNr formulas.

- ATNr matrices use rows as resources and columns as consumers. In
  `ATNModel._allometric_rate()` calls, `M_prey` and `M_pred` are broadcast in
  the opposite row/column orientation before being multiplied by `adj_mat`.
- `ATNModel._L_matrix()` also appears transposed relative to ATNr's
  `create_Lmatrix()`: ATNr zeroes basal columns because basal species cannot be
  consumers, while the Python code zeroes basal rows.
- Basal loss to consumers should use predator biomass times the consumer's
  feeding rate on basal resource `i`. The current basal loop uses consumer
  biomass inside the resource-biomass power term.
- Consumer loss to predators should sum over all predators feeding on consumer
  `j`. The current consumer loop reuses the focal consumer's own functional
  response.
- The default Python configuration differs strongly from ATNr defaults. Most
  visibly, `CONFIG` has `r0 == X0 == 0.5`; for basal species this makes low-
  density growth and metabolic loss cancel before density dependence and
  herbivory are applied.

These points explain why the local model runs but does not reproduce the ATNr
Figure 3 and Figure 4 dynamics closely yet.
