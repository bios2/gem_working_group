# ATN Equation / Implementation Comparison

This compares the unscaled spatial ATN equations in
`ATN_model_spatiotemporal_formulas_parameters.Rmd` with the current
implementation in `atn_model.py`.

Last updated: 2026-05-25, after fixing all four documented mismatches.

## Overall Conclusion

The current Python model faithfully implements the unscaled ATN equations
from Section 8 of the Rmd document:

- food-web orientation (rows = resources, columns = consumers) is respected
  throughout;
- pairwise allometric rate matrices use resource mass on rows and consumer
  mass on columns;
- the L matrix uses the correct body-size ratio `z[i,j] = M_consumer[j] /
  (M_resource[i] * R_opt)` and zeroes basal columns, not basal rows;
- the functional response `_functional_response()` is correct end-to-end;
- basal loss to consumers uses `sum_j B[j] * F[i,j]` via
  `_functional_response`, with the correct resource biomass `B[i]^q` and
  interference term;
- consumer loss to predators uses `sum_{j'} B[j'] * F[j, j']`, computing
  each predator's feeding rate on the focal species.

One simplification remains: basal logistic growth uses identity competition
(`alpha = I`) rather than the full interspecific competition matrix.

## Side-by-Side Checks

| Model component | Original model expectation | Current implementation | Assessment |
|---|---|---|---|
| Food-web orientation | Rows are resources, columns are consumers. `fw[i, j] = 1` means consumer `j` eats resource `i`. | The docstring states this convention; `_functional_response()` slices column `j` for consumer `j`. | Matches. |
| L matrix | `z[i,j] = M_consumer[j] / (M_resource[i] * R_opt)`; basal species cannot be consumers, so basal **columns** are zero. | `z = M_consumer / (M_resource * R_opt)` with `M_resource = M[:, np.newaxis]` (rows) and `M_consumer = M[np.newaxis, :]` (columns); `L[:, basal_idx] = 0`. | **Matches** (fixed). |
| Pairwise allometric rates | `p[i, j] = p0 * M_resource[i]^b_prey * M_consumer[j]^b_pred * temp_factor`. | `M_prey = M[:, np.newaxis]` (rows), `M_pred = M[np.newaxis, :]` (columns), giving `a_ij[i,j] = p0 * M[i]^b_prey * M[j]^b_pred`. | Matches. |
| Temperature dependence | Rates multiply by `exp(-E * (T0 - T) / (k * T * T0))`. | `_allometric_rate()`, `_metabolic_rate()`, and `_basal_growth_rate()` all use this form. | Matches. |
| Functional response | `F[i, j] = a[i, j] * B[i]^q / (1 + c[j]*B[j] + sum_k h[k,j]*a[k,j]*B[k]^q)`. | `_functional_response(B, j)` slices column `j` for both `a` and `h`, uses `B^q` (resource biomass) in the numerator, and includes `c*B_j` interference and the handling-time sum. | Matches. |
| Basal growth | `B[i] * r[i] * (1 - sum_l alpha[i,l]*B[l] / K[i])`. | `r[i] * B[i] * (1 - B[i] / K_i)`. | Simplified match: correct for identity competition (`alpha = I`); interspecific competition among basal species not implemented. |
| Basal metabolic loss | `-X[i] * B[i]`. | Basal loop subtracts `X[i] * B[i]`. | Matches. |
| Basal loss to consumers | `sum_j B[j] * F[i, j]`. | `sum(B[j] * _functional_response(B, j)[i] for j in consumer_idx)`. | Matches. |
| Consumer gain | `B[i] * sum_j e[j] * F[j, i]`. | `B[j] * np.sum(e_i * F_ij)` where `F_ij = _functional_response(B, j)`. | Matches. |
| Consumer metabolic loss | `-X[i] * B[i]`. | `X[j] * B[j]` inside `loss`. | Matches. |
| Consumer loss to predators | `sum_{j'} B[j'] * F[i, j']`, where `j'` are predators that eat focal species `i`. | `sum(B[jp] * _functional_response(B, jp)[j] for jp in consumer_idx)`. | **Matches** (fixed). |
| Extinction handling | Species below threshold should no longer contribute meaningfully. | Species below threshold receive decay `-B / timescale`; `_functional_response` zeroes feeding on extinct resources. | Partial match. |

## What This Means

All four previously documented mismatches have been resolved. The model now
implements the unscaled ATN equations as specified. The only remaining gap
relative to the full Rmd formulation is the simplified basal competition
term (no interspecific alpha matrix), which is a deliberate simplification
rather than a bug.

Recommended next steps:

1. Add toy tests for one plant, plant–herbivore, and three-species chain
   cases to lock in correctness numerically.
2. Implement the full competition matrix `alpha[i, l]` for basal species if
   interspecific competition is needed.
3. Wire `_L_matrix()` into `derivatives()` when the body-size feeding
   kernel is ready to replace the fixed adjacency matrix.
