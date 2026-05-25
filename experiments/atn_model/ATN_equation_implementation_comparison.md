# ATN Equation / Implementation Comparison

This compares the unscaled spatial ATN equations in
`ATN_model_spatiotemporal_formulas_parameters.Rmd` with the current
implementation in `atn_model.py`.

## Overall Conclusion

The current Python model follows the broad structure of the unscaled ATN:

- one local ODE system per grid cell;
- basal logistic growth;
- metabolic losses;
- consumer gains from feeding;
- consumer-resource functional response;
- optional temperature dependence.

However, the implementation does **not yet faithfully implement the original
resource-by-consumer equations**. The main problems are matrix orientation and
the way predation losses are calculated.

## Side-by-Side Checks

| Model component | Original model expectation | Current implementation | Assessment |
|---|---|---|---|
| Food-web orientation | Rows are resources, columns are consumers. `fw[i, j] = 1` means consumer `j` eats resource `i`. | The docstring states this convention, and `_functional_response()` reads consumer `j` from column `j`. | Mostly matches. |
| L matrix | `L[i, j]` is resource `i`, consumer `j`; basal species cannot be consumers, so basal **columns** should be zero. | `_L_matrix()` builds predator masses on rows and prey masses on columns, then zeroes basal **rows**. | Mismatch if `_L_matrix()` is used. |
| Pairwise allometric rates | `p[i, j] = p0 * M_resource[i]^b_prey * M_consumer[j]^b_pred * temp_factor`. | In `derivatives()`, `M_prey = M[np.newaxis, :]` and `M_pred = M[:, np.newaxis]`, so rows receive consumer/predator masses and columns receive resource/prey masses before multiplying by `adj_mat`. | Mismatch / transposed. |
| Temperature dependence | Rates multiply by `exp(-E * (T0 - T) / (k * T * T0))`. | `_allometric_rate()`, `_metabolic_rate()`, and `_basal_growth_rate()` use this same form. | Matches structurally. |
| Functional response | `F[i, j] = a[i, j] * B[i]^q / (1 + c[j]B[j] + sum_k h[k, j]a[k, j]B[k]^q)`. | `_functional_response(B, j)` uses column `j`, `a_ij * B^q`, and the correct denominator structure. | Matches structurally, assuming rate matrices have the correct orientation. |
| Basal growth | `B[i] * r[i] * (1 - sum(alpha[i,l]B[l]) / K[i])`. | Basal loop uses `r[i] * B[i] * (1 - B[i] / K_i)`. | Simplified match: correct for identity competition only. |
| Basal metabolic loss | `-X[i] * B[i]`. | Basal loop subtracts `X[i] * B[i]`. | Matches. |
| Basal loss to consumers | `-sum_j B[j] * F[i, j]`; consumer biomass multiplies feeding rate on resource `i`. | Basal loop multiplies by `B[consumer]` but also uses `B[consumer]^q` in the functional-response numerator instead of resource biomass `B[i]^q`. It also omits the interference term. | Mismatch. |
| Consumer gain | `B[i] * sum_j e[j] * F[j, i]`, where focal consumer `i` gains from its resources `j`. | Consumer loop computes `B[j] * sum_i e_i * F_i,j` for focal consumer `j`. | Matches structurally. |
| Consumer metabolic loss | `-X[i] * B[i]`. | Consumer loop includes `X[j] * B[j]`. | Matches. |
| Consumer loss to predators | `-sum_j B[j] * F[i, j]`, where predators `j` consume focal species `i`. | Consumer loop reuses the focal consumer's own `F_ij[j]` instead of computing every predator's feeding rate on focal species `j`. | Mismatch. |
| Extinction handling | Species below threshold should no longer contribute meaningfully to dynamics. | Species below threshold are given decay `-B / timescale`; functional response zeroes feeding on extinct resources. | Partial match. |

## Most Important Equation-Level Mismatches

### 1. Pairwise Rate Orientation

The original equation uses rows as resources and columns as consumers:

```text
p[i, j] = p0 * M_resource[i]^b_prey * M_consumer[j]^b_pred * temperature
```

The current code constructs:

```python
M_prey = M[np.newaxis, :]
M_pred = M[:, np.newaxis]
self.a_ij = self._allometric_rate('attack', M_prey, M_pred, T_K) * self.adj_mat
```

This makes row variation come from `M_pred` and column variation come from
`M_prey`, which is the opposite of the resource-by-consumer convention.

### 2. Basal Loss Uses Consumer Biomass Where Resource Biomass Is Needed

The original basal loss term is:

```text
sum_j B[j] * F[i, j]
```

and `F[i, j]` contains `B[i]^q` in its numerator.

The current basal loop uses:

```python
B[self.consumer_idx] * self.a_ij[i, self.consumer_idx] *
np.power(B[self.consumer_idx], self.q_hill)
```

The multiplier `B[self.consumer_idx]` is correct because consumers remove
resource biomass, but the `B^q` term should be resource biomass `B[i]^q`, not
consumer biomass raised to `q`.

### 3. Consumer Predation Loss Reuses The Wrong Functional Response

The original consumer loss term for focal species `i` is:

```text
sum_j B[j] * F[i, j]
```

That requires every predator `j`'s feeding rate on focal species `i`.

The current consumer loop computes `F_ij` only for the focal consumer `j`, then
uses `F_ij[j]` in its own loss term. That is the focal consumer feeding on
itself, not predators feeding on the focal consumer.

## What This Means

The current code is a useful scaffold, and some pieces are aligned, but it is
not yet a faithful implementation of the original unscaled ATN equations. The
first fixes should be:

1. Correct pairwise rate matrix orientation.
2. Make all loss-to-consumers terms use a full feeding matrix `F[resource, consumer]`.
3. Recompute basal and consumer losses from that feeding matrix.
4. Add toy tests for one plant, plant-herbivore, and three-species chain cases.
