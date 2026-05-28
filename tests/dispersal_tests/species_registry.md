# SpeciesRegistry — attribute synthesis

Source: [src/gem/engine/species_registry.py](../../src/gem/engine/species_registry.py)

`SpeciesRegistry` is the "species book" of the engine: it holds the species list, what groups each species belongs to, their per-species parameters, and the food-web adjacency matrix. It is built once at the start of a simulation and then passed (by reference) into `EcosystemGridState`.

Species are referenced by **integer index** (0 to `num_species - 1`), never by name. Groups and parameters are layered on top of that integer indexing.

## Attributes

All four attributes are created by `__init__(num_species)`.

| Attribute | Type | Initial value (for `num_species = N`) | Purpose |
|---|---|---|---|
| `num_species` | `int` | `N` | Total species in the model. Used everywhere downstream to size `(S,)` parameter arrays and the species axis of grid layers. |
| `groups` | `dict[str, list[int]]` | `{"all": [0, 1, …, N-1]}` | Named groups → list of species indices. The default group `"all"` already contains every species. New groups are added via `add_species_to_group`. |
| `params` | `dict[str, np.ndarray]` | `{}` | Named per-species parameters. Each value is a flat `(N,)` numpy array; entry `i` is the parameter value for species `i`. Empty until populated. |
| `adj_mat` | `np.ndarray`, shape `(N, N)`, dtype `int8` | all zeros | Food-web adjacency. `adj_mat[i, j] == 1` means **species i is eaten by species j** (i.e. i is a resource for consumer j). Dispersal does not use this; ATN will. |
