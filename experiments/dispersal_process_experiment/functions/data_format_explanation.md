# Numpy array shapes

## 2D array — one species on a 3×3 grid

Shape: `(3, 3)`. Each element is the biomass at a given row and column.

```
biomass[:, :]

  col 0  col 1  col 2
  -----  -----  -----
   1.2    0.0    0.5   ← row 0
   0.8    2.1    0.3   ← row 1
   0.0    0.4    1.7   ← row 2
```

To read a single cell: `biomass[row, col]`. For example `biomass[1, 1]` gives `2.1`.

## 3D array — three species on a 3×3 grid

Shape: `(3, 3, 3)` — rows × columns × species. Think of it as a stack of 2D grids,
one per species, placed behind each other.

```
biomass[:, :, 0]   ← species 0 (e.g. grass)
  1.2    0.0    0.5
  0.8    2.1    0.3
  0.0    0.4    1.7

biomass[:, :, 1]   ← species 1 (e.g. rabbit)
  0.0    0.3    0.0
  0.1    0.0    0.2
  0.0    0.0    0.1

biomass[:, :, 2]   ← species 2 (e.g. fox)
  0.0    0.0    0.1
  0.0    0.0    0.0
  0.0    0.1    0.0
```

To read a single cell: `biomass[row, col, species]`. For example `biomass[1, 1, 0]`
gives `2.1` — the grass biomass at row 1, column 1.

## Why species goes on the trailing (last) axis

With species on the trailing axis, a per-species parameter like growth rate `r` of
shape `(3,)` broadcasts cleanly against `biomass` of shape `(3, 3, 3)`. Numpy lines
up axes from the right, so `r[s]` automatically applies to `biomass[:, :, s]` for
every cell:

```python
r = np.array([0.5, 0.2, 0.1])   # shape (3,)  — one rate per species
biomass * r                       # shape (3, 3, 3) — r broadcasts over all cells
```

If species were on the first axis instead, you would need to reshape `r` manually
before every operation. Putting species last avoids that.
