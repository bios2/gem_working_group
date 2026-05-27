"""
Plant vegetation model for the spatial ATN.

This module owns the NPP-driven basal growth calculation and the plant
functional-type bookkeeping used by the ATN model.
"""
from pathlib import Path
from typing import Dict

import numpy as np
<<<<<<< HEAD
from numpy.typing import NDArray
=======
>>>>>>> refs/remotes/origin/orchestrated_vegatn
import pandas as pd


class PlantVegetationModel:
    """
    NPP-driven plant biomass growth model.

<<<<<<< HEAD
    growth()          accepts B: (S,)         — used by the single-cell interface
    growth_all_cells() accepts B: (n_cells, S) — used by derivatives() and save_output()
    """

    def __init__(self, traits_df: pd.DataFrame, env_df: pd.DataFrame, config: Dict):
=======
    The ATN passes the current biomass vector and cell index to this class.
    The vegetation model returns a growth-rate vector G with one entry per
    species; non-basal species receive zero growth.
    """

    def __init__(self, traits_df: pd.DataFrame, env_df: pd.DataFrame, config: Dict):
        """
        Initialize plant vegetation state and parameters.

        Parameters:
            traits_df: DataFrame with species traits, including is_basal and vegetation_type
            env_df: DataFrame with per-cell environmental inputs, including NPP
            config: Dict with vegetation parameters
        """
>>>>>>> refs/remotes/origin/orchestrated_vegatn
        self.traits = traits_df
        self.env = env_df
        self.cfg = config
        self.n_species = len(traits_df)

        self.basal_idx = np.where(traits_df['is_basal'] == 1)[0]

        herb_mask = (traits_df['is_basal'] == 1) & (traits_df['vegetation_type'] == 'herb')
        tree_mask = (traits_df['is_basal'] == 1) & (traits_df['vegetation_type'] == 'tree')
        self.herb_idx = np.where(herb_mask)[0]
        self.tree_idx = np.where(tree_mask)[0]

        f_default = config.get('f_struct_default', 0.3)
        if 'f_struct' in traits_df.columns:
            self.f_struct = traits_df['f_struct'].fillna(f_default).values.astype(float)
        else:
            self.f_struct = np.full(self.n_species, f_default)

        self.alpha_herbs = config['alpha_herbs_default']
        self.psi = config['psi']

<<<<<<< HEAD
        # Pre-extract NPP for vectorised multi-cell calls — avoids repeated DataFrame lookup
        self.npp = env_df['NPP'].values.astype(float)  # (n_cells,)

    def growth(self, B: NDArray[np.float64], cell_idx: int,
               t: float = None) -> NDArray[np.float64]:
        """
        NPP-driven leaf biomass growth rate G_i for each species in one cell.

        B: (S,) — current biomass for this cell
        Returns G: (S,) — growth rate per species (zero for consumers)
        """
        del t  # reserved for future time-varying vegetation drivers

        G = np.zeros(self.n_species)
        NPP = self.npp[cell_idx]
        B_trees = B[self.tree_idx].sum() if len(self.tree_idx) > 0 else 0.0

        # All herbs share the same C_i (Michaelis-Menten on total tree biomass)
        G[self.herb_idx] = (NPP * self.psi
                            * (1.0 - self.f_struct[self.herb_idx])
                            * self.alpha_herbs / (self.alpha_herbs + B_trees))

        if len(self.tree_idx) > 0:
            C_tree = B[self.tree_idx] / (self.alpha_herbs + B_trees)
            G[self.tree_idx] = (NPP * self.psi
                                * (1.0 - self.f_struct[self.tree_idx])
                                * C_tree)
        return G

    def growth_all_cells(self, B: NDArray[np.float64]) -> NDArray[np.float64]:
        """
        NPP-driven leaf biomass growth rate G_i for all cells simultaneously.

        B: (n_cells, S) — must have the same shape
        Returns G: (n_cells, S) — growth rate per cell per species
        """
        assert B.ndim == 2 and B.shape[1] == self.n_species

        n_cells = B.shape[0]
        G = np.zeros((n_cells, self.n_species))

        B_trees = (B[:, self.tree_idx].sum(axis=1)
                   if len(self.tree_idx) > 0
                   else np.zeros(n_cells))          # (n_cells,)

        # Herbs: C_i is the same for every herb in a given cell
        # self.npp: (n_cells,); f_struct[herb_idx]: (n_herbs,)
        C_herb = self.alpha_herbs / (self.alpha_herbs + B_trees)   # (n_cells,)
        G[:, self.herb_idx] = (
            (self.npp * self.psi)[:, np.newaxis]
            * (1.0 - self.f_struct[self.herb_idx])[np.newaxis, :]
            * C_herb[:, np.newaxis]
        )

        # Trees: C_i = B[i] / (alpha + B_trees) differs per tree species
        if len(self.tree_idx) > 0:
            C_tree = B[:, self.tree_idx] / (self.alpha_herbs + B_trees[:, np.newaxis])  # (n_cells, n_trees)
            G[:, self.tree_idx] = (
                (self.npp * self.psi)[:, np.newaxis]
                * (1.0 - self.f_struct[self.tree_idx])[np.newaxis, :]
                * C_tree
            )

        return G

    def save_output(self, B_traj: NDArray[np.float64],
                    t_eval: NDArray[np.float64], output_dir) -> None:
        """
        Save instantaneous vegetation growth rates for all basal species to vegetation.txt.

        Evaluates growth_all_cells() at every recorded time point — no inner cell loop.
=======
    def growth(self, B: np.ndarray, cell_idx: int, t: float = None) -> np.ndarray:
        """
        Compute NPP-driven leaf biomass growth rate G_i for each basal species.

        Implements the vegetation.md equation:
            G_i = NPP * psi * (1 - f_struct_i) * C_i

        Competitive partition C_i:
            Herb i:  C_i = alpha / (alpha + B_trees)
            Tree i:  C_i = B[i]  / (alpha + B_trees)

        The optional time argument is accepted so this interface can later
        support seasonal or time-varying vegetation drivers.
        """
        del t  # Reserved for future time-dependent vegetation dynamics.

        G = np.zeros(self.n_species)
        NPP = self.env.loc[cell_idx, 'NPP']
        B_trees = float(np.sum(B[self.tree_idx])) if len(self.tree_idx) > 0 else 0.0

        for i in self.herb_idx:
            C_i = self.alpha_herbs / (self.alpha_herbs + B_trees)
            G[i] = NPP * self.psi * (1.0 - self.f_struct[i]) * C_i

        for i in self.tree_idx:
            C_i = B[i] / (self.alpha_herbs + B_trees)
            G[i] = NPP * self.psi * (1.0 - self.f_struct[i]) * C_i

        return G

    def save_output(self, B_traj: np.ndarray, t_eval: np.ndarray, output_dir) -> None:
        """
        Save instantaneous vegetation growth rates for all basal species to vegetation.txt.

        Evaluates growth(B, cell) at every recorded time point and cell using the
        saved biomass trajectory, then writes a long-format table.
>>>>>>> refs/remotes/origin/orchestrated_vegatn

        Columns: pixel_id, x, y, time, species, delta_biomass
          delta_biomass = G_i (g/m²/day) — the NPP-driven growth contribution
        """
        n_tp = len(t_eval)
        n_cells = len(self.env)
        n_basal = len(self.basal_idx)

        G_arr = np.zeros((n_tp, n_cells, n_basal))
        for t_idx in range(n_tp):
<<<<<<< HEAD
            G = self.growth_all_cells(B_traj[t_idx])      # (n_cells, S) — no inner cell loop
            G_arr[t_idx] = G[:, self.basal_idx]
=======
            for cell_idx in range(n_cells):
                G = self.growth(B_traj[t_idx, cell_idx, :], cell_idx)
                G_arr[t_idx, cell_idx, :] = G[self.basal_idx]
>>>>>>> refs/remotes/origin/orchestrated_vegatn

        cell_x = self.env['x'].values.astype(int)
        cell_y = self.env['y'].values.astype(int)

        t_rep    = np.repeat(t_eval, n_cells * n_basal)
        cell_rep = np.tile(np.repeat(np.arange(n_cells), n_basal), n_tp)
        x_rep    = np.tile(np.repeat(cell_x, n_basal), n_tp)
        y_rep    = np.tile(np.repeat(cell_y, n_basal), n_tp)
        sp_rep   = np.tile(self.basal_idx, n_tp * n_cells)
        g_rep    = G_arr.ravel()

        table = np.column_stack([cell_rep, x_rep, y_rep, t_rep, sp_rep, g_rep])
        np.savetxt(
            Path(output_dir) / 'vegetation.txt', table,
            fmt=['%d', '%d', '%d', '%.4f', '%d', '%.6e'],
            header='pixel_id x y time species delta_biomass',
            comments=''
        )
