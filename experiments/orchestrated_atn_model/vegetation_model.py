"""
Plant vegetation model for the spatial ATN.

This module owns the NPP-driven basal growth calculation and the plant
functional-type bookkeeping used by the ATN model.
"""
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd


class PlantVegetationModel:
    """
    NPP-driven plant biomass growth model.

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

        Columns: pixel_id, x, y, time, species, delta_biomass
          delta_biomass = G_i (g/m²/day) — the NPP-driven growth contribution
        """
        n_tp = len(t_eval)
        n_cells = len(self.env)
        n_basal = len(self.basal_idx)

        G_arr = np.zeros((n_tp, n_cells, n_basal))
        for t_idx in range(n_tp):
            for cell_idx in range(n_cells):
                G = self.growth(B_traj[t_idx, cell_idx, :], cell_idx)
                G_arr[t_idx, cell_idx, :] = G[self.basal_idx]

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
