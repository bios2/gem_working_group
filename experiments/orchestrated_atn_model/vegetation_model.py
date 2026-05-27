"""
Plant vegetation model for the spatial ATN.

This module owns the NPP-driven basal growth calculation and the plant
functional-type bookkeeping used by the ATN model.
"""
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
