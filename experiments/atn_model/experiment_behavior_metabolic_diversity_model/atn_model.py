"""
Spatially explicit unscaled ATN model with metabolic diversity.

This version preserves the behavior-model trophic equations and adds an
optional thermal-group metabolism model for consumers:
  - ectotherms use environmental temperature;
  - endotherms use fixed bird or mammal body temperature.
"""

from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd
from scipy.integrate import odeint

from metabolism import calculate_metabolic_loss_per_day


class ATNModel:
    """Unscaled ATN model with optional endotherm/ectotherm metabolism."""

    def __init__(self, adj_mat: np.ndarray, traits_df: pd.DataFrame, env_df: pd.DataFrame, config: Dict):
        self.adj_mat = adj_mat
        self.traits = traits_df.copy()
        self.env = env_df.copy()
        self.cfg = config

        self.n_species = len(traits_df)
        self.n_cells = len(env_df)
        self.basal_idx = np.where(self.traits["is_basal"].values == 1)[0]
        self.consumer_idx = np.where(self.traits["is_basal"].values == 0)[0]

        self.r0 = config["r0"]
        self.br = config["b_r"]
        self.X0 = config["X0"]
        self.bX = config["b_X"]
        self.a0 = config["a0"]
        self.ba_prey = config["b_a_prey"]
        self.ba_pred = config["b_a_pred"]
        self.h0 = config["h0"]
        self.bh_prey = config["b_h_prey"]
        self.bh_pred = config["b_h_pred"]
        self.q_hill = config["q_hill"]
        self.Ropt = config["R_opt"]
        self.gamma = config["gamma"]
        self.link_threshold = config["link_threshold"]
        self.T0_K = config["T0_K"]
        self.k_B = config["k_B"]
        self.E_a = config["E_a"]
        self.ext_threshold = config["ext_threshold"]

    def _L_matrix(self, cell_idx: int) -> np.ndarray:
        M = self.traits["body_mass_g"].values
        z = M[np.newaxis, :] / (M[:, np.newaxis] * self.Ropt)
        L = np.power(z * np.exp(1.0 - z), self.gamma)
        L[L < self.link_threshold] = 0
        L[:, self.basal_idx] = 0
        return L * self.adj_mat

    def _temperature_factor(self, T_K: float) -> float:
        return np.exp(-self.E_a * (self.T0_K - T_K) / (self.k_B * T_K * self.T0_K))

    def _allometric_rate(self, rate_type: str, M_prey: np.ndarray, M_pred: np.ndarray, T_K: float) -> np.ndarray:
        if rate_type == "attack":
            p0, b_prey, b_pred = self.a0, self.ba_prey, self.ba_pred
        elif rate_type == "handling":
            p0, b_prey, b_pred = self.h0, self.bh_prey, self.bh_pred
        else:
            raise ValueError(f"Unknown rate type: {rate_type}")

        p_allom = p0 * np.power(M_prey, b_prey) * np.power(M_pred, b_pred)
        if self.cfg["use_temperature"]:
            p_allom *= self._temperature_factor(T_K)
        return p_allom

    def _legacy_metabolic_rate(self, M: np.ndarray, T_K: float) -> np.ndarray:
        X = self.X0 * np.power(M, self.bX)
        if self.cfg["use_temperature"]:
            X *= self._temperature_factor(T_K)
        return X

    def _trait_value(self, species_idx: int, column: str, default):
        if column not in self.traits.columns:
            return default
        value = self.traits.iloc[species_idx][column]
        if pd.isna(value):
            return default
        return value

    def _metabolic_rate(self, M: np.ndarray, T_K: float) -> np.ndarray:
        """
        Return X_i in day^-1 for every species.

        In thermal_group mode, basal species keep the legacy ATN metabolism
        unless apply_thermal_metabolism_to_basal is True. Consumers use the
        Blyth/Gillooly metabolism implemented in metabolism.py.
        """
        if self.cfg.get("metabolism_model", "atn") != "thermal_group":
            return self._legacy_metabolic_rate(M, T_K)

        X = self._legacy_metabolic_rate(M, T_K)
        target_idx = np.arange(self.n_species)
        if not self.cfg.get("apply_thermal_metabolism_to_basal", False):
            target_idx = self.consumer_idx

        for idx in target_idx:
            default_group = (
                "ectotherm"
                if idx in self.basal_idx
                else self.cfg.get("default_consumer_thermal_group", "ectotherm")
            )
            thermal_group = str(self._trait_value(idx, "thermal_group", default_group)).lower()
            if thermal_group in {"basal", "plant", "producer", "legacy", "atn"}:
                continue
            endotherm_group = self._trait_value(
                idx,
                "endotherm_group",
                self.cfg.get("default_endotherm_group", "mammal"),
            )
            if thermal_group == "ectotherm":
                endotherm_group = None
            X[idx] = calculate_metabolic_loss_per_day(
                M[idx],
                T_K,
                thermal_group,
                endotherm_group=endotherm_group,
                field_metabolic_multiplier=self.cfg.get("field_metabolic_multiplier", 3.0),
                seconds_per_day=self.cfg.get("seconds_per_day", 86400.0),
                joules_per_g_wet_biomass=self.cfg.get("joules_per_g_wet_biomass", 7000.0),
            )

        return X * self.cfg.get("metabolic_rate_multiplier", 1.0)

    def _basal_growth_rate(self, M: np.ndarray, T_K: float) -> np.ndarray:
        r = self.r0 * np.power(M, self.br)
        if self.cfg["use_temperature"]:
            r *= self._temperature_factor(T_K)
        return r

    def _functional_response(self, B: np.ndarray, j_idx: int) -> np.ndarray:
        B_j = B[j_idx]
        a_ij = self.a_ij[:, j_idx]
        h_kj = self.h_ij[:, j_idx]
        B_q = np.power(B, self.q_hill)
        numerator = a_ij * B_q
        denominator = 1.0 + self.cfg["interference"] * B_j + np.sum(h_kj * a_ij * B_q)
        F_ij = numerator / denominator
        F_ij[B < self.ext_threshold] = 0
        return F_ij

    def derivatives(self, y: np.ndarray, t: float, cell_idx: int) -> np.ndarray:
        B = np.maximum(y, 0)
        dydt = np.zeros_like(B)
        env_row = self.env.iloc[cell_idx]
        T_K = env_row["temperature_K"]
        M = self.traits["body_mass_g"].values

        M_prey = M[:, np.newaxis]
        M_pred = M[np.newaxis, :]
        self.a_ij = self._allometric_rate("attack", M_prey, M_pred, T_K) * self.adj_mat
        self.h_ij = self._allometric_rate("handling", M_prey, M_pred, T_K) * self.adj_mat
        X = self._metabolic_rate(M, T_K)
        r = self._basal_growth_rate(M, T_K)

        for i in self.basal_idx:
            k_col = f"K_plant_{i}"
            K_i = env_row[k_col] if k_col in self.env.columns else self.cfg["K_default"]
            growth = r[i] * B[i] * (1.0 - B[i] / K_i) if K_i > 0 else 0
            loss_to_consumers = sum(B[j] * self._functional_response(B, j)[i] for j in self.consumer_idx)
            dydt[i] = growth - X[i] * B[i] - loss_to_consumers

        e_i = np.array(
            [self.cfg["e_plant"] if self.traits.iloc[i]["is_basal"] else self.cfg["e_animal"] for i in range(self.n_species)]
        )
        for j in self.consumer_idx:
            F_ij = self._functional_response(B, j)
            gain = B[j] * np.sum(e_i * F_ij)
            loss_to_predators = sum(B[jp] * self._functional_response(B, jp)[j] for jp in self.consumer_idx)
            dydt[j] = gain - X[j] * B[j] - loss_to_predators

        dydt[B < self.ext_threshold] = -B[B < self.ext_threshold] / self.cfg["extinction_timescale"]
        return dydt

    def run_cell(self, B_initial: np.ndarray, cell_idx: int, t_eval: np.ndarray) -> np.ndarray:
        def deriv(y, t):
            return self.derivatives(y, t, cell_idx)

        return odeint(deriv, B_initial, t_eval, full_output=False, rtol=1e-6, atol=1e-8)

    def run_all_cells(self, B_initial: np.ndarray, t_eval: np.ndarray) -> np.ndarray:
        B_traj = np.zeros((len(t_eval), self.n_cells, self.n_species))
        for cell_idx in range(self.n_cells):
            B_traj[:, cell_idx, :] = self.run_cell(B_initial[cell_idx, :], cell_idx, t_eval)
        return B_traj
