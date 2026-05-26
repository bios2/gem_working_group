"""
Reproduce the ATNr paper Figure 3 and Figure 4 scenarios with this Python ATN.

This script intentionally does not modify the model equations in atn_model.py.
It uses the installed R package ATNr only to generate the paper's input
scenarios and reference outputs, then runs the same species/body-mass/food-web
inputs through the local Python ATNModel where the current model has an
equivalent concept.

Important model-scope notes:
  * Figure 3 in Gauzens et al. uses the unscaled ATN with explicit nutrients.
    The local Python model is unscaled without nutrients, so the Python run is
    a closest-available analogue rather than an exact reproduction.
  * Figure 4 is described in the paper text as using the scaled model. The
    local Python model is unscaled, so K=1 and K=10 are tested with the local
    logistic basal-growth implementation.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/matplotlib-cache")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from atn_model import ATNModel
from config import CONFIG


R_SCENARIO_SCRIPT = r"""
suppressPackageStartupMessages(library(ATNr))
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[[1]]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

write_csv <- function(x, name, row.names = FALSE) {
  utils::write.csv(x, file.path(out_dir, name), row.names = row.names)
}

# ---------------- Figure 3: temperature versus extinctions ----------------
set.seed(12)
n_species <- 50
n_basal <- 20
n_nut <- 2
masses <- 10 ^ c(sort(runif(n_basal, 1, 3)),
                 sort(runif(n_species - n_basal, 2, 9)))
L <- create_Lmatrix(masses, n_basal, Ropt = 50, gamma = 2, th = 0.01)
fw <- L
fw[fw > 0] <- 1
model <- create_model_Unscaled_nuts(n_species, n_basal, n_nut, masses, fw)
temperatures <- seq(4, 22, by = 2)
biomasses <- runif(n_species + n_nut, 2, 3)
times <- seq(0, 100000, 100)
extinctions <- rep(NA, length(temperatures))

for (i in seq_along(temperatures)) {
  temp <- temperatures[[i]]
  model <- initialise_default_Unscaled_nuts(model, L, temperature = temp)
  model$q <- rep(1.4, n_species - n_basal)
  model$S <- rep(10, n_nut)
  sol <- lsoda_wrapper(times, biomasses, model, verbose = FALSE)
  extinctions[[i]] <- sum(sol[nrow(sol), 4:ncol(sol)] < 1e-6)
}

write_csv(data.frame(species = seq_len(n_species) - 1,
                     body_mass_g = masses,
                     is_basal = as.integer(seq_len(n_species) <= n_basal)),
          "figure3_species.csv")
write_csv(fw, "figure3_food_web.csv", row.names = FALSE)
write_csv(data.frame(state = c(paste0("nutrient_", seq_len(n_nut) - 1),
                               paste0("species_", seq_len(n_species) - 1)),
                     initial = biomasses),
          "figure3_initial_state.csv")
write_csv(data.frame(temperature_C = temperatures,
                     extinctions_atnr = extinctions),
          "figure3_atnr_extinctions.csv")

# ---------------- Figure 4: K=1 versus K=10 enrichment scenario ------------
set.seed(1234)
S <- 10
fw4 <- create_niche_model(S, C = .15)
TL <- TroLev(fw4)
masses4 <- as.numeric(0.01 * 100 ^ (TL - 1))
n_basal4 <- sum(colSums(fw4) == 0)
mod <- create_model_Scaled(nb_s = S, nb_b = n_basal4, BM = masses4, fw = fw4)
mod <- initialise_default_Scaled(mod)
times4 <- seq(0, 300, by = 2)
biomasses4 <- runif(S, 2, 3)

mod$K <- 1
sol1 <- lsoda_wrapper(times4, biomasses4, mod, verbose = FALSE)
mod$K <- 10
sol10 <- lsoda_wrapper(times4, biomasses4, mod, verbose = FALSE)

write_csv(data.frame(species = seq_len(S) - 1,
                     body_mass_g = masses4,
                     trophic_level = as.numeric(TL),
                     is_basal = as.integer(seq_len(S) <= n_basal4)),
          "figure4_species.csv")
write_csv(fw4, "figure4_food_web.csv", row.names = FALSE)
write_csv(data.frame(species = seq_len(S) - 1,
                     initial_biomass = biomasses4),
          "figure4_initial_biomass.csv")
write_csv(as.data.frame(sol1), "figure4_atnr_K1_timeseries.csv")
write_csv(as.data.frame(sol10), "figure4_atnr_K10_timeseries.csv")
"""


def run_r_reference(output_dir: Path) -> None:
    """Generate ATNr reference inputs and outputs from the paper code."""
    r_file = output_dir / "_generate_atnr_reference.R"
    r_file.write_text(R_SCENARIO_SCRIPT)
    subprocess.run(
        ["Rscript", str(r_file), str(output_dir)],
        check=True,
    )


def read_matrix(path: Path) -> np.ndarray:
    return pd.read_csv(path).to_numpy(dtype=float)


def make_traits(species_df: pd.DataFrame, initial_biomass: np.ndarray) -> pd.DataFrame:
    traits = pd.DataFrame(
        {
            "species_id": species_df["species"].astype(int),
            "body_mass_g": species_df["body_mass_g"].astype(float),
            "is_basal": species_df["is_basal"].astype(int),
            "initial_biomass_g_per_m2": initial_biomass.astype(float),
        }
    )
    return traits.set_index("species_id")


def make_env(temperature_c: float, basal_indices: np.ndarray, k_value: float) -> pd.DataFrame:
    row = {"cell_id": 0, "temperature_K": temperature_c + 273.15}
    for i in basal_indices:
        row[f"K_plant_{int(i)}"] = k_value
    return pd.DataFrame([row]).set_index("cell_id")


def run_python_single_cell(
    adj_mat: np.ndarray,
    traits: pd.DataFrame,
    temperature_c: float,
    t_eval: np.ndarray,
    k_value: float,
    config_updates: dict | None = None,
    integrator: str = "euler",
    euler_dt: float = 1.0,
) -> np.ndarray:
    cfg = CONFIG.copy()
    if config_updates:
        cfg.update(config_updates)
    basal_indices = np.where(traits["is_basal"].to_numpy() == 1)[0]
    env = make_env(temperature_c, basal_indices, k_value)
    model = ATNModel(adj_mat, traits, env, cfg)
    y0 = traits["initial_biomass_g_per_m2"].to_numpy(dtype=float)
    if integrator == "odeint":
        traj = model.run_cell(y0, 0, t_eval)
    elif integrator == "euler":
        traj = run_cell_euler(model, y0, 0, t_eval, euler_dt)
    else:
        raise ValueError(f"Unknown integrator: {integrator}")
    traj = np.maximum(traj, 0)
    return traj


def run_cell_euler(
    model: ATNModel,
    y0: np.ndarray,
    cell_idx: int,
    t_eval: np.ndarray,
    max_dt: float,
) -> np.ndarray:
    """Small explicit integrator for bounded reproduction sweeps."""
    if max_dt <= 0:
        raise ValueError("Euler dt must be positive")
    y = np.maximum(y0.astype(float), 0)
    traj = np.zeros((len(t_eval), len(y)))
    traj[0] = y
    t = float(t_eval[0])

    for out_idx in range(1, len(t_eval)):
        target = float(t_eval[out_idx])
        while t < target:
            dt = min(max_dt, target - t)
            dydt = model.derivatives(y, t, cell_idx)
            y = y + dt * dydt
            y = np.nan_to_num(y, nan=0.0, posinf=1e12, neginf=0.0)
            y = np.clip(y, 0.0, 1e12)
            t += dt
        traj[out_idx] = y
    return traj


def reproduce_figure3(
    output_dir: Path,
    python_t_max: float,
    python_dt: float,
    integrator: str,
    euler_dt: float,
) -> pd.DataFrame:
    species = pd.read_csv(output_dir / "figure3_species.csv")
    fw = read_matrix(output_dir / "figure3_food_web.csv")
    init_state = pd.read_csv(output_dir / "figure3_initial_state.csv")
    atnr = pd.read_csv(output_dir / "figure3_atnr_extinctions.csv")

    n_species = len(species)
    species_initial = init_state["initial"].to_numpy(dtype=float)[-n_species:]
    traits = make_traits(species, species_initial)
    t_eval = np.arange(0, python_t_max + python_dt, python_dt, dtype=float)

    python_ext = []
    for temp in atnr["temperature_C"].to_numpy(dtype=float):
        print(f"Running Python Figure 3 analogue at {temp:g} deg C "
              f"to t={python_t_max:g} ...", flush=True)
        traj = run_python_single_cell(
            fw,
            traits,
            temperature_c=temp,
            t_eval=t_eval,
            k_value=CONFIG["K_default"],
            config_updates={"q_hill": 1.4},
            integrator=integrator,
            euler_dt=euler_dt,
        )
        python_ext.append(int(np.sum(traj[-1] < CONFIG["ext_threshold"])))

    out = atnr.copy()
    out["python_t_max"] = python_t_max
    out["extinctions_python_unscaled_no_nutrients"] = python_ext
    out.to_csv(output_dir / "figure3_python_vs_atnr_extinctions.csv", index=False)

    fig, ax = plt.subplots(figsize=(5.5, 3.6), dpi=160)
    ax.plot(out["temperature_C"], out["extinctions_atnr"], "o-", color="#1f77b4",
            label="ATNr unscaled + nutrients")
    ax.plot(out["temperature_C"], out["extinctions_python_unscaled_no_nutrients"],
            "s--", color="#d62728", label="Python unscaled, no nutrients")
    ax.set_xlabel("Temperature (deg C)")
    ax.set_ylabel("Number of extinctions")
    ax.set_ylim(0, 50)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "figure3_temperature_extinctions_python_vs_atnr.png")
    plt.close(fig)
    return out


def reproduce_figure4(output_dir: Path, integrator: str, euler_dt: float) -> tuple[np.ndarray, np.ndarray]:
    species = pd.read_csv(output_dir / "figure4_species.csv")
    fw = read_matrix(output_dir / "figure4_food_web.csv")
    initial = pd.read_csv(output_dir / "figure4_initial_biomass.csv")[
        "initial_biomass"
    ].to_numpy(dtype=float)
    traits = make_traits(species, initial)
    t_eval = np.arange(0, 300 + 2, 2, dtype=float)

    traj_k1 = run_python_single_cell(
        fw, traits, 20.0, t_eval, k_value=1.0,
        integrator=integrator, euler_dt=euler_dt
    )
    traj_k10 = run_python_single_cell(
        fw, traits, 20.0, t_eval, k_value=10.0,
        integrator=integrator, euler_dt=euler_dt
    )

    pd.DataFrame(np.column_stack([t_eval, traj_k1])).to_csv(
        output_dir / "figure4_python_K1_timeseries.csv", index=False
    )
    pd.DataFrame(np.column_stack([t_eval, traj_k10])).to_csv(
        output_dir / "figure4_python_K10_timeseries.csv", index=False
    )

    colors = plt.cm.coolwarm(np.linspace(0.05, 0.95, traj_k1.shape[1]))
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.0), dpi=160, sharex="col")
    atnr_k1 = pd.read_csv(output_dir / "figure4_atnr_K1_timeseries.csv")
    atnr_k10 = pd.read_csv(output_dir / "figure4_atnr_K10_timeseries.csv")

    for i, color in enumerate(colors):
        axes[0, 0].plot(atnr_k1.iloc[:, 0], atnr_k1.iloc[:, i + 1], color=color, lw=1)
        axes[1, 0].plot(atnr_k10.iloc[:, 0], atnr_k10.iloc[:, i + 1], color=color, lw=1)
        axes[0, 1].plot(t_eval, traj_k1[:, i], color=color, lw=1)
        axes[1, 1].plot(t_eval, traj_k10[:, i], color=color, lw=1)

    for ax in axes.ravel():
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_ylabel("Biomass")
    axes[1, 0].set_xlabel("Time")
    axes[1, 1].set_xlabel("Time")
    axes[0, 0].set_title("ATNr scaled, K=1")
    axes[1, 0].set_title("ATNr scaled, K=10")
    axes[0, 1].set_title("Python unscaled, K=1")
    axes[1, 1].set_title("Python unscaled, K=10")
    fig.tight_layout()
    fig.savefig(output_dir / "figure4_enrichment_python_vs_atnr.png")
    plt.close(fig)
    return traj_k1, traj_k10


def write_notes(output_dir: Path, fig3: pd.DataFrame, traj_k1: np.ndarray, traj_k10: np.ndarray) -> None:
    notes = {
        "scope": [
            "ATNr Figure 3 uses the unscaled-with-nutrients model; this Python model currently has no nutrient state variables.",
            "ATNr Figure 4 code uses create_model_Scaled; this Python model currently runs the unscaled logistic basal model.",
            "The Python scripts therefore test whether the local model can run the paper scenarios, not exact numerical equality.",
        ],
        "figure3_python_extinctions": fig3.to_dict(orient="records"),
        "figure4_python_final_biomass": {
            "K1": traj_k1[-1].tolist(),
            "K10": traj_k10[-1].tolist(),
        },
    }
    (output_dir / "reproduction_notes.json").write_text(json.dumps(notes, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        default="atnr_reproduction_output",
        help="Directory for generated inputs, CSVs, and figures.",
    )
    parser.add_argument(
        "--skip-r-reference",
        action="store_true",
        help="Reuse existing ATNr reference CSVs instead of regenerating them.",
    )
    parser.add_argument(
        "--fig3-python-t-max",
        type=float,
        default=10000.0,
        help=(
            "Time horizon for the Python no-nutrient Figure 3 analogue. "
            "ATNr reference still uses the paper value 100000."
        ),
    )
    parser.add_argument(
        "--fig3-python-dt",
        type=float,
        default=100.0,
        help="Output time step for the Python no-nutrient Figure 3 analogue.",
    )
    parser.add_argument(
        "--python-integrator",
        choices=["euler", "odeint"],
        default="euler",
        help="Integrator used for Python-model reproduction runs.",
    )
    parser.add_argument(
        "--euler-dt",
        type=float,
        default=1.0,
        help="Maximum internal step for the explicit Euler reproduction integrator.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not args.skip_r_reference:
        run_r_reference(output_dir)

    fig3 = reproduce_figure3(
        output_dir,
        args.fig3_python_t_max,
        args.fig3_python_dt,
        args.python_integrator,
        args.euler_dt,
    )
    traj_k1, traj_k10 = reproduce_figure4(
        output_dir,
        args.python_integrator,
        args.euler_dt,
    )
    write_notes(output_dir, fig3, traj_k1, traj_k10)

    print(f"Wrote ATNr reproduction outputs to: {output_dir.resolve()}")
    print("Key figures:")
    print(f"  {output_dir / 'figure3_temperature_extinctions_python_vs_atnr.png'}")
    print(f"  {output_dir / 'figure4_enrichment_python_vs_atnr.png'}")


if __name__ == "__main__":
    main()
