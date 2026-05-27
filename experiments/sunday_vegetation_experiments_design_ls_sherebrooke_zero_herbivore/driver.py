"""
driver.py -- Python driver for the Madingley C++ executable.

Drives the precompiled mac_exec/madingley binary via subprocess, with no
R in the loop. Writes all the CSV inputs the binary expects, runs spin,
forces zero biomass starting conditions, runs the simulation, and parses
outputs.

The argv layout, CSV formats, and quirks are documented in
notes/cpp_interface_report.md and notes/r_defaults_report.md.

Usage from a Python shell or notebook:
    import driver
    driver.full_run(years=3, output_dir="run_outputs")
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import rasterio


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

EXPT_DIR = Path(__file__).resolve().parent
MADINGLEY_PKG = EXPT_DIR / "context" / "MadingleyR" / "Package"
BINARY = MADINGLEY_PKG / "inst" / "mac_exec" / "madingley"
RASTERS_DIR = MADINGLEY_PKG / "inst" / "spatial_input_rasters"


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class RunConfig:
    """A single Madingley run configuration."""

    # Sherbrooke window -- 1 cell at centre (-72.5, 45.5).
    # The CSV crop uses strict bounds (x > xmin AND x < xmax), so a
    # 1-degree window picks the single cell whose centre is the half-degree
    # interior of that window.
    xmin: int = -73
    xmax: int = -72
    ymin: int = 45
    ymax: int = 46

    years: int = 3
    max_cohort: int = 500
    grid_size: int = 1

    output_dir: Path = field(default_factory=lambda: EXPT_DIR / "run_outputs")

    @property
    def input_dir(self) -> Path:
        return self.output_dir / "input"

    @property
    def spatial_dir(self) -> Path:
        return self.output_dir / "spatial_inputs" / "1deg"


# ---------------------------------------------------------------------------
# Spatial raster -> CSV conversion
# ---------------------------------------------------------------------------

# Single-layer rasters: tif filename -> CSV name (the C++ binary reads
# the CSV name, not the tif name). The hanpp tif is hanpp_2005.tif but
# is written as hanpp.csv.
SINGLE_LAYER_RASTERS: dict[str, str] = {
    "realm_classification.tif": "realm_classification",
    "land_mask.tif": "land_mask",
    "hanpp_2005.tif": "hanpp",
    "available_water_capacity.tif": "available_water_capacity",
    "Ecto_max.tif": "Ecto_max",
    "Endo_C_max.tif": "Endo_C_max",
    "Endo_H_max.tif": "Endo_H_max",
    "Endo_O_max.tif": "Endo_O_max",
}

# Multi-layer rasters: 12 monthly files each. Different tif naming
# conventions (some zero-padded, some not). Output CSVs are unpadded _1..12.
MULTI_LAYER_RASTERS: dict[str, str] = {
    "terrestrial_net_primary_productivity": "0padded",
    "near-surface_temperature": "0padded",
    "precipitation": "0padded",
    "ground_frost_frequency": "0padded",
    "diurnal_temperature_range": "unpadded",
}


def _raster_to_csv(tif_path: Path, csv_path: Path, cfg: RunConfig) -> int:
    """Convert one tif to var,x,y CSV cropped to the spatial window.

    Returns the number of rows written.

    The C++ binary uses strict inequalities (x > xmin AND x < xmax+1,
    y > ymin AND y < ymax+1). NA values are written as -999.
    Rows are sorted (y ascending, x ascending).
    """
    with rasterio.open(tif_path) as src:
        arr = src.read(1)  # first band
        nodata = src.nodata
        # cell centres from the affine transform
        rows, cols = np.indices(arr.shape)
        xs, ys = rasterio.transform.xy(src.transform, rows.ravel(), cols.ravel())
        xs = np.array(xs)
        ys = np.array(ys)
        vals = arr.ravel().astype(float)

    # mask nodata to NA, then write -999
    if nodata is not None:
        vals = np.where(vals == nodata, np.nan, vals)
    vals = np.where(np.isnan(vals), -999.0, vals)

    df = pd.DataFrame({"var": vals, "x": xs, "y": ys})

    # crop to window using R's strict bounds (x > xmin AND x < xmax).
    # The C++ binary uses a wider filter (x < xmax+1) but cropping the CSV
    # tighter here matches what MadingleyR would produce and keeps the
    # 1-cell Sherbrooke run as the plan intends.
    mask = (df["x"] > cfg.xmin) & (df["x"] < cfg.xmax) & \
           (df["y"] > cfg.ymin) & (df["y"] < cfg.ymax)
    df = df.loc[mask].sort_values(["y", "x"]).reset_index(drop=True)

    df.to_csv(csv_path, index=False)
    return len(df)


def write_spatial_inputs(cfg: RunConfig) -> dict[str, int]:
    """Write all 68 spatial CSV files to cfg.spatial_dir.

    Returns a {filename: nrows} dict for verification.
    """
    cfg.spatial_dir.mkdir(parents=True, exist_ok=True)
    out: dict[str, int] = {}

    for tif_name, csv_base in SINGLE_LAYER_RASTERS.items():
        tif = RASTERS_DIR / tif_name
        csv = cfg.spatial_dir / f"{csv_base}.csv"
        out[csv.name] = _raster_to_csv(tif, csv, cfg)

    for base, padding in MULTI_LAYER_RASTERS.items():
        for month in range(1, 13):
            if padding == "0padded":
                tif_name = f"{base}_{month:02d}.tif"
            else:
                tif_name = f"{base}_{month}.tif"
            tif = RASTERS_DIR / tif_name
            csv = cfg.spatial_dir / f"{base}_{month}.csv"
            out[csv.name] = _raster_to_csv(tif, csv, cfg)

    return out


# ---------------------------------------------------------------------------
# Control CSV writers (input/)
# ---------------------------------------------------------------------------

COHORT_DEF_HEADER = (
    "DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,"
    "DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,"
    "DEFINITION_Reproductive strategy,DEFINITION_Endo/Ectotherm,"
    "PROPERTY_Herbivory assimilation,PROPERTY_Carnivory assimilation,"
    "PROPERTY_Proportion suitable time active,PROPERTY_Minimum mass,"
    "PROPERTY_Maximum mass,PROPERTY_Initial number of GridCellCohorts,"
    "NOTES_group description"
)

# 9 terrestrial rows from get_default_cohort_def.R, with +1 applied to
# PROPERTY_Maximum.mass per write_cohort_def.R. We will replace the
# Initial_number column with floor(max_cohort/9) at write time.
COHORT_DEF_ROWS = [
    ("Heterotroph", "Herbivore", "All", "Terrestrial", "Mobile", "iteroparity", "Endotherm", 0.5, 0, 0.5, 1, 7000001),
    ("Heterotroph", "Carnivore", "All", "Terrestrial", "Mobile", "iteroparity", "Endotherm", 0, 0.8, 0.5, 5, 800001),
    ("Heterotroph", "Omnivore", "All", "Terrestrial", "Mobile", "iteroparity", "Endotherm", 0.38, 0.64, 0.5, 5, 150001),
    ("Heterotroph", "Herbivore", "All", "Terrestrial", "Mobile", "semelparity", "Ectotherm", 0.5, 0, 0.5, 0.04, 501),
    ("Heterotroph", "Carnivore", "All", "Terrestrial", "Mobile", "semelparity", "Ectotherm", 0, 0.8, 0.5, 0.08, 2001),
    ("Heterotroph", "Omnivore", "All", "Terrestrial", "Mobile", "semelparity", "Ectotherm", 0.36, 0.64, 0.5, 0.04, 2001),
    ("Heterotroph", "Herbivore", "All", "Terrestrial", "Mobile", "iteroparity", "Ectotherm", 0.5, 0, 0.5, 1, 100001),
    ("Heterotroph", "Carnivore", "All", "Terrestrial", "Mobile", "iteroparity", "Ectotherm", 0, 0.8, 0.5, 1.5, 100001),
    ("Heterotroph", "Omnivore", "All", "Terrestrial", "Mobile", "iteroparity", "Ectotherm", 0.36, 0.64, 0.5, 1.5, 55001),
]


def write_cohort_def(input_dir: Path, max_cohort: int = 500) -> None:
    n_per_fg = max_cohort // len(COHORT_DEF_ROWS)
    path = input_dir / "CohortFunctionalGroupDefinitions.csv"
    with path.open("w") as f:
        f.write(COHORT_DEF_HEADER + "\n")
        for row in COHORT_DEF_ROWS:
            cells = [str(c) for c in row] + [str(n_per_fg), "None"]
            f.write(",".join(cells) + "\n")


STOCK_DEF_HEADER = (
    "DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,"
    "DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,"
    "DEFINITION_Leaf strategy,PROPERTY_Herbivory assimilation,"
    "PROPERTY_Carnivory assimilation,PROPERTY_Proportion herbivory,"
    "PROPERTY_Individual mass"
)


def write_stock_def(input_dir: Path) -> None:
    path = input_dir / "StockFunctionalGroupDefinitions.csv"
    with path.open("w") as f:
        f.write(STOCK_DEF_HEADER + "\n")
        # Row 0: Deciduous, Row 1: Evergreen. NA fields blank.
        f.write("Autotroph,Photosynthesis,,Terrestrial,Sessile,Deciduous,,,,0\n")
        f.write("Autotroph,Photosynthesis,,Terrestrial,Sessile,Evergreen,,,,0\n")


SIM_PARAMS_ROWS = [
    ("RootDataDirectory", "../MadingleyData-master/NETCDF/"),
    ("RunParallel", "1"),
    ("ThreadNumber", "24"),
    ("TimeStepUnits", "month"),
    ("LengthOfSimulationInYears", "2"),  # ignored at runtime
    ("RunUntilStable", "1"),
    ("GridCellSize", "1"),
    ("ExtinctionThreshold", "1"),
    ("MaximumNumberOfCohorts", "1000"),
    ("DrawRandomly", "yes"),
    ("HumanNPPScenarioType", "none"),
    ("HumanNPPExtractionScale", "0"),
    ("HumanNPPScenarioDuration", "0"),
    ("BurninSteps", "0"),
    ("ImpactSteps", "0"),
    ("RecoverySteps", "0"),
    ("TimeStepStartExtinction", "1200"),
    ("StartBodyMass", "25000"),
    ("EndBodyMass", "21000"),
    ("StepBodyMass", "2000"),
    ("SelectCarnivores", "1"),
    ("SelectOmnivores", "0"),
    ("SelectHerbivores", "0"),
]


def write_sim_params(input_dir: Path) -> None:
    path = input_dir / "SimulationControlParameters.csv"
    with path.open("w") as f:
        f.write("Parameter,Value\n")
        for k, v in SIM_PARAMS_ROWS:
            f.write(f"{k},{v}\n")


# 78 mass bin values (matches write_mass_bin_def.R)
def _mass_bins() -> list[float]:
    bins = [1_000_000.0]
    # outer product 9..1 x 10^5..10^-2, column-major
    for power in range(5, -3, -1):
        for n in range(9, 0, -1):
            bins.append(n * 10.0 ** power)
    bins.extend([0.001, 0.0001, 0.00001, 0.000001, 0.0])
    return bins


def write_mass_bin_def(input_dir: Path) -> None:
    path = input_dir / "MassBinDefinitions.csv"
    with path.open("w") as f:
        f.write("Mass bin lower bound\n")
        for b in _mass_bins():
            # avoid scientific notation
            if b == 0:
                f.write("0\n")
            elif b >= 1:
                f.write(f"{b:.0f}\n")
            else:
                # plain decimal for small values
                f.write(f"{b:.10f}".rstrip("0").rstrip(".") + "\n")


def write_control_csvs(cfg: RunConfig) -> None:
    cfg.input_dir.mkdir(parents=True, exist_ok=True)
    write_cohort_def(cfg.input_dir, cfg.max_cohort)
    write_stock_def(cfg.input_dir)
    write_sim_params(cfg.input_dir)
    write_mass_bin_def(cfg.input_dir)


# ---------------------------------------------------------------------------
# Subprocess wrappers
# ---------------------------------------------------------------------------

def _common_argv(cfg: RunConfig,
                 output_ts_months: tuple[int, int, int, int],
                 gridout_bool: int,
                 cohort_csv: str,
                 stock_csv: str) -> list[str]:
    """Return argv[3..23] (the 21 positional args after the subcommand)."""
    # paths end with /
    out_dir_s = str(cfg.output_dir) + "/"
    input_dir_s = str(cfg.input_dir) + "/"
    spatial_dir_s = str(cfg.spatial_dir) + "/"

    return [
        str(cfg.xmin), str(cfg.xmax), str(cfg.ymin), str(cfg.ymax),  # argv 3-6
        out_dir_s,                                                  # argv 7
        str(output_ts_months[0]), str(output_ts_months[1]),         # argv 8-9
        str(output_ts_months[2]), str(output_ts_months[3]),         # argv 10-11
        str(gridout_bool),                                          # argv 12
        input_dir_s,                                                # argv 13
        str(cfg.max_cohort),                                        # argv 14
        cohort_csv,                                                 # argv 15
        stock_csv,                                                  # argv 16
        "0",                                                        # argv 17 RestartedFromTimeStepMonth
        spatial_dir_s,                                              # argv 18
        str(cfg.grid_size),                                         # argv 19
        "0",                                                        # argv 20 USE_HANPP
        "1",                                                        # argv 21 RunWithoutDispersal
        "1",                                                        # argv 22 RunInParallel
        "0",                                                        # argv 23 UseNonDefaultModelParameters
    ]


def run_spin(cfg: RunConfig, log_path: Path | None = None) -> subprocess.CompletedProcess:
    """Run the init/spin step. Writes FullCohortProperties_99999.csv and
    StockProperties_99999.csv into the output directory's subfolders."""
    argv = ["spin", "0"] + _common_argv(
        cfg,
        output_ts_months=(11988, 11988, 11988, 11988),  # never write intermediate
        gridout_bool=0,
        cohort_csv="none",
        stock_csv="none",
    )
    return _invoke(argv, log_path)


def run_simulation(cfg: RunConfig,
                   cohort_csv: Path,
                   stock_csv: Path,
                   log_path: Path | None = None) -> subprocess.CompletedProcess:
    """Run the actual simulation. Writes monthly outputs."""
    # write every month: thresholds = 0 for stock and grid; keep cohort outputs sparse
    argv = ["run", str(cfg.years)] + _common_argv(
        cfg,
        # bin cohort, full cohort, foodweb, stock
        # Write stock every month (0). Keep others sparse to limit disk.
        output_ts_months=(11988, 11988, 11988, 0),
        gridout_bool=1,
        cohort_csv=str(cohort_csv),
        stock_csv=str(stock_csv),
    )
    return _invoke(argv, log_path)


def _invoke(argv: list[str], log_path: Path | None) -> subprocess.CompletedProcess:
    cmd = [str(BINARY)] + argv
    res = subprocess.run(cmd, capture_output=True, text=True)
    if log_path is not None:
        log_path.write_text(
            "$ " + " ".join(cmd) + "\n"
            "--- STDOUT ---\n" + res.stdout +
            "\n--- STDERR ---\n" + res.stderr +
            f"\n--- exit={res.returncode} ---\n"
        )
    return res


# ---------------------------------------------------------------------------
# Zero biomass enforcement (between spin and run)
# ---------------------------------------------------------------------------

def write_zero_biomass_restart(cfg: RunConfig) -> tuple[Path, Path]:
    """Build a C.csv (empty, headerless) and S.csv (zero biomass, headerless)
    for the run step. Returns (cohort_path, stock_path)."""
    init_stock = cfg.output_dir / "stock_properties" / "StockProperties_99999.csv"
    if not init_stock.exists():
        raise FileNotFoundError(f"init step did not produce {init_stock}")

    # parse init stock file to learn how many cells / stocks there are
    df = pd.read_csv(init_stock)
    df["TotalBiomass"] = 0.0

    c_csv = cfg.input_dir / "C.csv"
    s_csv = cfg.input_dir / "S.csv"

    # empty cohort restart: zero rows, no header
    c_csv.write_text("")
    # stock restart: 3 cols, no header, sorted by GridcellIndex then FG
    df = df.sort_values(["GridcellIndex", "FunctionalGroupIndex"])
    df.to_csv(s_csv, index=False, header=False)
    return c_csv, s_csv


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

def parse_monthly_stocks(cfg: RunConfig) -> pd.DataFrame:
    """Read all stock_properties/StockProperties_NNNNN.csv files into one
    tidy frame: columns [month, GridcellIndex, FunctionalGroupIndex,
    TotalBiomass]."""
    stock_dir = cfg.output_dir / "stock_properties"
    frames = []
    for p in sorted(stock_dir.glob("StockProperties_*.csv")):
        # filename: StockProperties_NNNNN.csv
        stem = p.stem  # StockProperties_00000
        try:
            month = int(stem.split("_")[-1])
        except ValueError:
            continue
        if month == 99999:
            month_label = "final"
            continue  # also written at end; will be the last regular month
        d = pd.read_csv(p)
        d["month"] = month
        frames.append(d)
    if not frames:
        return pd.DataFrame(columns=["month", "GridcellIndex",
                                     "FunctionalGroupIndex", "TotalBiomass"])
    return pd.concat(frames, ignore_index=True).sort_values(["month", "GridcellIndex",
                                                              "FunctionalGroupIndex"])


def parse_monthly_stock_timeline(cfg: RunConfig) -> pd.DataFrame:
    """Read timelines/MonthlyStockBiomass.csv -- the aggregate total stock
    biomass per month (sum over all stocks and all cells)."""
    p = cfg.output_dir / "timelines" / "MonthlyStockBiomass.csv"
    if not p.exists():
        return pd.DataFrame(columns=["month", "TotalStockBiomass_kg"])
    vals = pd.read_csv(p, header=None).iloc[:, 0].astype(float)
    return pd.DataFrame({"month": np.arange(1, len(vals) + 1),
                         "TotalStockBiomass_kg": vals.values})


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------

def setup(cfg: RunConfig, clean: bool = True) -> None:
    """Create output_dir, write inputs."""
    if clean and cfg.output_dir.exists():
        shutil.rmtree(cfg.output_dir)
    cfg.output_dir.mkdir(parents=True, exist_ok=True)
    write_control_csvs(cfg)
    write_spatial_inputs(cfg)


def full_run(years: int = 3, output_dir: str | Path = "run_outputs",
             clean: bool = True) -> dict:
    """End-to-end: setup, spin, zero biomass, run. Returns summary dict."""
    out_dir = (EXPT_DIR / output_dir) if not Path(output_dir).is_absolute() else Path(output_dir)
    cfg = RunConfig(years=years, output_dir=out_dir)

    setup(cfg, clean=clean)

    spin_log = cfg.output_dir / "spin.log"
    spin = run_spin(cfg, spin_log)
    if spin.returncode != 0:
        raise RuntimeError(f"spin failed: see {spin_log}\nstderr: {spin.stderr[:2000]}")

    c_csv, s_csv = write_zero_biomass_restart(cfg)

    run_log = cfg.output_dir / "run.log"
    runp = run_simulation(cfg, c_csv, s_csv, run_log)
    if runp.returncode != 0:
        raise RuntimeError(f"run failed: see {run_log}\nstderr: {runp.stderr[:2000]}")

    return {
        "cfg": cfg,
        "spin_exit": spin.returncode,
        "run_exit": runp.returncode,
        "stocks": parse_monthly_stocks(cfg),
        "timeline": parse_monthly_stock_timeline(cfg),
    }


if __name__ == "__main__":
    import sys
    cmd = sys.argv[1] if len(sys.argv) > 1 else "all"
    cfg = RunConfig()
    if cmd == "setup":
        setup(cfg)
        print("setup done")
    elif cmd == "spin":
        spin = run_spin(cfg, cfg.output_dir / "spin.log")
        print(f"spin exit={spin.returncode}")
        print(spin.stderr[-1000:] if spin.returncode != 0 else "ok")
    elif cmd == "zero":
        c, s = write_zero_biomass_restart(cfg)
        print(f"wrote {c} and {s}")
    elif cmd == "run":
        c_csv = cfg.input_dir / "C.csv"
        s_csv = cfg.input_dir / "S.csv"
        runp = run_simulation(cfg, c_csv, s_csv, cfg.output_dir / "run.log")
        print(f"run exit={runp.returncode}")
    elif cmd == "all":
        result = full_run()
        print(f"spin exit={result['spin_exit']} run exit={result['run_exit']}")
        print(f"stocks: {len(result['stocks'])} rows; timeline: {len(result['timeline'])} months")
    else:
        print(f"unknown cmd: {cmd}")
