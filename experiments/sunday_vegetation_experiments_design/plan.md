PLAN -- VEGETATION-ONLY MADINGLEY RUN FROM PYTHON
Experiment: sunday_vegetation_experiments_design
Target: Sherbrooke (Quebec) grid cell, zero biomass start, short run
Date: 2026-05-25


GOAL

    Drive the MadingleyR C++ executable directly from a Python Jupyter
    notebook, with no R in the loop. The notebook prepares all input
    CSVs, calls the precompiled binary via subprocess, parses outputs,
    and plots autotroph (leaf) biomass dynamics.

    Scientific scope: a single grid cell over Sherbrooke, Quebec, with
    zero starting autotroph (stock) biomass AND zero heterotroph
    cohorts. We want to observe vegetation regrowth from bare ground
    under the local climate forcings (NPP, temperature, precipitation,
    frost frequency).


KEY DESIGN DECISIONS (from clarifying Q&A)

    1. C++ binding strategy:    subprocess to precompiled mac_exec
                                binary. No source recompile, no
                                pybind11 wrapping for v1.
    2. Zero biomass meaning:    stocks AND cohort table both emptied.
                                Pure autotroph-only spin-up.
    3. Simulation length:       short test, 1-5 years. We pick 3 years
                                for v1 to keep wall time tiny while
                                exercising at least three Miami-NPP
                                annual cycles.


SCIENTIFIC TARGET REGION

    Sherbrooke city centroid:   45.4042 N, -71.8929 W
    Spatial window for Madingley (xmin,xmax,ymin,ymax):
                                c(-73, -71, 45, 46)

    Rationale: default Madingley spatial inputs are 1 degree grids
    (~111 km / cell). The window above guarantees that the cell
    containing Sherbrooke (lon -72.5, lat 45.5 at 1 deg centres) is
    inside the cropped region. The C++ code uses strict > / < on the
    window, so we need at least a 2x2 deg box to keep a cell on
    either side of the centroid -- using a 2-deg window gives us a
    handful of cells with Sherbrooke at or near the centre.

    Note: 1 degree is coarse for "Sherbrooke city" specifically. The
    plan reports this honestly; if a finer resolution is needed we
    must rework the spatial input pipeline (out of scope here).


PIPELINE OVERVIEW (NOTEBOOK STRUCTURE)

    The notebook will have the following sections, each with a
    markdown cell explaining purpose, equations or file format,
    and a code cell that does the work. Cell order is the run order.

    SECTION 0   Title, scientific context, references (markdown only)
    SECTION 1   Environment setup (system + Python deps, paths)
    SECTION 2   Verify the precompiled MadingleyR binary
    SECTION 3   Define the spatial window for Sherbrooke
    SECTION 4   Read default spatial input rasters (terra-equivalent
                via rasterio); crop to Sherbrooke window
    SECTION 5   Write the 13 spatial input CSVs in the format the
                C++ code expects
    SECTION 6   Write the four control CSVs:
                  - CohortFunctionalGroupDefinitions.csv
                  - StockFunctionalGroupDefinitions.csv
                  - SimulationControlParameters.csv
                  - MassBinDefinitions.csv
    SECTION 7   Run the init step (the C++ "spin 0" command).
                This produces FullCohortProperties_99999.csv and
                StockProperties_99999.csv as starting state.
    SECTION 8   Force zero biomass starting conditions:
                  - empty the cohort CSV (header only)
                  - set stock TotalBiomass values to 0
    SECTION 9   Run the C++ "run" command for 3 years with the
                zeroed input state.
    SECTION 10  Parse outputs: stocks, cohorts, time series.
    SECTION 11  Plot terrestrial leaf biomass over time (evergreen
                and deciduous separately), plus NPP from the model
                state file for sanity checking.
    SECTION 12  Notes, caveats, follow-up questions.


C++ EXECUTABLE INTERFACE (FROM SOURCE INSPECTION)

    Binary path (precompiled, x86_64, runs under Rosetta on arm64):
      context/MadingleyR/Package/inst/mac_exec/madingley

    Two subcommands:
      ./madingley spin 0 <args>     init / spin-up
      ./madingley run   <years> <args>

    Argument order (matches madingley_init.R and madingley_run.R):
      1.  xmin xmax ymin ymax              (4 space-separated floats)
      2.  output_dir (quoted)              ends with /
      3.  output timesteps in months       (4 space-separated ints:
                                            bin cohort, full cohort,
                                            bin food-web, full stock)
      4.  gridout_bool                     0 for init, 1 for run
      5.  input_dir (quoted)               ends with /
      6.  max_cohort                       e.g. 500
      7.  cohort_csv path or "none"        full path to C.csv
      8.  stock_csv  path or "none"        full path to S.csv
      9.  start_t                          0
      10. spatial_inputs_dir (quoted)      ends with /1deg/
      11. grid_size                        1
      12. hanpp                            0 (off) for v1
      13. NoDispersal                      1 (single cell, no need)
      14. RunInParallel                    1
      15. model_params                     0 (use defaults)

    The binary writes results into:
      <output_dir>/cohort_properties/
      <output_dir>/stock_properties/
      <output_dir>/timeline_*.csv
      <output_dir>/foodweb_*.csv (if requested)


INPUT FILE FORMATS

    All under <out_dir>/input/:

      CohortFunctionalGroupDefinitions.csv
          header line + rows from get_default_cohort_def()
          for vegetation-only run: keep terrestrial definitions but
          we will provide an empty C.csv so no cohorts initialise

      StockFunctionalGroupDefinitions.csv
          two rows: Terrestrial Evergreen, Terrestrial Deciduous
          (the marine row is dropped in get_default_stock_def())

      SimulationControlParameters.csv
          23 rows of model parameters (see get_simulation_parameters)
          we set LengthOfSimulationInYears = 3, TimeStepUnits = month

      MassBinDefinitions.csv
          standard log-spaced body mass bins

    Under <out_dir>/spatial_inputs/1deg/:
      realm_classification.csv
      land_mask.csv
      hanpp.csv
      available_water_capacity.csv
      Ecto_max.csv
      Endo_C_max.csv
      Endo_H_max.csv
      Endo_O_max.csv
      terrestrial_net_primary_productivity_1.csv ... _12.csv
      near-surface_temperature_1.csv ... _12.csv
      precipitation_1.csv ... _12.csv
      ground_frost_frequency_1.csv ... _12.csv
      diurnal_temperature_range_1.csv ... _12.csv

      Each CSV has columns var,x,y (NA written as -999), one row per
      grid cell intersecting the spatial window, sorted by (y, x).


ZERO BIOMASS IMPLEMENTATION

    After SECTION 7 the init step produces:
      cohort_properties/FullCohortProperties_99999.csv
      stock_properties/StockProperties_99999.csv

    We use these as the C.csv and S.csv inputs to the run step, but
    modified:

      C.csv:  write the same header as FullCohortProperties_99999.csv
              with zero data rows (the C++ reader treats empty body
              as no cohorts; we will validate this in SECTION 8).
              Fallback if empty body is not accepted: write one row
              with CohortAbundance = 0 and IndividualBodyMass at the
              functional group default. A zero-abundance cohort is
              functionally absent.

      S.csv:  copy StockProperties_99999.csv but set the
              TotalBiomass column (column index per the stocks
              schema) to 0.0 for every row. Leaves the two terrestrial
              stocks defined (evergreen, deciduous) but with no
              standing biomass at t = 0.

    This is the cleanest realisation of "zero biomass starting
    conditions": vegetation pools exist as accounting categories but
    contain no mass; consumers do not exist.


DEPENDENCIES TO INSTALL

    System / pre-existing on this machine:
      macOS Darwin 25.2.0, Apple Silicon (arm64)
      Python 3.14.5 already present
      Rosetta 2 (needed to run the x86_64 madingley binary).
        Install with:
          softwareupdate --install-rosetta --agree-to-license

    Python (install into a project venv at experiments/.../.venv):
      jupyter            -- to run the notebook
      ipykernel          -- register the venv as a kernel
      numpy
      pandas             -- CSV IO, modifying init outputs
      rasterio           -- read the .tif spatial input rasters
                            (replaces R's terra)
      matplotlib         -- plotting biomass time series
      geopandas          (optional) -- only if we want to overlay
                            the Sherbrooke admin boundary on a map

    None of these need compilation; rasterio wheels exist for
    arm64 macOS on PyPI.

    No R is required at any point.


STEP-BY-STEP TASK LIST

    T1  Confirm context/MadingleyR is populated (init.sh has run).
        If not, run ./init.sh from the repo root.

    T2  Verify mac_exec/madingley runs:
          file context/.../mac_exec/madingley     (expect x86_64)
          chmod u+x context/.../mac_exec/madingley
          context/.../mac_exec/madingley           (should print help
                                                    or arg-count error)
        If "bad CPU type", install Rosetta 2.

    T3  Create a Python venv inside the experiment folder and install
        the dependency list above. Register it with Jupyter.

    T4  Create vegetation_sherbrooke.ipynb with the 12 sections from
        the pipeline overview. Each section starts with a markdown
        cell, then has its code cell.

    T5  Run the notebook end-to-end on a 3-year simulation with the
        zero-biomass C.csv and S.csv. Confirm:
          - init step exits 0 and writes the expected files
          - run step exits 0 and writes timeline_*.csv
          - parsed stocks show leaf biomass starting at 0 and growing

    T6  Write a short follow-up note in notes/ summarising the run.


OUTPUT METRICS TO REPORT IN THE NOTEBOOK

    Primary:
      - Terrestrial evergreen leaf biomass over time
      - Terrestrial deciduous leaf biomass over time
      - Total terrestrial autotroph biomass over time

    Secondary (sanity checks):
      - Monthly NPP at the Sherbrooke cell from input rasters
      - Evergreen fraction f_ever from frost frequency at this cell
      - Annual NPP from Miami-model variables (T and P inputs)


OPEN QUESTIONS / CAVEATS

    Q1  Sherbrooke is at ~45 N which is mid-latitude with strong
        winter frost. The Miami model and evergreen-fraction
        equations may give a deciduous-dominated cell, which means
        winter leaf biomass should drop sharply. Worth checking
        Ffrost values in the input rasters at this cell first.

    Q2  Default spatial rasters are 1 degree, which is a much larger
        area than "Sherbrooke city". For a city-scale run we would
        need higher-resolution climate inputs (out of scope for v1).

    Q3  The x86_64 binary will be slower under Rosetta than a native
        arm64 build. For a 3-year single-cell run this is irrelevant;
        if we move to longer runs we should compile from source.

    Q4  Cohort CSV with zero data rows: must be validated by reading
        Main.cpp / LoadCSVEnvironment.cpp once. The fallback
        (CohortAbundance=0) is safe regardless.

    Q5  Whether to additionally turn off the heterotroph functional
        groups via cohort_def (remove rows) instead of relying on an
        empty C.csv. Either should give a vegetation-only run;
        choosing one is a follow-up decision.


REFERENCES

    Harfoot et al. (2014) PLoS Biology 12(4): e1001841
    Hoeks et al. (2021) Global Ecology and Biogeography 30: 1922-1933
    notes/madingley_session_synthesis.md (prior session equations)
    CLAUDE.md (format and scientific context for this experiment)
