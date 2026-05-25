R DEFAULTS EXTRACTION REPORT
MadingleyR Package -- Terrestrial-Only Inputs for Python Notebook
Source: Package/R/ directory, read as documentation (no R executed)
Date: 2026-05-25


===========================================================================
SECTION 1  COHORT FUNCTIONAL GROUP DEFINITION TABLE (TERRESTRIAL ONLY)
===========================================================================

Source file: get_default_cohort_def.R
Write logic: write_cohort_def.R
Header source: get_input_header.R (index 1)

The original data frame has 19 rows. Rows 1-10 have DEFINITION_Realm =
"Marine", rows 11-19 have DEFINITION_Realm = "Terrestrial". The function
filters with df[df$DEFINITION_Realm=="Terrestrial",] and resets row names,
yielding 9 rows.

IMPORTANT: write_cohort_def.R applies a +1 offset to PROPERTY_Maximum.mass
before writing. The table below shows post-offset values (as written to the
CSV that the C++ binary reads).

The CSV file is written to:
  <out_dir>/input/CohortFunctionalGroupDefinitions.csv

The header line (one line, no quotes around values in the data rows because
fwrite uses plain CSV):

DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,DEFINITION_Reproductive strategy,DEFINITION_Endo/Ectotherm,PROPERTY_Herbivory assimilation,PROPERTY_Carnivory assimilation,PROPERTY_Proportion suitable time active,PROPERTY_Minimum mass,PROPERTY_Maximum mass,PROPERTY_Initial number of GridCellCohorts,NOTES_group description

Data rows (9 rows, terrestrial only, maximum mass has +1 applied):

Heterotroph,Herbivore,All,Terrestrial,Mobile,iteroparity,Endotherm,0.5,0,0.5,1,7000001,50,None
Heterotroph,Carnivore,All,Terrestrial,Mobile,iteroparity,Endotherm,0,0.8,0.5,5,800001,50,None
Heterotroph,Omnivore,All,Terrestrial,Mobile,iteroparity,Endotherm,0.38,0.64,0.5,5,150001,50,None
Heterotroph,Herbivore,All,Terrestrial,Mobile,semelparity,Ectotherm,0.5,0,0.5,0.04,501,50,None
Heterotroph,Carnivore,All,Terrestrial,Mobile,semelparity,Ectotherm,0,0.8,0.5,0.08,2001,50,None
Heterotroph,Omnivore,All,Terrestrial,Mobile,semelparity,Ectotherm,0.36,0.64,0.5,0.04,2001,50,None
Heterotroph,Herbivore,All,Terrestrial,Mobile,iteroparity,Ectotherm,0.5,0,0.5,1,100001,50,None
Heterotroph,Carnivore,All,Terrestrial,Mobile,iteroparity,Ectotherm,0,0.8,0.5,1.5,100001,50,None
Heterotroph,Omnivore,All,Terrestrial,Mobile,iteroparity,Ectotherm,0.36,0.64,0.5,1.5,55001,50,None

Dimensions: 9 rows x 14 columns.

Factor decoding notes (for reproducibility):
  DEFINITION_Nutrition source labels: ["Carnivore","Herbivore","Omnivore"] (1-indexed)
  DEFINITION_Diet labels:             ["All","AllSpecial","Planktivore"] (1-indexed)
  DEFINITION_Realm labels:            ["Marine","Terrestrial"] (1-indexed)
  DEFINITION_Mobility labels:         ["Mobile","Planktonic"] (1-indexed)
  DEFINITION_Reproductive strategy:   ["iteroparity","semelparity"] (1-indexed)
  DEFINITION_Endo/Ectotherm labels:   ["Ectotherm","Endotherm"] (1-indexed)

The PROPERTY_Initial.number.of.GridCellCohorts column is overwritten by
madingley_init.R before writing. The value used is:
  n_cohorts_per_fg = floor(max_cohort / nrow(terrestrial_rows))
For max_cohort=500 and 9 terrestrial rows: floor(500/9) = 55.
So the actual written value for all 9 rows will be 55, not 50.
The 50 in the source is only the default before the init override.

For the Python notebook's "run" step (Section 9 in the plan), if calling
madingley_run.R logic with no new cohorts, this column value does not matter
because cohort state is read from C.csv, not reinitialized.


===========================================================================
SECTION 2  STOCK FUNCTIONAL GROUP DEFINITION TABLE (TERRESTRIAL ONLY)
===========================================================================

Source file: get_default_stock_def.R
Write logic: write_stock_def.R
Header source: get_input_header.R (index 4)

The raw structure has 3 rows:
  Row 1: Marine,   Planktonic, Leaf strategy = "na",       mass = 1e-4
  Row 2: Terrestrial, Sessile, Leaf strategy = "Deciduous", mass = 0
  Row 3: Terrestrial, Sessile, Leaf strategy = "Evergreen", mass = 0

The filter is def[2:3,] (R 1-indexed), which keeps rows 2 and 3 only.
Leaf strategy factor labels: ["Deciduous","Evergreen","na"] (1-indexed).
Codes [3,1,2] map to: row1="na", row2="Deciduous", row3="Evergreen".

The CSV file is written to:
  <out_dir>/input/StockFunctionalGroupDefinitions.csv

Header line:

DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,DEFINITION_Leaf strategy,PROPERTY_Herbivory assimilation,PROPERTY_Carnivory assimilation,PROPERTY_Proportion herbivory,PROPERTY_Individual mass

Data rows (2 rows, both terrestrial):

Autotroph,Photosynthesis,NA,Terrestrial,Sessile,Deciduous,NA,NA,NA,0
Autotroph,Photosynthesis,NA,Terrestrial,Sessile,Evergreen,NA,NA,NA,0

Dimensions: 2 rows x 10 columns.

Notes:
  - DEFINITION_Diet, PROPERTY_Herbivory assimilation, PROPERTY_Carnivory
    assimilation, and PROPERTY_Proportion herbivory are all NA for stocks.
  - PROPERTY_Individual mass is 0 for both terrestrial stocks.
  - data.table::fwrite writes R NA as empty string in CSV by default. The
    Python notebook should write empty fields (not the string "NA") for
    these columns unless testing confirms the C++ reads "NA" as NA. The R
    source writes the data.frame directly with fwrite so those NA fields
    become empty in the CSV. Worth confirming with a hex dump of a real
    MadingleyR-generated file.


===========================================================================
SECTION 3  SIMULATION CONTROL PARAMETERS TABLE
===========================================================================

Source file: get_simulation_parameters.R
Write logic: write_simulation_parameters.R
Header source: get_input_header.R (index 3)

The factor structure has 23 rows. The factor levels are decoded below in
the same row order as the data frame (which is controlled by the integer
codes for the Parameter factor).

The CSV file is written to:
  <out_dir>/input/SimulationControlParameters.csv

Header line:

Parameter,Value

All 23 data rows (default values as stored in the R package):

RootDataDirectory,../MadingleyData-master/NETCDF/
RunParallel,1
ThreadNumber,24
TimeStepUnits,month
LengthOfSimulationInYears,2
RunUntilStable,1
GridCellSize,1
ExtinctionThreshold,1
MaximumNumberOfCohorts,1000
DrawRandomly,yes
HumanNPPScenarioType,none
HumanNPPExtractionScale,0
HumanNPPScenarioDuration,0
BurninSteps,0
ImpactSteps,0
RecoverySteps,0
TimeStepStartExtinction,1200
StartBodyMass,25000
EndBodyMass,21000
StepBodyMass,2000
SelectCarnivores,1
SelectOmnivores,0
SelectHerbivores,0

KEY PARAMETER NOTES:

  LengthOfSimulationInYears  Default is 2. Set to 3 for the Sherbrooke
                             experiment.

  TimeStepUnits              "month" -- the model runs in monthly steps.
                             A 3-year run = 36 time steps.

  MaximumNumberOfCohorts     1000 (default). madingley_init.R overrides
                             PROPERTY_Initial.number.of.GridCellCohorts
                             in the cohort table, but this parameter
                             caps total cohorts in the simulation.

  DrawRandomly               "yes" -- stochastic initialization.

  Output sample frequency:   There is no explicit output-sample-frequency
                             parameter in SimulationControlParameters.csv.
                             Output timestep control is passed as a
                             command-line argument to the C++ binary
                             (output_ts_months, described in Section 6).
                             For madingley_init, those are hardcoded to
                             999 years (effectively no intermediate output).
                             For madingley_run with output_timestep=c(0,0,0,0),
                             the code sets all four to (years-1)*12 months.
                             The four values are:
                               bin cohort, full cohort, bin food-web, full stock.

  TimeStepStartExtinction    1200 (months = 100 years). Not relevant for
                             a 3-year run.

  RootDataDirectory          Path is relative and points to a NETCDF
                             subdirectory, used internally by the C++
                             binary. For runs using CSV spatial inputs
                             this may be ignored, but include it anyway.


===========================================================================
SECTION 4  MASS BIN DEFINITIONS TABLE
===========================================================================

Source file: write_mass_bin_def.R

Construction logic:
  mass_bins = c(10^6,
                c(9:1 %o% 10^(5:-2)),   <- outer product, column-major order
                0.001, 0.0001, 0.00001, 0.000001,
                0)

The outer product 9:1 %o% 10^(5:-2) in R:
  Row vector: [9,8,7,6,5,4,3,2,1]
  Column vector (powers): [1e5,1e4,1e3,1e2,1e1,1e0,1e-1,1e-2]
  R stores matrices column-major, then c() reads them column by column.
  Column 1 (power 1e5): 9e5, 8e5, 7e5, ..., 1e5
  Column 2 (power 1e4): 9e4, 8e4, ..., 1e4
  ... and so on through power 1e-2.

Total bin count: 1 + 72 + 4 + 1 = 78 values.

The CSV file is written to:
  <out_dir>/input/MassBinDefinitions.csv

Header line:

Mass bin lower bound

All 78 data rows (written with scipen=999, so no scientific notation):

1000000
900000
800000
700000
600000
500000
400000
300000
200000
100000
90000
80000
70000
60000
50000
40000
30000
20000
10000
9000
8000
7000
6000
5000
4000
3000
2000
1000
900
800
700
600
500
400
300
200
100
90
80
70
60
50
40
30
20
10
9
8
7
6
5
4
3
2
1
0.9
0.8
0.7
0.6
0.5
0.4
0.3
0.2
0.1
0.09
0.08
0.07
0.06
0.05
0.04
0.03
0.02
0.01
0.001
0.0001
0.00001
0.000001
0

Note: The file uses data.table::fwrite with append=TRUE. The header is
written first via writeLines, then the data column is appended with no
column name and no row names. Python should write this as a plain list of
numbers after the header, one per line. For values <= 1e-4 write as
decimal, not scientific notation (scipen=999 suppresses sci notation).


===========================================================================
SECTION 5  MODEL PARAMETERS TABLE (DEFAULT VALUES)
===========================================================================

Source file: get_default_model_parameters.R

This table is NOT written to a CSV by default. The C++ binary receives
model_params=0 in the command-line arguments when using defaults, which
tells it to use its own internal defaults. The R data frame is only
written to the CSV when the user supplies a custom model_parameters
argument to madingley_run.

The format when passed as a non-default argument:
  model_params = paste(c(1, model_parameters$values), collapse=" ")
  The leading "1" tells C++ that custom params follow.
  The values are space-separated floats in the exact order below.

Total parameter rows: 85

Group ordering in the data frame (and if passed to C++):
  Activity_Parameters       (4 values)
  Dispersal_Parameters      (6 values)
  EatingCarnivory_Parameters (8 values) + EatingOmnivory_Parameters (1)
  EatingHerbivory_Parameters (9 values)
  MetabolismEctotherm_Parameters (7 values)
  MetabolismEndotherm_Parameters (6 values)
  MetabolismHeterotroph_Parameters (3 values)
  Mortality_Parameters      (5 values)
  Reproduction_Parameters   (4 values)
  VegetationModel_Parameters (32 values)

Full table (params,values,notes):

Activity_Parameters,6.61,Terrestrial Warming Tolerance Intercept
Activity_Parameters,1.6,Terrestrial Warming Tolerance Slope
Activity_Parameters,1.51,Terrestrial TSM Intercept
Activity_Parameters,1.53,Terrestrial TSM Slope
Dispersal_Parameters,0.0278,Diffusive Speed Body Mass Scalar
Dispersal_Parameters,0.48,Diffusive Speed Body Mass Exponent
Dispersal_Parameters,50000,Responsive Density Threshold Scaling
Dispersal_Parameters,0.0278,Responsive Speed Body Mass Scalar
Dispersal_Parameters,0.48,Responsive Speed Body Mass Exponent
Dispersal_Parameters,0.8,Starvation Dispersal Body Mass Threshold
EatingCarnivory_Parameters,0.5,Handling Time Scalar Terrestrial
EatingCarnivory_Parameters,0.7,Handling Time Exponent Terrestrial
EatingCarnivory_Parameters,0.5,Handling Time Scalar Marine
EatingCarnivory_Parameters,0.7,Handling Time Exponent Marine
EatingCarnivory_Parameters,1.0,Reference Mass
EatingCarnivory_Parameters,0.000006,Kill Rate Constant
EatingCarnivory_Parameters,1.0,Kill Rate Constant Mass Exponent
EatingCarnivory_Parameters,0.7,Feeding Preference Standard Deviation
EatingOmnivory_Parameters,0.1,Max Allowed Prey Ratio Omnivores
EatingHerbivory_Parameters,0.7,Handling Time Scalar Terrestrial
EatingHerbivory_Parameters,0.7,Handling Time Scalar Marine
EatingHerbivory_Parameters,0.7,Handling Time Exponent Terrestrial
EatingHerbivory_Parameters,0.7,Handling Time Exponent Marine
EatingHerbivory_Parameters,1.0,Reference Mass
EatingHerbivory_Parameters,1.0E-11,Herbivory Rate Constant
EatingHerbivory_Parameters,1.0,Herbivory Rate Mass Exponent
EatingHerbivory_Parameters,2.1,Attack Rate Exponent Terrestrial
EatingHerbivory_Parameters,0.1,Fraction Edible Stock Mass
MetabolismEctotherm_Parameters,0.88,Metabolism Mass Exponent
MetabolismEctotherm_Parameters,148984000000,Normalization Constant
MetabolismEctotherm_Parameters,0.69,Activation Energy
MetabolismEctotherm_Parameters,8.617e-5,Boltzmann Constant
MetabolismEctotherm_Parameters,41918272883,Normalization Constant BMR
MetabolismEctotherm_Parameters,0.69,Basal Metabolism Mass Exponent
MetabolismEctotherm_Parameters,0.036697248,Energy Scalar
MetabolismEndotherm_Parameters,0.7,Metabolism Mass Exponent
MetabolismEndotherm_Parameters,9.0809083973E+11,Normalization Constant
MetabolismEndotherm_Parameters,0.69,Activation Energy
MetabolismEndotherm_Parameters,8.617e-5,Boltzmann Constant
MetabolismEndotherm_Parameters,0.0366972,Energy Scalar
MetabolismEndotherm_Parameters,37.00,Endotherm Body Temperature
MetabolismHeterotroph_Parameters,0.71,Metabolism Mass Exponent
MetabolismHeterotroph_Parameters,0.69,Activation Energy
MetabolismHeterotroph_Parameters,8.617e-5,Boltzmann Constant
Mortality_Parameters,0.001,Background Mortality Rate
Mortality_Parameters,0.003,Senescence Mortality Rate
Mortality_Parameters,0.6,Starvation Logistic Inflection Point
Mortality_Parameters,0.05,Starvation Logistic Scaling Parameter
Mortality_Parameters,1.0,Starvation Maximum Starvation Rate
Reproduction_Parameters,1.5,Mass Ratio Threshold
Reproduction_Parameters,0.25,Mass Evolution Probability Threshold
Reproduction_Parameters,0.05,Mass Evolution Standard Deviation
Reproduction_Parameters,0.5,Semelparity Adult Mass Allocation
VegetationModel_Parameters,0.961644704,Max NPP (Miami model NPPmax)
VegetationModel_Parameters,0.237468183,T1NPP (Miami temperature param 1)
VegetationModel_Parameters,0.100597089,T2NPP (Miami temperature param 2)
VegetationModel_Parameters,0.001184101,PNPP (Miami precipitation param)
VegetationModel_Parameters,7.154615416,Fraction Structure Scalar
VegetationModel_Parameters,1.270782192,Fraction Evergreen A (a_fever)
VegetationModel_Parameters,-1.828591558,Fraction Evergreen B (b_fever)
VegetationModel_Parameters,0.844864063,Fraction Evergreen C (c_fever)
VegetationModel_Parameters,0.040273936,Evergreen Annual Leaf Mortality Slope
VegetationModel_Parameters,1.013070062,Evergreen Annual Leaf Mortality Intercept
VegetationModel_Parameters,0.020575964,Deciduous Annual Leaf Mortality Slope
VegetationModel_Parameters,-1.195235464,Deciduous Annual Leaf Mortality Intercept
VegetationModel_Parameters,0.04309283,Fine Root Mortality Rate Slope
VegetationModel_Parameters,-1.478393163,Fine Root Mortality Rate Intercept
VegetationModel_Parameters,0.139462774,Structural Mortality P2
VegetationModel_Parameters,-4.395910091,Structural Mortality P1
VegetationModel_Parameters,0.362742634,Leaf Carbon Fixation MaxFracStruct
VegetationModel_Parameters,0.388125108,Half Saturation Fire Mortality Rate
VegetationModel_Parameters,19.98393943,Scalar Fire Mortality Rate
VegetationModel_Parameters,1.148698636,NPP Half Saturation Fire Mortality Rate
VegetationModel_Parameters,8.419032427,NPP Scalar Fire Mortality Rate
VegetationModel_Parameters,0.01,Min Evergreen Annual Leaf Mortality
VegetationModel_Parameters,24.0,Max Evergreen Annual Leaf Mortality
VegetationModel_Parameters,0.01,Min Deciduous Annual Leaf Mortality
VegetationModel_Parameters,24.0,Max Deciduous Annual Leaf Mortality
VegetationModel_Parameters,0.01,Min Fine Root Mortality Rate
VegetationModel_Parameters,12.0,Max Fine Root Mortality Rate
VegetationModel_Parameters,1,Max Structural Mortality
VegetationModel_Parameters,0.001,Min Structural Mortality
VegetationModel_Parameters,2.0,Base Scalar Fire
VegetationModel_Parameters,0.00000226032940698105,Min Return Interval
VegetationModel_Parameters,0.476,Mass Carbon Per Mass Leaf Dry Matter

Cross-reference to key autotroph equations (from synthesis document
Section 3.9 and CLAUDE.md quick reference):
  Parameters 54-57 (VegetationModel indices 1-4):
    NPPmax = 0.961644704 (kg C / m2 / year)
    cP     = 0.237468183  (Miami T param, denominator intercept)
    mP     = 0.100597089  (Miami T param, denominator slope)
    rho    = 0.001184101  (Miami P param, exp decay rate)
  Parameters 58-61 (indices 5-8):
    a_fever = 1.270782192
    b_fever = -1.828591558
    c_fever = 0.844864063
    (fever = a_fever * Ffrost^2 + b_fever * Ffrost + c_fever)
  Parameter 57 (index 5): fStruct scalar = 7.154615416


===========================================================================
SECTION 6  SHELL COMMAND CONSTRUCTION (madingley_init.R AND madingley_run.R)
===========================================================================

HELPERS CALLED BY BOTH FUNCTIONS

  madingley_inputs()                  loads default spatial rasters from
                                      R package inst/ directory if
                                      spatial_inputs argument is not a list

  write_cohort_def(out_dir, cohort_def)
    -> writes input/CohortFunctionalGroupDefinitions.csv

  write_stock_def(out_dir, stock_def)
    -> writes input/StockFunctionalGroupDefinitions.csv

  write_simulation_parameters(out_dir)
    -> writes input/SimulationControlParameters.csv

  write_mass_bin_def(out_dir)
    -> writes input/MassBinDefinitions.csv

  write_spatial_inputs_to_temp_dir(spatial_inputs, XY_window, crop=TRUE,
                                   input_dir=out_dir, silenced)
    -> writes spatial_inputs/1deg/<name>.csv for single-layer rasters
    -> writes spatial_inputs/1deg/<name>_1.csv ... _12.csv for brick rasters

madingley_init ALSO calls:
  return_output_list_ini(cohort_def, stock_def, out_dir)
    -> reads init output files back into an R list

madingley_run ALSO calls:
  write_madingley_data_cohorts_stocks_to_temp_dir_fast(input_dir, madingley_data)
    -> writes input/C.csv and input/S.csv from the madingley_data object
  return_output_list_run(cohort_def, stock_def, out_dir, out_dir_name)


DIRECTORY CONVENTIONS

  madingley_init always writes to tempdir() regardless of any out_dir arg.
  The subdir created is:
    tempdir()/madingley_outs_DD_MM_YY_HH_MM_SS/

  madingley_run writes to the user-supplied out_dir (default tempdir()).
  The subdir created is:
    <out_dir>/madingley_outs_DD_MM_YY_HH_MM_SS/

  Within the run directory:
    input/                     all four control CSVs, C.csv, S.csv
    spatial_inputs/1deg/       all spatial input CSVs
    cohort_properties/         output: FullCohortProperties_*.csv
    stock_properties/          output: StockProperties_*.csv
    timeline_*.csv             timeseries output


INIT COMMAND (Mac/Linux, generalised form)

  ./madingley spin 0
    <xmin> <xmax> <ymin> <ymax>
    "<out_dir>/"
    999 999 999 999
    0
    "<input_dir>/"
    <max_cohort>
    none
    none
    0
    "<out_dir>/spatial_inputs/1deg/"
    <grid_size>
    0.0
    0
    1
    0

  Positional arguments for spin:
    arg 1-4:   spatial window (xmin xmax ymin ymax, space-separated)
    arg 5:     out_dir (quoted, trailing slash)
    arg 6-9:   output timesteps in months (bin cohort, full cohort,
               bin food-web, full stock) -- hardcoded 999 for init
    arg 10:    gridout_bool = 0 (no spatial output during init)
    arg 11:    input_dir (quoted, trailing slash)
    arg 12:    max_cohort
    arg 13:    cohort_csv = "none" (model generates its own)
    arg 14:    stock_csv = "none"
    arg 15:    start_t = 0
    arg 16:    spatial_inputs_location (quoted, ends in /1deg/)
    arg 17:    grid_size (1 for 1-degree rasters)
    arg 18:    NotYetAssignedVariableToPassToCPP = 0.0
    arg 19:    NoDispersal = 0
    arg 20:    RunInParallel = 1
    arg 21:    use_non_default_mp = 0


RUN COMMAND (Mac/Linux, generalised form)

  ./madingley run <years>
    <xmin> <xmax> <ymin> <ymax>
    "<out_dir>/"
    <output_ts_months[0]> <output_ts_months[1]> <output_ts_months[2]> <output_ts_months[3]>
    1
    "<input_dir>/"
    <max_cohort>
    "<input_dir>/C.csv"
    "<input_dir>/S.csv"
    0
    "<out_dir>/spatial_inputs/1deg/"
    <grid_size>
    <hanpp>
    <NoDispersal>
    <RunInParallel>
    <model_params>

  Differences from init:
    - "run" subcommand, not "spin 0"
    - <years> follows the subcommand before the spatial window
    - gridout_bool = 1 (enable spatial output)
    - cohort_csv and stock_csv are real file paths (C.csv and S.csv)
    - hanpp = 0 (off), 1 (fractional), or 2 (absolute)
    - model_params = 0 for defaults, or "1 v1 v2 ... v85" for custom

  output_ts_months derivation:
    If output_timestep=c(0,0,0,0) (default in madingley_run):
      output_ts_years = c(years-1, years-1, years-1, years-1)
      output_ts_months = output_ts_years * 12
    For years=3: output_ts_months = "24 24 24 24"
    This means output is written at month 24 (once, near the end).
    To get output at every time step set output_timestep to c(0,0,0,0)
    and understand the derivation above, or supply explicit values.

  For the Sherbrooke experiment with NoDispersal=1, years=3, max_cohort=500:
    ./madingley run 3
      -73 -71 45 46
      "<out_dir>/"
      24 24 24 24
      1
      "<input_dir>/"
      500
      "<input_dir>/C.csv"
      "<input_dir>/S.csv"
      0
      "<out_dir>/spatial_inputs/1deg/"
      1
      0
      1
      1
      0


===========================================================================
SECTION 7  SPATIAL RASTER TO CSV CONVERSION (write_spatial_inputs_to_temp_dir.R)
===========================================================================

Source files: write_spatial_inputs_to_temp_dir.R,
              check_and_rewrite_spatial_inputs.R

OUTPUT DIRECTORY

  All spatial CSVs are written to:
    <input_dir>/spatial_inputs/1deg/

  The function first calls unlink() on that directory recursively, then
  recreates it. Any previous contents are deleted every time.

SINGLE-LAYER RASTERS (class SpatRaster, nlyr == 1)

  The 8 single-layer inputs are:
    realm_classification, land_mask, hanpp, available_water_capacity,
    Ecto_max, Endo_C_max, Endo_H_max, Endo_O_max

  Conversion steps (directly from the R code):

    1. Extract all cell values as a vector: r[]
    2. Extract x,y coordinates for all cells: terra::xyFromCell(r, 1:ncell(r))
    3. Combine into a data.frame with columns (var, x, y)
    4. Replace NA values: df$var = ifelse(is.na(df$var), -999, df$var)
       NA encoding is -999 (integer, no decimal).
    5. Sort rows: df = df[with(df, order(y, x)),]
       Sort order is ascending y first (latitude south to north), then
       ascending x (longitude west to east) within each y value.
    6. Crop to spatial window using STRICT inequality:
       keep rows where x > XY_window[1] AND x < XY_window[2]
                   AND y > XY_window[3] AND y < XY_window[4]
       Note: strictly greater/less than, NOT >= or <=.
       Cell centres exactly on the window boundary are EXCLUDED.
    7. Write with write.csv(), which:
       - writes a header line: var,x,y
       - writes row indices as the first column (row.names = FALSE
         is NOT set, so default write.csv adds row numbers)

  IMPORTANT: The R call is write.csv(df, filename, row.names=FALSE).
  row.names=FALSE is explicitly set, so no row index column is written.
  The output CSV has exactly 3 columns: var, x, y.

  For the Sherbrooke spatial window c(-73, -71, 45, 46) with 1-degree
  rasters, cell centres are at half-degree offsets. With centres at lon
  -72.5 and lat 45.5, the window check is:
    x = -72.5 > -73? YES; x = -72.5 < -71? YES
    y = 45.5  > 45?  YES; y = 45.5  < 46?  YES
  So the single cell covering Sherbrooke is included.
  The 2-deg window (-73 to -71, 45 to 46) gives exactly 2 x 1 = 2 cells
  at 1 degree resolution (centres at -72.5 and -71.5, lat 45.5).
  Wait: lon range -73 to -71 spans 2 degrees -> centres at -72.5 only
  (since -71.5 < -71 is false). Lat range 45 to 46 spans 1 degree ->
  centre at 45.5 only. So exactly 1 cell is selected by that window.
  To get 2x2 cells use a window like c(-74,-71,44,47) or similar.

MULTI-LAYER RASTERS (class SpatRaster, nlyr == 12)

  The 5 multi-layer inputs are:
    terrestrial_net_primary_productivity,
    near-surface_temperature,
    precipitation,
    ground_frost_frequency,
    diurnal_temperature_range

  These are processed identically to single-layer rasters but once per
  layer (j = 1 to 12). Output filenames are:
    <name>_1.csv, <name>_2.csv, ..., <name>_12.csv

  Each file has the same 3-column structure: var,x,y.
  Same NA encoding (-999), same sort order (y asc, x asc), same strict
  crop, same row.names=FALSE.

COMPLETE FILE LIST UNDER spatial_inputs/1deg/

  realm_classification.csv
  land_mask.csv
  hanpp.csv
  available_water_capacity.csv
  Ecto_max.csv
  Endo_C_max.csv
  Endo_H_max.csv
  Endo_O_max.csv
  terrestrial_net_primary_productivity_1.csv
  terrestrial_net_primary_productivity_2.csv
  terrestrial_net_primary_productivity_3.csv
  terrestrial_net_primary_productivity_4.csv
  terrestrial_net_primary_productivity_5.csv
  terrestrial_net_primary_productivity_6.csv
  terrestrial_net_primary_productivity_7.csv
  terrestrial_net_primary_productivity_8.csv
  terrestrial_net_primary_productivity_9.csv
  terrestrial_net_primary_productivity_10.csv
  terrestrial_net_primary_productivity_11.csv
  terrestrial_net_primary_productivity_12.csv
  near-surface_temperature_1.csv
  near-surface_temperature_2.csv
  ... (through _12.csv)
  precipitation_1.csv ... precipitation_12.csv
  ground_frost_frequency_1.csv ... ground_frost_frequency_12.csv
  diurnal_temperature_range_1.csv ... diurnal_temperature_range_12.csv

  Total: 8 + 5*12 = 68 CSV files.

NOTE ON check_and_rewrite_spatial_inputs.R

  This function is a GUARD that checks whether the spatial rasters have
  changed since the last write. It computes sd and mean statistics over
  all raster layers and compares them against a cached file at:
    <tempdir>/spatial_inputs/spatial_inputs.csv

  It returns TRUE if the spatial inputs need to be rewritten, FALSE if
  they are unchanged. The calling code (check_temp_dir.R presumably) uses
  this return value to skip the expensive raster conversion step if nothing
  has changed.

  There is a bug in the source as committed: line 11 references
  "spatial_inputs2" (with a trailing "2") rather than "spatial_inputs".
  This means the function would error at runtime in R. It appears to be
  a dead-code path or an internal debugging artifact; in the main workflow
  madingley_init.R and madingley_run.R call write_spatial_inputs_to_temp_dir
  directly and do not call this check function.

  The Python notebook does not need to implement this caching guard at all.
  It can write the spatial CSVs unconditionally each run.


===========================================================================
SECTION 8  RESTART HEADERS (cohort and stock CSV formats for C.csv / S.csv)
===========================================================================

Source files: get_cohort_restart_header.R, get_stock_restart_header.R

COHORT RESTART CSV (C.csv) HEADER

GridcellIndex,FunctionalGroupIndex,JuvenileMass,AdultMass,IndividualBodyMass,CohortAbundance,LogOptimalPreyBodySizeRatio,BirthTimeStep,ProportionTimeActive,TrophicIndex,IndividualReproductivePotentialMass,MaturityTimeStep

12 columns. For the zero-cohort run, write only this header line and no
data rows, or write data rows with CohortAbundance = 0.

STOCK RESTART CSV (S.csv) HEADER

GridcellIndex,FunctionalGroupIndex,TotalBiomass

3 columns. For the zero-biomass start, write this header plus one row per
stock per grid cell with TotalBiomass = 0.

Example S.csv content for 1 grid cell and 2 stock functional groups:
  GridcellIndex,FunctionalGroupIndex,TotalBiomass
  0,0,0.0
  0,1,0.0

FunctionalGroupIndex 0 = Deciduous (first row in the filtered stock def),
FunctionalGroupIndex 1 = Evergreen (second row).


===========================================================================
SECTION 9  SPATIAL INPUT ORDER AND RASTER NAMES (madingley_inputs.R)
===========================================================================

The R package reads rasters from its installed inst/spatial_input_rasters/
directory. The Python notebook must supply equivalent data from the same
.tif files. The canonical names and their single vs. multi-layer status:

Single-layer (1 CSV each):
  1.  realm_classification
  2.  land_mask
  3.  hanpp
  4.  available_water_capacity
  5.  Ecto_max
  6.  Endo_C_max
  7.  Endo_H_max
  8.  Endo_O_max

Multi-layer 12-monthly (12 CSVs each, suffix _1 through _12):
  9.  terrestrial_net_primary_productivity
  10. near-surface_temperature
  11. precipitation
  12. ground_frost_frequency
  13. diurnal_temperature_range

The R code re-orders spatial_inputs by this list before writing. The
Python code must write files with these exact names.

The .tif files are in the R package installed directory at:
  <R_lib>/MadingleyR/spatial_input_rasters/

For the Python notebook they are accessible at the same location via the
MadingleyR package clone:
  context/MadingleyR/Package/inst/spatial_input_rasters/

(Verify the path with: ls context/MadingleyR/Package/inst/)


===========================================================================
SECTION 10  SUMMARY OF FILE LOCATIONS AND NAMES AS EXPECTED BY C++
===========================================================================

Given:
  out_dir    = /path/to/working/directory/madingley_outs_*/
  input_dir  = out_dir + "input/"
  sp_dir     = out_dir + "spatial_inputs/1deg/"

Files the C++ binary reads from input_dir:
  CohortFunctionalGroupDefinitions.csv
  StockFunctionalGroupDefinitions.csv
  SimulationControlParameters.csv
  MassBinDefinitions.csv
  C.csv   (cohort restart, only for "run" subcommand)
  S.csv   (stock restart,  only for "run" subcommand)

Files the C++ binary reads from sp_dir:
  68 CSV files listed in Section 7.

Files the C++ binary writes into out_dir:
  cohort_properties/FullCohortProperties_<timestep>.csv
  stock_properties/StockProperties_<timestep>.csv
  timeline_*.csv
  foodweb_*.csv (if gridout_bool = 1)
  BasicOutputs.csv (summary statistics)

The timestep suffix in output filenames is the month number (e.g., 99999
for the init spin-up output, or the month index for regular run outputs).


===========================================================================
END OF REPORT
===========================================================================
