CPP INTERFACE REPORT -- MADINGLEY EXECUTABLE
Prepared from source inspection of SourceCode/madingley/src/
Source version: CPP_VERSION = "2.02"  (Main.cpp line 8)
Date: 2026-05-25
Author: automated source read, no R execution

All line references are to files under
  experiments/sunday_vegetation_experiments_design/context/MadingleyR/SourceCode/madingley/src/


=============================================================
SECTION 1  EXACT ARGV PARSING IN MAIN.CPP
=============================================================

SUBCOMMAND POSITION

    argv[0]  binary path  (OS-supplied)
    argv[1]  subcommand string:  "spin", "run", or "version"
    argv[2]  years / spinup argument  (integer or "0")

    The subcommand "version" exits immediately after printing
    CPP_VERSION ("2.02") and does not read any further arguments.

SPIN VS RUN MEANING

    argv[1] = "spin"
        NumberOfSimulationYears = stoi(argv[2])  -- must be 0 for init
        OverwriteApplySpinUp is NOT set (remains false)
        TypeOfRun == "spin" branch in Madingley.cpp constructor
          --> WriteCohortProperties(..., 99999)
          --> WriteStockProperties(..., 99999)
        The model runs through MadingleyInitialisation to seed cohorts
        and stocks from the functional group definition CSVs (NOT from
        C.csv / S.csv; those are ignored in the spin branch).
        Outputs:
          cohort_properties/FullCohortProperties_99999.csv
          stock_properties/StockProperties_99999.csv

    argv[1] = "run"
        NumberOfSimulationYears = stoi(argv[2])  -- the actual year count
        OverwriteApplySpinUp = true  (Main.cpp line 174)
        TypeOfRun == "run" branch loads cohorts and stocks from
        C.csv and S.csv via readSpinUpStateCohort / readSpinUpStateStock.
        The model runs for (years * 12) months.
        Final state is also written as FullCohortProperties_99999.csv
        and StockProperties_99999.csv (Madingley.cpp lines 199-200).

    IMPORTANT NOTE ABOUT "spin 0":
        The R function madingley_init.R calls the binary as:
          ./madingley spin 0 <args>
        So argv[2] = "0" --> NumberOfSimulationYears = 0.
        Parameters::SetLengthOfSimulationInMonths receives 0.
        mLengthOfSimulationInMonths = 0 * 12 = 0.
        The run loop in Madingley::Run() is not entered (TypeOfRun=="spin"
        never calls Run()). The constructor does the init and writes the
        99999 files. This is correct and intended.

THE 15 ARGUMENTS AFTER THE SUBCOMMAND

    Below, "argv[N]" counts from argv[0] = binary name.
    The two arguments argv[1] and argv[2] consume the subcommand and
    years, so the "15 args" begin at argv[3].

    argv[3]   int       Xmin -- minimum longitude (integer degrees)
                        Source: Xmin = stoi(argv[3])  Main.cpp line 89
                        For Sherbrooke window: -73

    argv[4]   int       Xmax -- maximum longitude (integer degrees)
                        Source: Xmax = stoi(argv[4])  Main.cpp line 91
                        For Sherbrooke window: -71

    argv[5]   int       Ymin -- minimum latitude (integer degrees)
                        Source: Ymin = stoi(argv[5])  Main.cpp line 93
                        For Sherbrooke window: 45

    argv[6]   int       Ymax -- maximum latitude (integer degrees)
                        Source: Ymax = stoi(argv[6])  Main.cpp line 95
                        For Sherbrooke window: 46

    COORDINATE CONVENTION:
        The C++ code stores these as int and passes them to
        ExtractVarGrid which uses strict inequalities:
          test_env_vec1[i][1] > Xmin  AND  test_env_vec1[i][1] < NewXmax
          test_env_vec1[i][2] > Ymin  AND  test_env_vec1[i][2] < NewYmax
        where NewXmax = Xmax + 1, NewYmax = Ymax + 1.
        (LoadCSVEnvironment.cpp lines 105-106, 118-119)
        So a window of Xmin=-73, Xmax=-71, Ymin=45, Ymax=46 will
        include cells with longitude centres in (-73, -72) and latitude
        centres in (45, 47).  For 1-degree cells with centres at
        half-degree offsets (e.g. -72.5, 45.5), the window -73,-71,45,46
        includes exactly one cell: (-72.5, 45.5).  Use a wider window
        to include more cells.

    argv[7]   string    output_dir -- full path ending with /
                        Quotes are stripped by C++.
                        Example: "/tmp/madingley_outs_25_05_26_10_30_00/"
                        The directory must already exist before running.

    argv[8]   int       TimestepWritingBinnedCohortStatistics
                        Threshold month at which binned cohort stats
                        start being written. Checks: timeStep >= this.
                        R sets this to 999 * 12 = 11988 for spin.
                        For a 3-year run, set to 35 to write every
                        month from month 35 onward, or 0 to write always.

    argv[9]   int       TimestepWritingFullCohortProperties
                        Threshold month for writing FullCohortProperties_N.csv.
                        R sets this to 999 * 12 = 11988 for spin (never
                        written during run loop; only 99999 is written
                        at end).

    argv[10]  int       TimestepWritingPreyBinnedFoodwebConnections
                        Threshold month for foodweb and consumption files.
                        Also sets TimeStep_Months_Calc_FoodWeb.
                        Set to 999 * 12 = 11988 to disable.

    argv[11]  int       TimestepWritingStockProperties
                        Threshold month for writing StockProperties_N.csv.
                        Set to 0 to write every month.
                        The R init call uses 999 * 12 = 11988.

    NOTE ON OUTPUT TIMESTEP ARGS:
        The R code computes output_ts_months = output_ts_years * 12.
        For a 3-year run with output_timestep = c(0,0,0,0), the R code
        sets output_ts_years = c(2,2,2,2) (years - 1), so the months
        thresholds become c(24,24,24,24).  All four output types are
        written from month 24 onward (i.e. last year of the 3-year run).
        The stock timeline (MonthlyStockBiomass.csv) is written EVERY
        month regardless of this threshold.

    argv[12]  int       WriteGridProperties (0 or 1)
                        If 1, writes GridProperties_00000.csv through
                        GridProperties_00011.csv to grid_properties/.
                        R uses gridout_bool = 0 for init, 1 for run.

    argv[13]  string    ConfigurationDirectory -- path to input/ dir
                        Full path, quoted, ending with /.
                        Files read from this directory:
                          CohortFunctionalGroupDefinitions.csv
                          StockFunctionalGroupDefinitions.csv
                          SimulationControlParameters.csv  (read but
                            most parameters are ignored; see Section 3)
                          MassBinDefinitions.csv
                        Example: "/tmp/madingley_outs_.../input/"

    argv[14]  int       MaxCohortNumber
                        Maximum number of cohorts per grid cell before
                        merging is triggered.
                        R default: 500

    argv[15]  string    CohortCSVLocation
                        Full path to the cohort restart CSV (C.csv).
                        For spin: set to "none" (literal string).
                        For run: full path, e.g. input/C.csv.
                        The string "none" causes the run to crash on
                        the "run" subcommand if C.csv is not valid.

    argv[16]  string    StockCSVLocation
                        Full path to the stock restart CSV (S.csv).
                        For spin: "none".
                        For run: full path, e.g. input/S.csv.
                        NOTE: There is a dead substr() call on line 134
                        of Main.cpp that does nothing:
                          StockCSVLocation.substr(1, StockCSVLocation.size() - 2)
                        The return value is discarded. This does NOT
                        modify StockCSVLocation. The variable set is
                        the one passed to InputParameters.

    argv[17]  int       RestartedFromTimeStepMonth
                        Stored but appears unused in all run logic.
                        Set to 0.

    argv[18]  string    SpatialInputLocation
                        Full path to the spatial_inputs/1deg/ directory.
                        Must end with /.
                        Quotes are stripped.
                        Example: "/tmp/madingley_outs_.../spatial_inputs/1deg/"

    argv[19]  double    GridcellSize
                        Spatial resolution in degrees. Only 1 is
                        supported (the R code enforces grid_size <= 1).
                        Passed to Parameters::SetGridCellSize.

    argv[20]  double    USE_HANPP
                        0.0   -- no HANPP applied (plan.md says 0)
                        ~1.0  -- fractional raster (0.7 < val < 1.4)
                        ~2.0  -- absolute gC/m^2/year raster (val > 1.7)
                        The R variable is named "hanpp" in madingley_run.R.

    argv[21]  int       RunWithoutDispersal
                        1 = dispersal disabled (single cell, no movement)
                        0 = dispersal enabled
                        Plan.md correctly sets this to 1 for single cell.

    argv[22]  int       RunInParallel
                        1 = use OpenMP parallel loop
                        0 = serial execution
                        Plan.md sets this to 1.

    argv[23]  int       UseNonDefaultModelParameters
                        0 = use hardcoded defaults in TerrestrialCarbon.cpp
                            and all other process files
                        1 = read the additional parameter block from
                            argv[24] onward (see Section 6 below)
                        Plan.md sets this to 0, which is correct.

    TOTAL POSITIONAL ARGUMENTS (including binary and subcommand):
        argv[0..23] = 24 arguments for UseNonDefaultModelParameters = 0.
        For UseNonDefaultModelParameters = 1, additional parameters
        follow starting at argv[24] (see Section 6).

FULL COMMAND LINE TEMPLATE (spin, Mac)

    ./madingley spin 0 \
      <Xmin> <Xmax> <Ymin> <Ymax> \
      "<output_dir>/" \
      999 999 999 999 \
      0 \
      "<input_dir>/" \
      500 \
      none none \
      0 \
      "<spatial_inputs_dir>/1deg/" \
      1 \
      0 1 1 0

FULL COMMAND LINE TEMPLATE (run, Mac)

    ./madingley run <years> \
      <Xmin> <Xmax> <Ymin> <Ymax> \
      "<output_dir>/" \
      <ts_bin_cohort> <ts_full_cohort> <ts_foodweb> <ts_full_stock> \
      1 \
      "<input_dir>/" \
      500 \
      "<input_dir>/C.csv" "<input_dir>/S.csv" \
      0 \
      "<spatial_inputs_dir>/1deg/" \
      1 \
      0 1 1 0

    Note: The output_ts_* values are in months. R uses (years-1)*12
    for all four if output_timestep argument is c(0,0,0,0).

VALIDATION / DEFAULTS

    The binary does NO validation of argument count. Fewer than 24
    arguments causes undefined behaviour (stoi/stod on uninitialized
    argv slots). There is no argc check or error message for wrong
    argument count.

    Xmin, Xmax, Ymin, Ymax are stored as int (not float). The R code
    rounds them with as.integer() implicitly via paste. Pass integer
    degrees only.


=============================================================
SECTION 2  SPATIAL INPUT CSV FORMAT
=============================================================

DIRECTORY

    All spatial input CSVs live in:
      <spatial_inputs_dir>/  (argv[18])
    which the R package always sets to
      <output_dir>/spatial_inputs/1deg/

    The C++ code appends filenames directly to the SpatialInputLocation
    string without adding a separator, so SpatialInputLocation MUST end
    with /.

EXACT FILENAMES READ BY MAIN.CPP

    Single-layer files (read once, no month suffix):
      realm_classification.csv         (Main.cpp line 328)
      available_water_capacity.csv      (Main.cpp line 329)
      land_mask.csv                     (Main.cpp line 330)
      hanpp.csv                         (Main.cpp line 331)
      Ecto_max.csv                      (Main.cpp line 354)
      Endo_C_max.csv                    (Main.cpp line 355)
      Endo_H_max.csv                    (Main.cpp line 356)
      Endo_O_max.csv                    (Main.cpp line 357)

    Monthly files (read 12 times each, suffix _1 through _12):
      diurnal_temperature_range_1.csv  ...  diurnal_temperature_range_12.csv
      ground_frost_frequency_1.csv     ...  ground_frost_frequency_12.csv
      near-surface_temperature_1.csv   ...  near-surface_temperature_12.csv
      precipitation_1.csv              ...  precipitation_12.csv
      terrestrial_net_primary_productivity_1.csv  ...  _12.csv

    Total: 8 + 5*12 = 68 CSV files.
    All filenames are EXACT -- case-sensitive on Linux and Mac.
    Note the hyphen in "near-surface_temperature" (not underscore).

FILE FORMAT (from LoadCSVEnvironment.cpp)

    Delimiter:    comma (",")  -- line 49: getline(row_stream, column, ',')
    Header:       YES, exactly one header line, skipped by the reader
                  (loop starts at i=1:  line 65  "start from 1 not 0")
                  The header content is not inspected -- any string works.
                  The R code writes:  "var,x,y"  (write_spatial_inputs_to_temp_dir.R line 53)
                  Use exactly that header.
    Column order: col[0] = var (the environmental value)
                  col[1] = x   (longitude, decimal degrees)
                  col[2] = y   (latitude, decimal degrees)
                  This is CONFIRMED by the comment in ExtractVarGrid:
                  "col 0: var, col 1: long, col 2: lat"  (line 98-100)
    Row ordering: sorted by (y ascending, x ascending) -- the R code
                  does df = df[with(df, order(y, x)),]  before writing.
                  The C++ reader does NOT enforce ordering; it reads rows
                  sequentially and filters by the spatial window.
                  However, the index tracking code in ExtractVarGrid
                  assumes sorted order for long_index_counter and
                  lat_index_counter (lines 126-137). Use (y, x) sort.
    NA encoding:  -999  (the R code writes: ifelse(is.na(df$var),-999,df$var))
                  The C++ code converts -999 to Constants::cMissingValue
                  (-9999 as int) at line 69-73:
                    if(stod(items[i][0]) == -999)
                      items_doubles[i-1][0] = Constants::cMissingValue;
    Extension:    .csv  (confirmed: all filenames above end in .csv)

COORDINATE COVERAGE

    The file does NOT need to contain the full global grid. It only needs
    to contain cells that survive the window filter:
      cell.x > Xmin  AND  cell.x < Xmax+1
      cell.y > Ymin  AND  cell.y < Ymax+1
    Cells outside this range are silently ignored.
    The file CAN contain the full global grid -- the reader reads all rows
    and filters. For efficiency, crop to the window before writing.

MISSING VALUE NOTE

    cMissingValue in Constants.h is -9999 (int), not -999.
    The input NA marker -999 is converted to -9999 inside Read_Env_CSV.
    The x and y columns are NEVER checked against -999 (lines 75-76
    use stod without checking), so do not use -999 for coordinates.


=============================================================
SECTION 3  CONTROL CSV FORMATS UNDER <input_dir>/input/
=============================================================

The four control files are read from ConfigurationDirectory (argv[13]),
which is the <output_dir>/input/ directory.

FILE PATHS (from Constants.h and FunctionalGroupDefinitions.cpp):
  CohortFunctionalGroupDefinitions.csv  -- cCohortDefinitionsFileName
  StockFunctionalGroupDefinitions.csv   -- cStockDefinitionsFileName
  MassBinDefinitions.csv                -- cMassBinDefinitionsFileName
  SimulationControlParameters.csv       -- cInputParametersFileName (read but mostly ignored)

----------------------------------------------------------------
3A. CohortFunctionalGroupDefinitions.csv
----------------------------------------------------------------

EXACT HEADER LINE (from get_input_header.R, index 1):
  DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,DEFINITION_Reproductive strategy,DEFINITION_Endo/Ectotherm,PROPERTY_Herbivory assimilation,PROPERTY_Carnivory assimilation,PROPERTY_Proportion suitable time active,PROPERTY_Minimum mass,PROPERTY_Maximum mass,PROPERTY_Initial number of GridCellCohorts,NOTES_group description

The FunctionalGroupDefinitions reader (FunctionalGroupDefinitions.cpp)
splits headers by comma, then by underscore to get category and property.
The two recognised categories are "definition" and "property"; "notes"
columns are parsed but not stored in any lookup used by the model.
All values are lowercased before storage (line 65: tolower).

COLUMN MEANINGS:

  DEFINITION_Heterotroph/Autotroph
      "Heterotroph" for all cohort rows.
      The slash is part of the trait name; it is split on underscore at
      the first underscore, so the category is "definition" and the
      key is "heterotroph/autotroph". Values are lowercased to
      "heterotroph".

  DEFINITION_Nutrition source
      "Carnivore", "Herbivore", or "Omnivore"  (lowercased in storage).

  DEFINITION_Diet
      "All", "AllSpecial", or "Planktivore" (lowercased).

  DEFINITION_Realm
      "Terrestrial" or "Marine" (lowercased).
      Cohorts with Realm != current grid cell realm are not seeded.

  DEFINITION_Mobility
      "Mobile" or "Planktonic" (lowercased).

  DEFINITION_Reproductive strategy
      "iteroparity" or "semelparity" (lowercased).
      Affects which spatial body mass maximum layer is used during init.

  DEFINITION_Endo/Ectotherm
      "Endotherm" or "Ectotherm" (lowercased).

  PROPERTY_Herbivory assimilation    (float)
  PROPERTY_Carnivory assimilation    (float)
  PROPERTY_Proportion suitable time active  (float)
  PROPERTY_Minimum mass              (float, grams)
  PROPERTY_Maximum mass              (float, grams)
      NOTE: write_cohort_def.R adds 1 to Maximum mass before writing:
        cohort_def$PROPERTY_Maximum.mass = cohort_def$PROPERTY_Maximum.mass + 1
      This is a quirk in the R package -- the Python notebook must
      replicate this if it writes cohort_def directly.

  PROPERTY_Initial number of GridCellCohorts  (int)
      Controls how many cohorts of this functional group to seed per
      cell during spin. Set to 0 to seed no cohorts.
      For the vegetation-only run, setting ALL rows to 0 means no
      cohorts are seeded during spin, so the C.csv is empty.

  NOTES_group description
      Free text, not used in model calculations.

DEFAULT VALUES (from get_default_cohort_def.R filtered to Terrestrial):
  9 terrestrial functional groups (rows 11-19 of the full 19-row table).
  The default PROPERTY_Initial.number.of.GridCellCohorts for Terrestrial
  rows is 50 each. The R init function overrides this based on max_cohort:
    n_cohorts_per_fg = floor(max_cohort / nrow(terrestrial_cohort_rows))
  With max_cohort=500 and 9 terrestrial FGs: 500/9 = 55 each.

DELIMITER: comma. No quoting expected.

----------------------------------------------------------------
3B. StockFunctionalGroupDefinitions.csv
----------------------------------------------------------------

EXACT HEADER LINE (from get_input_header.R, index 4):
  DEFINITION_Heterotroph/Autotroph,DEFINITION_Nutrition source,DEFINITION_Diet,DEFINITION_Realm,DEFINITION_Mobility,DEFINITION_Leaf strategy,PROPERTY_Herbivory assimilation,PROPERTY_Carnivory assimilation,PROPERTY_Proportion herbivory,PROPERTY_Individual mass

COLUMN MEANINGS:

  DEFINITION_Heterotroph/Autotroph
      "Autotroph" for all stock rows.

  DEFINITION_Nutrition source
      "Photosynthesis" for all stock rows.

  DEFINITION_Diet
      NA (empty) -- the R code writes NA as blank.
      Leave the column blank or write empty string.

  DEFINITION_Realm
      "Terrestrial" or "Marine".
      Only Terrestrial stocks are seeded in non-marine cells (Stock.cpp
      line 20: checks stockDefinitions.GetTraitNames("Realm",fg)=="terrestrial")

  DEFINITION_Mobility
      "Sessile" for terrestrial; "Planktonic" for marine.

  DEFINITION_Leaf strategy
      "Deciduous" or "Evergreen" (lowercased in storage).
      This is the key trait that determines which TerrestrialCarbon
      function branch is called. The two terrestrial stocks should be:
        Row 1: Realm=Terrestrial, Leaf strategy=Deciduous
        Row 2: Realm=Terrestrial, Leaf strategy=Evergreen

  PROPERTY_Herbivory assimilation    NA (blank)
  PROPERTY_Carnivory assimilation    NA (blank)
  PROPERTY_Proportion herbivory      NA (blank)
  PROPERTY_Individual mass           float (0 for terrestrial stocks)

DEFAULT VALUES (from get_default_stock_def.R, rows 2-3 of 3-row table):
  Row 1: Autotroph, Photosynthesis, NA, Terrestrial, Sessile, Deciduous, NA, NA, NA, 0
  Row 2: Autotroph, Photosynthesis, NA, Terrestrial, Sessile, Evergreen, NA, NA, NA, 0

  The marine row (row 1 of the full table) is DROPPED by the R package
  (def = def[2:3,]). If you include the marine row the model will read
  it but not seed it in terrestrial cells.

DELIMITER: comma.

----------------------------------------------------------------
3C. SimulationControlParameters.csv
----------------------------------------------------------------

EXACT HEADER LINE (from get_input_header.R, index 3):
  Parameter,Value

This is a two-column CSV with 23 data rows. However, the current
C++ code has effectively DISABLED reading of this file. In
FileReader::ReadInputParameters (FileReader.cpp lines 62-64):
  bool success = true;  // ReadTextFile call is COMMENTED OUT
  if(success == true) success = Parameters::Get()->Initialise(mMetadata);
  return success;

The mMetadata is always empty, so Parameters::Initialise is called
with empty data and sets hardcoded defaults in Parameters.cpp lines 41-67:
  TimeStepUnits = "months"
  PlanktonSizeThreshold = 0.01
  ApplyModelSpinup = 10
  ExtinctionThreshold = 1
  DrawRandomly = "yes"
  etc.

The critical parameter LengthOfSimulationInYears is set from:
  InputParameters::Get()->GetNumberOfSimulationYears()
which comes from argv[2], NOT from SimulationControlParameters.csv.

CONCLUSION: SimulationControlParameters.csv is READ (the file must
exist and must have the correct header) but its values are IGNORED.
The binary reads LengthOfSimulationInYears from the command line.
You still MUST write this file because MadingleyInitialisation calls
ReadInitialisationFiles() which opens it. If it is missing, the model
may print an error but does not exit (the ReadTextFile call is commented
out, so the file is never actually opened in the current code path).

NONETHELESS: Write the file with the correct 23 rows to be safe.
The R function get_simulation_parameters() produces the default content.

KNOWN ROWS (from get_simulation_parameters.R):
  RootDataDirectory, ../MadingleyData-master/NETCDF/
  DrawRandomly, yes
  LengthOfSimulationInYears, 1       <-- IGNORED; use argv[2]
  TimeStepUnits, month               <-- IGNORED
  GridCellSize, 2                    <-- IGNORED; use argv[19]
  RunUntilStable, 1                  <-- IGNORED
  ExtinctionThreshold, 1             <-- set in Parameters.cpp
  MaximumNumberOfCohorts, 1000       <-- IGNORED; use argv[14]
  ImpactSteps, 1000                  <-- stored but unused
  RecoverySteps, 1                   <-- stored but unused
  BurninSteps, 1                     <-- stored but unused
  HumanNPPScenarioType, none         <-- stored but unused in run
  HumanNPPExtractionScale, 0         <-- stored but unused in run
  HumanNPPScenarioDuration, 0        <-- stored but unused in run
  ThreadNumber, 1                    <-- stored
  RunParallel, 1                     <-- IGNORED; use argv[22]
  TimeStepStartExtinction, 1200      <-- stored but unused
  StartBodyMass, 25000               <-- stored but unused
  StepBodyMass, 21000                <-- stored but unused
  EndBodyMass, 24                    <-- stored but unused
  SelectCarnivores, 1                <-- stored but unused
  SelectHerbivores, 1                <-- stored but unused
  SelectOmnivores, 1                 <-- stored but unused

DELIMITER: comma.

----------------------------------------------------------------
3D. MassBinDefinitions.csv
----------------------------------------------------------------

EXACT HEADER LINE (from get_input_header.R, index 2):
  Mass bin lower bound

Single column: one float per row, sorted ascending.
The MassBinsHandler reads with operator>> which splits on whitespace.

DEFAULT VALUES (from write_mass_bin_def.R):
  10^6, then 9:1 x 10^5:-2, then 0.001, 0.0001, 0.00001, 0.000001, 0

In practice the R code writes:
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
  1
  0.1
  0.01
  0.001
  0.0001
  0.00001
  0.000001
  0

This file is used only for binned output statistics. It is required
but the mass values it contains do not affect core vegetation dynamics.


=============================================================
SECTION 4  COHORT AND STOCK RESTART FILES
=============================================================

----------------------------------------------------------------
4A. FullCohortProperties_<NNNNN>.csv
----------------------------------------------------------------

EXACT HEADER LINE (from WriteModelState.cpp lines 79-95):
  GridcellIndex,FunctionalGroupIndex,JuvenileMass,AdultMass,IndividualBodyMass,CohortAbundance,LogOptimalPreyBodySizeRatio,BirthTimeStep,ProportionTimeActive,TrophicIndex,IndividualReproductivePotentialMass,MaturityTimeStep,IsAdult,AgeMonths,TimeStepsJuviline,TimeStepsAdult

16 columns. Delimiter: comma.

NOTE: The R helper function get_cohort_restart_header() returns only
12 columns (up to MaturityTimeStep). The ACTUAL file written by the
C++ binary has 16 columns. The last 4 are:
  IsAdult            (int, 0 or 1)
  AgeMonths          (int)
  TimeStepsJuviline  (int)
  TimeStepsAdult     (int)

The R restart reader (write_madingley_data_cohorts_stocks_to_temp_dir_fast.R)
writes the full cohorts data.frame including all 16 columns using
fwrite(..., col.names=FALSE), so no header is written to C.csv.

C.csv WRITTEN BY THE R PACKAGE:
  NO HEADER.  16 data columns.
  Format of C.csv passed to the C++ run command: headerless, 16-col CSV.

C.csv READER IN C++ (readSpinUpStateCohort):
  Allocates v with 16 columns: std::vector<std::vector<std::string>> v(16,...)
  Reads by scanning comma-delimited tokens. Does NOT skip a header.
  Therefore C.csv MUST NOT have a header row.

COLUMN MAPPING IN readSpinUpStateCohort / SeedCohortsApplySpinUpFast:
  CohortData[0]   GridcellIndex             (stoul)
  CohortData[1]   FunctionalGroupIndex      (stoi)
  CohortData[2]   JuvenileMass              (stod)
  CohortData[3]   AdultMass                 (stod)
  CohortData[4]   IndividualBodyMass        (stod)
  CohortData[5]   CohortAbundance           (stod)
  CohortData[6]   LogOptimalPreyBodySizeRatio (stod, then exp() applied!)
  CohortData[7]   BirthTimeStep             (stod cast to unsigned)
  CohortData[8]   ProportionTimeActive      (stod)
  CohortData[9]   TrophicIndex              (stod)
  CohortData[10]  IndividualReproductivePotentialMass (stod)
  CohortData[11]  MaturityTimeStep          (stod cast to unsigned)
  CohortData[12]  IsAdult                   (stoi)
  CohortData[13]  AgeMonths                 (stoi)
  CohortData[14]  TimeStepsJuviline         (stoi)
  CohortData[15]  TimeStepsAdult            (stoi)

CRITICAL NOTE ON LogOptimalPreyBodySizeRatio:
  The writer (WriteModelState.cpp line 108) writes c->mLogOptimalPreyBodySizeRatio
  directly. The reader (SeedCohortsApplySpinUpFast line 562) applies exp():
    optimalPreyBodySizeRatio = exp(stod(CohortData[6][jj]))
  So the file stores the LOG of the optimal prey body size ratio, and exp()
  is applied on read. This is consistent: write log, read and apply exp.

EMPTY COHORT FILE BEHAVIOUR:
  If C.csv has zero data rows (file exists, no content, or header only),
  readSpinUpStateCohort returns a vector of 16 empty vectors.
  SeedCohortsApplySpinUpFast: nrows = CohortData[0].size() = 0.
  The for loop does not execute.
  HOWEVER: The code then loops over all cells and inserts "ghost" cohorts
  for any cell not represented in C.csv (lines 625-643):
    Cohort* NewCohort = new Cohort(cell, 1, 1, 1, 1, 1, 1, 0, 0, mNextCohortID, 2, 0, UINT_MAX, 0, 0, 0, 0)
  So a completely empty C.csv causes one ghost cohort to be seeded per
  terrestrial grid cell. These ghosts have:
    functionalGroup = 1, JuvenileMass = 1g, AdultMass = 1g,
    IndividualBodyMass = 1g, Abundance = 1, FG=1 (index depends on your
    cohort def ordering).
  They will die in the first timestep (abundance = 1 falls below
  ExtinctionThreshold = 1 immediately, or starvation kills them).

  IF C.csv FILE IS COMPLETELY ABSENT:
  The std::ifstream constructor does not throw on missing file. The file
  stream will simply be empty, equivalent to an empty file. Behaviour
  is the same as zero-row file.

  RECOMMENDED APPROACH FOR ZERO COHORTS:
  Use an empty C.csv (no header, zero data rows). The ghost seeding
  introduces one negligible cohort per cell that dies immediately.
  This is the cleanest approach. Alternatively set all
  PROPERTY_Initial.number.of.GridCellCohorts to 0 in cohort_def to
  prevent cohort seeding during spin, then use the resulting empty
  (or trivially-populated) C.csv.

WHAT COLUMNS MUST BE PRESERVED FOR A VALID RESTART:
  All 16 columns. The reader parses exactly 16 columns by column index.
  A file with fewer columns causes stod/stoi calls on empty strings and
  will throw an exception (caught by the try/catch, incrementing
  totalCohorts_init_failed). Missing values result in skipped cohorts,
  not a crash, because of the try/catch wrapper.

----------------------------------------------------------------
4B. StockProperties_<NNNNN>.csv
----------------------------------------------------------------

EXACT HEADER LINE (from WriteModelState.cpp lines 138-141):
  GridcellIndex,FunctionalGroupIndex,TotalBiomass

3 columns. Delimiter: comma.

Confirmed by get_stock_restart_header.R:
  "GridcellIndex,FunctionalGroupIndex,TotalBiomass"

S.csv WRITTEN BY THE R PACKAGE:
  NO HEADER.  3 data columns.
  (write_madingley_data_cohorts_stocks_to_temp_dir_fast.R:
   fwrite(madingley_data$stocks, sRdir, row.names=FALSE, col.names=FALSE))

S.csv READER IN C++ (readSpinUpStateStock):
  Allocates v with 3 columns: std::vector<std::vector<std::string>> v(3,...)
  Does NOT skip a header. S.csv MUST NOT have a header row.

COLUMN MAPPING in SeedStocksApplySpinUp:
  StockData[0]   GridcellIndex    (stoi for binary search)
  StockData[1]   FunctionalGroupIndex  (not used directly; iteration order used)
  StockData[2]   TotalBiomass     (stod -- line 449: s.mTotalBiomass = stod(StockData[2][StockCounter]))

TotalBiomass COLUMN:
  Column index 2 (0-based). Third column.
  Units: grams (wet matter).
  The writer divides by 1000 only for the timeline output
  (WriteStockTimeline line 459: StockBiomass/1000). The state files
  store raw grams.

HOW SeedStocksApplySpinUp USES StockData:
  It searches GCindics (GridcellIndex column) using binary_search.
  For each cell that has an entry in StockData, it iterates over all
  stocks in the cell (ApplyFunctionToAllStocks) and assigns
  s.mTotalBiomass = stod(StockData[2][StockCounter]) in order.
  StockCounter starts at firstIndex (the first row for this gridcell)
  and increments for each stock in the cell.
  The functional group index stored in column 1 is NOT used to match
  stocks -- only row order within the cell's block matters.

COLUMN ORDER MATCHES STOCK ORDER:
  The rows in S.csv for a given GridcellIndex must be in the same
  order as the stocks were seeded (iterating over
  mStockFunctionalGroupDefinitions.mAllFunctinoalGroupsIndex).
  For the two terrestrial stocks (Deciduous = fg 0, Evergreen = fg 1
  assuming 0-indexed default order), the rows must be:
    GridcellIndex, 0, <deciduous_biomass>
    GridcellIndex, 1, <evergreen_biomass>

  The plan.md wants TotalBiomass = 0.0 for both. A valid zeroed S.csv
  for a single cell (index 0) looks like:
    0,0,0.0
    0,1,0.0

  IMPORTANT: The GridcellIndex in S.csv must match the actual grid cell
  indices produced by the run. For a window of -73,-71,45,46 with 1-deg
  resolution the grid is constructed with cell index 0 for the first
  cell (Parameters::CalculateParameters iterates latitude outer, longitude
  inner: for lat in range for lon in range cellIndex++). For the 1-cell
  window, the only cell has index 0.

WHAT COLUMNS MUST BE PRESERVED FOR A VALID RESTART:
  All 3 columns. The StockCounter iterates sequentially so column 2
  (TotalBiomass) is the critical value to modify for the zero-biomass start.
  Column 1 (FunctionalGroupIndex) is stored in StockData but not used
  by the loader -- only the row order within a gridcell block matters.


=============================================================
SECTION 5  OUTPUT FILES
=============================================================

OUTPUT DIRECTORY STRUCTURE

    The C++ creates these subdirectories of output_dir on the SPIN run:
      cohort_properties/
      stock_properties/

    On the RUN run, additional directories are created:
      cohort_properties/
      stock_properties/
      grid_properties/         (only if WriteGridProperties == 1)
      consumption_statistics/  (only if TypeOfRun == "run")
      timelines/               (only if TypeOfRun == "run")

    Source: WriteModelState::CreateOutputSubfolders  (WriteModelState.cpp
    lines 394-435, platform-specific blocks).

FILENAME PATTERNS

    cohort_properties/FullCohortProperties_NNNNN.csv
        NNNNN is zero-padded to 5 digits using a "00000" template.
        Written at each timestep >= TimestepWritingFullCohortProperties.
        Always written at end of run as FullCohortProperties_99999.csv.
        Also written at end of spin as FullCohortProperties_99999.csv.

    cohort_properties/BinnedCohortStatistics_NNNNN.csv
        Written at each timestep >= TimestepWritingBinnedCohortStatistics.

    stock_properties/StockProperties_NNNNN.csv
        Written at each timestep >= TimestepWritingStockProperties.
        Always written at end of run/spin as StockProperties_99999.csv.

    grid_properties/GridProperties_NNNNN.csv
        Written for months 00000 through 00011 (0-11) once at start of
        run if WriteGridProperties == 1.

    consumption_statistics/PreyBinnedFoodwebConnections_NNNNN.csv
        Written at each timestep >= TimestepWritingPreyBinnedFoodwebConnections.

    consumption_statistics/CohortConsumptionSummary_NNNNN.csv
        Written at each timestep >= TimestepWritingPreyBinnedFoodwebConnections.

    timelines/MontlyCohortBiomass.csv  (NOTE: TYPO -- "Montly" not "Monthly")
        OVERWRITTEN every timestep (open with ios::out, not ios::app).
        Contains cumulative data from month 0 to current month.
        Columns: Month, Year, Biomass_FG_0, Biomass_FG_1, ...
        One row per month elapsed, with Biomass_FG_N being total kg
        biomass of cohorts in functional group N summed over all cells.
        Month is 1-indexed (ii+1). Year is 1-indexed (increments every 12).
        NOTE: the R reader reads this as MontlyCohortBiomass.csv --
        same typo. Use that exact filename.

    timelines/MonthlyStockBiomass.csv
        APPENDED every timestep (open with ios::app).
        One value per line: total stock biomass in kg summed over ALL
        stock functional groups and all grid cells.
        NO header. The R reader reads as a headerless single-column CSV.

LEAF BIOMASS FOR EVERGREEN AND DECIDUOUS STOCKS

    The timelines/MonthlyStockBiomass.csv contains ONLY the sum of all
    stock biomass. It does NOT separate evergreen from deciduous.

    To get per-stock per-cell time series, read StockProperties_NNNNN.csv
    files. Each file has columns GridcellIndex, FunctionalGroupIndex,
    TotalBiomass. Row with FunctionalGroupIndex matching the deciduous
    or evergreen stock FG index gives that stock's biomass at that timestep.
    To get these every month, set TimestepWritingStockProperties = 0
    when calling the run command.

    The stock functional group index is determined by row order in
    StockFunctionalGroupDefinitions.csv. With the default 2-row stock
    def (Deciduous first, Evergreen second) the indices are:
      FunctionalGroupIndex 0 = Deciduous
      FunctionalGroupIndex 1 = Evergreen
    (Indices are assigned in the order rows appear in the definitions file.)

TIMESTEP INDEXING

    The main run loop uses: timeStep = 0 ... (LengthOfSimulationInMonths - 1)
    In the output files, the month index N in the filename equals the
    zero-based timestep. Month 0 is the first month of simulation.
    The timeline file uses 1-based month (ii+1) for display.
    LengthOfSimulationInMonths = NumberOfSimulationYears * 12
    (Parameters::CalculateParameters line 86).
    For a 3-year run: months 0 through 35, files _00000 through _00035.
    Final state always written as _99999.

GRID PROPERTIES FILE

    Written once at start of run if WriteGridProperties == 1.
    Contains per-cell environmental values for each of the 12 calendar
    months (GridProperties_00000.csv = month 0 = January, etc.).
    Includes TerrestrialNPP column which is the monthly NPP value from
    the input raster. This file is useful for sanity-checking that the
    spatial inputs were read correctly.


=============================================================
SECTION 6  THE model_params ARGUMENT (argv[23])
=============================================================

argv[23] = UseNonDefaultModelParameters (int, stored in that variable)

VALUE = 0 (plan.md sets this):
    The model uses all hardcoded default parameter values defined in:
      TerrestrialCarbon.cpp (InitialisePlantModelParameters)
      Activity.cpp, DispersalDiffusive.cpp, EatingCarnivory.cpp, etc.
    No additional argv values are consumed.
    The R code sets model_params = 0 when model_parameters is not a
    data.frame (the default condition).

VALUE = 1:
    The model reads 78 additional parameter values from argv[24] onward.
    The R code constructs this as:
      model_params = paste(c(1, model_parameters$values), collapse=" ")
    where model_parameters is a data.frame from get_default_model_parameters().

    When UseNonDefaultModelParameters = 1, argv is consumed starting
    at var_counter = 24 in this exact order (Main.cpp lines 191-291):

      argv[24..27]   Activity_Parameters (4 values)
      argv[28..33]   Dispersal_Parameters (6 values)
      argv[34..42]   EatingCarnivory_Parameters (9 values)
      argv[43..51]   EatingHerbivory_Parameters (9 values)
      argv[52..58]   MetabolismEctotherm_Parameters (7 values)
      argv[59..64]   MetabolismEndotherm_Parameters (6 values)
      argv[65..67]   MetabolismHeterotroph_Parameters (3 values)
      argv[68..72]   Mortality_Parameters (5 values)
      argv[73..76]   Reproduction_Parameters (4 values)
      argv[77..108]  VegetationModel_Parameters (32 values)

      Total additional args: 4+6+9+9+7+6+3+5+4+32 = 85 values.
      With argv[23]=1 prepended by R, the total string passed as
      model_params is "1 <v0> <v1> ... <v84>".
      So argv[23] = "1" and argv[24..108] = the 85 parameter values.

    NOTE: The C++ code declares VegetationModel_Parameters with size 32
    (Main.cpp line 65). TerrestrialCarbon.cpp defines 32 named parameters
    (indices [0] through [31]). The 32nd is mMassCarbonPerMassLeafDryMatter.
    The R get_default_model_parameters() lists 32 VegetationModel values.

THERE IS NO SEPARATE PARAMS CSV FILE:
    When UseNonDefaultModelParameters = 1, the parameters are passed
    on the command line as separate argv tokens, not as a file path.
    The plan.md comment "path to a custom params CSV" is INCORRECT.
    The arg is always an integer (0 or 1), and when it is 1 the 85
    parameter values follow as additional positional argv tokens.

DEFAULT VEGETATION PARAMETERS (from TerrestrialCarbon.cpp):
  [0]  mMaxNPP                               = 0.961644704
  [1]  mT1NPP                                = 0.237468183
  [2]  mT2NPP                                = 0.100597089
  [3]  mPNPP                                 = 0.001184101
  [4]  mFracStructScalar                     = 7.154615419
  [5]  mAFracEvergreen                       = 1.270782192
  [6]  mBFracEvergreen                       = -1.828591558
  [7]  mCFracEvergreen                       = 0.844864063
  [8]  mMEGLeafMortality                     = 0.040273936
  [9]  mCEGLeafMortality                     = 1.013070062
  [10] mMDLeafMortality                      = 0.020575964
  [11] mCDLeafMortality                      = -1.195235464
  [12] mMFRootMort                           = 0.04309283
  [13] mCFRootMort                           = -1.478393163
  [14] mP2StMort                             = 0.139462774
  [15] mP1StMort                             = -4.395910091
  [16] mMaxFracStruct                        = 0.362742634
  [17] mLFSHalfSaturationFire               = 0.388125108
  [18] mLFSScalarFire                        = 19.98393943
  [19] mNPPHalfSaturationFire               = 1.148698636
  [20] mNPPScalarFire                        = 8.419032427
  [21] mErMin (min evergreen leaf mort)      = 0.01
  [22] mErMax (max evergreen leaf mort)      = 24.0
  [23] mDrMin (min deciduous leaf mort)      = 0.01
  [24] mDrMax (max deciduous leaf mort)      = 24.0
  [25] mFrmMin (min fine root mort)          = 0.01
  [26] mFrmMax (max fine root mort)          = 12.0
  [27] mStmMax (max structural mort)         = 1
  [28] mStmMin (min structural mort)         = 0.001
  [29] mBaseScalarFire                       = 2.0
  [30] mMinReturnInterval                    = 0.00000226032940698105
  [31] mMassCarbonPerMassLeafDryMatter       = 0.476
  [32] mMassLeafDryMatterPerMassLeafWetMatter = 0.213  (hardcoded, NOT
       in VegetationModel_Parameters array even when non-default params used;
       see TerrestrialCarbon.cpp line 89)


=============================================================
SECTION 7  ADDITIONAL TECHNICAL NOTES
=============================================================

BUG IN UpdateLeafStock (TerrestrialCarbon.cpp line 207):
    CalculateFracEvergreen is called with GetGridcellFractionYearFire()
    instead of GetGridcellFractionYearFrost(). This is a copy-paste
    error (the comment says "Frost" but the call is "Fire"). The spin-up
    function CalculateEquilibriumLeafMass correctly uses
    GetGridcellFractionYearFrost(). This means the dynamic leaf stock
    update uses Fire fraction as the evergreen fraction driver instead
    of Frost fraction. This is a confirmed bug in the existing source.
    It affects how fever is computed each timestep but not during init.

STOCK LOADING DURING RUN:
    SeedStocksApplySpinUp uses binary_search on GCindics. For this to
    work correctly, the GridcellIndex column in S.csv must be SORTED
    ASCENDING. The R function writes stocks in the order returned by
    the model (which follows grid cell traversal order -- ascending
    cell index). Ensure S.csv rows are sorted by GridcellIndex.

GHOST COHORT SEEDING:
    Even if C.csv is empty, SeedCohortsApplySpinUpFast seeds one ghost
    cohort (functionalGroup=1, all masses=1g, abundance=1) in every
    terrestrial cell not represented in C.csv. These ghosts will be
    killed in timestep 1 by the extinction check. They do not affect
    vegetation dynamics but they do cause the model to initialise the
    cohort data structures correctly. The plan.md notes this as Q4 --
    the answer is: empty C.csv is accepted but ghosts are seeded.

PARALLEL EXECUTION NOTE:
    With RunInParallel=1, OpenMP is used for grid cell ecology. The
    stock ecology (RunWithinCellStockEcology) is called inside the
    parallel loop. For a single-cell run the parallelism is irrelevant
    but does not cause errors.

SIMULATION LENGTH FROM argv[2]:
    Parameters::SetLengthOfSimulationInMonths receives
    InputParameters::Get()->GetNumberOfSimulationYears() which is
    NumberOfSimulationYears = stoi(argv[2]).
    The function name is misleading -- it actually sets
    mLengthOfSimulationInYears = the passed value.
    Then CalculateParameters multiplies by 12 to get months.
    So argv[2] = 3 --> 36 months of simulation. This is correct.

COORDINATE PARSING:
    Xmin, Xmax, Ymin, Ymax are parsed as int (stoi). Fractional degree
    windows cannot be specified via command line. Only integer degree
    boundaries.

OUTPUT DIRECTORY MUST EXIST:
    The binary does not create output_dir. It only creates subdirectories
    within it. Create output_dir before calling the binary.
    The R package does: dir.create(out_dir, showWarnings = F)

INPUT DIRECTORY LOCATION:
    argv[13] (ConfigurationDirectory) must point to a directory that
    contains the four control CSVs. The FunctionalGroupDefinitions
    constructor prepends ConfigurationDirectory to the filename:
      fileName = ConfigurationDirectory + fileName
    So if ConfigurationDirectory = "/path/to/input/" and the file is
    CohortFunctionalGroupDefinitions.csv, the binary opens
    "/path/to/input/CohortFunctionalGroupDefinitions.csv".
    The trailing slash on ConfigurationDirectory is required.

C.csv AND S.csv FORMAT SUMMARY FOR THE PYTHON NOTEBOOK:

    C.csv for zero cohorts (minimum valid file):
      No header line.
      Zero data rows.
      (One ghost cohort will be seeded per cell but will die immediately.)

    S.csv for zero stock biomass (minimum valid file for 1 cell, 2 stocks):
      No header line.
      Row 1: 0,0,0.0
      Row 2: 0,1,0.0
      (GridcellIndex=0, FunctionalGroupIndex=0 (Deciduous), TotalBiomass=0)
      (GridcellIndex=0, FunctionalGroupIndex=1 (Evergreen), TotalBiomass=0)

    For multi-cell windows, repeat the pattern for each cell index.
    Cell indices are assigned row-major (latitude outer, longitude inner)
    by Parameters::CalculateParameters.

=============================================================
END OF REPORT
=============================================================
