RUN 1 SYNTHESIS -- VEGETATION-ONLY MADINGLEY AT SHERBROOKE
Experiment: sunday_vegetation_experiments_design
Date: 2026-05-25
Author: Claude Code session (Opus 4.7, 1M context)


WHAT WAS DONE

    Drove the precompiled Madingley C++ binary (mac_exec/madingley,
    Mach-O x86_64 v2.02) directly from Python via subprocess. No R was
    invoked. Produced a Jupyter notebook (vegetation_sherbrooke.ipynb)
    backed by a Python module (driver.py) that:

      - Writes the four control CSVs (cohort_def, stock_def,
        SimulationControlParameters, MassBinDefinitions) by hard-coding
        the R-package default tables extracted from the R source.
      - Converts the 68 input rasters from inst/spatial_input_rasters
        to the var,x,y CSV format the C++ binary expects, cropped to a
        1-degree window over Sherbrooke (-72.5 lon, 45.5 lat).
      - Calls the spin (init) step.
      - Overwrites the cohort and stock restart CSVs to enforce
        zero starting biomass for both autotroph stocks and empty
        heterotroph cohorts.
      - Runs a 3-year (36-month) simulation with monthly stock output.
      - Parses StockProperties_NNNNN.csv files and plots deciduous,
        evergreen, and total leaf biomass.

    Wall time: spin 0.016 sec, run 0.038 sec. Trivially fast for
    a 1-cell vegetation-only run.


KEY FINDINGS

    1. THE RUN WORKED CLEANLY. From a zero-biomass start the model
       grows leaf biomass over month 1-8 to a year-1 peak of
       ~2.0e13 g, drops sharply over autumn-winter (deciduous loses
       ~90%, evergreen loses ~12%), then repeats with higher amplitude
       in year 2 (peak 3.07e13 g) and year 3 (peak 3.71e13 g),
       approaching a quasi-equilibrium.

    2. SIMULATED EVERGREEN FRACTION DOES NOT MATCH THE PUBLISHED
       EQUATION. At the year-3 growing-season peak the simulated
       evergreen share is 0.848. The synthesis equation
       f_ever = a*F_frost^2 + b*F_frost + c with the published
       coefficients (a=1.27, b=-1.83, c=0.85) and the model-computed
       FractionYearFrost = 0.589 at Sherbrooke predicts 0.21.

    3. THE DISCREPANCY EXACTLY MATCHES THE KNOWN SOURCE BUG.
       TerrestrialCarbon.cpp line 207 calls
         CalculateFracEvergreen( gcl.GetGridcellFractionYearFire() )
       in UpdateLeafStock, where the spin-up path (line 111) and the
       supplementary mathematics both require FractionYearFrost.
       At Sherbrooke FractionYearFire = 0, so the buggy expression
       evaluates to c_fever = 0.845 -- nearly identical to the
       simulated 0.848. Confirmed via the model's own diagnostic
       output (GridProperties_00000.csv).

    4. SIMULATIONCONTROLPARAMETERS.CSV IS IGNORED AT RUNTIME.
       The C++ FileReader::ReadInputParameters has the ReadTextFile
       call commented out. LengthOfSimulationInYears, RunInParallel,
       etc. are read from command-line arguments instead. The CSV
       must still be present (header check) but its values are inert.

    5. C.CSV AND S.CSV MUST BE HEADERLESS. The C++ readers parse
       from byte 0 without skipping a header line. The R helper
       get_cohort_restart_header() that returns 12 columns is
       misleading -- the actual restart writer produces 16 columns,
       and the writer/reader pair uses no header.

    6. EMPTY C.CSV IS ACCEPTED. The reader seeds one "ghost" cohort
       (FG 1, mass 1g, abundance 1) per terrestrial cell, which dies
       in timestep 1 below the extinction threshold. The run log
       confirms n_cohorts = 0 from month 1 onward, so vegetation
       dynamics are isolated.

    7. SPATIAL WINDOW SUBTLETY. The C++ binary uses
         x > Xmin AND x < Xmax + 1
       (NewXmax = Xmax + 1) but the R workflow CSV-crops with
         x > Xmin AND x < Xmax
       For the Sherbrooke window c(-73, -71, 45, 46) both bounds
       admit two cells (-72.5 and -71.5 at lat 45.5), contradicting
       the plan's "exactly one cell". Used xmax = -72 in driver.py
       to deliver exactly one cell as intended.

    8. PER-STOCK MONTHLY BIOMASS REQUIRES StockProperties_NNNNN.csv.
       The timelines/MonthlyStockBiomass.csv file is a single column
       summed over all stocks and cells. To separate deciduous from
       evergreen we set TimestepWritingStockProperties = 0
       (argv[11] = 0) so the binary writes a state file every month.


CAVEATS

    Q1 RESOLVED. The plan asked whether deciduous biomass would drop
       sharply in winter given Sherbrooke's frost regime. Answer: yes,
       the deciduous stock drops ~90% over October-January each year.
       However, due to finding 3 above, the evergreen-deciduous
       *split* is unreliable in this version of the source code. The
       *seasonal cycle* of each pool is still physically reasonable.

    Q2 OPEN. 1-degree spatial resolution remains too coarse for
       city-scale studies. Higher-resolution forcing rasters would
       require regenerating all 60 monthly + 8 static inputs.

    Q4 RESOLVED. Empty headerless C.csv is the cleanest realisation
       of "no heterotrophs". The ghost-cohort kludge does not affect
       vegetation.

    Q5 NOT NEEDED. Empty C.csv already suffices; no need to also
       zero out cohort_def Initial-number-of-cohorts.


SUGGESTED NEXT STEPS

    1. PATCH THE FROST/FIRE BUG locally and rerun. Replace
       GetGridcellFractionYearFire() with
       GetGridcellFractionYearFrost() at TerrestrialCarbon.cpp:207
       (and likely :220), recompile, and confirm that f_ever_sim
       drops from ~0.85 to ~0.21 at Sherbrooke.

    2. EXTEND TO LONGER RUNS (10-50 years) at Sherbrooke to test
       equilibration time and inter-annual variability.

    3. LATITUDINAL TRANSECT (say 30 N to 60 N at fixed longitude) to
       map the frost-driven f_ever gradient predicted by the synthesis
       equation. This needs the patched binary to be meaningful.

    4. COMPARE PEAK NPP to MODIS MOD17 annual NPP at this cell.

    5. RE-ENABLE HETEROTROPHS (do not zero C.csv) and check the
       terrestrial herbivore:autotroph biomass ratio against the 0.93%
       benchmark from the synthesis document.


FILES PRODUCED

    driver.py                                  Python module: CSV writers,
                                               subprocess wrappers, output
                                               parsers.

    vegetation_sherbrooke.ipynb                The narrative notebook with
                                               climate forcings, biomass
                                               trajectory, f_ever
                                               comparison, and caveats.

    run_outputs/                               Full set of input CSVs,
                                               spin and run logs, monthly
                                               stock/cohort/grid CSVs.

    notes/cpp_interface_report.md              Source-derived report on
                                               the C++ executable
                                               interface and CSV formats.

    notes/r_defaults_report.md                 Source-derived report on
                                               the R default tables and
                                               the R-package shell call
                                               construction.

    notes/run1_sherbrooke_synthesis.md         This document.


REFERENCES

    Harfoot et al. (2014) PLoS Biology 12(4): e1001841
    Hoeks et al. (2021) Global Ecology and Biogeography 30: 1922-1933
    notes/madingley_session_synthesis.md (prior session)
