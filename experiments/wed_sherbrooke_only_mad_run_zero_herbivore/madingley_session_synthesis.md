# Madingley Model Working Group: Session Synthesis

## Purpose of This Document

This document synthesizes a working session conducted with a scientific
working group focused on testing code for models of biodiversity and
climate change. It is intended to allow a Claude Code session (or any
follow-up session) to continue the work without loss of context. All
responses in the original session were in plain ASCII format.

---

## 1. Papers and Documents Read

### Document 1: Supplementary Material (Text S1)
- File: `pbio_1001841_s018.docx`
- Content: Full technical and mathematical details of the Madingley
  General Ecosystem Model (GEM)
- Authors: Harfoot et al.
- Affiliation: UNEP-WCMC and Microsoft Research Computational Science
  Laboratory

### Document 2: Main Paper
- File: `Harfoot_et_al_2014.pdf`
- Full title: "Emergent Global Patterns of Ecosystem Structure and
  Function from a Mechanistic General Ecosystem Model"
- Journal: PLoS Biology 12(4): e1001841 (2014)
- DOI: 10.1371/journal.pbio.1001841

### Document 3: R Package Paper
- File: `Global_Ecology_and_Biogeography_-_2021_-_Hoeks_-_MadingleyR__An_R_package_for_mechanistic_ecosystem_modelling.pdf`
- Full title: "MadingleyR: An R package for mechanistic ecosystem
  modelling"
- Journal: Global Ecology and Biogeography 30: 1922-1933 (2021)
- DOI: 10.1111/geb.13354
- Authors: Hoeks S, Tucker MA, Huijbregts MAJ, Harfoot MBJ,
  Bithell M, Santini L

---

## 2. The Madingley Model: Summary

The Madingley Model is the first mechanistic General Ecosystem Model
(GEM) that is both global in scope and applies to all terrestrial and
marine environments. It is also known as the GEM.

### Key Characteristics
- Simulates organisms with body masses from 10 mg to 150,000 kg
  (14 orders of magnitude)
- Applies to both terrestrial and marine realms (freshwater excluded)
- Spatially explicit: 2D grid, typically 1x1 or 2x2 degree cells,
  65N to 65S latitude
- Time step: monthly
- Derives ecosystem-level properties from individual-level biology
  without imposing top-down constraints
- Originally written in C#, later translated to C++
- Now available as an R package: MadingleyR (terrestrial realm only)

### Model Components
- Autotrophs represented as stocks (total biomasses)
- Heterotrophs represented as cohorts (groups of organisms with
  identical functional traits in the same grid cell)
- Each cohort tracks: individual body mass, abundance, reproductive
  potential mass

### Ecological Processes Modelled
- Autotrophs: growth (NPP-driven), mortality (herbivory)
- Heterotrophs: eating, metabolism, reproduction, mortality
  (background/starvation/senescence), dispersal

### Environmental Drivers (inputs per grid cell)
- Terrestrial: air temperature, precipitation, soil water
  availability, frost days, NPP seasonality (MODIS)
- Marine: sea-surface temperature, NPP (VGPM), ocean current
  velocity

---

## 3. Autotroph Ecology: Equations and Variables

This was the primary mathematical focus of the session. Full
equation-by-variable breakdowns were produced for all autotroph
ecology equations.

### 3.1 Biomass Dynamics

**Terrestrial leaf biomass:**
```
Bl,(t+dt) = Bl,t + dBl,Growth - dBl,Mort
```

**Marine phytoplankton biomass:**
```
Bp,(t+dt) = Bp,(t) + dBp,Growth - dBp,Mort
```

### 3.2 Phytoplankton Growth and Mortality
```
dBp,Growth = NPPmM * xi * Acell * dtNPP
dBp,Mort   = Lp,Herbivory
```

### 3.3 Terrestrial Leaf Growth
```
dBl,Growth = Gever + Gdecid

Gever  = NPPmT * Acell * psi * dtNPP * (1-fStruct) * fLeafMort * fever
Gdecid = NPPmT * Acell * psi * dtNPP * (1-fStruct) * fLeafMort * (1-fever)
```

### 3.4 Terrestrial Leaf Mortality
```
dBl,Mort = dtmu,l * (mu_ever * Bl,t * fever
           + mu_decid * Bl,t * (1-fever)) + Ll,Herbivory
```

### 3.5 Structural Allocation
```
fStruct = min( fStructMin * exp(phi_fstruct * NPPTerr)
               / (1 + fStructMin * (exp(phi_fstruct * NPPTerr) - 1)),
               0.99 * fStructMax )
```

### 3.6 Leaf and Root Mortality
```
fLeafMort   = mu_Leaf / (mu_Leaf + mu_FineRoot)
mu_Leaf     = exp( fever * ln(mu_ever) + (1-fever) * ln(mu_decid) )
mu_ever     = exp( me * T(t)C - ce )
mu_decid    = exp( -(md * T(t)C + cd) )
mu_FineRoot = exp( mf * T(t)C + cf )
```

### 3.7 Evergreen Fraction
```
fever = a_fever * Ffrost^2 + b_fever * Ffrost + c_fever
```

### 3.8 Annual Terrestrial NPP (Miami Model)
```
NPPyT = min(NPPT, NPPP)
NPPT  = NPPmax / (1 + exp(cP - mP * T(t)C))
NPPP  = NPPmax * (1 - exp(-rho * P))
NPPmT = NPPyT * omega_cell,m
```

### 3.9 Key Variable Definitions

| Symbol        | Definition                                               |
|---------------|----------------------------------------------------------|
| Bl,t          | Leaf biomass at time t (g wet biomass)                   |
| Bp,(t)        | Phytoplankton biomass at time t (g wet biomass)          |
| NPPmM         | Monthly marine NPP from VGPM satellite model             |
| NPPmT         | Monthly terrestrial NPP (kg C m-2 month-1)               |
| NPPyT         | Annual terrestrial NPP (Miami model)                     |
| NPPmax        | Maximum possible NPP                                     |
| xi            | Carbon to wet matter conversion factor (marine)          |
| psi           | Carbon to wet matter conversion factor (terrestrial)     |
| Acell         | Grid cell area                                           |
| dtNPP         | Scalar: monthly NPP to model timestep                    |
| dtmu,l        | Scalar: annual leaf mortality to model timestep          |
| omega_cell,m  | Fraction of yearly NPP in month m (from MODIS)           |
| fStruct       | Fraction of NPP allocated to structural tissue           |
| fLeafMort     | Fraction of total mortality that is leaf mortality       |
| mu_Leaf       | Mean leaf mortality rate (per year)                      |
| mu_ever       | Evergreen leaf mortality rate (per year)                 |
| mu_decid      | Deciduous leaf mortality rate (per year)                 |
| mu_FineRoot   | Fine root mortality rate (per year)                      |
| fever         | Proportion of NPP from evergreen leaves                  |
| Ffrost        | Frost frequency at grid cell                             |
| T(t)C         | Monthly average temperature (degrees C)                  |
| P             | Total annual precipitation                               |
| me, ce        | Slope and intercept for evergreen mortality vs temp      |
| md, cd        | Slope and intercept for deciduous mortality vs temp      |
| mf, cf        | Slope and intercept for fine root mortality vs temp      |
| cP, mP        | Coefficients for temperature-NPP logistic function       |
| rho           | Coefficient relating NPP to precipitation                |
| Lp,Herbivory  | Phytoplankton biomass consumed by herbivores/omnivores   |
| Ll,Herbivory  | Leaf biomass consumed by herbivores/omnivores            |

---

## 4. Recommended Steps to Model Autotrophic Ecology

### Terrestrial
1. Obtain monthly terrestrial NPP from MODIS Terra
2. Obtain monthly climate data: temperature, precipitation,
   frost days, diurnal temperature range
3. Calculate annual NPP using Miami model equations, then
   disaggregate to monthly using MODIS seasonal weights
4. Calculate fever from frost frequency (quadratic relationship)
5. Calculate fStruct from annual NPP (sigmoid function)
6. Calculate fLeafMort from temperature-dependent mortality rates
7. Apply Gever and Gdecid growth equations each timestep
8. Apply leaf mortality each timestep using temperature-dependent
   rates scaled to model timestep
9. Subtract herbivory losses after heterotroph cohorts have acted

### Marine
1. Obtain monthly marine NPP from VGPM satellite data
2. Apply conversion factor xi (carbon to wet matter)
3. Scale monthly NPP to model timestep using dtNPP, multiply by Acell
4. Subtract herbivory losses after heterotroph cohorts have acted
5. Note: background phytoplankton mortality is assumed negligible

### General
1. Initialise biomass stocks at simulation start
2. Run autotroph equations first at each timestep, before
   heterotroph cohorts act
3. Track biomass through time; check stocks persist
4. Validate against empirical benchmarks:
   - Terrestrial herbivore:autotroph ratio ~ 0.93%
   - Marine herbivore:autotroph ratio ~ 52%
5. Note: herbivory does not affect plant allocation strategies
   in current model formulation

---

## 5. MadingleyR R Package

### Installation (requires internet access on user machine)
```r
# Step 1: Install remotes
install.packages('remotes')

# Step 2: Install MadingleyR from GitHub
library('remotes')
install_github('MadingleyR/MadingleyR',
               subdir='Package',
               build_vignettes = TRUE)

# Step 3: Load and verify
library('MadingleyR')
madingley_version()
vignette('MadingleyR')
```

Dependencies installed automatically: terra, sf, data.table
Also downloaded automatically: pre-compiled C++ executable,
default spatio-temporal input layers, default model parameters

### Key Functions
```r
madingley_inputs()   # load default inputs
madingley_init()     # initialise model (generates cohorts and stocks)
madingley_run()      # run simulation for specified years
madingley_plot()     # produce output plots
```

### Workflow
```r
# Step 1: Load inputs
sptl_inp  <- madingley_inputs('spatial inputs')
chrt_def  <- madingley_inputs('cohort definition')
stck_def  <- madingley_inputs('stock definition')
mdl_pars  <- madingley_inputs('model parameters')

# Step 2: Modify inputs if needed (optional)

# Step 3: Initialise
m_data <- madingley_init(spatial_window = c(31,35,-5,-1),
                         cohort_def     = chrt_def,
                         stock_def      = stck_def,
                         spatial_inputs = sptl_inp)

# Step 4: Spin-up (100-1000 years recommended)
m_data2 <- madingley_run(m_data, years = 100)

# Step 5: Modify output object if running a scenario
# e.g. remove large herbivores, reduce autotroph production

# Step 6: Run scenario simulation
m_data3 <- madingley_run(m_data2, years = 50)

# Step 7: Plot outputs
madingley_plot(m_data3)
```

### Three Case Studies from Hoeks et al. (2021)
1. Large herbivore removal (>100 kg endothermic herbivores removed
   from Serengeti after 100-year spin-up; 50-year follow-up)
2. Small-scale land-use intensity (autotroph production reduced in
   10% increments; nonlinear response in endotherm biomass)
3. Continental-scale land-use (95% autotroph reduction across Africa;
   herbivores persisted longer in high-productivity areas)

### Limitations
- Terrestrial realm only in current R package version
- Marine realm not yet included
- Fundamental process equations require C++ edits to change
- Computationally intensive for large spatial domains

### Code and Data Availability
- R package:   https://github.com/MadingleyR/MadingleyR
- C++ source:  https://github.com/MadingleyR/MadingleyR/tree/master/SourceCode
- Docs:        https://madingleyr.github.io/MadingleyR/
- Zenodo:      https://doi.org/10.5281/zenodo.4790806
- GitHub release: https://github.com/MadingleyR/MadingleyR/releases/tag/GEB

---

## 6. Installation Attempt in This Session

R version 4.3.3 was successfully installed in the compute
environment. However, CRAN repositories are not accessible
from this environment due to network restrictions. The
remotes package could not be downloaded. Installation must
be completed on the working group's own machines.

Error encountered:
```
Warning: unable to access index for repository
https://cloud.r-project.org/src/contrib
Warning: package 'remotes' is not available for this version of R
```

---

## 7. Context Needed to Continue in Claude Code

To continue this work in a Claude Code session, the following
context and files will be needed:

### Files to Carry Forward
- `pbio_1001841_s018.docx` -- supplementary mathematical details
- `Harfoot_et_al_2014.pdf` -- main Madingley model paper
- `Hoeks_et_al_2021_MadingleyR.pdf` -- R package paper

### Information Needed from the Working Group
1. What operating system will be used for running MadingleyR?
   (Windows, Mac, or Linux -- affects installation steps)
2. Is R already installed on the target machine? What version?
3. What is the primary scientific goal of the working group?
   Options include:
   - Testing/validating existing model code
   - Implementing the autotroph equations from scratch
   - Running MadingleyR simulations for specific scenarios
   - Modifying model parameters or functional forms
   - Comparing model output to empirical data
4. Is the focus terrestrial only, marine only, or both?
   (Note: MadingleyR currently covers terrestrial only)
5. What spatial domain or region is of interest?
6. Are there specific biodiversity or climate change scenarios
   to be tested?
7. Is there existing code (R, C++, Python, or other) already
   written by the group that needs to be reviewed or debugged?
8. What empirical datasets are available for validation?
9. What is the computational resource available?
   (local laptop vs. HPC cluster)
10. What output metrics are most important to the group?
    (biomass, trophic structure, species distributions, etc.)

### Suggested Next Steps in Claude Code
1. Install MadingleyR on the working machine
2. Run the tutorial vignette to confirm installation
3. Run a basic spin-up simulation for a focal grid cell
4. Implement and test the autotroph equations as standalone
   R or Python functions for unit testing
5. Compare model autotroph outputs to MODIS NPP data

---

## 8. Format Note

All responses in the original session were delivered in plain
ASCII format as requested by the working group. This should
be maintained in any follow-up session.

---

*Document generated from working session with scientific
working group on biodiversity and climate change modelling.*
*Date: May 25, 2026*
