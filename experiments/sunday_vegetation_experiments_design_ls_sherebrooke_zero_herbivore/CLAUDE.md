# CLAUDE.md — Vegetation Experiment Design
# Madingley GEM Working Group
# Experiment: sunday_vegetation_experiments_design

---

## FORMAT AND TONE

Use plain ASCII text in all responses. No markdown formatting symbols
(no #, *, **, __, >, -, backtick code fences, etc.). Use indentation
and blank lines for structure. Use ALL CAPS for section headings.
Use plain hyphens for list items where needed.

This is a scientific working group with strong domain knowledge in
ecology and biodiversity modelling. Be direct and technically precise.
Do not over-explain basic concepts. Use correct scientific notation
and terminology. When presenting equations, use plain-text notation
as used in the synthesis document and in the Harfoot et al. 2014
supplementary material.

---

## SCIENTIFIC CONTEXT

This experiment is part of a working group studying the Madingley
General Ecosystem Model (GEM), a mechanistic model of global
ecosystem structure and function. The primary scientific papers are:

  Harfoot et al. (2014). Emergent Global Patterns of Ecosystem
  Structure and Function from a Mechanistic General Ecosystem Model.
  PLoS Biology 12(4): e1001841. DOI: 10.1371/journal.pbio.1001841

  Hoeks et al. (2021). MadingleyR: An R package for mechanistic
  ecosystem modelling. Global Ecology and Biogeography 30: 1922-1933.
  DOI: 10.1111/geb.13354

The current experiment focuses on the VEGETATION (AUTOTROPH) component
of the Madingley model — specifically the terrestrial leaf biomass
dynamics, NPP calculations, and the equations governing evergreen vs.
deciduous partitioning.

---

## PRIOR SESSION WORK

A prior Claude session (see notes/madingley_session_synthesis.md)
covered the following:

  - Full read-through of the Harfoot et al. (2014) main paper and
    supplementary mathematical details (Text S1)
  - Complete breakdown of all autotroph ecology equations and
    variable definitions (see Section 3 and 3.9 of synthesis)
  - Step-by-step recommended workflow for implementing terrestrial
    and marine autotroph ecology (see Section 4 of synthesis)
  - Full MadingleyR R package workflow documentation (see Section 5)
  - Attempted MadingleyR installation (failed due to network
    restrictions in the prior compute environment)

All equation derivations and variable tables from that session are
preserved in the synthesis document. Do not repeat or re-derive
them unless asked. Reference the synthesis document when relevant.

---

## SETUP FOR NEW CONTRIBUTORS

  After cloning the repository, run the init script from the repo root:

    ./init.sh

  This clones MadingleyR from GitHub into the context folder. If the
  folder already exists, it pulls the latest changes instead. The
  context/MadingleyR/ directory is git-ignored and must be populated
  this way -- it is not committed to the repository.

  Then install the R package from the local clone (run inside R):

    remotes::install_local(
      'experiments/sunday_vegetation_experiments_design/context/MadingleyR/Package'
    )

  Or install directly from GitHub (requires internet):

    remotes::install_github('MadingleyR/MadingleyR',
                            subdir = 'Package',
                            build_vignettes = TRUE)

  Confirm installation:

    library('MadingleyR')
    madingley_version()
    vignette('MadingleyR')

---

## EXPERIMENT STRUCTURE AND CONVENTIONS

  Each experiment lives in experiments/<experiment_name>/ and should
  contain the following:

    CLAUDE.md          Instructions for Claude Code in this experiment.
                       Include: format/tone, scientific context, prior
                       session work, file map, open questions, next steps.
                       Keep it updated as the experiment progresses.

    context/           External packages, datasets, and reference
                       materials. Large or cloned external dependencies
                       go here and must be added to the root .gitignore.
                       Populate them with init.sh entries.

    notes/             Session synthesis documents and working notes.
                       After each substantive Claude Code session, save
                       a synthesis document here so context is not lost
                       across sessions. Use plain ASCII.

  When starting a new experiment:
    1. Copy this CLAUDE.md as a template and update all sections.
    2. Add any new external dependencies to init.sh and .gitignore.
    3. Write a brief synthesis note after the first session.

---

## FILES IN THIS EXPERIMENT

  context/MadingleyR/              Full GitHub clone of MadingleyR
                                   (MadingleyR/MadingleyR, master branch)
                                   NOT in git -- populated by init.sh.
    Package/                       R package source (DESCRIPTION,
                                   NAMESPACE, R/, man/, vignettes/)
    SourceCode/                    C++ Madingley model source code
    CaseStudies/                   Three worked case studies from
                                   Hoeks et al. (2021)
    Documentation/                 PDFs: R functions, spatial input
                                   units, model parameters
    Tests/                         Workflow test scripts v1-v4

  notes/madingley_session_synthesis.md
                                   Full synthesis of the prior session
                                   including all equations, variables,
                                   workflow steps, and open questions

---

## KEY AUTOTROPH EQUATIONS (QUICK REFERENCE)

Terrestrial leaf biomass update:
  Bl,(t+dt) = Bl,t + dBl,Growth - dBl,Mort

Growth contributions:
  dBl,Growth = Gever + Gdecid
  Gever  = NPPmT * Acell * psi * dtNPP * (1-fStruct) * fLeafMort * fever
  Gdecid = NPPmT * Acell * psi * dtNPP * (1-fStruct) * fLeafMort * (1-fever)

Leaf mortality:
  dBl,Mort = dtmu,l * (mu_ever * Bl,t * fever
             + mu_decid * Bl,t * (1-fever)) + Ll,Herbivory

Annual terrestrial NPP (Miami model):
  NPPyT = min(NPPT, NPPP)
  NPPT  = NPPmax / (1 + exp(cP - mP * T(t)C))
  NPPP  = NPPmax * (1 - exp(-rho * P))
  NPPmT = NPPyT * omega_cell,m

Evergreen fraction from frost frequency:
  fever = a_fever * Ffrost^2 + b_fever * Ffrost + c_fever

Full variable definitions are in the synthesis document Section 3.9.

---

## OPEN QUESTIONS (FROM SYNTHESIS, SECTION 7)

These questions have not yet been answered by the working group.
Prompt for them if they are needed to proceed:

  1. Target operating system (Windows / Mac / Linux)
  2. R version installed on the target machine
  3. Primary scientific goal:
       - testing/validating existing model code
       - implementing autotroph equations from scratch
       - running MadingleyR scenarios
       - modifying parameters or functional forms
       - comparing output to empirical data
  4. Spatial domain / region of interest
  5. Terrestrial only, marine only, or both
     (note: MadingleyR R package covers terrestrial only)
  6. Specific biodiversity or climate change scenarios
  7. Existing code to review or debug
  8. Empirical datasets available for validation
  9. Computational resources (local laptop vs. HPC)
  10. Priority output metrics (biomass, trophic structure,
      species distributions, etc.)

---

## SUGGESTED NEXT STEPS

  1. Install MadingleyR from local clone:
       library('remotes')
       install_local('context/MadingleyR/Package')

  2. Run the tutorial vignette to confirm installation:
       vignette('MadingleyR')

  3. Run a basic spin-up for a focal grid cell

  4. Implement and unit-test autotroph equations as standalone
     R or Python functions

  5. Compare model autotroph output to MODIS NPP data

---

## PACKAGE VERSION NOTE

MadingleyR v1.0.6 (source code v2.02, June 2024).
Uses terra, sf, data.table. Does not use deprecated raster/rgdal.
Marine realm is NOT included in the current R package version.
Fundamental process equations require C++ edits to modify.

---

## VALIDATION BENCHMARKS

  Terrestrial herbivore:autotroph ratio ~ 0.93%
  Marine herbivore:autotroph ratio ~ 52%

---
