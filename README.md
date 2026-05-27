# GEM Working Group

This repository supports the GEM working group focused on designing and prototyping a new General Ecosystem Model.

The project is inspired by Madingley, but the goal is not to reproduce it directly. The working group is exploring a new model structure that can better represent biodiversity change, species-level dynamics, trophic interactions, spatial structure, dispersal, vegetation, coexistence, and ecosystem processes.

A second goal is to document how AI coding agents can support collaborative ecological modelling, while keeping scientific assumptions, model structure, tests, and validation explicit.

## Objectives

- Design a new process-based ecosystem model.
- Compare the effects of biodiversity change and climate change on ecosystem functioning.
- Explore species-level or hybrid species/guild representations.
- Develop modular model components for trophic dynamics, vegetation, dispersal, mortality, reproduction, and spatial processes.
- Define ecological tests, diagnostics, and validation targets.
- Document lessons from AI-assisted collaborative model development.

## Resources and documentation

[📂 SharePoint Documents](https://usherbrooke.sharepoint.com/sites/ielabworkinggroup)

[💬 Teams Chat](https://teams.microsoft.com/l/team/19%3A5UwBzgYI52ESK9znTZVksNqlEEprrzf8AeMWIcRYWpU1%40thread.tacv2/conversations?groupId=f8341ee7-d260-4f7a-9890-6f13bbc5ec80&tenantId=3a5a8744-5935-45f9-9423-b32c3a5de082)

[📅 Logistics](https://usherbrooke.sharepoint.com/:f:/r/sites/ielabworkinggroup/Documents%20partages/logistics?csf=1&web=1&e=nsdYNf)

## General Collaboration guidelines

- This repo is where all notes, papers, and documents for the duration of the row group will live
- The `main` branch is protected, so please create a new branch for any changes and submit a pull request.

## Project structure



## Python packaging and environment

<!-- Dependancies -->

<!-- .toml -->

<!-- .venv and and getting started with venv and pip install -e . -->


## Style guide and naming conventions

<!-- Style guide, casing for modules, functions and naming -->


## Experiments and prototyping

The `experiments` folder is where we will develop and share code for model prototyping, testing, and experiments. Each experiment should have its own subfolder with a README describing the purpose, methods, and results.

Naming convention for experiment folders: `DAY_GROUPNAME_experimentNAME`. Example - `sunday_atn_bylot_experiment1`.

We recommend using Jupyter notebooks for prototyping and documentation, but feel free to use other formats as needed. The key is to keep everything organized and well-documented for future reference.

## Model development

<!-- Implementation details for processes into module from processes contract -->

- <!-- Big idea -->
- Numpy structures ...
- ...

## Input data files


## Geographic grid


## Simulation engine

<!-- Describing Alex's engine state management species registry. There should be a doc describing that  -->

- State management
- <!-- Processes adapters  -->
- ...Broadcasting
- Initialization


## Running simulations

Inside experiments folder. Store relevant data.

---
<!-- 
- Dependancies modules (ex. metabolism) that do not return biomass but are reused by multiple processes. Should be specified in processes contract. Processes should return biomass delta, or biomass delta and other outputs (e.g. fluxes, rates, etc.) that are relevant to the process and can be used by other processes or for analysis. Their outputs should be stored in broadcasting-friendly data structures that can be easily accessed by other processes and engine.
- How to handle and store simulation runs (notebooks ?,  saved outputs ?) and how to make them accessible to the team. Make minimal requirements for reproducibility of runs (e.g. saving the random seed, saving the configuration file, etc.). Naming convention for runs and outputs with date and time and description.
- Input data. Script to download and preprocess input data (e.g. environmental data, species traits, etc.) and store them in a standardized format that can be easily accessed by the engine and processes. We recommend using geotiff to store spatial data and csv or json for tabular data. We also recommend using a standardized directory structure for input data (e.g. data/raw, data/processed, etc.) and a naming convention for files (e.g. data/raw/environmental_data_2024-06-01.tif). Local gitignores for large .tiff datasets that are not stored in the repository but can be downloaded and processed by the script.
- Initialization scripts for the engine (e.g. to set up the grid, load initial conditions, species list and traits, load the input data). To be described in engine section. Should be modular and reusable for different runs and configurations. Should also include error handling and logging to facilitate debugging and tracking of runs.

 -->