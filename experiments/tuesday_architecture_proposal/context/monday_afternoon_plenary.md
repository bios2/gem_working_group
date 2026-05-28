# Meeting notes synthesis — Monday afternoon GEM working session

**Meeting:** Monday afternoon working session  
**Date in transcript:** 2026-05-25  
**Transcript duration:** 1 h 51 min 29 s  
**Main groups discussed:** vegetation, ATN, distribution/dispersal, software architecture.

**Sources used:** Monday afternoon transcript, distribution/dispersal whiteboard, Vincent’s afternoon notes, and prior Harfoot/Madingley context.

**Scope note:** These notes synthesize only what appears in the provided documents. Where the transcript is unclear, I mark the point as uncertain instead of filling gaps.

**Transcript correction note:** the transcript contains several automatic transcription errors. I interpret “Madding / Madding Lane / Manning Lee / admittingly” as **Madingley**. I interpret the transcript’s “ATM” transcription errors as **ATN**, meaning **Allometric Trophic Network**.

---

## 1. Overall purpose of the afternoon session

The afternoon session focused on turning the GEM working group ideas into more concrete modelling and implementation pieces. Three scientific/model components were discussed in turn:

1. vegetation dynamics;
2. distribution/dispersal dynamics;
3. ATN dynamics.

The session also included a software architecture discussion, mainly about how to organize model processes, states, data, and simulations so that the model remains modifiable and testable.

---

## 2. Vegetation model

### 2.1 Functional group decisions

The vegetation discussion defines a slightly more diversified set of plant functional groups than the current Madingley-style vegetation representation.

The functional groups retained for the current formulation are:

- **herbaceous plants / herbs**;
- **evergreen trees**, described as conifer trees;
- **deciduous trees**, described as trees that lose their leaves.

Vincent’s personal notes record the functional groups as:

- deciduous trees: **D**;
- herbaceous: **H**;
- evergreens: **E**.

A fourth group, **shrubs**, was mentioned as a possible later addition. It remains an open question because it was described as potentially more complicated and was not part of the immediate formulation.

### 2.2 Motivation

The motivation is that the current Madingley-style plant model has evergreens and deciduous plants but lacks herbs, which makes it unable to represent grasslands properly.

A target behaviour is that different ecosystem types should emerge across environmental gradients such as precipitation and temperature. The transcript explicitly says the goal is for ecosystem differences in space to emerge from the model rather than being directly imposed by tunable parameters.

The personal notes also mention a **Whittaker-style community function of climate** as an expected dynamic or output to reproduce.

### 2.3 Hierarchical regional/local structure

The proposed vegetation model is hierarchical, mixing regional and local dynamics.

At the **regional** level:

- the region is the grid cell used by the simulation;
- occupancy is tracked as a proportion of the landscape;
- disturbance, especially fire, is introduced at this level.

At the **local** level:

- each grid cell is considered to contain an infinite number of very small patches;
- local patches are small enough that each patch follows only one tree type/species trajectory at a time;
- a patch can be on an evergreen trajectory or a deciduous trajectory;
- herbs are always present, but herbaceous biomass is dynamic and can be high or low depending on the vegetation dynamics.

The transcript clarifies that “local” could be thought of as a very small area such as a square metre or 100 square metres, but the exact size is not fixed. The purpose is to work with proportions rather than stochastic fluctuations.

### 2.4 Fire and succession

Fire is treated as essential for North American vegetation dynamics, especially in a climate-change context.

The proposed succession logic is:

- after fire, the community is dominated by herbaceous plants;
- herbs grow quickly and dominate early succession;
- tree biomass then grows;
- eventually forest dominates and herbs decline;
- herbivores can delay tree takeover;
- if herbivory on trees is very strong, trees may be prevented from taking over and the system may remain herb-dominated.

Vincent’s notes list expected succession directionality as:

- **E → H**;
- **H → E**;
- **H → D**;
- **D → H**.

### 2.5 Occupancy and time since disturbance

A new regional variable is introduced:

\[
P_{i,a,t}
\]

Interpreted in the transcript and notes as the proportion of the landscape occupied by functional group \(i\), with \(a\) representing time since fire/disturbance and \(t\) representing model time.

The notes define:

- \(P_{i,a}\): occupancy, where \(a\) is time since disturbance/fire;
- \(e_i\): probability of observing fire for functional group \(i\);
- \(G_{it}\): growth;
- \(M_{it}\): mortality.

The transcript takes priority for the boundary condition: the written boundary-condition equation was explicitly described as wrong or incomplete and needing to be rewritten. The conceptual idea was nevertheless stated: the forest area disturbed by fire starts again at age 0.

### 2.6 Boundary condition after fire

The boundary condition after disturbance was discussed conceptually.

After a fire:

- patches restart at age 0;
- the patches are initially herb-dominated;
- tree seeds or very small tree biomass may already be present but are considered negligible in biomass;
- the future evergreen/deciduous trajectory of disturbed patches is assigned according to the relative biomass of evergreen and deciduous trees in the landscape.

Example given in the transcript:

- if evergreens represent 80% of the tree biomass in the landscape, then 80% of disturbed patches would follow the evergreen trajectory and 20% the deciduous trajectory.

The transcript clarifies that disturbed patches “all go to herb” in terms of initial dominance, but are split between evergreen and deciduous trajectories based on tree biomass proportions.

### 2.7 Growth term modification from Madingley

The discussion refers to the original Madingley vegetation growth term as depending on monthly net primary productivity, constants, non-structural tissue, and a functional-group fraction.

The proposed modification is:

- drop the functional-group fraction \(F_i\) from the growth term;
- introduce a competition factor \(C_i\).

Vincent’s personal notes record this explicitly as: “Drop the Functional group factor \(F_i\) from the Growth term and introduce \(C_i\) for competition.”

### 2.8 Competition factor and herb/tree succession

The proposed competition factor is used to model succession between herbs and trees.

For herbs:

- when there is no tree biomass, the factor is near 1;
- as tree biomass increases, the herb factor declines;
- the proposed shape was described as a decreasing curve;
- the transcript links this to Beer’s law / light availability: more foliage reduces light at the ground.

For trees:

- the tree competition factor follows the opposite pattern;
- trees increase as succession proceeds and eventually exclude herbs.

The model keeps one tree type/species in a local patch, while herbs can coexist with trees within the patch.

### 2.9 Herb dispersal vs tree dispersal

The transcript states that herbs are assumed to disperse everywhere, so their dispersal does not need to be tracked explicitly.

Tree dispersal is more restricted than herb dispersal and matters because it determines the relative abundance of evergreen and deciduous trajectories after fire.

### 2.10 Herbivory and mortality

Herbivory is expected to slow vegetation growth by increasing mortality.

The transcript states that high mortality or intense herbivory can reduce biomass growth and may even make biomass decline. This was discussed as one mechanism by which herbivores could delay or prevent tree takeover.

A question was raised about whether mortality should also be dynamic in a way that allows succession. The answer given was that mortality affects the rate at which biomass grows.

### 2.11 Link to ATN

A question was raised about how the ATN group would receive vegetation information.

The answer in the transcript was that the model would follow the equation in the ATN to determine the herbivory rate, following the same principle as in the original model. The transcript states that this is compatible.

### 2.12 Climate and vegetation dominance

The vegetation dynamics need to change with climate so that some regions become dominated by evergreens and others by deciduous trees.

Parameters mentioned as relevant are:

- mortality rate;
- fraction of structural tissue;
- disturbance rate.

Disturbance rate was also discussed as depending on temperature.

The original Madingley model was described as using a polynomial depending on a parameter transcribed as **F-Frost**. The exact meaning of this parameter was not known in the meeting, but it was presumed to relate to days with negative temperature.

### 2.13 Vegetation checks and experiments

Points identified for follow-up:

- check assumptions with plant knowledge/literature;
- check whether the vegetation model produces the intended emergence;
- decide whether time since disturbance is continuous or updated yearly;
- reproduce Dominique’s local-growth graphs;
- use a Whittaker-style plot as a separate vegetation experiment;
- potentially use the herbaceous group in that experiment.

The transcript mentions that vegetation needed to choose among three options for zero biomass before running an experiment. The details of those options are not clear in the transcript.

---

## 3. Distribution / dispersal model

### 3.1 Current progress

The distribution/dispersal group had a simple diffusion function working for one species.

The current implementation uses **reflective boundaries**: at the border, species “hit a wall.”

### 3.2 Model reference discussed

The group looked into a paper name that is unclear in the transcript, transcribed approximately as “right surprise” or “riser.” The transcript says the paper was more complicated than expected, but the core idea was not too complicated.

The core dispersal structure discussed was:

- emigration from cells;
- immigration from adjacent cells;
- a dispersal/emigration rate that may depend on conditions.

### 3.3 Emigration function

The whiteboard shows an **emigration** section with a sigmoid/logistic curve and the equation:

\[
E_{ij} = d B_{ij}
\]

The transcript describes emigration as initially considered as a constant outflow, but then discusses a more complex version where dispersal changes when conditions degrade.

The condition-dependence is described as:

- when conditions are bad, dispersal approaches a maximal dispersal rate;
- “bad conditions” are represented by the difference between metabolism and growth rate within a cell, or more generally by whether the population is struggling for survival within the cell;
- if metabolism is higher than population growth rate, the population is only “fighting for survival” and cannot reproduce;
- if the inverse holds, dispersal can be near zero.

A parameter controls the steepness of the logistic function.

### 3.4 Immigration function

The whiteboard shows an **immigration** formula resembling:

\[
I = \sum E_n \; f(\text{survival movement}) \; \left(\frac{1}{4}\right)
\]

The exact notation is partly unclear, but the visible components are:

- immigration is a sum over neighbouring cells;
- the formula includes emigration from neighbouring cells;
- there is a survival-during-movement term;
- there is a division by 4.

The transcript states that immigration is computed from adjacent cells. Since four neighbours are used, the amount is divided by 4.

The whiteboard also shows a small grid/neighbourhood sketch and the notation **k = 10**, but the transcript does not clearly define what \(k = 10\) means.

### 3.5 Survival during movement

The transcript mentions a survival function during movement.

For now, survival can be set to 1. Later, it could depend on landscape type or vegetation type. A question was raised about whether moving through forest is harder than moving through grass, but no firm decision was made.

### 3.6 State-dependence and density-dependence

The transcript distinguishes the discussed dispersal from standard diffusion because the dispersal rate depends on the state of the source cell.

The discussion clarified:

- there is no preferential direction for now;
- direction could potentially be added later;
- dispersal depends only on the state of the cell the organism is in, not on knowledge of neighbouring-cell quality;
- the process was described as density-dependent or condition-dependent;
- immigration/emigration is proportional to biomass or density.

A concern was raised that dispersal may also occur when conditions are good, through spillover or reproduction. Foxes were mentioned as an example where movement may occur when breeding, and breeding happens when conditions are good.

A possible response in the discussion was that the model may need flexibility in the functional form:

- constant;
- density-independent;
- density-dependent;
- state-dependent;
- possibly allowing different shapes depending on assumptions.

No final functional form was fixed.

### 3.7 Temperature and mortality

A question was raised about whether dispersal is sensitive to temperature, for example under heat-dome conditions.

The response was that dispersal may not be directly temperature-sensitive, but temperature could affect ABM/biomass/metabolism through the ATN, which would then influence dispersal or survival.

This point was discussed but not resolved as an implementation decision.

### 3.8 Parameterization

For distribution/dispersal, the group expected to decide how to parameterize the distribution model and what is most important relative to what the modellers are doing.

The transcript identifies movement rate by species as a critical parameterization issue.

Potential empirical relationships mentioned:

- range size and body size;
- maximum speed and body size.

The transcript says a distinction will be needed between:

- birds;
- organisms that walk or crawl.

Freshwater and semi-aquatic species were mentioned briefly, but the discussion returned to terrestrial first.

### 3.9 Ocean dispersal

Ocean dispersal was briefly discussed.

The transcript says that in the ocean, currents make movement non-random and can push individuals in the direction of the current. Coastal/upwelling/open-ocean differences were mentioned. A “fall off the end of the earth” type model was mentioned as a possible way of losing individuals if they are swept into some ocean region.

No implementation decision was made.

### 3.10 Distribution/dispersal experiments

The group expected that by the end of the week, it should be possible to have a function that updates distribution at each time step according to a diffusion model.

Vincent challenged the group to cut the task into smaller testable pieces and show plots using fake data, for example focusing only on immigration or migration terms before integrating the final model.

Later discussion also says calibration with real species is a later problem. The immediate priority is to have functions that work, even with two random species.

---

## 4. ATN model

### 4.1 Current progress

The ATN group was checking equations, experimenting with prompts, and editing a README.

The workflow described was:

- provide equations to an AI tool;
- ask it to turn them into Python scripts;
- in a separate session, ask another AI/tool to compare the original paper against the scripts;
- list differences side by side;
- check those differences manually.

Judicaël was mentioned as checking the equation/script differences.

### 4.2 Outputs

The ATN outputs were initially in `.npy` format.

The group wanted to convert outputs into a more useful text format, with biomass per species through time.

A question was raised about whether plots would be produced. The response was that equations should be checked before focusing on plots.

### 4.3 Parameters and temperature dependence

The transcript says that parameters/constants come from the ATN papers.

The group stated that the model has temperature dependence for the functional response. Temperature values or parameters were described as random for now, used only to see what happens, with the intention to tune later once the team is more comfortable.

### 4.4 Generality test

A proposed eventual test is to check generality by changing:

- the parameter matrix;
- the species list;
- associated traits.

The purpose is to make sure the model still runs when those inputs change.

### 4.5 Global-variable / session-state risk

A risk was raised from previous coding experience: a model can appear to work because hidden global variables exist in the active session, then fail after the session is closed and reopened.

This was identified as something to watch for in the ATN code or configuration.

### 4.6 Species-composition variability

Dominique proposed an important next step: introduce variability in species composition so that the same food web is not used everywhere.

This was described as an obvious next step after the ATN starts working.

### 4.7 Numerical discretization

Later in the transcript, a question was raised about how the ATN equations are solved numerically.

The concern was that continuous equations must be discretized for simulation, and the discretization method can change model behaviour. It was stated that some parameter values may cause the model to fail to converge or explode.

This point was raised as important but not resolved in the transcript.

---

## 5. Software architecture

### 5.1 Purpose of the architecture discussion

Alexis introduced the architecture discussion as a short software-engineering course. The purpose was to think differently about the project so that model components can be modified and updated more easily.

The stated motivation was to avoid reproducing the problems encountered with Madingley: hard-to-change code where processes are intertwined and difficult to isolate.

### 5.2 Janitor/store/cart metaphor

Alexis used a janitor metaphor.

In the first scenario, the janitor calls the store for specific cleaning supplies every time he enters a new room. He knows the building from experience, but the workflow is hard to transfer, especially to an intern who does not share the same language or implicit knowledge.

In the second scenario, the janitor carries a cart containing commonly used tools and products. The cart can be updated when a new need appears, and useful tools can be reused across rooms.

The analogy was then mapped to the model:

- the old janitor represents Madingley-style architecture;
- the “store” represents data sources in the project;
- the “cart” represents a shared object or shared structure carrying commonly needed information across processes.

### 5.3 Process registration and modularity

The proposed architecture uses processes that can be registered and chained.

Examples of processes mentioned:

- dispersal;
- vegetation;
- reproduction;
- fire or other future processes.

The transcript describes a code structure where new processes are registered by passing functions. Processes can be removed by removing the corresponding line from the registry, meaning the process no longer acts in the simulation.

The advantage discussed is that processes become easier to change, replace, and test.

### 5.4 Common inputs and outputs

The architecture discussion proposed that processes should share common input/output structures.

A process may need access to:

- the list of individuals or species/functional groups being modelled;
- available resources;
- time;
- the state of the environment;
- coordinates or spatial information, at least in the older example shown.

The process should return or update a standard state so the simulation can continue.

The transcript says passing or returning extra information can feel redundant, but makes the system easier to modify.

### 5.5 Definition of architecture

In the discussion, architecture was defined as the structure of information and processes, including how they are organized through:

- files;
- modules;
- functions;
- classes;
- guidelines for how components interact.

The transcript also emphasizes dependency direction: inner-core changes should not break the outer core, and rules should be established about how information is exchanged.

### 5.6 Modularity, atomicity, hierarchy

The discussion distinguishes modularity from simply nesting functions.

Points made:

- modularity helps maintainability and reworking;
- isolated modules are easier to test;
- sequential algorithms are easier to modularize than nested or hierarchical structures;
- hierarchy is acceptable, but the points where information enters and leaves must be clear;
- the call sequence is not necessarily the same as the designed system.

A large-matrix / nested-blocks example was discussed as a separate issue. The general response was that there is no single answer; maintainability depends on deciding which parts should be kept together and which parts may need to be replaced later.

### 5.7 “Donut” / layered architecture

Vincent described a layered or concentric architecture.

Components mentioned:

- state variables or state classes at the centre;
- processes that consume a previous state and return a new state;
- simulation layer handling time;
- spatial layer, with uncertainty about where exactly distribution/spatial processes should fit;
- data layer;
- file/folder organization.

The transcript says processes should focus on “where it is happening right here right now” and not concern themselves with time or space unless that is their role.

### 5.8 State classes and data contracts

The group needs to decide on state classes.

State classes may include:

- biomass;
- species lists;
- variables needed by processes at a given time and place.

These were described as **data contracts**: the way processes talk to one another.

The transcript says the state classes should be computationally efficient. NumPy arrays were mentioned as a possible representation, but exchanging arrays repeatedly could become computationally expensive. This was identified as a future decision, not something to optimize immediately.

### 5.9 Tests

The architecture discussion returned to testing.

The desired codebase should ideally include:

- tests for individual processes;
- integration tests that combine multiple processes;
- tests to check whether process behaviour matches expectations.

The transcript says this is not necessarily for immediate implementation, but it should be considered later.

---

## 6. AI workflow reflection

The transcript proposes a plenary discussion at the end of the day to discuss how AI was used, how it helped, what challenges were encountered, and which workflows were good or bad.

Participants were asked to write this down inside their experiment folder so they could later show it to the group.

---

## 7. End-of-afternoon expectations by group

### Vegetation

Expected direction for the afternoon:

- finish or clarify the formulation;
- choose among three zero-biomass options mentioned in the transcript;
- run and explore first results if possible;
- possibly start the next vegetation experiment, linked to the Whittaker plot and herbaceous dynamics.

### ATN

Expected direction:

- continue checking equation-to-code differences;
- ensure equations are correct before focusing on plots;
- convert outputs into useful formats;
- later test generality by changing parameter matrices, species lists, and traits.

### Distribution/dispersal

Expected direction:

- make decisions about parameterization;
- cut the problem into smaller experiments;
- show small plots with fake data if possible;
- eventually build a function that updates distributions at each time step according to a diffusion model.

---

## 8. Open questions explicitly raised

### Vegetation

- Should shrubs be added as a fourth functional group beyond the three retained groups?
- Should time since fire/disturbance be continuous or updated yearly?
- What is the correct boundary-condition equation after fire?
- How should zero biomass be handled?
- How exactly should the original Madingley frost-related parameter be interpreted?
- How should vegetation assumptions be checked with plant-specific knowledge?
- Does the model actually produce the intended emergent vegetation patterns?

### Distribution/dispersal

- Should dispersal be constant, density-independent, density-dependent, or state-dependent?
- Should dispersal increase when conditions are bad, when conditions are good, or both?
- Should dispersal be directly temperature-sensitive?
- Should survival during movement depend on landscape or vegetation type?
- How should movement rates be parameterized for each species?
- How should birds differ from walking/crawling organisms?
- How should ocean currents affect dispersal?
- What does **k = 10** on the whiteboard refer to?

### ATN

- Are the AI-generated Python scripts faithful to the original ATN paper equations?
- Which numerical integration/discretization method is used?
- Could discretization cause non-convergence or exploding dynamics?
- How should temperature parameters be set after the current random exploratory values?
- How should species composition vary across space?

### Architecture

- What are the standard state classes?
- What are the process data contracts?
- Which layer is responsible for time?
- Which layer is responsible for space?
- Where does distribution/dispersal belong in the architecture?
- How should fire be organized, and which component is responsible for it?
- What should be unit-tested now, and what should wait for integration tests later?

---

## 9. Harfoot / Madingley context used only for orientation

The afternoon discussion repeatedly refers to Madingley as the inherited model structure being modified or replaced. The relevant Harfoot/Madingley context is:

- Madingley is a mechanistic General Ecosystem Model;
- it represents autotrophs as stocks and heterotrophs as cohorts;
- core processes include primary production, eating, metabolism, growth, reproduction, dispersal, and mortality;
- dynamics are spatially explicit on grid cells;
- dispersal exists in the original model, including diffusive and active/responsive forms, and marine advective dispersal;
- Harfoot et al. also identify future development needs around detrital loops, directed dispersal, hibernation/stasis, complex plant–herbivore and predator–prey interactions, data-constrained parameterization, rigorous evaluation, numerical methods, and flexible model infrastructure.

This context supports why the afternoon focused on vegetation, dispersal, ATN, and architecture, but it does not add decisions beyond the transcript and notes.

