# External data needs for calibrating the maximum dispersal rate (*a*)

This document lists external data needed to calibrate the maximum per-capita dispersal
rate parameter *a* in `diffuse_density_dependent.r`, in order of priority. For each
data need, potential literature and database sources are listed.

> **Note:** Citations are drawn from knowledge available at the time of writing (May 2026).
> Verify journal titles, volumes, and page numbers before citing in a manuscript.

---

## Priority 1 — Observed dispersal distances

**Why:** Direct empirical measurements of natal dispersal distance are the most
defensible calibration target for *a*. All traits listed below are proxies; this is the
response variable.

**What we need:** Mean or median natal dispersal distance (km or grid-cells) per
species, from mark-recapture or telemetry studies.

**Potential sources:**
- Clobert, J. et al. (2012) *Dispersal Ecology and Evolution.* Oxford University Press.
  Broad cross-taxa synthesis with dispersal distance compilations and drivers.
- Sutherland, G.D. et al. (2000) A comparative study of vertebrate dispersal and
  colonization rates. *Journal of Animal Ecology* 69:564–575.
- Paradis, E. et al. (1998) Patterns of natal and breeding dispersal in birds.
  *Journal of Animal Ecology* 67:518–536. Bird-specific natal dispersal distances.
- Trochet, A. et al. (2016) Population-level variation in dispersal ability: a
  standardised assessment across 21 amphibian species. *Oikos* 125:1661–1671.
  Dispersal distances for amphibians and reptiles.

---

## Priority 2 — Home range size

**Why:** Home range area is the closest single proxy for how much space an individual
occupies per unit time, and scales with *a* more directly than any other indirect trait.

**What we need:** Mean home range area (km²) per species.

**Potential sources:**
- Jones, K.E. et al. (2009) PanTHERIA: a species-level database of life history,
  ecology, and geography of extant and recently extinct mammals. *Ecology* 90:2648.
  Comprehensive mammal home range data.
- Jetz, W. et al. (2004) The scaling of animal space use. *Science* 306:266–268.
  Cross-taxa allometric scaling of home range with body mass — useful for deriving
  estimates where direct measurements are absent.
- Myhrvold, N.P. et al. (2015) An amniote life-history database to perform comparative
  analyses with birds, mammals, and reptiles. *Ecology* 96:3109. Includes home range
  estimates for reptiles and birds.

---

## Priority 3 — Movement mode (locomotion type)

**Why:** The largest gap in the existing dataset (`SpeciesTraitsFull.csv`). Flying,
running, swimming, and fossorial species differ by orders of magnitude in *a* even at
similar body masses. Taxonomic Class only partially captures this.

**What we need:** A categorical variable: `{flyer, swimmer, runner, fossorial}` or
similar.

**Potential sources:**
- Wilman, H. et al. (2014) EltonTraits 1.0: species-level foraging attributes of the
  world's birds and mammals. *Ecology* 95:2027. The `ForStrat` (foraging stratum)
  variable partially proxies locomotion mode for birds and mammals.
- Myhrvold, N.P. et al. (2015) (see above). Includes locomotion and activity mode for
  amniotes.
- Manual coding from IUCN Red List species accounts is often necessary for reptiles
  and amphibians, where coverage in trait databases is thinner.

---

## Priority 4 — Morphological dispersal proxy for birds (Hand-Wing Index)

**Why:** For Aves — likely the highest-dispersal class in the dataset — wing morphology
is a well-validated proxy for natal dispersal ability and substantially improves
prediction of *a* beyond body mass and Class alone.

**What we need:** Hand-Wing Index (HWI) per bird species.

**Potential sources:**
- Tobias, J.A. et al. (2022) AVONET: morphological, ecological and geographical data
  for all birds. *Ecology Letters* 25:581–597. Contains HWI, wing length, and other
  morphological traits for all ~11,000 bird species.
- Sheard, C. et al. (2020) Ecological drivers of global gradients in avian dispersal
  inferred from wing morphology. *Nature Communications* 11:2463. Validated HWI as
  a cross-species predictor of natal dispersal distance.

---

## Priority 5 — Mass-specific metabolic rate

**Why:** In `diffuse_density_dependent.r`, the sigmoid inflection point is the species
metabolic rate *x_i*. Using the same empirical metabolic rate to set both *x_i* and
to inform the calibration of *a* ensures internal consistency between the two
parameters. Without this, *a* and *x_i* may be estimated independently and contradict
each other.

**What we need:** Basal or field metabolic rate (W or J/s) per species, or taxon-level
allometric intercepts.

**Potential sources:**
- Brown, J.H. et al. (2004) Toward a metabolic theory of ecology. *Ecology*
  85:1771–1789. Allometric scaling equations for metabolic rate — can generate estimates
  for any species from body mass and temperature alone.
- de Magalhães, J.P. & Costa, J. (2009) A database of vertebrate longevity records and
  their relation to other life-history traits. *Journal of Evolutionary Biology*
  22:1770–1774. AnAge database; includes metabolic rates for ~4,000 vertebrate species.
- Brose, U. et al. (2006) Consumer–resource body-size relationships in natural food
  webs. *Ecology* 87:2411–2417. ATN-specific allometric scaling parameters for
  metabolic rate, directly compatible with the model structure used here.

---

## Priority 6 — Generation time

**Why:** *a* is a per-time-step rate. If the model time step represents different
biological durations for different species, *a* values are not comparable across taxa.
Generation time allows rescaling so that *a* has consistent biological meaning across,
for example, an annual amphibian and a long-lived bird of prey.

**What we need:** Mean generation time (years) or age at first reproduction per species.

**Potential sources:**
- Myhrvold, N.P. et al. (2015) (see above). Age at first reproduction, gestation
  length, and other life-history timings for ~4,900 amniote species.
- Salguero-Gómez, R. et al. (2016) COMADRE: a global database of animal demography.
  *Journal of Animal Ecology* 85:371–384. Population matrix models from which
  generation time can be derived.

---

## Priority 7 — Thermal tolerance (ectotherms only)

**Why:** For Reptilia and Amphibia, movement capacity is temperature-limited. Thermal
tolerance range constrains how many hours per day an animal can be active and
dispersing, and is not predictable from body mass alone.

**What we need:** Critical thermal maximum and minimum (CTmax, CTmin) per species.

**Potential sources:**
- Bennett, J.M. et al. (2018) GlobTherm: a global database on thermal tolerances for
  aquatic and terrestrial organisms. *Scientific Data* 5:180022. Covers ~2,000
  vertebrate and invertebrate species.

---

## Priority 8 — Social behavior / territorial system

**Why:** Territorial species disperse to establish or defend new territories (higher *a*);
colonial or philopatric species show strong site fidelity (lower *a*). After controlling
for body mass and movement mode, social system adds independent predictive power for
dispersal propensity.

**What we need:** Categorical variable `{solitary-territorial, group-territorial,
colonial, philopatric}` or continuous group size.

**Potential sources:**
- Jones, K.E. et al. (2009) PanTHERIA (see above). Includes `SocialGrpSize` for
  mammals.
- Wilman, H. et al. (2014) EltonTraits 1.0 (see above). Partial coverage for birds
  and mammals.
- Coverage for reptiles and amphibians largely requires manual extraction from
  species accounts.

---

## Summary table

| Priority | Trait | Class coverage | Best source |
|---|---|---|---|
| 1 | Observed dispersal distance | All | Clobert et al. 2012; Sutherland et al. 2000 |
| 2 | Home range size | Mammals, birds | PanTHERIA; Jetz et al. 2004 |
| 3 | Movement mode | All (partial) | EltonTraits; Myhrvold et al. 2015 |
| 4 | Hand-Wing Index | Aves only | AVONET; Sheard et al. 2020 |
| 5 | Metabolic rate | All | Brown et al. 2004; Brose et al. 2006 |
| 6 | Generation time | All | Myhrvold et al. 2015; COMADRE |
| 7 | Thermal tolerance | Reptilia, Amphibia | GlobTherm; Bennett et al. 2018 |
| 8 | Social behavior | Mammals (best) | PanTHERIA; EltonTraits |
