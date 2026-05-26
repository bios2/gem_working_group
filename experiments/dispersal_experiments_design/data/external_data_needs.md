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
