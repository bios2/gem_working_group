# How more sophisticated leaf biomass simulations can increase the realism of modelled animal populations

**Jens Krause**$^{a,*}$, **Mike Harfoot**$^{b}$, **Selwyn Hoeks**$^{c}$, **Peter Anthoni**$^{a}$, **Calum Brown**$^{a}$, **Mark Rounsevell**$^{a}$, **Almut Arneth**$^{a}$

$^{a}$ KIT-Campus Alpin, Institute of Meteorology and Climate Research (IMK-IFU), Garmisch-Partenkirchen, Germany  
$^{b}$ UN Environment Programme World Conservation Monitoring Center, Cambridge, United Kingdom  
$^{c}$ Department of Environmental Science, Radboud University Nijmegen, Netherlands  

*Ecological Modelling* 471 (2022) 110061  
Available online 2 July 2022  
https://doi.org/10.1016/j.ecolmodel.2022.110061  
Received 25 October 2021; Received in revised form 19 June 2022; Accepted 21 June 2022

> An update to this article is included at the end.

---

**Keywords:** Modelling animal populations in terrestrial ecosystem; Advances/refinement in methods for ecological modelling

---

## Abstract

Animal biodiversity, and its key roles in ecosystem state and functioning, is facing critical challenges in the wake of anthropogenic activities. It is urgently necessary to improve understanding of the interconnections between animals and the vegetation within ecosystems. Process-based modelling has shown to be a mighty tool in making assessments on ecological processes. We assess the effect of different vegetation models on simulated animal biodiversity by replacing the vegetation module within Madingley, a multi-trophic model of functional diversity with LPJ-GUESS, a dynamic global vegetation model. We compare the output metrics of the model system to Madingley's default version for four ecosystem types around the globe and analyse whether the realism of the simulation results increased as a result of the coupling between Madingley and LPJ-GUESS. Simulated animal populations react to the coupling by shifting towards smaller individuals with a higher abundance. General shifts in body mass and animal distributions can be traced back to ecological processes, allowing in-depth analysis of heterotrophic responses to changes in leaf biomass. We also derive power-law relationships for herbivory to NPP and herbivore biomass to NPP and conclude that the coupled model system simulates animal populations that follow reasonable power-laws which are similar to power-laws derived from empirical data. Our results indicate that developing process-based model systems is a viable way to assess multi-trophic interconnections between animal populations and the ecosystems vegetation.

---

## 1. Introduction

Animals can play a key role in controlling the state and function of all terrestrial ecosystems (Schmitz et al., 2014, 2018). By consuming autotrophic biomass, herbivory enhances light transfer into plant canopies and affects net carbon assimilation. Cascading trophic effects triggered by top predators or the largest herbivores propagate through food webs, regulating levels of herbivory and affecting soil carbon and nitrogen cycling through excreta and dead bodies (Schmitz et al., 2014). However, anthropogenic activities such as habitat modification, harvesting and, increasingly, anthropogenic climate change are driving large declines in biodiversity (Cardinale et al., 2012; Pacifici et al., 2015; Arneth et al., 2020). Therefore, a better understanding of how changes and losses of plant and animal functional diversity will affect ecosystem functioning is needed from the perspective of both an ecosystems role for both climate change mitigation and climate change adaptation. However, the diverse mixture of species and the interconnectivity in food webs makes quantifying the link between ecosystem biogeochemical cycling and an ecosystems functional diversity a challenge (Schmitz et al., 2018).

Animals comprise lower total biomass overall compared to plants (Bar-On et al., 2018) which could easily result in the assumption that animals play a relatively minor role in terrestrial ecosystem carbon cycling. But this view is being challenged by an increasing body of literature. For instance, for eastern Africa (Holdo et al., 2009) extrapolated impacts of strongly reduced wildebeest numbers due to the rinderpest into the future, and hypothesised that major future wildebeest die-offs might lead to a significant increase in grass biomass, leading to changes in fire regimes and a potential shift of ecosystems into a net carbon source. Hooper et al. (2012) suggested that the absence or presence of herbivores and carnivores within an ecosystem can lead to similar magnitude impacts as caused by other environmental change drivers such as global warming, enhanced CO$_2$ or nitrogen pollution. Wilmers and Schmitz (2016) estimated, based on measurements, that wolves preying upon moose in Northern American boreal forests enhance net primary productivity (NPP) and net ecosystem carbon uptake by up to 30% by decreasing moose CO$_2$ respiration and increasing growth of deciduous trees, due to releasing herbivory pressure. Sobral et al. (2017) carried out an analysis based on observations in the Amazonian region which showed animal diversity influencing the carbon cycle directly via metabolism and indirectly via seed dispersal and altering the vegetations species composition. Nevertheless, so far studies and extrapolations of their results that seek to link animal impacts to ecosystem biogeochemical cycles are typically limited by their temporal or spatial scale (Holdo et al., 2009; Staver and Bond, 2014; Sobral et al., 2017) or focus on a few distinct animal species (Schmitz et al., 2018; Nichols et al., 2016) so limiting the general applicability of their conclusion.

Process-based models can help to overcome observation limitations and improve our knowledge about the complex interactions within ecosystems by expanding the scope to the whole ecosystem. A few model studies have started to explore the role of animals. Berzaghi et al. (2019) included elephant disturbances into an ecosystem demography model and showed that the resulting reduction in tree stem density altered competition for light and water and actually enhanced above-ground biomass. Pachzelt et al. (2015) coupled a physiological grazer population model with a dynamic vegetation model and identified NPP and the length of the dry season as the main determinants of grazer population size in African savannah ecosystems. They found that the presence or absence of grazers substantially impacts the standing grass and tree vegetation biomass, as well as the area vulnerable to burning. Dangal et al. (2017) coupled a herbivore population dynamics model with a dynamic land ecosystem model and were able to preproduce observed herbivore populations for livestock species. The presence of these herbivores reduced the NPP in all nine modelled regions while revealing a vulnerability of herbivores to climate extremes.

While these studies contribute important information to the understanding of animals' roles in ecosystem functioning, they also focus on distinct animal species or species classes and neglect complexities such as the top-down control predators impose on herbivores and omnivores and thus food-webs in ecosystems. The Madingley model (Harfoot et al., 2014) adopts a different approach, linking fluxes of biomasses between above-ground autotrophs and heterotrophs by modelling the fundamental processes affecting the ecology of these organisms. Trophic structures in the modelled ecosystems emerge from the representation of individual and community level biological processes, and interactions between organisms and their environment. Madingley has been shown to reproduce patterns and functional processes in animal populations and has been used to model terrestrial and marine ecosystems (Harfoot et al., 2014; Enquist et al., 2020) and to analyse complex food chains (Hoeks et al., 2020). At present, the vegetation, which drives herbivory, is modelled by an empirical carbon cycling model (Smith et al., 2013), in which primary production is determined by the Miami Model (Lieth, 1975). However, the Miami Model does not capture important dynamics such as interannual variability and trends in vegetation composition and productivity introduced by weather, climate change, carbon–nitrogen interactions and changes in atmospheric CO$_2$, which modern dynamic vegetation models do. This limits Madingleys applicability to explore future global environmental change impacts on ecosystem functional diversity.

Here, we couple Madingley to the advanced dynamic global vegetation model (DGVM) LPJ-GUESS to test how more complex vegetation dynamics affects characteristics of simulated animal populations. The work presented here is the first step towards creating a fully coupled multi-trophic model of functional diversity which combines the advantages of complex process-based vegetation and animal modelling, allowing comprehensive assessment of animal-vegetation feedbacks in ecosystems.

---

## 2. Methods

### 2.1. The Madingley model

Madingley represents herbivores, omnivores and carnivores of all size classes, from very large to very small. The heterotrophs are organised into functional groups according to a set of categorical functional traits related to feeding modes, metabolic pathways, reproduction and movement. The model is defined at the level of individual organisms. To avoid impractical computational requirements, within each functional group, individuals are grouped together into cohorts, in which all individuals possess the same set of categorical traits plus identical continuous traits defining their juvenile, adult and current body mass and age. These continuous traits vary across cohorts which share a functional group. Cohorts of individuals undergo a set of ecological processes: metabolism, predation, feeding, reproduction and mortality (Harfoot et al., 2014). Dividing the model domain into distinct grid cells makes it possible for cohorts to disperse between those grid cells.

The lowest levels of the heterotrophic pyramid, the herbivores, are feeding upon evergreen and deciduous leaf biomass stocks simulated by the carbon cycle model of Smith et al. (2013). While the deciduous stock is showing seasonal availability, herbivores assimilate relatively more biomass when compared to the same amount of evergreen leaf biomass being consumed. Evergreen and deciduous stock therefore represent a trade-off between food availability and food quality. At higher trophic levels, carnivores prey upon other animal cohorts, with a prey size preference. Omnivores feed from the leaf biomass stocks and heterotroph cohorts alike but are less effective doing so, due to a lower assimilation rate for both herbivory and carnivory. These feeding processes are limited to the stocks and cohorts of the respective grid cell. Herbivory is not uniform (e.g, different types of herbivores access different amounts and types of leaf biomass in an ecosystem), processes that are currently not yet captured in Madingley. As a conservative assumption (implemented in Madingley as $\phi_{herb,f} = 0.1$, Table 6 in Harfoot et al. (2014)), the total biomass consumed is thus limited to 10% of the gridcells leaf-biomass stock.

The growth of individuals is handled as the difference between their nutrition uptake and metabolic cost. Once an individual reaches its maturity, it starts accumulating its nutrition uptake towards reproduction. For this study, we adopted and modified the latest version of the model described by Hoeks et al. (2020).

So far, the leaf biomass stocks, which apply a bottom-up control on the simulated trophic pyramid, are modelled by the carbon cycle model of Smith et al. (2013), which includes the Miami Model by Lieth (1975). The Miami Model predicts the total annual net primary production (NPP) based on the yearly mean temperature and total precipitation sum in the following way:

$$NPP(T, p) = \min\!\left(\frac{0.916645}{1 + \exp(0.23747 - 0.006 \cdot T)},\; \frac{0.916645}{1 - \exp(0.00118 \cdot p)}\right) \tag{1}$$

In Madingley, the annual simulated NPP is converted to wet matter, which animals consume, via an estimate of 9.813 g$_{wet matter}$/g$_C$ (Kattge et al., 2011). The annual values are then split into monthly fractions using a seasonality factor. Leaf mortality and carbon partitioning to leaves is also parameterised through temperature related factors. The resulting monthly wet matter leaf increment is divided between evergreen and deciduous vegetation according to the fraction of the year prone to frost in the specific grid cell.

### 2.2. Forcing Madingley with LPJ-GUESS vegetation data

LPJ-GUESS is a biogeochemical dynamic global vegetation model that combines the advantages of a global macroscopic nutrient cycle and carbon assimilation model with the advantages of an individual-based growth model (Smith et al., 2014; Wårlind et al., 2014). Plant individuals are represented by age-cohorts, which share key ecological properties and traits, categorising them into plant functional types (PFTs) distinguished by their bioclimatic preferences, photosynthesis pathways and growth strategies. Plant cohorts compete for nitrogen, light, water and space. Annual processes, like biomass allocation, leaf, root and sapwood turnover, disturbances or mortality are simulated at the beginning of a year. Short-term processes like soil hydrology, stomata regulation, photosynthesis, plant respiration, decomposition and phenology are simulated on a daily basis. The overall accuracy of the modelled processes and interactions has been evaluated extensively (Hickler et al., 2006; Wramneby et al., 2008; Lindeskog et al., 2021; Smith et al., 2014) and has shown LPJ-GUESSs skills in capturing large-scale vegetation patterns and vegetation dynamics. LPJ-GUESSs physiological vegetation process representation goes far beyond Madingleys current approach to estimate the biomass of deciduous and evergreen leave.

LPJ-GUESS requires multiple years of climate data to ensure sufficient variability for computing the models stochastics, interannual disturbances and for the model spin-up. For our purposes here, we chose to cycle the CRUJRA v2.1 climatology (Harris, 2020) based on the years 1950 to 1979. This approach allows us to compare vegetation simulated with LPJ-GUESS and Miami, since Miami was parameterised based on empirical data available during the 1960s. We cycled atmospheric CO$_2$ concentrations ranging from 311 to 335 ppm during the LPJ-GUESS simulations. Since the CRUJRA dataset has a 6-h resolution, we used daily postprocessing methods (supplementary information S1).

LPJ-GUESS simulates 12 woody and 2 grassy plant functional types, which were converted, based on their phenology, into the evergreen and deciduous leaf-biomass stocks that are available for herbivory in Madingley (Table 1). Since LPJ-GUESS is running on a daily timestep, we account leaf biomass from all PFT cohorts during every timestep and produce a monthly average of the total evergreen and deciduous leaf biomass per grid cell. We then pass this information to Madingley using a file-transfer method and replace Madingleys parameterised autotroph stocks during every timestep. Fig. 1 shows a simplified technical blueprint to illustrate how we realised the coupling of the two models.

The latest version of Madingley (Hoeks et al., 2020) cycles a single year of environmental input data throughout the simulation. Climate does not only affect the vegetation, but also the animals metabolism due to its ecto- and endothermic thermoregulation method. To ensure that the animals in Madingley are exposed to the same climate used for the vegetation prediction in LPJ-GUESS, we expanded the length of Madingleys environmental forcing from the single-year cycle to a 30-year time series. We used the same CRUJRA v2.1 climate data but due to Madingleys monthly timestep, we needed to process the data monthly (supplementary information S1).

Both model setups represent a solid and stable version of the model. The computational time necessary to carry out a single simulation did not change noticeably when coupling Madingley with LPJ-GUESS. The additional time needed to precompile the LPJ-GUESS is small in comparison to the average Madingley runtime. Thus, we did not increase the computational cost for a single run significantly. However, the Madingley simulations with the C++ version are quite time and memory-consuming so we decided to carry out smaller scale simulations at different locations instead of a continental or global run. Other recent studies with Madingley followed a similar approach (Newbold et al., 2020) or used a coarser grid for their simulations (Hoeks et al., 2020).

**Table 1.** Plant functional types in LPJ-GUESS with evergreen or deciduous leaf biomass. Summergreen, raingreen and grassy PFTs, are accounted towards Madingley's deciduous stock.

| PFT-ID | Climate type | Lifeform   | Leaf physiognomy | Phenology              | Light behaviour  |
|--------|-------------|------------|-----------------|------------------------|-----------------|
| BNE    | Boreal      | Tree       | Needle-leaved   | Evergreen              | Shade tolerant  |
| BINE   | Boreal      | Tree       | Needle-leaved   | Evergreen              | Shade intolerant|
| TeNE   | Temperate   | Tree       | Needle-leaved   | Evergreen              | Shade tolerant  |
| TeBE   | Temperate   | Tree       | Broadleaved     | Evergreen              | Shade tolerant  |
| TrBE   | Tropical    | Tree       | Broadleaved     | Evergreen              | Shade tolerant  |
| TrIBE  | Tropical    | Tree       | Broadleaved     | Evergreen              | Shade intolerant|
| BLSE   | Boreal      | Low Shrub  | Needle-leaved   | Evergreen              | Shade intolerant|
| BNS    | Boreal      | Tree       | Needle-leaved   | Summergreen            | Shade intolerant|
| IBS    | Boreal      | Tree       | Broadleaved     | Summergreen            | Shade intolerant|
| TeBS   | Temperate   | Tree       | Broadleaved     | Summergreen            | Shade tolerant  |
| TrBR   | Tropical    | Tree       | Broadleaved     | Raingreen              | Shade intolerant|
| BLSS   | Boreal      | Low Shrub  | Broadleaved     | Summergreen            | Shade intolerant|
| C3G    | All         | Grass      | C3 grass        | Cold season or yearlong|                 |
| C4G    | All         | Grass      | C4 grass        | Warm season            |                 |

### 2.3. Madingley simulation setup

We compared Madingley output simulated with its standard Miami-based vegetation with simulations in which, at the beginning of each timestep, the leaf biomass was replaced by the modelled leaf biomass from LPJ-GUESS, summed up for deciduous and evergreen pfts. Four experiments are presented here to test how LPJ-GUESS vs. Miami vegetation affects model stability, community structures and individual traits. To investigate how ecosystem productivity and climatic conditions affect the coupling responses, we performed simulations for four locations representing three different types of forest ecosystems and one savanna ecosystem (Table 2). Each location was resolved by a grid of three-by-three cells, spanning a 1.5° latitude by 1.5° longitude model domain with a resolution of 0.5° in each direction. Both Madingley and LPJ-GUESS run on a 0.5° resolution, so we did not need to up- or downscale any of the exchanged data.

**Table 2.** Description coordinates of the four experiment locations.

| Simulation tag              | Vegetation model | Longitude range | Latitude range | Ecosystem type              |
|-----------------------------|-----------------|-----------------|----------------|-----------------------------|
| Hyytiälä (Finland) FIN M-M  | Miami Model     | 24°E–25°E       | 61°N–62°N      | Boreal Coniferous Forest    |
| FIN M-LPJG                  | LPJ-GUESS       | (0.5° resolution)|               |                             |
| Waldstein (Germany) GER M-M | Miami Model     | 11°E–12°E       | 50°N–51°N      | Temperate Mixed Forest      |
| GER M-LPJG                  | LPJ-GUESS       | (0.5° resolution)|               |                             |
| Lake Mburo Nat. Park (Southern Uganda) UGA M-M | Miami Model | 28°E–29°E | 0°N–1°N  | Tropical Rainforest         |
| UGA M-LPJG                  | LPJ-GUESS       | (0.5° resolution)|               |                             |
| Pretoria (South Africa) SAF M-M | Miami Model | 28°E–29°E       | 26°S–25°S      | Subtropical Savanna         |
| SAF M-LPJG                  | LPJ-GUESS       | (0.5° resolution)|               |                             |

In all four studies, we compared both setups, Madingley+Miami-Model (M-M) and Madingley+LPJ-GUESS (M-LPJG) and check what effects our coupling implies on the model system.

All simulations include the same animal functional group definitions. These groups are distinguished by their feeding mode, their reproductive strategy and their thermoregulation strategy (Table 3). Each grid cell was initialised with 100 cohorts of each functional group with initial body masses that are drawn randomly from the cohorts juvenile and adult body mass. During a simulation, cohorts die out or are newly created in response to predation, mortality and reproduction. The total number of cohorts that are allowed to coexist in one grid cell is limited to 1000.

**Table 3.** Animal functional groups in Madingley and their corresponding key ecological traits.

| Feeding mode | Reproductive strategy | Thermoregulation strategy | Min. body mass [g] | Max. body mass [g] | Herbivory assimilation efficiency | Carnivory assimilation efficiency |
|-------------|----------------------|--------------------------|-------------------|-------------------|----------------------------------|----------------------------------|
| Herbivore   | Iteroparity          | Endotherm                | 1.5               | 5,000,000         | 0.8                              | 0.0                              |
| Omnivore    | Iteroparity          | Endotherm                | 3                 | 200,000           | 0.65                             | 0.65                             |
| Carnivore   | Iteroparity          | Endotherm                | 3                 | 400,000           | 0.0                              | 0.8                              |
| Herbivore   | Iteroparity          | Ectotherm                | 0.0004            | 10                | 0.8                              | 0.0                              |
| Omnivore    | Iteroparity          | Ectotherm                | 0.0004            | 20                | 0.65                             | 0.65                             |
| Carnivore   | Iteroparity          | Ectotherm                | 0.0008            | 20                | 0.0                              | 0.8                              |
| Herbivore   | Semelparity          | Ectotherm                | 0.0004            | 1000              | 0.8                              | 0.0                              |
| Omnivore    | Semelparity          | Ectotherm                | 0.0004            | 2000              | 0.65                             | 0.8                              |
| Carnivore   | Semelparity          | Ectotherm                | 0.0008            | 2000              | 0.0                              | 0.65                             |

For all experiments, we carried out runs at each of the four locations from Table 2. All the data extracted from the model were collected at the end of each timestep and thus represent the state of the ecosystem after all ecological processes have been executed. This is true for all animal functional groups and autotroph stocks alike.

In the first experiment, we ran an ensemble of ten long-term simulations to investigate the development of leaf biomass and heterotrophic functional-group dynamics. Each simulation covered 1000 years of simulation time. The objective here was to test whether Madingleys model dynamics reaches equilibrium conditions with LPJ-GUESSs, as well as with Miamis vegetation input. When a simulation reaches a state where aggregated biomass densities of both autotroph stocks and heterotroph cohorts are not trending up or downwards anymore, we assume the model system has reached its dynamic equilibrium stage.

For the second experiment, we compared the canopy compositions of the simulated vegetation from both LPJ-GUESS and the Miami Model. The objective here was to visualise the differences between the two vegetation models. Since these differences are the key drivers for all coupling-related changes, understanding them is essential for the interpretation of all other experiments.

In our third experiment, we focussed on how the coupling alters the basic shape of the trophic pyramid and how it affects the animal population on a community level as a whole. To do so, we investigated the biomass pools of the trophic pyramid and the fluxes between those pools. We used an ensemble of ten simulations to capture stochastic variation within the pools. Shifts in leaf biomass density are likely to affect the size distribution of animals as a result of the bottom-up regulation of autotrophs in the trophic pyramid, so we also compared the emergent size distribution of cohorts driven by the two vegetation models.

After investigating the effects of our coupling on the community level, our fourth experiment focussed on the effects on an individual level by taking a closer look at power-law relationships between individual body mass and growth rate per timestep, days needed to reach maturity, lifespan and lifetime reproduction successes.

The simulations for the third and fourth experiments differ in terms of simulation length due to the longer time needed for the UGA setup to reach its dynamic equilibrium (as found in experiment 1). The UGA simulations covered 600 years of simulation time, while all other setups covered 250 years of simulation time. We extracted detailed information on cohort characteristics during the last year of every simulation. A table with all setup descriptions and the most relevant output for each experimental setup can be found in the supporting information S4.

To set the effects of our coupling into context, we compared both M-LPJGs and M-Ms simulated NPP for the four experiment locations. Since the NPP parameterised by the Miami Model forms the foundation of the M-M simulations, we wanted to test if the NPP simulated by LPJ-GUESS is more realistic. For the comparison we used NPP data derived from the MODIS 17ASHGFv061 dataset (Running and Zhao, 2021) as the average NPP of the area corresponding to the model domains. The MODIS data layers are derived from the RUE (radiation use efficiency) concept (Hatfield and Dold, 2019) and therefore are produced by a model themselves.

We also compared the results of our simulations with power-law relationships (Cebrian, 2004), derived from a huge body of literature on net primary production and herbivore biomass and consumption by primary consumers. Most authors estimate empirical data for herbivore biomass from a mean herbivore individual weight multiplied by the herbivore abundance, which is similar to the total herbivore body mass calculation in Madingley. Empirical data for consumption by primary consumers is commonly recreated based on animal evacuation, exclosures/enclosures, parameterisation of herbivore metabolism, reconstruction methods based on bite marks and leaf growth rates. Studies using these techniques derive their herbivory rate reconstruction by comparing ecosystem states with and without reduced herbivory stress. It is worth noting that Cebrian (2004) combined herbivore consumption and detritus consumption in one pool as primary consumption, while Madingley only includes herbivores. We derived power-laws for the same functional relationships from the last ten years of each simulation.

Since long-term measurement data is available for the GER and FIN sites, we included comparisons between Madingley and LPJ-GUESS model runs and data from Cebrian (2004), Rebmann et al. (2010) and Launiainen et al. (2022) in the supplementary information (S5).

---

## 3. Results

### 3.1. Experiment 1: Long-term analysis

The left side of Fig. 2 shows the simulated biomasses resulting from the M-M simulations, while the right side shows those from M-LPJG. In M-M, the simulations reach a state of dynamic equilibrium after about 100 years. The time needed for a simulation to reach equilibrium in M-LPJG ranges from about 100 years (locations FIN, GER, SAF) to about 500 years (UGA).

Simulated leaf biomass density in the two vegetation models differs markedly. In the M-LPJG simulations, the amount of evergreen leaf biomass barely changes over the course of a year, in contrast to the seasonal oscillation seen in deciduous leaf biomass. This is clearly visible, for example, in the boreal and temperate climate regions with a high seasonality – represented by FIN and GER M-LPJG. In the FIN, GER and SAF M-M setups – both evergreens and deciduous stocks fluctuate strongly over the course of a year, while also having a higher amplitude and lower winter minima when compared to M-LPJG simulations. In both, UGA M-M and UGA M-LPJG, the leaf biomass stocks fluctuate only little. Leaf biomass stocks in SAF M-LPJG are more variable than in SAF M-M.

In the FIN M-LPJG simulation, simulated carnivores and omnivores show similar biomass densities, while omnivore biomass density exceeds carnivore biomass density in FIN M-M. All other simulations show similar hierarchies in heterotroph biomass densities between the M-M and M-LPJG simulations. In all simulations, herbivores dominate heterotroph biomass densities.

### 3.2. Experiment 2: Canopy composition

As indicated already in Fig. 2, significant differences between LPJ-GUESS and the Miami Model emerge from the simulated vegetation composition. Fig. 2 shows simulated autotrophic biomass split into the evergreen and deciduous stocks for the M-LPJG and M-M setups. A more detailed pft composition for each location can be found in the supplementary information S2. In the FIN and GER location, the Miami Model simulates large seasonal fluctuations throughout a year for both evergreen and deciduous stocks. LPJ-GUESSs simulated vegetation stocks are much more stable throughout the year for evergreen leaf biomass and show the expected seasonal fluctuations in deciduous leaf biomass. Under sub-tropical and tropical climate conditions, M-LPJG simulates substantially less deciduous leaf biomass than M-M.

In SAF, variability in evergreen leaf biomass is higher in the M-LPJG simulations and less regular, likely due to fire impact in these fire-prone environments. While M-M is parameterised to burn a fixed percentage of leaf biomass annually depending on climatic conditions, M-LPJG uses the statistically-based fire model SIMFIRE-BLAZE, which simulates fire-frequencies, fire-intensities, fire-related fluxes and responses in vegetation (Launiainen et al., 2022; Rabin et al., 2017).

### 3.3. Experiment 3: Community level analysis

All community level simulations show differences between the M-M and M-LPJG simulations (Figs. 4 and 5). The average leaf biomass is smaller in FIN M-LPJG (−16%), GER M-LPJG (−5%) and SAF M-LPJG (−15%) compared to the M-M results. In UGA, M-LPJG simulates 66% more average leaf biomass than M-M. The overall herbivore biomass of the M-LPJG simulations exceed those of the M-M simulations in every location (Fig. 3), most dominantly in UGA (+170%). Average carnivore biomass in the M-LPJG simulations is greater in FIN (+93%), SAF (+34%) and again most dominantly in UGA (+542%). Average omnivore biomass in the M-LPJG simulations is decreased in FIN (−39%) and GER (−51%), and is increased in UGA (+453%), compared to M-M.

In FIN and GER, all biomass fluxes besides L→H, H→C and C→O show similar medians in both M-M and M-LPJG simulations. In SAF all fluxes show a slight increase in M-LPJG up to +20%. In UGA, all fluxes show a significant increase in M-LPJG up to +500% (Fig. 5).

In FIN, GER and SAF with the M-LPJG coupling, the size distribution spectrum shows a tendency towards individuals with a higher abundance and lower biomasses (Fig. 6). In the M-LPJG runs at the UGA location, the high and low end of the body mass spectrum of endotherms shows higher individual numbers when compared to the M-M simulations. Especially ectotherms show larger individual numbers at the lower end of the body mass spectrum.

Under boreal and temperate climate conditions, the M-LPJG coupling showed the largest increases in the number of endotherm cohorts. In contrast, at the tropical and sub-tropical sites, ectotherm cohorts increase more in the runs using M-LPJG vegetation (Fig. 6). At all four sites, the M-LPJG coupling consistently results in larger cohort numbers of carnivore endotherms. In FIN M-M, the whole functional group of carnivore endotherms disappear in the first years of the simulation without any reestablishment afterwards. Herbivore endotherms are concentrated at the higher end of the body mass spectrum. In FIN M-LPJG, carnivore endotherms are present throughout the simulation and herbivore endotherms populate a wider range of the body mass spectrum.

### 3.4. Experiment 4: Individual level analysis

The individual level analysis shows general similarities between the M-M and M-LPJG setups for each location (Fig. 7A). Growth rates are slightly larger in the M-LPJG runs. Most notable is an increased growth rate of the smallest individuals in SAF and UGA. The lifespan of bigger individuals (individuals that have body mass >100 g) decreases in the M-LPJG simulations for all locations (Fig. 7B). Consequently, the mortality rate for cohorts containing bigger individuals increases. This effect is most dominant in UGA. In FIN and SAF, the time endotherm cohorts need to reach maturity state decreases noticeably in M-LPJG (Fig. 7C). Besides an increased growth rate and a shorter time to reach maturity, the shorter lifespan seems to be the dominant effect, leading to the cohorts having a reduced lifetime reproduction success rate in all M-LPJG simulations (Fig. 7D).

### 3.5. Comparison to external sources

At all experiment locations, M-LPJG consistently predicts the NPP better than M-M when compared to the MODIS dataset (Fig. 8). The best LPJ-GUESS predictions were found in the FIN (−9%), GER (−5%) and SAF (−3%) setup, while the only larger difference was found in UGA (+17%). The M-M simulations show a larger gap between Miamis predicted NPP and the data derived from the MODIS dataset (FIN −21%, GER −20%, SAF +45% and UGA −25%).

Beyond comparing the vegetation model outputs to the MODIS V6.1 dataset, we also made an assessment of whether the realism of the entire model system output was increased. According to Cebrian (2004), terrestrial ecosystems NPP is related to herbivore biomass and herbivore consumption by the following logarithmic power laws:

$$\log(biomass) = -3.43 + 1.30 \cdot \log(NPP) \tag{2}$$

$$\log(primary\,consumption) = 0.24 + 0.82 \cdot \log(NPP) \tag{3}$$

We were able to determine similar shaped power-law relations from the M-LPJG simulations by applying a logarithmic fit for the last ten annual sums of herbivore biomass and consumption. We derived the following power-law relationships:

$$\log(biomass) = -1.75 + 1.03 \cdot \log(NPP) \tag{4}$$

$$\log(primary\,consumption) = -1.94 + 1.55 \cdot \log(NPP) \tag{5}$$

These power-law relationships are shown in Fig. 9. Notable is the high herbivore biomass and herbivory consumption in SAF M-LPJG, which is the ecosystem with the lowest simulated productivity.

Since the Miami model predicts a smaller range of NPP throughout diverse ecosystems when compared to M-LPJG, all data points for the M-M simulations are clustered in a narrow range. This renders a logarithmic fit for herbivore biomass and herbivory consumption unreasonable, so these are not provided here. However, it is worth noting, that in the M-M simulations both herbivore biomass and herbivory consumption consistently show lower values in ecosystems with slightly higher NPP, providing weak support for an inverse relationship in that model.

---

## 4. Discussion

The large differences between M-M and M-LPJG simulations show a complex interconnectivity of coupling-related effects throughout the experiments. To analyse their validity, we first want to explain the observed changes and the underlying processes. Secondly, we want to discuss the differences to the MODIS dataset and the power laws derived by Cebrian (2004).

### 4.1. Development of the animal population

The model reaches a dynamic steady state in all M-LPJG and M-M setups. The longer time needed for the UGA M-LPJG simulation is most likely caused by the much larger amount of edible leaf biomass simulated for this location in LPJ-GUESS, compared with Miami. Cohorts of heterotrophs were initialised in the same way for all four locations. Their subsequent growth is therefore not only limited by food availability, but also by climatic conditions and predation. UGAs warm climate lengthen the active time of ectotherms and thus increases their food uptake and metabolic cost. Therefore, UGA M-LPJGs larger vegetation biomass supply funnels through to the higher trophic levels which increases predation stress for the lower trophic levels. These interactions increase the growth limit of the cohorts, and thus result in a longer period necessary to reach the model systems dynamic steady state. This effect appears to be most visible in enhancing herbivore and carnivore biomasses. Omnivores also responded to the changes in vegetation stocks, but not as much as the other groups, possibly because while omnivores have the ability to feed both from plants and other cohorts, they do so with a less efficient assimilation strategy.

In contrast to UGA, leaf biomass was lower at the other three locations in the M-LPJG simulations compared to M-M. The larger total herbivore biomass in the M-LPJG runs was thus unexpected. Since the reported leaf biomass data represents the stock state after the cohorts feed from it, one might argue that the reduced leaf biomass stock simply reflects a larger degree of herbivory. However, this seemingly simple explanation does not capture the complexity of the underlying processes that affect the development of the emergent animal populations. In the following, we want to discuss the ecological responses to the coupling-related changes in total leaf biomass and varying climate conditions.

**Total leaf biomass.** is the chief alteration to the original Madingley model. All differences that can be observed between the M-M and M-LPJG setups originate from shifts in the simulated vegetation. In M-LPJG, more leaf biomass is available as a food source during the cold season compared to M-M. Thus, M-LPJG supports a larger biomass flux from autotrophs to herbivores over winter. The larger total cohort biomasses that correspond to the lower average leaf biomasses but higher minimum leaf biomasses (as computed in M-LPJG compared to M-M), indicates that total cohort biomasses are more sensitive to the minimum leaf biomass. This is plausible, since at least at locations with a temperate or boreal climate — the cohorts growth limitations are more likely to be reached during the cold season. In FIN M-LPJG, endotherms need significantly less time to reach maturity than in FIN M-M. This indicates that an increased minimum leaf biomass also increases the leaf biomass flux to herbivores on an individual level. Thus, the herbivores assimilate more biomass throughout the M-LPJG simulation and reach their adult body mass more quickly. A similar effect can be found in SAF M-LPJG, where a higher minimum in the leaf biomass pool results in a higher flux from autotrophic biomass to herbivores. This is quite remarkable since the average leaf biomass in SAF M-LPJG is lower than in SAF M-M and explains how despite the lower leaf biomass pool, the herbivore biomass pool can be higher in M-LPJG.

In three out of four locations, M-LPJG predicts a lower average annual leaf biomass when compared to M-M. This seems to shift the size distributions towards smaller individuals, most visible in the increased abundance of individuals lighter than 1 g (Fig. 6). This result is consistent with the increased number of cohorts going extinct in the M-LPJG simulations and their lower lifespans, which in turn reflects increased predation stress on the smaller individuals.

The larger total annual leaf biomass supports cohorts of heavier individual body mass in the M-M simulations. While this observation is consistent with the megafauna theory (Evans et al., 2012), the default Madingley version nonetheless overestimates the size of large herbivores by four to six times compared to field studies (Harfoot et al., 2014). It was also observed that the smaller the gap between the predicted M-M and the M-LPJG vegetation, the higher the herbivore biomass in the M-LPJG run (compare GER M-LPJG and SAF M-LPJG). This indicates that a uniform increase in vegetation biomass without changing its evergreen/deciduous composition still enhances herbivore biomass.

**Climatic conditions.** Under colder conditions, animals are less active and thus have less time to fulfil their metabolic cost (Harfoot et al., 2014). While ectotherms become inactive below a certain temperature threshold, endotherms constantly need to be active under tropical and boreal conditions alike. Climatic conditions thus contribute to the cohorts stress by enhancing the effect of calorific shortages.

Fig. 6 shows that the coupling leads to an increase in ectotherm herbivore abundance under warm climate (UGA and SAF). In these locations, the smallest individuals, which are typically ectotherm herbivores, show overall higher growth rates in M-LPJG, compared to M-M. Fig. 7B shows this effect also on an individual level. In M-M, the smallest individuals, which are typically ectotherms, are showing a wide pattern of growth rates in similar weight classes. This is an indicator that most individuals are not reaching their maximum potential growth rate, therefore the amount of leaf-biomass seems to be the limiting factor. In contrast, in M-LPJG, the smallest individuals all have similar growth rates under similar climatic conditions. This again indicates that the higher leaf-biomass no longer limits the growth rate.

Under colder climate conditions (FIN and GER) the ectotherm herbivore abundance did not show significant changes due to the coupling. In these locations, the growth rates of similar individuals are not enhanced by LPJ-GUESS vegetation input significantly. We expect the climatic conditions to be the limiting factors in both M-M and M-LPJG simulations.

These two examples have shown that climatic conditions can weaken or enhance the effect of changes in leaf-biomass. As the most noticeable example, we want to highlight the extinction of carnivore endotherms in FIN M-M simulations. During the first few timesteps, the endotherm herbivore population grows faster in FIN M-LPJG as a result of a higher minimum in leaf biomass during the cold season than in FIN M-M. This keeps the total biomass of herbivores in FIN M-LPJG large enough to feed endotherm carnivores in contrast to the FIN M-M simulations. The reduced growth of endotherm herbivores in FIN M-M leads to the total herbivore biomass being too low to support a carnivore endotherm population. As a result, endothermic carnivores in FIN M-M go extinct during the first few timesteps. Since there are no mechanisms that would allow these extinct functional group to enter the simulation, besides the initialisation, they stay extinct for the rest of the simulation. After the first decades, the established herbivore population grows without the top-down control of predation and thus their body mass size distribution shows fewer, heavier individuals than in the FIN M-LPJG simulation. This effect of the absence of large carnivores in Madingley, leading to the non-regulated growth of the lower trophic levels, is consistent with observations (Elmhagen et al., 2010; Brose et al., 2019) and highlights the importance of representing all trophic levels in a simulation.

Overall, the individual level characteristics only show minor differences when comparing the M-LPJG simulations to the M-M simulations, but, at the same time, we observe huge variations in the traits on a community level. This leads us to the conclusion that relatively small changes in individual processes can lead to substantial changes in the dynamics of a whole animal community. Minor differences on the individual level also indicate that Madingley still successfully predicts ecological processes for individuals and thus still produces reasonable assumptions for animal communities.

### 4.2. Comparison to external sources

M-LPJG consistently predicts NPP that is closer to empirical observations than the NPP predicted by M-M. Madingleys chosen NPP parameterisation (Miami) does not seem to fit well for temperature or precipitation extremes (very hot, very arid or very wet). In arid ecosystems such as SAF, M-LPJG also includes wildfires in its simulations. While the Miami model predicts a lower NPP than LPJ-GUESS in FIN and GER, we see higher leaf biomass stocks in M-M rather than in M-LPJG. This is caused by Madingley assuming that 64% (FIN) and 68% (GER) of the monthly available NPP is allocated towards leaves (see also supplementary information S3). LPJ-GUESS is typically accounting about 30% of the allocation towards leaves (De Kauwe et al., 2014).

The simulation derived power-law for herbivore biomass to NPP is more consistent with empirical relationships in the M-LPJG setup. Still, while the power-laws slope is similar, the magnitude of the simulated herbivore biomass is higher than seen in the empirical data from Cebrian (2004). There is an ongoing discussion in the scientific community about whether a higher biomass density of herbivores increases the ecosystems productivity through accelerating nutrient cycles (Enquist et al., 2020), significantly reduces plant biomass through damaging individuals (Jia et al., 2018) or shifts ecosystems plant species distributions (Schmitz et al., 2014). The strength of these effects likely are highly site and ecosystem specific. We aim to include the effects described by the scientific community in our future model development, such as implementing C:N stoichiometry in Madingley, implementing a feedback loop from Madingley to LPJ-GUESS, affecting leaf area and photosynthesis, and installing a litter pool in Madingley to track animal faeces. Such additional model development will allow to explore further which process underpins the overestimation of herbivore biomass in Fig. 9. Still, the Madingley model can, for the first time, predict similar kinds of power-laws as observable in nature, which is a major improvement of the model and lays the foundation for future model developments.

In SAF M-LPJG, both herbivore biomass and herbivory consumption are overestimated when compared to the power-law relationships. We expect this to result from every cohort having access to 10% of the leaf biomass stock in Madingley. While this value may be accurate for large insect swarms, we find it unlikely for a herd of mammals to access 10% of the leaf biomass in an area of ~2,500 km$^2$ (one grid cell in SAF) over the course of one month. While the cohorts are limited by the accessible leaf biomass, they do not necessarily consume all of it, but only try to fulfil their metabolic cost. In SAF, herbivores consume a higher portion of the standing stock due to the ecosystems lower productivity.

**Interpretation of previous Madingley publications.** We have shown that Madingleys build in NPP parameterisation by the Miami is only predicting realistic NPP in one out of the four simulated locations. Even in European temperate forest, NPP is underestimated by 30%. In ecosystems with extremely low or high productivity, it is very likely that Miami has made a completely unrealistic NPP prediction. We have also shown that not only the average NPP, but also the Miamis high yearly fluctuations of both evergreen and deciduous stocks more specifically the vegetation stocks winter minima and the length of the growing season- are very likely to affect the development of animal populations. In conclusion, we expect that previous publications have underestimated the abundance and growth rates of ectothermic animals in high or low-productivity ecosystems, while also overestimating the cohorts lifespan and reproduction cycles in all terrestrial ecosystems.

### 4.3. Limitations and future priorities

While we explicitly chose the four experiment locations to represent ecosystems under different climatic conditions, it remains to be tested how M-LPJGs leaf biomass simulation compares globally to that of M-M and whether this will align herbivore biomass better with observations. We assume that the improved modelling of the vegetation by LPJ-GUESS and thus shifts to more realistic canopy compositions enhance the realism of the ecological processes in Madingley. The next step will be the full coupling of both models and thus reduce the leaf mass from the individuals in LPJ-GUESS according to the herbivory activity in Madingley. During this process, we also want to make the leaf biomass accessible to the herbivores based on their traits, instead of giving every cohort access to 10% of the stocks. We expect this to also improve the representation of low-productivity ecosystems in out model system. In addition to the full coupling, we aim to implement a litter pool in Madingley and track animal faeces, which are currently escaping the model environment. There is an ongoing effort to include nitrogen in Madingley and to include this compartment of the code in the fully coupled version. This would enable further investigations if the model system can track an accelerated nutrient cycle and thus an enhanced ecosystem productivity.

Another plan we are pursuing is the implementation of different herbivory types, like grazing and browsing. A new version of Madingley is also being developed which would enable ground vs. canopy feeding.

### 4.4. Conclusion

The upcoming UN decade of restoration undoubtedly imposes huge challenges upon decision makers and emphasises the urgent need of understanding processes and interconnectivity in ecosystems. A large focus lies on nature-based solutions, and the increased interest in protecting animal population and biodiversity (European Commission, 2020). Modelling is a strong tool for making assessments on ecosystem functioning and interconnectivity of living organisms. With this study, we show that coupling the Madingley model with the process-based dynamic global vegetation model LPJ-GUESS is improving the realism of Madingleys modelled animal populations on multiple levels. Simulated animal cohorts were facing a wider range of net primary production and realistic fluctuations in seasonal vegetation biomass and composition. These changes lead to significant differences when compared to Madingleys default version, such as the persistence of whole ecological groups, which otherwise go extinct. We found general shifts in the animal populations towards smaller but more abundant individuals and demonstrated that our model system is consistently portraying and preserving all ecological groups, which is a major improvement. At the same time, the model still presents the underlying ecological processes, which is the foundation of the simulation of functional diversity. This ensures that that changes in the model output correspond to changes in the models boundary conditions and are in-fact not based on statistical disruptions in the food-chain or numeric instabilities.

Ultimately, we conclude that our implementations were in-fact improving the realism of the model systems prediction by comparing our results to empirical data. Madingley is a powerful tool in assessing an ecosystems functioning by modelling its whole trophic pyramid and we are keen to implement further processes and feedbacks in our model system in the future. With further development, out model system will help making informed decisions on ecosystem management in regards of animal biodiversity and overall ecosystem functioning.

### 4.5. Code availability

The LPJ-GUESS model code is managed and maintained by the Department of Physical Geography and Ecosystem Science at the Lund University, Sweden. The original Madingley model code is managed by the UN Environment Programme World Conservation Monitoring Centre (UNEP-WCMC) at the University of Cambridge, England. The intellectual property right for the translated version we used in this study is hold by the Radboud University in Nimwegen, Netherlands. Therefore, a DOI for both the LPJ-GUESS, as well as the Madingley model code cannot be provided. The source can be made available under a collaboration agreement under the acceptance of certain conditions.

---

## CRediT authorship contribution statement

Jens Krause: Conceptualization, Methodology, Writing – original draft, Visualization. Mike Harfoot: Conceptualization, Methodology, Software, Writing – review & editing. Selwyn Hoeks: Software, Writing – review & editing. Peter Anthoni: Conceptualization, Software, Validation, Writing – review & editing. Calum Brown: Writing – review & editing. Mark Rounsevell: Writing – review & editing. Almut Arneth: Validation, Writing – review & editing, Supervision, Funding acquisition.

## Declaration of competing interest

The authors declare that they have no known competing financial interests or personal relationships that could have appeared to influence the work reported in this paper.

## Data availability

The authors do not have permission to share data.

## Acknowledgements

Almut Arneth, Peter Anthoni and Jens Krause acknowledge funding via the Helmholz Foundation Impulse and Networking, Germany fund and the Helmholz ATMO programme, Germany. We also thank the reviewers for their constructive suggestion and the time they invested in supporting the review process of this study.

## Appendix A. Supplementary data

Supplementary material related to this article can be found online at https://doi.org/10.1016/j.ecolmodel.2022.110061. The supplementary material shows climate input, model setup, further simulation results and validation comparisons.

---

## References

Arneth, A., Shin, Y.-J., Leadley, P., Rondinini, C., Bukvareva, E., Kolb, M., Midgley, G.F., Oberdorff, T., Palomo, I., Saito, O., 2020. Post-2020 biodiversity targets need to embrace climate change. *Proc. Natl. Acad. Sci. USA* 117 (49), 30882–30891. http://dx.doi.org/10.1073/pnas.2009584117.

Bar-On, Y.M., Phillips, R., Milo, R., 2018. The biomass distribution on earth. *Proc. Natl. Acad. Sci. USA* 115 (25), 6506. http://dx.doi.org/10.1073/pnas.1711842115.

Berzaghi, F., Longo, M., Ciais, P., Blake, S., Bretagnolle, F., Vieira, S., Scaranello, M., Scarascia-Mugnozza, G., Doughty, C., 2019. Carbon stocks in central african forests enhanced by elephant disturbance. *Nature Geosci.* 12. http://dx.doi.org/10.1038/s41561-019-0395-6.

Brose, U., Archambault, P., Barnes, A.D., Bersier, L.-F., Boy, T., Canning-Clode, J., Conti, E., Dias, M., Digel, C., Dissanayake, A., Flores, A.A.V., Fussmann, K., Gauzens, B., Gray, C., Häussler, J., Hirt, M.R., Jacob, U., Jochum, M., Kéfi, S., McLaughlin, O., MacPherson, M.M., Latz, E., Layer-Dobra, K., Legagneux, P., Li, Y., Madeira, C., Martinez, N.D., Mendonça, V., Mulder, C., Navarrete, S.A., O'Gorman, E.J., Ott, D., Paula, J., Perkins, D., Piechnik, D., Pokrovsky, I., Raffaelli, D., Rall, B.C., Rosenbaum, B., Ryser, R., Silva, A., Sohlström, E.H., Sokolova, N., Thompson, M.S.A., Thompson, R.M., Vermandele, F., Vinagre, C., Wang, S., Wefer, J.M., Williams, R.J., Wieters, E., Woodward, G., Iles, A.C., 2019. Predator traits determine food-web architecture across ecosystems. *Nature Ecol. Evol.* 3 (6), 919–927. http://dx.doi.org/10.1038/s41559-019-0899-x.

Cardinale, B.J., Duffy, J.E., Gonzalez, A., Hooper, D.U., Perrings, C., Venail, P., Narwani, A., Mace, G.M., Tilman, D., Wardle, D.A., Kinzig, A.P., Daily, G.C., Loreau, M., Grace, J.B., Larigauderie, A., Srivastava, D.S., Naeem, S., 2012. Biodiversity loss and its impact on humanity. *Nature* 486 (7401), 59–67. http://dx.doi.org/10.1038/nature11148.

Cebrian, J., 2004. Role of first-order consumers in ecosystem carbon flow: Carbon flow through first-order consumers. *Ecol. Lett.* 7 (3), 232–240. http://dx.doi.org/10.1111/j.1461-0248.2004.00574.x.

Dangal, S.R.S., Tian, H., Lu, C., Ren, W., Pan, S., Yang, J., Di Cosmo, N., Hessl, A., 2017. Integrating herbivore population dynamics into a global land biosphere model: Plugging animals into the earth system. *J. Adv. Model. Earth Syst.* 9 (8), 2920–2945. http://dx.doi.org/10.1002/2016MS000904.

De Kauwe, M., Medlyn, B., Zaehle, S., Walker, A., Dietze, M., Wang, Y., Luo, Y., Jain, A., El Masri, B., Hickler, T., Warlind, D., Weng, E., Parton, W., Thornton, P., Wang, S., Prentice, I., Asao, S., Smith, B., Mccarthy, H., Norby, R., 2014. Where does the carbon go? A model-data intercomparison of vegetation carbon allocation and turnover processes at two temperate forest free-air CO$_2$ enrichment sites. *New Phytol.* 203. http://dx.doi.org/10.1111/nph.12847.

Elmhagen, B., Ludwig, G., Rushton, S.P., Helle, P., Lindén, H., 2010. Top predators, mesopredators and their prey: interference ecosystems along bioclimatic productivity gradients. *J. Anim. Ecol.* 79 (4), 785–794. http://dx.doi.org/10.1111/j.1365-2656.2010.01678.x.

Enquist, B.J., Abraham, A.J., Harfoot, M.B.J., Malhi, Y., Doughty, C.E., 2020. The megabiota are disproportionately important for biosphere functioning. *Nature Commun.* 11 (1), 699. http://dx.doi.org/10.1038/s41467-020-14369-y.

European Commission and Directorate-General for Research and Innovation, 2020. Biodiversity and Nature-Based Solutions: Analysis of EU-Funded Projects. Publications Office. http://dx.doi.org/10.2777/183298.

Evans, A.R., Jones, D., Boyer, A.G., Brown, J.H., Costa, D.P., Ernest, S.K.M., Fitzgerald, E.M.G., Fortelius, M., Gittleman, J.L., Hamilton, M.J., Harding, L.E., Lintulaakso, K., Lyons, S.K., Okie, J.G., Saarinen, J.J., Sibly, R.M., Smith, F.A., Stephens, P.R., Theodor, J.M., Uhen, M.D., 2012. The maximum rate of mammal evolution. *Proc. Natl. Acad. Sci. USA* 109 (11), 4187–4190. http://dx.doi.org/10.1073/pnas.1120774109.

Harfoot, M., Newbold, T., Tittensor, D.P., Emmott, S., Hutton, J., Lyutsarev, V., Smith, M.J., Scharlemann, J.P., Purves, D.W., 2014. Emergent global patterns of ecosystem structure and function from a mechanistic general ecosystem model. *PLoS Biol.* 12 (4), e1001841. http://dx.doi.org/10.1371/journal.pbio.1001841.

Harris, I.C., 2020. CRU JRA v2.1: a forcings dataset of gridded land surface blend of climatic research unit (CRU) and japanese reanalysis (JRA) data. Jan.1901 - Dec.2019. URL https://catalogue.ceda.ac.uk/uuid/10d2c73e5a7d46f4ada08b0a26302ef7.

Hatfield, J.L., Dold, C., 2019. Chapter 1 - photosynthesis in the solar corridor system. In: Deichman, C.L., Kremer, R.J. (Eds.), *The Solar Corridor Crop System*. Academic Press, pp. 1–33. http://dx.doi.org/10.1016/B978-0-12-814792-4.00001-2.

Hickler, T., Prentice, I., Smith, B., Sykes, M., Zaehle, S., 2006. Implementing plant hydraulic architecture within the LPJ dynamic global vegetation model. *Glob. Ecol. Biogeography* 15, 567–577. http://dx.doi.org/10.1111/j.1466-8238.2006.00254.x.

Hoeks, S., Huijbregts, M.A.J., Busana, M., Harfoot, M.B.J., Svenning, J.-C., Santini, L., 2020. Mechanistic insights into the role of large carnivores for ecosystem structure and functioning. *Ecography* 43 (12), 1752–1763. http://dx.doi.org/10.1111/ecog.05191.

Holdo, R.M., Sinclair, A.R.E., Dobson, A.P., Metzger, K.L., Bolker, B.M., Ritchie, M.E., Holt, R.D., 2009. A disease-mediated trophic cascade in the serengeti and its implications for ecosystem c. *PLoS Biol.* 7 (9), e1000210. http://dx.doi.org/10.1371/journal.pbio.1000210.

Hooper, D.U., Adair, E.C., Cardinale, B.J., Byrnes, J.E.K., Hungate, B.A., Matulich, K.L., Gonzalez, A., Duffy, J.E., Gamfeldt, L., O'connor, M.I., 2012. A global synthesis reveals biodiversity loss as a major driver of ecosystem change. *Nature* 486, 105–109. http://dx.doi.org/10.1038/nature11118.

Jia, S., Wang, X., Yuan, Z., Lin, F., Ye, J., Hao, Z., Luskin, M.S., 2018. Global signal of top-down control of terrestrial plant communities by herbivores. *Proc. Natl. Acad. Sci. USA* 115 (24), 6237. http://dx.doi.org/10.1073/pnas.1707984115.

Kattge, J., Díaz, S., Lavorel, S., Prentice, I.C., Leadley, P., Bönisch, G., Garnier, E., Westoby, M., Reich, P.B., Wright, I.J., Cornelissen, J.H.C., Violle, C., Harrison, S.P., Van Bodegom, P.M., Reichstein, M., Enquist, B.J., Soudzilovskaia, N.A., Ackerly, D.D., Anand, M., Atkin, O., Bahn, M., Baker, T.R., Baldocchi, D., Bekker, R., Blanco, C.C., Blonder, B., Bond, W.J., Bradstock, R., Bunker, D.E., Casanoves, F., Cavender-Bares, J., Chambers, J.Q., Chapin III, F.S., Chave, J., Coomes, D., Cornwell, W.K., Craine, J.M., Dobrin, B.H., Duarte, L., Durka, W., Elser, J., Esser, G., Estiarte, M., Fagan, W.F., Fang, J., Fernández-Méndez, F., Fidelis, A., Finegan, B., Flores, O., Ford, H., Frank, D., Freschet, G.T., Fyllas, N.M., Gallagher, R.V., Green, W.A., Gutierrez, A.G., Hickler, T., Higgins, S.I., Hodgson, J.G., Jalili, A., Jansen, S., Joly, C.A., Kerkhoff, A.J., Kirkup, D., Kitajima, K., Kleyer, M., Klotz, S., Knops, J.M.H., Kramer, K., Kühn, I., Kurokawa, H., Laughlin, D., Lee, T.D., Leishman, M., Lens, F., Lenz, T., Lewis, S.L., Lloyd, J., Llusià, J., Louault, F., Ma, S., Mahecha, M.D., Manning, P., Massad, T., Medlyn, B.E., Messier, J., Moles, A.T., Müller, S.C., Nadrowski, K., Naeem, S., Niinemets, U., Nöllert, S., Nüske, A., Ogaya, R., Oleksyn, J., Onipchenko, V.G., Onoda, Y., Ordoñez, J., Overbeck, G., Ozinga, W.A., Patiño, S., Paula, S., Pausas, J.G., Peñuelas, J., Phillips, O.L., Pillar, V., Poorter, H., Poorter, L., Poschlod, P., Prinzing, A., Proulx, R., Rammig, A., Reinsch, S., Reu, B., Sack, L., Salgado-Negret, B., Sardans, J., Shiodera, S., Shipley, B., Siefert, A., Sosinski, E., Soussana, J.-F., Swaine, E., Swenson, N., Thompson, K., Thornton, P., Waldram, M., Weiher, E., White, M., White, S., Wright, S.J., Yguel, B., Zaehle, S., Zanne, A.E., Wirth, C., 2011. TRY - a global database of plant traits. *Glob. Change Biol.* 17 (9), 2905–2935. http://dx.doi.org/10.1111/j.1365-2486.2011.02451.x.

Launiainen, S., Katul, G.G., Leppä, K., Kolari, P., Aslan, T., Grönholm, T., Korhonen, L., Mammarella, I., Vesala, T., 2022. Does growing atmospheric CO$_2$ explain increasing carbon sink in a boreal coniferous forest? *Glob. Change Biol.* n/a. http://dx.doi.org/10.1111/gcb.16117.

Lieth, H., 1975. Modeling the primary productivity of the world. In: Lieth, H., Whittaker, R.H. (Eds.), *Primary Productivity of the Biosphere*. Springer Berlin Heidelberg, pp. 237–263. http://dx.doi.org/10.1007/978-3-642-80913-2_12.

Lindeskog, M., Lagergren, F., Smith, B., Rammig, A., 2021. Accounting for forest management in the estimation of forest carbon balance using the dynamic vegetation model LPJ-GUESS (v4.0, r9333): Implementation and evaluation of simulations for Europe. *Geosci. Model Dev. Discuss.* 2021, 1–42. http://dx.doi.org/10.5194/gmd-2020-440.

Newbold, T., Tittensor, D.P., Harfoot, M.B.J., Scharlemann, J.P.W., Purves, D.W., 2020. Non-linear changes in modelled terrestrial ecosystems subjected to perturbations. *Sci. Rep.* 10 (1), 14051. http://dx.doi.org/10.1038/s41598-020-70960-9.

Nichols, E., Peres, C.A., Hawes, J.E., Naeem, S., 2016. Multitrophic diversity effects of network degradation. *Ecol. Evol.* 6 (14), 4936–4946. http://dx.doi.org/10.1002/ece3.2253.

Pachzelt, A., Forrest, M., Rammig, A., Higgins, S., Hickler, T., 2015. Potential impact of large ungulate grazers on african vegetation, carbon storage and fire regimes: Grazer impacts on african savannas. *Glob. Ecol. Biogeography* 24. http://dx.doi.org/10.1111/geb.12313.

Pacifici, M., Foden, W.B., Visconti, P., Watson, J.E., Butchart, S.H., Kovacs, K.M., Scheffers, B.R., Hole, D.G., Martin, T.G., Akçakaya, H.R., Corlett, R.T., Huntley, B., Bickford, D., Carr, J.A., Hoffmann, A.A., Midgley, G.F., Pearce-Kelly, P., Pearson, R.G., Williams, S.E., Willis, S.G., Young, B., Rondinini, C., 2015. Assessing species vulnerability to climate change. *Nature Clim. Change* 5 (3), 215–225. http://dx.doi.org/10.1038/nclimate2448.

Rabin, S.S., Melton, J.R., Lasslop, G., Bachelet, D., Forrest, M., Hantson, S., Kaplan, J.O., Li, F., Mangeon, S., Ward, D.S., Yue, C., Arora, V.K., Hickler, T., Kloster, S., Knorr, W., Nieradzik, L., Spessa, A., Folberth, G.A., Sheehan, T., Voulgarakis, A., Kelley, D.I., Prentice, I.C., Sitch, S., Harrison, S., Arneth, A., 2017. The fire modeling intercomparison project (firemip), phase 1: experimental and analytical protocols with detailed model descriptions. *Geosci. Model Dev.* 10 (3), 1175–1197. http://dx.doi.org/10.5194/gmd-10-1175-2017.

Rebmann, C., Zeri, M., Lasslop, G., Mund, M., Kolle, O., Schulze, E.-D., Feigenwinter, C., 2010. Treatment and assessment of the CO$_2$-exchange at a complex forest site in Thuringia, Germany. *Agric. Forest Meteorol.* 150 (5), 684–691. http://dx.doi.org/10.1016/j.agrformet.2009.11.001.

Running, S., Zhao, M., 2021. MODIS/Terra net primary production gap-filled yearly L4 global 500 m SIN grid V061. 2021. URL https://doi.org/10.5067/MODIS/MOD17A3HGF.061.

Schmitz, O.J., Raymond, P.A., Estes, J.A., Kurz, W.A., Holtgrieve, G.W., Ritchie, M.E., Schindler, D.E., Spivak, A.C., Wilson, R.W., Bradford, M.A., Christensen, V., Deegan, L., Smetacek, V., Vanni, M.J., Wilmers, C.C., 2014. Animating the carbon cycle. *Ecosystems* 17 (2), 344–359. http://dx.doi.org/10.1007/s10021-013-9715-7.

Schmitz, O.J., Wilmers, C.C., Leroux, S.J., Doughty, C.E., Atwood, T.B., Galetti, M., Davies, A.B., Goetz, S.J., 2018. Animals and the zoogeochemistry of the carbon cycle. *Science* 362 (6419), eaar3213. http://dx.doi.org/10.1126/science.aar3213.

Smith, M.J., Purves, D.W., Vanderwel, M.C., Lyutsarev, V., Emmott, S., 2013. The climate dependence of the terrestrial carbon cycle, including parameter and structural uncertainties. *Biogeosciences* 10 (1), 583–606. http://dx.doi.org/10.5194/bg-10-583-2013.

Smith, B., Wärlind, D., Arneth, A., Hickler, T., Leadley, P., Siltberg, J., Zaehle, S., 2014. Implications of incorporating n cycling and n limitations on primary production in an individual-based dynamic vegetation model. *Biogeosciences* 11 (7), 2027–2054. http://dx.doi.org/10.5194/bg-11-2027-2014.

Sobral, M., Silvius, K.M., Overman, H., Oliveira, L.F.B., Raab, T.K., Fragoso, J.M.V., 2017. Mammal diversity influences the carbon cycle through trophic interactions in the amazon. *Nature Ecol. Evol.* 1 (11), 1670–1676. http://dx.doi.org/10.1038/s41559-017-0334-0.

Staver, A.C., Bond, W.J., 2014. Is there a 'browse trap'? Dynamics of herbivore impacts on trees and grasses in an African savanna. *J. Ecol.* 102, 595–602. http://dx.doi.org/10.1111/1365-2745.12230.

Wilmers, C.C., Schmitz, O.J., 2016. Effects of gray wolf-induced trophic cascades on ecosystem carbon cycling. *Ecosphere* 7 (10), e01501. http://dx.doi.org/10.1002/ecs2.1501.

Wramneby, A., Smith, B., Zaehle, S., Sykes, M., 2008. Parameter uncertainties in the modelling of vegetation dynamics—Effects on tree community structure and ecosystem functioning in European forest biomes. *Ecol. Modell.* 216, 277–290. http://dx.doi.org/10.1016/j.ecolmodel.2008.04.013.

Wårlind, D., Smith, B., Hickler, T., Arneth, A., 2014. Nitrogen feedbacks increase future terrestrial ecosystem carbon uptake in an individual-based dynamic vegetation model. *Biogeosciences* 11 (21), 6131–6146. http://dx.doi.org/10.5194/bg-11-6131-2014.

---

## Corrigendum

*Ecological Modelling* 492 (2024) 110706  
Available online 12 April 2024  
https://doi.org/10.1016/j.ecolmodel.2024.110706

**Corrigendum to "How more sophisticated leaf biomass simulations can increase the realism of modelled animal populations" [Ecological Modelling 471 (2022) 110061]**

Jens Krause$^{a,*}$, Mike Harfoot$^{b}$, Selwyn Hoeks$^{c}$, Peter Anthoni$^{a}$, Calum Brown$^{a}$, Mark Rounsevell$^{a}$, Almut Arneth$^{a}$

$^{a}$ KIT-Campus Alpin, Institute of Meteorology and Climate Research (IMK-IFU), Garmisch-Partenkirchen, Germany  
$^{b}$ UN Environment Programme World Conservation Monitoring Center, Cambridge, United Kingdom  
$^{c}$ Department of Environmental Science, Radboud University Nijmegen, Netherlands

The authors regret that the data presented in Fig. 9 were incorrect. We created an algorithm to estimate each grid cell's area. This algorithm contained an error, which led to an overestimation of the calculated area and thus an underestimation of the herbivore biomass density and the herbivory consumption density. We corrected the algorithm and used the same model state and setup as in our original manuscript and repeated the simulations used to determine the power-law relationships. We expanded the analysis time-frame from 10 to 30 years to ensure that low-productivity ecosystems are represented well. The corrected power-law relationships are:

$$\log(biomass) = 0.9 + 0.31 \cdot \log(NPP)$$

$$\log(primary\,consumption) = 0.26 + 0.66 \cdot \log(NPP)$$

The corrected Fig. 9 is the following: [Fig. 1 in the corrigendum — Power-law relationships. Blue triangles represent an annual sum during a M-LPJG simulation. Green dots represent an annual sum during a M-M simulation. For the M-LPJG setups, herbivore biomass and herbivory consumption are related to NPP. There are no reasonable corresponding logarithmic fits for the M-M setups.]

Our original claims still uphold with the corrected Fig. 9. Our model system shows a positive response of herbivore biomass density to the ecosystem's NPP, while the default Madingley model version is showing a narrow range of NPP data points, and lacks a logarithmic relationship. With the corrected grid cell area algorithm, the derived power-law relationship for herbivory consumption is showing a very close fit to (Cebrian, 2004).

The authors would like to apologise for any inconvenience caused.

### Reference

Cebrian, J., 2004. Role of first-order consumers in ecosystem carbon flow: Carbon flow through first-order consumers. *Ecology Letters* 7, 232–240. https://doi.org/10.1111/j.1461-0248.2004.00574.x.
