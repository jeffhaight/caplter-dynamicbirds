# caplter-dynamicbirds

# A repository for:
Haight, Jeffrey D.; Bateman, Heather L.; Frazier, Amy E.; Larson, Kelli L.; Susannah B. Lerman; Paige S. Warren; Albuquerque, Fabio S. Hot in the city: impacts of climate variability and urban landscape change on birds

This repository contains code and data used for fitting a temporally autoregressive community occupancy model to bird occurrence data across multiple spring sampling seasons, analyzing model parameters, and producing visualizations based on these results. Occurrence data were collected at 51 long-term bird point-count censusing locations across central Arizona, USA, as part of the Central Arizona-Phoenix Long-term Ecological Research (CAP LTER) program.


This `README` file includes information on the various scripts and datasets used for this analysis, though not every data source is saved in this repository due to data size limitations (the fitted occupany model object; GIS data). The manuscript includes citation of sources of environmental data. This `README` file is organized into three sections corresponding to the three main folders into which the repository is organized
1. [code](#code)  
2. [data](#data)
3. [figures](#figures)

Please direct all questions to Jeffrey Haight jdhaight.eco(at)gmail.com



---
---
<div align="center"> <h3>code</h3> </div>

This folder contains all .R script and RMarkdown (.rmd) files used for conducting the analyses that produce the `model_outputs` and `model_summaries` data and the `figures`. There are a total of eight R scripts, numbered by the order in which they are to be run. Scripts 1_0 and 1_1 are for preparing data for modeling and examining environmental relationships, 2_0 through 2_4 are for conducting the analyses, and 3_0 through 3_2 are for examining model outputs and producing figures components.

| File  | Description  |
|---|---|
|**./code/1_DataPrep_SpeciesObservations_caplter-dynamicbirds.rmd**   | Code for cleaning and organizing [CAP LTER bird monitoring data](https://doi.org/10.6073/pasta/1d1badb48df8a5c5fcfafb3275713bc1), in preparation for the fitting of occupancy models |
|**./code/2_DataPrep_SpeciesTraits_caplter-dynamicbirds.rmd**   | Code for joining lists of Central Arizona bird species lists with species trait datasets: [EltonTraits 1.0 (Wilman et al. 2016)](https://doi.org/10.6084/m9.figshare.c.3306933), [AVONET (Tobias et al. 2022)](https://doi.org/10.1111/ele.13898), and State of Arizona's [State Wildlife Action Plan: 2012-2022](https://www.azgfd.com/wildlife-conservation/on-the-ground-conservation/state-wildlife-action-plan/) and its corresponding [Heritage Data Management System](https://www.azgfd.com/wildlife-conservation/on-the-ground-conservation/cooperative-programs/az-natural-heritage-program/) |
|**./code/3_ModelDataFormatting_caplter-dynamicbirds.rmd**   | Code for filtering species and environmental data down to the spatiotemporal scope of the analysis and examining relationships between and long-term trends in environmental conditions|
|**./code/4_CommunityOccupancyModel_caplter-dynamicbirds.R**   | Script for using JAGS to fit the global Bayesian autologistic community occupancy model containing scale-optimized covariates within R, using the R package `jagsUI`|
|**./code/5_OutputAnalysisandVisualization_caplter-dynamicbirds.R**  | Script for summarizing species-environment relationships and creating figures illustrating those relationships, as well as among-species variation in relationships paritally associated with functional traits |
|**./code/communityoccupancy_autologistic.R**   |  The final, published autologistic community occupancy model that we fit to multi-year species occurrence data   |
|**./code/communityoccupancy_autologistic_original.R**   |  The original (pre-peer review) version of the autologistic community occupancy model that we fit to multi-year species occurrence data   |

Code was run using the following R packages:  
| Package  | Version  |
|---|---|
| tidyverse | 2.0.0 |
| sf        | 1.1-0 |
| abind | 1.4-8 |
| lubridate | 1.9.5 |
| jagsUI | 1.6.3 |
| RColorBrewer| 1.1-3 |
| png | 0.1-9 |
| gghighlight | 0.4.1 |
| ggpubr| 0.6.3 |
| ggcorrplot | 0.1.4.1 |
| GGally | 2.2.1 |
| viridis| 0.6.5 |
| beepr | 2.0 |
| scales | 1.3.0 |
| ggridges | 0.5.7 |



---
<div align="center"> <h3>data</h3> </div>
This contains two subfolders with all the data files used for conducting the analyses and producing the `figures`:   
**`input`**  
and  
**`output`**    

 

Within **`input`**, there are several files:  

**./data/input/birdtraits_corebirds2025.csv**    
Contains data of species-level attributes, including each species' names, taxonomic groupings, and species traits. Created through the running of the **./code/2_DataPrep_SpeciesTraits_caplter-dynamicbirds.rmd** script.


**./data/input/core_birds_countmeanbyseason.csv** 
A table containing multi-year average counts (`abundance_mean`) and summarized detections of all identified and unidentified species 

**./data/input/core_birds_obs_countbysurvey.rds** 
By-survey bird count data; a five-dimensional (site X species X survey X year X season) array of species counts, including both identified/known and unidentified/unknown species  

**./data/input/core_birds_obs_surveycount.rds**
Dataset representing the number of days on which each site was surveyed by point count during each season  


**./data/input/indices_Landsat_corebirds_mean1000m_ISadjusted_pc1.csv**
A dataset of impervious surface values derived from [land cover-adjusted rasters of the Enhanced Normalized Difference in Impervious Surface Index (ENDISI)](https://doi.org/10.6073/PASTA/147EFCA87F731ED8F16433B5F692AF30) around each bird sampling site, following the methods of related [Haight et al. 2024 environmental summaries for the bird monitoring locations](https://doi.org/10.6073/PASTA/9D44CD85F881586D6D06E7A7293E833C), with which it was merged.


**./data/input/birdtraits_corebirds2025.csv** 
A list of bird species (and unidentified "species") observed across all bird monitoring locations and years (2000-2024)


**./data/input/modelinputs_dynamicbirds2025.RData**
R Studio environment containing all R objects for running the autologistic community occupancy analysis. These data were compiled using the **./code/3_ModelDataFormatting_caplter-dynamicbirds.rmd** script, and they include:  
 

| Object Name	| Description   |
|---------------------------|--------|
| a.cov.det| A three-dimensional (site X year X intercept/covariate) array storing standardized detection covariates, with the first (and only) set of values containing a '1' corresponding to the intercept|
| a.cov.occ| A three-dimensional (site X year X intercept/covariate) array storing standardized occupancy covariates, with the first set of values containing a '1' corresponding to the intercept|
| bird.counts| Bird detection/non-detection data stored in **./data/input/core_birds_obs_countbysurvey.rds** ; a five-dimensional (site X species X survey X year X season) array of species counts, including both identified/known and unidentified/unknown species|
| covariates.det | Names corresponding to the model intercept and covariates for detection probability|
| covariates.occ | Names corresponding to the model intercept and covariates for occupancy|
| data.spp		| A dataframe of traits for all observed bird species	|
| data.spp.model		| A dataframe of traits for the 139 modeled bird species	|
| env.model	| A dataframe containing environmental variables at 51 bird monitoring sites, summarized across five timeframes (four seasons and annually)|
| env.spr	| A dataframe containing environmental variables at 51 bird monitoring sites during 15 spring monitoring seasons (included in analysis)|
| env.win	| A dataframe containing environmental variables at 51 bird monitoring sites during 15 winter monitoring seasons (not included in analysis)|
| family_vec| A numeric index corresponding to the taxonomic family of each modeled bird species; not utilized in the final analysis|
| K 	         | A matrix including the number of days on which each site was surveyed by point count during each season; also stored in **./data/input/core_birds_obs_surveycount.rds** |
| K.model      | A matrix including the number of days on which each site was surveyed by point count during 15 spring seasons |
| n.season		| The number of sampling "seasons" in the occupancy model i.e., survey years|
| n.site		  | The number of sites surveyed in each year|
| n.spp		    | The number of modeled species observed across all years	|
| n.survey		| The number of survey occasions per year	|
| names.common		| A vector of common names for each modeled species	|
| obs.model| A dataframe containing multi-year average counts (`abundance_mean`) and summarized detections of all identified and unidentified species; a subset of **./data/input/core_birds_countmeanbyseason.csv** |
| order_vec| A numeric index corresponding to the taxonomic order of each modeled bird species; not utilized in the final analysis|
| spp.known| Four-letter alpha codes corresponding to all bird species observed and identified bird species across all surveys years and seasons (source: Institute for Bird Populations) | 
| spp.model| Four-letter alpha codes corresponding to all modeled bird species; identical to `spp.spring` (source: Institute for Bird Populations) |
| spp.spring| Four-letter alpha codes corresponding to all bird species observed and identified bird species across all surveys years in the spring season (source: Institute for Bird Populations)|
| spp.unknown| Four-letter alpha codes representing all observed but generically classified "Unknown" bird species (e.g., "UNHA" = "unknown hawk")|
| spp.winter| Four-letter alpha codes corresponding to all bird species observed and identified bird species across all surveys years in the winter season (source: Institute for Bird Populations)|
| spp.model		| A vector of shortened common names for each modeled species, the names by which species occurrence data was alphabetized	|
| t.families| a list of the taxonomic families of all modeled bird species; not utilized in the final analysis|
| t.orders| a list of the taxonomic orders of all modeled bird species; not utilized in the final analysis|
| traits.response | A matrix of standardized species functional traits used as covariates in the occupancy model|
| y	|A five-dimensional (site X species X survey X year X season) array of species detections (1 = detected, 0 = not detected), including only identifed/known species|
| year_vec| A numeric index of all survey years|
| ysum	|A four-dimensional (site X species X year X season) array of species detections summed across surveys (max value = n.survey = 6) |
| ysum.spr	|A three-dimensional (site X species X year) array of species detections summed across surveys during the spring season|
| ysum.win	|A three-dimensional (site X species X year) array of species detections summed across surveys during the winter season|
| ysum.model	| Identical to `ysum.spr`; the occurrence dataset used to fit the community occupancy model model |  


Within **`output`**, there are two file:
**"~/projects/DynamicBirds/data/output/MSOM_CAPbirds_spring_summary40ksamp.csv"**
Contains the JAGS summary of the model occupancy model parameters  

**"~/data/output/MSOMoutput_specieseffectsummary.csv"**
Contains the summary of the model occupancy model parameters, combined with species information  


Please note that file containing the full JAGS output with the fitted occupancy model **"~/projects/DynamicBirds/data/output/MSOM_CAPbirds_spring_40ksamp.rds"** is not included in this repository, as the file is too large. This file can be replicated by fitting the community occupancy model code using the provided R script (**./code/4_CommunityOccupancyModel_caplter-dynamicbirds.R**) and input files, and can also be shared by the lead author via an email request to jdhaight.eco(at)gmail.com.  
  

---
<div align="center"> <h3>figures</h3> </div>

This folder contains all images utilized in the production of the manuscript figures and tables. This folder includes visualizations produced using the above R scripts (see **code**). These visualizations were then combined within with one another within Inkscape (https://inkscape.org/) to create the figures in the published manuscript. 

The `figures` folder also includes the subfolder **./figures/bird_images** which contains image files used to key bird species in Figure 4 of the manuscript. All bird graphics were sourced from PhyloPic (https://www.phylopic.org/) and were utilized as part of the public domain (https://creativecommons.org/publicdomain/zero/1.0/). Bird graphics are accompanied by the text file **./figures/mammalgraphics/ImageCopyrights_BirdSilhouettes.txt**, which specifies each image's source and provides attribution for each image, as outlined below.
| File Name  | Creative Commons License | Attribution  |  Source  |
|---|---|---|---|
| Anas crecca  | CC0 1.0 Universal Public Domain Dedication   | Andy Wilson  |  https://www.phylopic.org/images/be6f4ddf-d1ae-4cfd-bc67-0fd1a10413ea/quiscalus-mexicanus  |
| Bald Eagle Haliaeetus leucocephalus  | CC0 1.0 Universal Public Domain Dedication   | Andy Wilson  |  https://www.phylopic.org/images/3002e4de-eb9b-48a4-8921-57dfe0d2e169/haliaeetus-leucocephalus  |
| Callipepla  | CC0 1.0 Universal Public Domain Dedication   | Ferran Sayol  |  https://www.phylopic.org/images/48053987-63c1-47bc-87c3-2abf84480a6f/callipepla-californica  |
| Campylorhyncus rufinucha  | CC0 1.0 Universal Public Domain Dedication   | Andy Wilson  |  https://www.phylopic.org/images/b75457e6-e19d-43e5-96c4-ba08404930c9/campylorhynchus-rufinucha  |
| hummingbird Trochilidae  | CC0 1.0 Universal Public Domain Dedication   | Ferran Sayol  |  https://www.phylopic.org/images/2bf1e800-5384-45cd-a533-ac940b8eadd6/trochilidae  |
| Mourning Dove Zenaida macroura  | CC0 1.0 Universal Public Domain Dedication   | Andy Wilson  |  https://www.phylopic.org/images/fc5ff059-bb2d-4e91-8422-9039ef106979/zenaida-macroura  |
| Quiscalus mexicanus  | CC0 1.0 Universal Public Domain Dedication   | Andy Wilson  |  https://www.phylopic.org/images/be6f4ddf-d1ae-4cfd-bc67-0fd1a10413ea/quiscalus-mexicanus  |

---
 