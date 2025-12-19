This repository contains data for the manuscript, **Abundance redistribution increases predator-prey interaction potentials among North American birds**, by Heng-Xing Zou, Chia Hsieh, Phoebe L. Zarnetske, Kai Zhu, and Brian C. Weeks. Heng-Xing Zou curated and analyzed the data using the code provided in the repository.

# Overview

This repository contains the following directories:

`Code` contains all code used for curating the data and the analyses.
`Data` contains all data used for and generated from the analyses.
`Figures` contains all figures generated from the analyses, including those in the manuscript.

**Currently, `Data` folder is empty because some datasets are large and cannot be uploaded to github. We can provide the data as needed before depositing it permanently upon publication.**

# Data Source

**Bird Species Abundances.** Species abundance data were obtained from the North American Breeding Bird Survey (publicly available at the [USGS website](https://www.usgs.gov/centers/eesc/science/north-american-breeding-bird-survey); also available to download with the package [`bbsBayes2`](https://bbsbayes.github.io/bbsBayes2/index.html)). The analysis uses the 2022 release of the data. See `Database/BBS_List_PIF.csv` for the full species list. The PIF adjustment factors correct species abundances based on their detectability under the NABBS protocol. See [Zou et al. 2025 bioRXiv](https://www.biorxiv.org/content/10.1101/2025.07.23.666479v2.abstract) for more information on the model fitting; see also the [github resository](https://github.com/hengxingzou/Zouetal2025bioRXiv) for that preprint.

**Functional Traits.** Functional trait data were obtained from various sources, including [AVONET](https://doi.org/10.1111/ele.13898), [EltonTraits 1.0](https://doi.org/10.1890/13-1917.1), [Myhrvold et al.](https://doi.org/10.1890/15-0846R.1), [Chia et al.](https://www.nature.com/articles/s41597-023-02837-1), and [Bird et al.](https://doi.org/10.1111/cobi.13486). We did not use all functional traits but kept them in the data file for code compatibility. See the manuscript for the traits used and appropriate citations.

**Bioclimatic Data.** Bioclimatic data were obtained from [WorldClim2](https://rmets.onlinelibrary.wiley.com/doi/10.1002/joc.5086) and processed to the resolution of one by one degree. 

**Human-induced Landscape Data.** Human-induced landscape data were obtained from [History Database of the Global Environment (HYDE)](https://www.pbl.nl/en/publications/new-anthropogenic-land-use-estimates-for-the-holocene-hyde-32) and processed to the resolution of one by one degree.

# Files in `Data`

`Relative_Abundance.csv` (large file) contains relative abundance data for all birds for all grids from 1970 to 2021, including all 470 species in the study. 

`Anthromes_Grid_Years_Prop.csv` is processed from the HYDE data on human-induced landscapes within each one by one grid. Each landscape type (column `Anthrome`) is calculated as proportions of the total land area of the grid (column `AREA_SQ_KM`), in the column `Coverage`. Data was processed from openly available data with the script `2_GenerateCommunityData/4_Process_HYDE_Data.R`.

`Bioclim.csv` is processed from the WorldClim2 data on bioclimatic variables within each one by one grid. All columns except `year` and `ST_12` (name of the grid) hold bioclimatic variables. Data was processed from openly available data with the script `2_GenerateCommunityData/5_Process_Climate_Data.R`.

`All_Funct_Data.csv` contains functional data from various sources. Columns `Beak.Length_Culmen` to `Primary.Lifestyle` are from [AVONET](https://doi.org/10.1111/ele.13898). Columns `Diet.Inv` to `Nocturnal` are from [EltonTraits 1.0](https://doi.org/10.1890/13-1917.1). Column `litter_or_clutch_size_n` is from [Myhrvold et al.](https://doi.org/10.1890/15-0846R.1). Columns `Parasite` to `NestAtt_pensile` are from [Chia et al.](https://www.nature.com/articles/s41597-023-02837-1). Column `GenLength` is from [Bird et al.](https://doi.org/10.1111/cobi.13486). See the original reference for definitions and details. See the manuscript for details on the rest of the columns.

`BBS_List_Num_Records.csv` is the complete species list from the BBS, with the common name (`english`), scientific name (`latin`), AOU number (`aou`), and number of records in BBS (`num_records`) for each species. 

`Filtered_Spp_List.csv` is the species list in the study filtered by data availability, with the common name (`english`), scientific name (`latin`), AOU number (`aou`), and the PIF adjustment factor (`PIF_adj_birds_km2`) for each species. 

`OTT_crosswalk_2023.csv` is the phylogeny obtained from [McTavish et al. 2025 PNAS](https://github.com/AvesTree/AvesData). Note that we are using the 2023 version of Clements taxonomy because it was the latest stable version at the time of the analyses.

# Code

`1_Process_Data_Hurlbert.R` processes predator-prey interaction data from the [Avian Diet Database](https://github.com/ahhurlbert/aviandietdb) by filtering and aligning species.

`2_Get_Similarities.R` calculates the functional simularities between recorded and unrecorded prey species. Calls `calculate_similarity.R`, which contains the function.

`3_Calculate_Intrinsic_IP.R` calculates the intrinsic interaction potential between a predator and prey pair.

`4_Calculate_IP.R` calculates the interaction potential, which also includes relative abundances of both predator and prey in addition to the intrinsic interaction potential.

`5_Analyze_Pair.R` analyzes predator-prey pairs.

`6_Analyze_by_Prey.R` analyzes results by prey functional groups.

`7_Analyze_Lat.R` analyzes latitudinal shifts of predators, prey, and their interactions.

